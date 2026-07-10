# generation.jl
#
# Spatially- and spectrally-resolved carrier generation rate G(x, λ) for
# a multilayer optical stack, computed by post-processing coh_tmm output.
#
# Direct translation of Steven Byrnes' tmm.tmm_core absorption machinery:
#   * absorp_analytic_fn          (Python class)
#   * position_resolved()         (point evaluation: poyn, absor, E-field)
#   * absorp_in_each_layer()      (integrated absorption per layer)
#
# Plus a Julia-native wrapper to assemble the [N_x × N_λ] generation matrix
# used by the EQE pipeline.
#
# Units convention (per Day 4 agreement):
#   x in nm, λ in nm, G in (carriers / nm³ / s) per (incident-photon-flux unit)
#   — TMM-native everywhere. We convert ONLY at the ChargeTransport interface.
#
# Physics convention:
#   For polarised light, McByrne's `absor` returns absorbed power per unit
#   LENGTH along x, as a fraction of incident power crossing the front face.
#   To get a per-volume rate you'd divide by the area-element; here, because
#   we incident a planewave (1D), `absor` IS already the per-depth absorption
#   rate of the incident intensity. See arXiv:1603.02720 §3 for the derivation.

# NOTE: tmm_core.jl (coh_tmm, list_snell, EPSILON, R_from_r, ...) is included
# by the top-level module PerovskiteTMMlite.jl BEFORE this file, so it is already
# in scope. If you `include("generation.jl")` directly (outside the module),
# make sure you `include("tmm_core.jl")` first.

# -----------------------------------------------------------------------------
# AbsorpAnalyticFn  -- analytic absorption-vs-depth inside one coherent layer
# -----------------------------------------------------------------------------
# McByrne shows (arXiv:1603.02720 §3) that absorption at depth z inside a
# coherent layer has the closed form
#
#     a(z) = A1 exp(+a1 z) + A2 exp(-a1 z)
#          + A3 exp(+i a3 z) + conj(A3) exp(-i a3 z)
#
# where A1, A2 are real, A3 is complex, and a1, a3 are real positive numbers.
# Storing (A1, A2, A3, a1, a3, d) is enough to evaluate a(z) at any depth, and
# to do algebra on multiple stacks (sum, flip, scale) before evaluating.
#
# Julia syntax notes:
#  * `mutable struct` — like `struct`, but fields can be reassigned. We need
#    this because `fill_in!`, `flip!`, `scale!`, `add!` rewrite the fields.
#    By convention, methods that mutate end with `!`.
#  * Field types are concrete (Float64 / ComplexF64) so the struct has a
#    known memory layout. No abstract `Number` field types.
mutable struct AbsorpAnalyticFn
    A1::Float64
    A2::Float64
    A3::ComplexF64
    a1::Float64
    a3::Float64
    d ::Float64        # layer thickness (nm)
end

# Default constructor — undefined state, to be filled by fill_in!.
# Julia note: `AbsorpAnalyticFn()` is a parameterless outer constructor. NaNs
# in a default-constructed instance make sure that using one without filling
# it in produces obviously-broken output rather than silent zeros.
AbsorpAnalyticFn() = AbsorpAnalyticFn(NaN, NaN, NaN + NaN*im, NaN, NaN, NaN)

"""
    fill_in!(f::AbsorpAnalyticFn, coh_data::Dict, layer::Int) -> f

Populate `f` from the output of `coh_tmm` (=`coh_data`) so it represents the
absorption in `layer` (a 1-based layer index, like everything else in Julia).

`layer` must be an INTERIOR coherent layer: 2 ≤ layer ≤ length(d_list)-1.

Returns `f` (so calls chain).
"""
function fill_in!(f::AbsorpAnalyticFn, coh_data::Dict, layer::Int)
    pol  = coh_data["pol"]
    v    = coh_data["vw_list"][layer, 1]    # forward amplitude in layer
    w    = coh_data["vw_list"][layer, 2]    # backward amplitude in layer
    kz   = coh_data["kz_list"][layer]
    n    = coh_data["n_list"][layer]
    n_0  = coh_data["n_list"][1]
    th_0 = coh_data["th_0"]
    th   = coh_data["th_list"][layer]
    f.d  = coh_data["d_list"][layer]

    # Real/imag parts of 2 kz. a1 is the evanescent decay rate of |Ef|², |Eb|².
    # a3 is the spatial-beating frequency between forward and backward waves.
    f.a1 = 2 * imag(kz)
    f.a3 = 2 * real(kz)

    if pol == "s"
        # Re-derived in arXiv:1603.02720 eq. (53)
        temp = imag(n * cos(th) * kz) / real(n_0 * cos(th_0))
        f.A1 = temp * abs2(w)
        f.A2 = temp * abs2(v)
        f.A3 = temp * v * conj(w)
    elseif pol == "p"
        # eq. (54) — A1, A2 still come from imag(kz), A3 from real(kz).
        temp = 2 * imag(kz) * real(n * cos(conj(th))) /
               real(n_0 * conj(cos(th_0)))
        f.A1 = temp * abs2(w)
        f.A2 = temp * abs2(v)
        # Note: in the p case, A3's prefactor has a DIFFERENT structure — the
        # cross-term picks up real(kz) and imag of the n cos term. Easy to
        # mistranscribe; matches Python tmm_core.py line ~412.
        f.A3 = v * conj(w) * (-2 * real(kz) * imag(n * cos(conj(th))) /
                              real(n_0 * conj(cos(th_0))))
    else
        error("Polarization must be \"s\" or \"p\"")
    end
    return f
end

"""
    evaluate(f::AbsorpAnalyticFn, z::Real) -> Float64

Absorbed power per unit length at depth `z` (z=0 = start of the layer),
as a fraction of incoming light power. The result is real-valued by
construction (A3 and its conjugate cancel imaginary parts).

I named this `evaluate` rather than `run` (McByrne's Python name) because
`run` collides with `Base.run` (used to launch subprocesses).
"""
function evaluate(f::AbsorpAnalyticFn, z::Real)
    return f.A1 * exp(+f.a1 * z) + f.A2 * exp(-f.a1 * z) +
           real(f.A3 * exp(+1im * f.a3 * z) +
                conj(f.A3) * exp(-1im * f.a3 * z))
end

# Make the struct callable: f(z) ≡ evaluate(f, z).
# Julia syntax: defining a method on the type itself overloads `f(args...)`.
# This is the same trick Python uses with `__call__`.
(f::AbsorpAnalyticFn)(z::Real) = evaluate(f, z)

"""
    flip!(f::AbsorpAnalyticFn) -> f

In-place: replace `a(z)` with `a(d - z)`, i.e. describe absorption from the
back face inward. Needed when combining forward and backward stacks in
inc_tmm. Returns `f` for chaining.
"""
function flip!(f::AbsorpAnalyticFn)
    newA1 = f.A2 * exp(-f.a1 * f.d)
    newA2 = f.A1 * exp(+f.a1 * f.d)
    f.A1, f.A2 = newA1, newA2
    f.A3 = conj(f.A3 * exp(1im * f.a3 * f.d))
    return f
end

"""
    scale!(f::AbsorpAnalyticFn, factor::Real) -> f

In-place: multiply absorption everywhere by `factor`.
"""
function scale!(f::AbsorpAnalyticFn, factor::Real)
    f.A1 *= factor
    f.A2 *= factor
    f.A3 *= factor
    return f
end

"""
    add!(f::AbsorpAnalyticFn, g::AbsorpAnalyticFn) -> f

In-place: add `g`'s absorption to `f`'s. Errors if `g` has different a1, a3
(i.e. describes a different physical layer at a different λ).
"""
function add!(f::AbsorpAnalyticFn, g::AbsorpAnalyticFn)
    if f.a1 != g.a1 || f.a3 != g.a3
        error("Incompatible absorption analytical functions " *
              "(different layers or wavelengths).")
    end
    f.A1 += g.A1
    f.A2 += g.A2
    f.A3 += g.A3
    return f
end

"""
    copy(f::AbsorpAnalyticFn) -> AbsorpAnalyticFn

Return an independent copy (so `flip!`/`scale!` don't affect the original).
"""
Base.copy(f::AbsorpAnalyticFn) =
    AbsorpAnalyticFn(f.A1, f.A2, f.A3, f.a1, f.a3, f.d)

# -----------------------------------------------------------------------------
# position_resolved  -- point evaluation of Poynting, absorption, E-field
# -----------------------------------------------------------------------------
"""
    position_resolved(layer::Int, distance::Real, coh_data::Dict) -> NamedTuple

At depth `distance` (nm) inside `layer` (1-based; layer 1 is the incidence
medium, in which only z ≤ 0 is meaningful), return a NamedTuple with fields

  poyn   :: Float64    — component of Poynting vector along z (normalised)
  absor  :: Float64    — absorbed power per unit length at this depth
  Ex, Ey, Ez :: ComplexF64 — electric-field amplitudes (|E_in| = 1)

Julia syntax: `(; a=1, b=2)` is a NamedTuple literal — accessed as `nt.a`.
Lighter than a Dict, immutable, fast.
"""
function position_resolved(layer::Int, distance::Real, coh_data::Dict)
    # Pull the things we need out of the coh_tmm output.
    if layer > 1
        v = coh_data["vw_list"][layer, 1]
        w = coh_data["vw_list"][layer, 2]
    else
        # Layer 1 (incidence medium) has no "vw" assigned; by convention the
        # forward amp is 1 (the incoming wave) and backward is r (reflected).
        v = ComplexF64(1)
        w = ComplexF64(coh_data["r"])
    end
    kz   = coh_data["kz_list"][layer]
    th   = coh_data["th_list"][layer]
    n    = coh_data["n_list"][layer]
    n_0  = coh_data["n_list"][1]
    th_0 = coh_data["th_0"]
    pol  = coh_data["pol"]

    # Sanity: distance should be inside the layer (or ≤ 0 if layer 1).
    # Note: @assert's message argument must be on the same logical line as the
    # condition. We pull it into a variable to avoid the multi-line parse error.
    _layer_err = "distance=$distance is outside layer $layer " *
                 "(thickness $(coh_data["d_list"][layer]))"
    @assert (layer >= 2 && 0 <= distance <= coh_data["d_list"][layer]) ||
            (layer == 1 && distance <= 0) _layer_err

    # Forward and backward field amplitudes at this depth.
    Ef = v * exp(+1im * kz * distance)
    Eb = w * exp(-1im * kz * distance)

    # Poynting vector (z-component, normalised to incoming = 1).
    poyn = if pol == "s"
        real(n * cos(th) * conj(Ef + Eb) * (Ef - Eb)) / real(n_0 * cos(th_0))
    else  # "p"
        real(n * conj(cos(th)) * (Ef + Eb) * conj(Ef - Eb)) /
            real(n_0 * conj(cos(th_0)))
    end

    # Absorbed power per unit length at this depth.
    absor = if pol == "s"
        imag(n * cos(th) * kz * abs2(Ef + Eb)) / real(n_0 * cos(th_0))
    else  # "p"
        imag(n * conj(cos(th)) *
             (kz * abs2(Ef - Eb) - conj(kz) * abs2(Ef + Eb))) /
            real(n_0 * conj(cos(th_0)))
    end

    # E-field components. z is normal to the interfaces; light rays in x,z plane.
    Ex, Ey, Ez = if pol == "s"
        (ComplexF64(0), Ef + Eb, ComplexF64(0))
    else  # "p"
        ((Ef - Eb) * cos(th), ComplexF64(0), (-Ef - Eb) * sin(th))
    end

    return (poyn = poyn, absor = absor, Ex = Ex, Ey = Ey, Ez = Ez)
end

# -----------------------------------------------------------------------------
# absorp_in_each_layer  -- integrated absorption per layer
# -----------------------------------------------------------------------------
"""
    absorp_in_each_layer(coh_data::Dict) -> Vector{Float64}

Return a length-`num_layers` vector whose entries sum to 1, giving the
fraction of incoming intensity absorbed in each layer. Entry 1 is the power
ultimately reflected (absorbed in the incidence medium by convention); the
last entry is `T` (absorbed in the exit medium).
"""
function absorp_in_each_layer(coh_data::Dict)
    num_layers = length(coh_data["d_list"])
    pe = zeros(Float64, num_layers)            # power entering each layer

    pe[1]   = 1.0                              # all incident light "enters" the 0th
    pe[2]   = coh_data["power_entering"]       # after first interface
    pe[end] = coh_data["T"]                    # power leaving the back

    for i in 3:(num_layers - 1)
        # Poynting at the start of layer i = power that passed into layer i.
        pe[i] = position_resolved(i, 0.0, coh_data).poyn
    end

    # Absorption in layer i = pe[i] - pe[i+1].
    # `diff` returns the successive differences; we negate to get
    # (entering - exiting) per layer.
    absorp = zeros(Float64, num_layers)
    absorp[1:end-1] .= -diff(pe)
    absorp[end] = pe[end]
    return absorp
end

# -----------------------------------------------------------------------------
# Building the G(x, λ) matrix for the EQE pipeline
# -----------------------------------------------------------------------------

"""
    Stack

Lightweight container describing the optical stack for the generation
calculation. `nk` is a vector of `NKData` (one per layer including the two
semi-infinite media). `d_nm` lists thicknesses with `Inf` first and last,
matching `coh_tmm`'s convention. `active_layer` is the 1-based index of the
photoactive layer (typically the perovskite) — flagged so downstream code
knows which slice of G to feed into ChargeTransport.

Julia note: this is a *struct*, not a Dict, so VS Code's Julia language
server gives you autocomplete on field names. Worth the few lines.
"""
struct Stack
    names::Vector{String}      # ["air", "MgF2", "ITO", "SnO2", "C60", "Pvk", "Ag"]
    nk::Vector                 # Vector{NKData} — left untyped to avoid
                               # tight-coupling to materials.jl's include order
    d_nm::Vector{Float64}      # [Inf, 100, 50, 15, 20, 800, Inf]
    active_layer::Int          # index into names/d_nm of the perovskite
end

"""
    layer_nk_at(stack::Stack, λ::Real) -> Vector{ComplexF64}

Assemble the `n_list` for `coh_tmm` at a single wavelength: one complex
index per layer. Semi-infinite media (air, Ag back contact) just sample
their NKData at λ; same for the finite layers. Falls back to `1.0+0im`
for any NKData entry that is `nothing` (used when we want "vacuum/air" as
a placeholder without a CSV).
"""
function layer_nk_at(stack::Stack, λ::Real)
    out = Vector{ComplexF64}(undef, length(stack.nk))
    for (i, mat) in enumerate(stack.nk)
        if mat === nothing
            out[i] = ComplexF64(1.0)            # vacuum
        else
            out[i] = interpolate_nk(mat, [λ])[1]
        end
    end
    return out
end

"""
    layer_starts_nm(d_nm::Vector{Float64}) -> Vector{Float64}

Cumulative distance (nm) from the front of the stack to the START of each
finite layer. Entries for the two semi-infinite layers are -Inf / +Inf
respectively. Mirrors McByrne's `layer_starts`.
"""
function layer_starts_nm(d_nm::Vector{Float64})
    out = zeros(Float64, length(d_nm))
    out[1] = -Inf
    if length(d_nm) >= 2
        out[2] = 0.0
        for i in 3:length(d_nm)
            out[i] = out[i-1] + d_nm[i-1]
        end
    end
    return out
end

"""
    compute_G_matrix(stack::Stack, λ_grid::AbstractVector, x_grid_nm::AbstractVector;
                     photon_flux = ones(length(λ_grid)),
                     pol = "s", th_0 = 0.0)
        -> (G::Matrix{Float64}, A_per_layer::Matrix{Float64}, R::Vector, T::Vector)

Spatially- and spectrally-resolved generation matrix.

  G[i, j] = generation rate at depth `x_grid_nm[i]` and wavelength
            `λ_grid[j]`, in (per nm) × (photon_flux[j]).

Specifically, for each λ:
  G(x, λ) = a(x, λ) * photon_flux(λ)
where a(x, λ) is the absorbed-power fraction per nm at depth x, computed
from `coh_tmm` + `position_resolved`. Note `a` already absorbs the α·|E|²
prefactor of your written formula — see McByrne arXiv:1603.02720 §3.

Crucially, when `photon_flux[j] = 1` (the default), the units of the
absorption fraction are (1/nm), so G has units of (photons/nm/incident
photon). To get a true volumetric rate in (carriers / nm³ / s) you can
either:
  (i) pass `photon_flux = AM1.5G_photon_flux_per_nm`  (Method 2 path), or
  (ii) leave `photon_flux = 1` and post-multiply outside (Method 1 path).

`x_grid_nm` is measured from the FRONT face of the multilayer stack
(i.e. just inside layer 2). Depths outside the active region are still
evaluated — useful for plotting parasitic absorption in C60/ITO — but
G is zero in the semi-infinite media.

Additional outputs (cheap to compute, useful for diagnostics):
  A_per_layer[i, j] — fraction of light absorbed in layer i at λ_j
  R[j], T[j]        — net reflectance / transmittance at λ_j
"""
function compute_G_matrix(stack::Stack,
                          λ_grid::AbstractVector{<:Real},
                          x_grid_nm::AbstractVector{<:Real};
                          photon_flux::AbstractVector{<:Real} = ones(length(λ_grid)),
                          pol::String = "s",
                          th_0::Real = 0.0)

    @assert length(photon_flux) == length(λ_grid) (
        "photon_flux must have the same length as λ_grid")

    N_x = length(x_grid_nm)
    N_λ = length(λ_grid)
    G   = zeros(Float64, N_x, N_λ)
    A_per_layer = zeros(Float64, length(stack.d_nm), N_λ)
    R_vec = zeros(Float64, N_λ)
    T_vec = zeros(Float64, N_λ)

    starts = layer_starts_nm(stack.d_nm)

    for (j, λ) in enumerate(λ_grid)
        n_list = layer_nk_at(stack, λ)
        coh    = coh_tmm(pol, n_list, stack.d_nm, th_0, λ)

        R_vec[j] = coh["R"]
        T_vec[j] = coh["T"]
        A_per_layer[:, j] = absorp_in_each_layer(coh)

        # Walk x_grid_nm and evaluate a(x, λ) layer-by-layer.
        # x positions in the two semi-infinite layers contribute 0 to G.
        # Cache one AbsorpAnalyticFn per interior layer at THIS wavelength so
        # we don't redo the fill_in work for every x in that layer.
        analytic = Vector{Union{AbsorpAnalyticFn, Nothing}}(
                       nothing, length(stack.d_nm))
        for ℓ in 2:(length(stack.d_nm) - 1)
            analytic[ℓ] = fill_in!(AbsorpAnalyticFn(), coh, ℓ)
        end

        for (i, x) in enumerate(x_grid_nm)
            ℓ = locate_layer(starts, stack.d_nm, x)
            if ℓ === nothing
                G[i, j] = 0.0                                # outside the stack
            else
                # depth into THIS layer:
                z = x - starts[ℓ]
                # Clamp to [0, d] to guard against floating-point drift on the
                # boundary.
                z = clamp(z, 0.0, stack.d_nm[ℓ])
                a = evaluate(analytic[ℓ], z)                 # 1/nm
                G[i, j] = a * photon_flux[j]
            end
        end
    end

    return (G = G, A_per_layer = A_per_layer, R = R_vec, T = T_vec)
end

"""
    locate_layer(starts, d_nm, x) -> Union{Int, Nothing}

Which 1-based layer index contains position `x` (nm from the front face)?
Returns `nothing` if x is in a semi-infinite layer (we don't compute G there).
"""
function locate_layer(starts::Vector{Float64}, d_nm::Vector{Float64}, x::Real)
    # Interior layers are 2 .. (end-1).
    for ℓ in 2:(length(d_nm) - 1)
        if starts[ℓ] <= x <= starts[ℓ] + d_nm[ℓ]
            return ℓ
        end
    end
    return nothing
end

# =============================================================================
# Day 7 patch for src/generation.jl
# -----------------------------------------------------------------------------
# REPLACE the stub `tmm_generation_slice(...)` (the block that ends with
#   error("tmm_generation_slice is a Day 7 task — not implemented yet")
# ) with the two functions below. Everything they depend on
# (coh_tmm, layer_nk_at, fill_in!, AbsorpAnalyticFn, evaluate) already lives
# in generation.jl / tmm_core.jl, so no new `include` is needed.
# =============================================================================

"""
    tmm_generation_slice(λ, stack, grid_m;
                         h_ndoping_m,
                         photon_flux = 1.0,
                         weight1 = 1.0, weight2 = 1.0,
                         pol = "s", th_0 = 0.0)
        -> Vector{Float64}

Carrier generation rate G [carriers / m³ / s] in the photoactive layer
(`stack.active_layer`, = perovskite = layer 6 for both Stack 1 and Stack 2/3),
evaluated at the ChargeTransport node positions `grid_m`.

Arguments
  λ            wavelength, nm (scalar)
  stack        the optical `Stack`; TMM is computed INSIDE this function
  grid_m       subg2[Coordinates][:] — intrinsic-region node positions, METRES
  h_ndoping_m  position of the perovskite FRONT face in CT coordinates, METRES
               (= p.h_ndoping; perovskite spans [h_ndoping, h_ndoping+h_intrinsic])

Keyword arguments
  photon_flux  incident photon flux, photons/m²/s (scalar at this λ)
  weight1      boundary weight applied to G[1]   (Day 8 control-volume term; 1.0 = off)
  weight2      boundary weight applied to G[end] (Day 8 control-volume term; 1.0 = off)
  pol, th_0    polarisation / incidence angle (project default: s, normal)

Returns a node-valued vector, `length == length(grid_m)`, units carriers/m³/s.

Unit chain (Day 7–10 doc, corrected):
  evaluate(fn,z) is a(z) [1/nm]  -- absorbed-power fraction per nm of depth
  G = a(z)[1/nm] * 1e9[nm/m] * photon_flux[1/m²/s]  ->  carriers/m³/s
The photon energy ħω does NOT appear: Byrnes' `absor` is already a *power*
fraction, and 1 absorbed photon → 1 e–h pair, so the ħω cancels (G = a·Φ).
Factor is 1e9, not 1e27.
"""
function tmm_generation_slice(λ, stack, grid_m;
                              h_ndoping_m,
                              photon_flux = 1.0,
                              weight1 = 1.0, weight2 = 1.0,
                              pol = "s", th_0 = 0.0)

    ℓ   = stack.active_layer            # 6 for both stacks; carried on the Stack
    d_ℓ = stack.d_nm[ℓ]                 # active-layer thickness, nm (= 800.0)

    # --- 1. coherent TMM at this λ -------------------------------------------
    n_list = layer_nk_at(stack, λ)
    coh    = coh_tmm(pol, n_list, stack.d_nm, th_0, λ)

    # --- 2. analytic absorption-vs-depth for the active layer ----------------
    fn = fill_in!(AbsorpAnalyticFn(), coh, ℓ)

    # --- 3. CT grid (m) → depth into the active layer (nm, z=0 at front) ------
    grid_nm = (grid_m .- h_ndoping_m) .* 1e9

    # Sanity: if h_ndoping_m or the grid is wrong, z won't start at 0 / span d.
    # Warn (don't error) so a slightly-off grid still plots while you debug.
    if abs(grid_nm[1]) > 1.0 || abs((grid_nm[end] - grid_nm[1]) - d_ℓ) > 1.0
        @warn "Active-layer grid looks off: front z=$(round(grid_nm[1],digits=3)) nm " *
              "(expect ~0), span=$(round(grid_nm[end]-grid_nm[1],digits=3)) nm " *
              "(expect ~$(d_ℓ)). Check h_ndoping_m and that grid_m is the " *
              "INTRINSIC subgrid."
    end

    # --- 4. evaluate a(z); clamp z to [0,d] to absorb floating-point drift ----
    #        evaluate() already returns Float64 (imag noise stripped internally).
    a_vec = [evaluate(fn, clamp(z, 0.0, d_ℓ)) for z in grid_nm]   # 1/nm

    # --- 5. per-nm → per-m, then × incident photon flux ----------------------
    G = a_vec .* (1e9 * photon_flux)        # carriers/m³/s

    # --- 6. boundary weights (no-op at default 1.0; real values supplied Day 8)
    G[1]   *= weight1
    G[end] *= weight2

    return G
end

"""
    beer_lambert_generation(grid_m, h_ndoping_m; α_per_m, photon_flux) -> Vector{Float64}

Ex104-style Beer–Lambert generation on the SAME `grid_m`, for the Day 7
order-of-magnitude / shape comparison:

    G_BL(z) = photon_flux * α * exp(-α * z),   z = depth into perovskite (m)

`α_per_m` is the absorption coefficient in 1/m (Ex104 placeholder: 1.3e7),
`photon_flux` in photons/m²/s (Ex104: 9.4e20). Use the SAME `photon_flux`
you pass to `tmm_generation_slice` so the two curves are directly comparable.
"""
function beer_lambert_generation(grid_m, h_ndoping_m; α_per_m, photon_flux)
    z_m = grid_m .- h_ndoping_m                      # depth into Pvk, metres
    return photon_flux .* α_per_m .* exp.(-α_per_m .* z_m)   # carriers/m³/s
end