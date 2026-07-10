# materials.jl
#
# Refractive-index data handling for PerovskiteTMMlite.
#
# Pipeline:
#   1. load_nk_data(path)  -> NKData      (raw, from CSV)
#   2. interpolate_nk(nk, lam_grid)  -> Vector{ComplexF64}  (ready for coh_tmm)
#
# Format expected: 3-column CSV `wl,n,k` with wavelengths in nm.
# This matches what refractiveindex.info exports for the Tejada perovskite
# dataset and the other layer materials (C60, ITO, MgF2, SnO2).

using CSV
using DataFrames

# -----------------------------------------------------------------------------
# NKData  -- compact, type-stable container for one material's (wl, n, k)
# -----------------------------------------------------------------------------
# Julia syntax notes:
#  * `struct` defines an immutable composite type. Immutable means you can't
#    reassign the fields after construction; you can still mutate the *contents*
#    of the vectors inside (`nk.wl[1] = ...` works, `nk.wl = ...` does not).
#  * Field types (`::Vector{Float64}`) make the struct CONCRETE. The compiler
#    knows the exact memory layout, so accesses are as fast as Python's numpy
#    arrays — no dictionary lookup per field, no boxing.
#  * `name::String` is a label we use in warnings and plots; carrying it on the
#    struct itself is cleaner than threading a separate string everywhere.
struct NKData
    name::String
    wl::Vector{Float64}      # wavelength grid, nm
    n::Vector{Float64}       # real part of refractive index
    k::Vector{Float64}       # imaginary part (absorption); k >= 0 for lossy
end

# Pretty printing for the REPL.  Julia's `show` is like Python's `__repr__`.
# Without this you'd see `NKData("perovskite", [250.0, 260.0, ...], ...)` dumped
# in full.
function Base.show(io::IO, nk::NKData)
    print(io, "NKData(\"$(nk.name)\", $(length(nk.wl)) points, ",
              "$(nk.wl[1])-$(nk.wl[end]) nm)")
end

# -----------------------------------------------------------------------------
# load_nk_data  -- parse a refractiveindex.info-style CSV
# -----------------------------------------------------------------------------
# `path`  : path to a `wl,n,k` CSV (header row required, nm)
# `name`  : human label; defaults to the filename stem
#
# Validation policy (per Day 3 agreed strategy): WARN on suspicious data
# but never throw. The TMM will still run; you'll just see the warning in
# the REPL and can decide whether to clean the file.
#
# Things we check:
#   * wl strictly increasing  (interpolation assumes this)
#   * no NaNs in any column
#   * k >= 0  (physics: passive media don't have gain)
#   * wl looks like nm, not um (sanity: any wl < 10 is almost certainly um)
function load_nk_data(path::AbstractString;
                     name::AbstractString = splitext(basename(path))[1])
    # CSV.read returns a DataFrame; we ask for Float64 columns explicitly so
    # any malformed entry trips an informative error early.
    df = CSV.read(path, DataFrame; types=Dict(:wl => Float64,
                                              :n  => Float64,
                                              :k  => Float64))

    # `df.wl` accesses a column by name. `copy()` because the DataFrame may
    # share memory with the file's parsed buffer; we want our own arrays.
    wl = copy(df.wl)
    n  = copy(df.n)
    k  = copy(df.k)

    # --- validations (warn-but-continue) ---
    # `@warn` is Julia's built-in logging macro. It prints to stderr with a
    # yellow box; doesn't halt. The `_id` keyword keeps each warning de-duped
    # if you call load_nk_data in a loop on bad files.
    if any(isnan, wl) || any(isnan, n) || any(isnan, k)
        @warn "NaN values found in $name; downstream results will be NaN"
    end

    # `diff` returns successive differences; `all(>(0), ...)` checks strictly
    # increasing. `>(0)` is shorthand for `x -> x > 0` (a partial-application
    # trick — useful Julia idiom).
    if !all(>(0), diff(wl))
        @warn "$name: wavelengths are not strictly increasing; " *
              "interpolation results will be unreliable"
    end

    if any(<(0), k)
        @warn "$name: k has negative values (gain medium? data error?); " *
              "the TMM will still run but is_forward_angle may complain"
    end

    if minimum(wl) < 10
        @warn "$name: smallest wavelength is $(minimum(wl)) — looks like " *
              "micrometres, not nm. Multiply by 1000 if so."
    end

    return NKData(name, wl, n, k)
end

# -----------------------------------------------------------------------------
# linear_interp  -- minimal 1D linear interpolation with clamping
# -----------------------------------------------------------------------------
# Why hand-roll instead of using Interpolations.jl?
#   * One fewer dependency for a 10-line function.
#   * Clamping behavior at the edges is exactly what we agreed (silent,
#     no extrapolation). Interpolations.jl can do this too but the API is
#     heavier than we need on Day 3.
#
# `xs` must be sorted ascending. `ys` is the y-values. `xq` is the query point.
# Returns ys[1] if xq <= xs[1], ys[end] if xq >= xs[end], else linear blend.
#
# Julia syntax notes:
#  * `searchsortedfirst(xs, xq)` returns the index of the first element >= xq.
#    Same as numpy's `np.searchsorted` with side="left". O(log N).
#  * `@inbounds` tells the compiler "I promise xs[i-1] and xs[i] are valid";
#    skips bounds checks for speed. Safe here because we just branched on i.
function linear_interp(xs::Vector{Float64}, ys::Vector{Float64}, xq::Real)
    if xq <= xs[1]
        return ys[1]
    elseif xq >= xs[end]
        return ys[end]
    end
    i = searchsortedfirst(xs, xq)   # xs[i-1] < xq <= xs[i]
    @inbounds begin
        x0, x1 = xs[i-1], xs[i]
        y0, y1 = ys[i-1], ys[i]
        t = (xq - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)
    end
end

# -----------------------------------------------------------------------------
# interpolate_nk  -- resample an NKData onto a target wavelength grid
# -----------------------------------------------------------------------------
# Returns Vector{ComplexF64} with element i = n(λᵢ) + i·k(λᵢ).
# This is exactly the format `coh_tmm` consumes for `n_list` (well, one slice
# of it — coh_tmm wants one ñ per layer at one wavelength; you'll splice these
# vectors together across layers at Day 4/5).
#
# `lam_grid` can be any AbstractVector (Vector, range, etc.). Common case:
#   lam_grid = 300:5:900   # range, 300 to 900 nm in 5 nm steps
#
# Julia syntax notes:
#  * `AbstractVector{<:Real}` accepts any vector of any real subtype
#    (Float64, Int, Float32, range...) without copying.
#  * Pre-allocate the output then fill it; faster than `[expr for x in xs]`
#    when the function call has any setup cost.
#  * `im` is Julia's imaginary unit (Python's `1j`).
function interpolate_nk(nk::NKData, lam_grid::AbstractVector{<:Real})
    out = Vector{ComplexF64}(undef, length(lam_grid))
    for (i, λ) in enumerate(lam_grid)
        n_i = linear_interp(nk.wl, nk.n, λ)
        k_i = linear_interp(nk.wl, nk.k, λ)
        out[i] = n_i + im * k_i
    end
    return out
end

# Convenience: load + interpolate in one call, when you don't need the raw data
function load_and_interp(path::AbstractString,
                         lam_grid::AbstractVector{<:Real};
                         name::AbstractString = splitext(basename(path))[1])
    return interpolate_nk(load_nk_data(path; name=name), lam_grid)
end
