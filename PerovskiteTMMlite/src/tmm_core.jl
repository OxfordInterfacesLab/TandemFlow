# tmm_core.jl
#
# Julia translation of the core transfer-matrix-method functions from
# Steven Byrnes' Python `tmm` library (tmm_core.py).
# Physics, conventions, derivations: https://arxiv.org/abs/1603.02720
#
# Ported 1:1 (including is_forward_angle gain/forward-wave checks).
# Indexing is 1-based throughout (Julia), so every loop bound is shifted
# by +1 relative to the Python original.
#
# Convention note: the whole pipeline is kept COMPLEX from the start.
# NumPy silently promotes reals to complex when needed; Julia does not,
# and `asin` of a real >1 throws a DomainError. Keeping indices/angles
# complex avoids that entire class of bug.

const EPSILON = eps(Float64)   # machine epsilon, ~2.22e-16
                               # Python used sys.float_info.epsilon (same value)

# -----------------------------------------------------------------------------
# make_2x2_array
# -----------------------------------------------------------------------------
# Python built [[a,b],[c,d]] via a hand-filled numpy array "for speed".
# In Julia the literal `[a b; c d]` IS the fast path — spaces separate
# columns, semicolons separate rows. We force the element type to ComplexF64
# so downstream matrix algebra never hits a type surprise.
#
# Julia syntax note: a function can be a one-liner with `=`. The `T` is a
# type parameter (like a generic); `::Type{T}` lets the caller pass a type
# as an argument, mirroring Python's `dtype=` keyword.
function make_2x2_array(a, b, c, d; dtype::Type{T}=ComplexF64) where {T}
    return T[a b; c d]      # T[...] makes an array with element type T
end

# -----------------------------------------------------------------------------
# interface_r  -- Fresnel reflection amplitude
# -----------------------------------------------------------------------------
# pol is "s" or "p". n_i,n_f are (complex) indices; th_i,th_f are (complex)
# angles in radians (0 = normal incidence).
#
# Julia syntax note: string comparison uses `==`. We `error(...)` instead of
# Python's `raise ValueError(...)`.
function interface_r(pol, n_i, n_f, th_i, th_f)
    if pol == "s"
        return (n_i * cos(th_i) - n_f * cos(th_f)) /
               (n_i * cos(th_i) + n_f * cos(th_f))
    elseif pol == "p"
        return (n_f * cos(th_i) - n_i * cos(th_f)) /
               (n_f * cos(th_i) + n_i * cos(th_f))
    else
        error("Polarization must be \"s\" or \"p\"")
    end
end

# -----------------------------------------------------------------------------
# interface_t  -- Fresnel transmission amplitude
# -----------------------------------------------------------------------------
function interface_t(pol, n_i, n_f, th_i, th_f)
    if pol == "s"
        return 2 * n_i * cos(th_i) / (n_i * cos(th_i) + n_f * cos(th_f))
    elseif pol == "p"
        return 2 * n_i * cos(th_i) / (n_f * cos(th_i) + n_i * cos(th_f))
    else
        error("Polarization must be \"s\" or \"p\"")
    end
end

# -----------------------------------------------------------------------------
# R_from_r and T_from_t  -- power from amplitudes
# -----------------------------------------------------------------------------
# Needed by coh_tmm at the end. abs2(z) == |z|^2 (faster than abs(z)^2).
R_from_r(r) = abs2(r)

function T_from_t(pol, t, n_i, n_f, th_i, th_f)
    if pol == "s"
        return abs2(t) * real(n_f * cos(th_f)) / real(n_i * cos(th_i))
    elseif pol == "p"
        return abs2(t) * real(n_f * conj(cos(th_f))) /
                         real(n_i * conj(cos(th_i)))
    else
        error("Polarization must be \"s\" or \"p\"")
    end
end

# power entering the first interface (≈ 1-R, but exact when n_i is complex)
function power_entering_from_r(pol, r, n_i, th_i)
    if pol == "s"
        return real(n_i * cos(th_i) * (1 + conj(r)) * (1 - r)) /
               real(n_i * cos(th_i))
    elseif pol == "p"
        return real(n_i * conj(cos(th_i)) * (1 + r) * (1 - conj(r))) /
               real(n_i * conj(cos(th_i)))
    else
        error("Polarization must be \"s\" or \"p\"")
    end
end

# -----------------------------------------------------------------------------
# is_forward_angle
# -----------------------------------------------------------------------------
# Given index n and angle theta, is this the forward-traveling wave?
# For real n & theta the criterion is -pi/2 < theta < pi/2, but for complex
# values it's subtler. See arXiv:1603.02720 Appendix D.
#
# Julia syntax notes:
#  * `n.real * n.imag` (Python) -> `real(n) * imag(n)`
#  * `assert cond, msg` (Python) -> `@assert cond msg`
#  * we `convert(Bool, ...)` to mirror the explicit numpy->python bool cast.
function is_forward_angle(n, theta)
    @assert real(n) * imag(n) >= 0 (
        "For materials with gain, it's ambiguous which beam is incoming vs " *
        "outgoing. See arXiv:1603.02720 Appendix C.\nn: $n   angle: $theta")
    ncostheta = n * cos(theta)
    answer = false
    if abs(imag(ncostheta)) > 100 * EPSILON
        # evanescent decay or lossy medium: the decaying one is forward
        answer = imag(ncostheta) > 0
    else
        # forward is the one with positive Poynting vector
        answer = real(ncostheta) > 0
    end
    answer = Bool(answer)

    error_string = ("It's not clear which beam is incoming vs outgoing. " *
                    "Weird index maybe?\nn: $n   angle: $theta")
    if answer
        @assert imag(ncostheta) > -100 * EPSILON error_string
        @assert real(ncostheta) > -100 * EPSILON error_string
        @assert real(n * cos(conj(theta))) > -100 * EPSILON error_string
    else
        @assert imag(ncostheta) < 100 * EPSILON error_string
        @assert real(ncostheta) < 100 * EPSILON error_string
        @assert real(n * cos(conj(theta))) < 100 * EPSILON error_string
    end
    return answer
end

# -----------------------------------------------------------------------------
# list_snell
# -----------------------------------------------------------------------------
# Return the angle in each layer given the incidence angle th_0 in layer 1.
# Angles may be complex.
#
# Julia syntax notes:
#  * `.` broadcasting: `asin.(x ./ n_list)` applies asin elementwise. This is
#    Julia's vectorization — the dot fuses the whole expression into one loop.
#  * `n_list[1]` is the FIRST layer (Python's n_list[0]); `n_list[end]` is last.
#  * We make the argument to asin explicitly Complex so asin never DomainErrors.
function list_snell(n_list, th_0)
    # n_list[1]*sin(th_0)/n_list, made complex, then elementwise asin
    arg = ComplexF64.(n_list[1] .* sin(th_0) ./ n_list)
    angles = asin.(arg)
    # Only the first and last entries must be the forward angle; the
    # intermediate layers don't matter (arXiv:1603.02720 Section 5).
    if !is_forward_angle(n_list[1], angles[1])
        angles[1] = pi - angles[1]
    end
    if !is_forward_angle(n_list[end], angles[end])
        angles[end] = pi - angles[end]
    end
    return angles
end

# -----------------------------------------------------------------------------
# coh_tmm  -- the main coherent transfer-matrix-method calculation
# -----------------------------------------------------------------------------
# pol      : "s" or "p"
# n_list   : refractive indices, light passes through in order. n_list[1] is
#            the semi-infinite incidence medium, n_list[end] the exit medium.
# d_list   : layer thicknesses; first and last MUST be Inf.
# th_0     : incidence angle (0 = normal). Complex if n_list[1] is complex.
# lam_vac  : vacuum wavelength (any length unit; be consistent with d_list).
#
# Returns a Dict (Julia's dictionary; Python returned a dict too) with keys
# "r","t","R","T","power_entering","vw_list","kz_list","th_list", plus echoes
# of the inputs.
function coh_tmm(pol, n_list, d_list, th_0, lam_vac)
    # Make sure we're working with complex indices and float thicknesses.
    n_list = ComplexF64.(n_list)
    d_list = Float64.(d_list)

    num_layers = length(n_list)

    # ----- input checks (mirror the Python asserts) -----
    if length(n_list) != length(d_list)
        error("Problem with n_list or d_list! Sizes differ.")
    end
    @assert d_list[1] == Inf && d_list[end] == Inf "d_list must start and end with Inf!"
    @assert abs(imag(n_list[1] * sin(th_0))) < 100 * EPSILON "Error in n0 or th0!"
    @assert is_forward_angle(n_list[1], th_0) "Error in n0 or th0!"

    # th_list: propagation angle in each layer (Snell). May be complex.
    th_list = list_snell(n_list, th_0)

    # kz: z-component of the (complex) wavevector for the forward wave.
    # Broadcasting again: elementwise over all layers.
    kz_list = 2 .* pi .* n_list .* cos.(th_list) ./ lam_vac

    # delta: total phase accrued traversing each layer = kz * thickness.
    # d_list has Inf at the ends, so this produces Inf/NaN there; that's fine,
    # those entries are never used. (NumPy warned about inf*; Julia is silent.)
    delta = kz_list .* d_list

    # For a very opaque layer, clamp imag(delta) to avoid overflow / divide-by-0.
    # imag(delta) > 35 corresponds to single-pass transmission < 1e-30.
    # Note 2:num_layers-1 in Julia == Python's range(1, num_layers-1).
    for i in 2:(num_layers - 1)
        if imag(delta[i]) > 35
            delta[i] = real(delta[i]) + 35im
        end
    end

    # t_list[i,j], r_list[i,j]: amplitudes going from layer i to layer j.
    # Only j = i+1 is ever needed. 2D arrays are overkill but match Python and
    # keep the indices unambiguous. zeros(ComplexF64, m, n) is an m×n complex 0 matrix.
    t_list = zeros(ComplexF64, num_layers, num_layers)
    r_list = zeros(ComplexF64, num_layers, num_layers)
    for i in 1:(num_layers - 1)
        t_list[i, i+1] = interface_t(pol, n_list[i], n_list[i+1],
                                     th_list[i], th_list[i+1])
        r_list[i, i+1] = interface_r(pol, n_list[i], n_list[i+1],
                                     th_list[i], th_list[i+1])
    end

    # M_list[i] is the 2x2 transfer matrix for layer i. Defined for the interior
    # layers only (i = 2 .. num_layers-1). We store them in a Vector of 2x2
    # matrices. Python used a (num_layers,2,2) array; a vector-of-matrices is
    # the idiomatic Julia equivalent and reads more naturally.
    #
    # Julia syntax note: `Matrix{ComplexF64}(undef, 2, 2)` is an uninitialized
    # 2x2; `[fill...]` comprehension below pre-allocates the vector.
    M_list = Vector{Matrix{ComplexF64}}(undef, num_layers)
    for i in 2:(num_layers - 1)
        phase = make_2x2_array(exp(-1im * delta[i]), 0, 0, exp(1im * delta[i]))
        fres  = make_2x2_array(1, r_list[i, i+1], r_list[i, i+1], 1)
        M_list[i] = (1 / t_list[i, i+1]) * (phase * fres)   # `*` = matrix product
    end

    # Mtilde = product of interior M's, then left-multiplied by the 0->1 interface.
    Mtilde = make_2x2_array(1, 0, 0, 1)   # 2x2 identity (complex)
    for i in 2:(num_layers - 1)
        Mtilde = Mtilde * M_list[i]
    end
    front = make_2x2_array(1, r_list[1, 2], r_list[1, 2], 1) / t_list[1, 2]
    Mtilde = front * Mtilde

    # Net complex reflection / transmission amplitudes.
    r = Mtilde[2, 1] / Mtilde[1, 1]    # Python Mtilde[1,0]/Mtilde[0,0]
    t = 1 / Mtilde[1, 1]               # Python 1/Mtilde[0,0]

    # vw_list[n,:] = [v_n, w_n], forward/backward amplitudes just past the
    # interface entering layer n. Layer 1 has no left interface (undefined).
    vw_list = zeros(ComplexF64, num_layers, 2)
    vw = ComplexF64[t; 0]              # a length-2 complex column vector
    vw_list[end, :] = vw
    # Python: range(num_layers-2, 0, -1)  ->  Julia: (num_layers-1):-1:2
    for i in (num_layers - 1):-1:2
        vw = M_list[i] * vw
        vw_list[i, :] = vw
    end

    # Net powers.
    R = R_from_r(r)
    T = T_from_t(pol, t, n_list[1], n_list[end], th_0, th_list[end])
    power_entering = power_entering_from_r(pol, r, n_list[1], th_0)

    return Dict(
        "r" => r, "t" => t, "R" => R, "T" => T,
        "power_entering" => power_entering,
        "vw_list" => vw_list, "kz_list" => kz_list, "th_list" => th_list,
        "pol" => pol, "n_list" => n_list, "d_list" => d_list,
        "th_0" => th_0, "lam_vac" => lam_vac)
end
