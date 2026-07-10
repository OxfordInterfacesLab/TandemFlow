# examples/02_layer_absorption.jl
# -----------------------------------------------------------------------------
# Where does each photon go? Per-layer absorptance A_i(λ) for Stack 1,
# separating USEFUL absorption (perovskite) from PARASITIC absorption
# (MgF2, ITO, SnO2, C60) and the light transmitted into the silicon.
#
# Demonstrates: absorp_in_each_layer via compute_G_matrix, and using the
# Stack's own layer names for labelling.
#
#     julia> include("examples/02_layer_absorption.jl")
# -----------------------------------------------------------------------------

using PerovskiteTMMlite
ENV["GKSwstype"] = "100"   # headless-safe PNG output (no display needed)
using Plots
using Printf
using DelimitedFiles

const OUTDIR = joinpath(@__DIR__, "outputs")
mkpath(OUTDIR)

λ     = collect(300.0:5.0:1000.0)
stack = build_stack1()

x_nm = collect(0.0:5.0:1000.0)
res  = compute_G_matrix(stack, λ, x_nm)

# Finite layers are indices 2 .. end-1 (index 1 = air/reflection, end = Si).
finite = 2:(length(stack.d_nm) - 1)
names  = stack.names

plt = plot(xlabel = "wavelength (nm)", ylabel = "absorptance",
           title = "Stack 1 — per-layer absorption", legend = :topright,
           ylims = (0, 1))
for ℓ in finite
    lw = ℓ == stack.active_layer ? 3 : 1.5      # emphasise the perovskite
    plot!(plt, λ, res.A_per_layer[ℓ, :]; label = names[ℓ], lw = lw)
end
# also show what leaves into Si
plot!(plt, λ, res.T; label = "→ Si (T)", lw = 1.5, ls = :dash, color = :black)
savefig(plt, joinpath(OUTDIR, "02_layer_absorption_stack1.png"))

# --- parasitic vs useful summary (integrated, flat weighting) ----------------
A_pvk   = sum(res.A_per_layer[stack.active_layer, :])
A_paras = 0.0
for ℓ in finite
    ℓ == stack.active_layer && continue
    global A_paras += sum(res.A_per_layer[ℓ, :])
end
@printf("Band-integrated (flat) — useful A_pvk : parasitic = %.1f : %.1f\n",
        A_pvk, A_paras)

# --- CSV: one column per finite layer ----------------------------------------
open(joinpath(OUTDIR, "02_layer_absorption_stack1.csv"), "w") do io
    hdr = "wavelength_nm," * join(names[finite], ",") * ",T_into_Si"
    println(io, hdr)
    M = hcat(λ, permutedims(res.A_per_layer[finite, :]), res.T)
    writedlm(io, M, ',')
end

println("wrote outputs/02_layer_absorption_stack1.png and .csv")
