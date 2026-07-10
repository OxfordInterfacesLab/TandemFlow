# examples/04_compare_stacks.jl
# -----------------------------------------------------------------------------
# Compare all three bundled stacks:
#
#   Stack 1     : ... / CsFAPbIBr / ITO_inter / Si        (paper optic, no HTL)
#   Stack 2     : ... / CsFAPbIBr / spiro-TTB / ITO_inter / Si   (with HTL)
#   Stack 1_Ag  : ... / CsFAPbIBr / Ag                    (single junction)
#
# Plots useful perovskite absorptance A_pvk(λ) and back-side transmittance T(λ).
# For the Ag stack T ≈ 0 (opaque back reflector); for the Si-terminated tandem
# stacks T(λ) is the light passed on to the silicon bottom cell.
#
#     julia> include("examples/04_compare_stacks.jl")
# -----------------------------------------------------------------------------

using PerovskiteTMMlite
ENV["GKSwstype"] = "100"   # headless-safe PNG output (no display needed)
using Plots
using Printf
using DelimitedFiles

const OUTDIR = joinpath(@__DIR__, "outputs")
mkpath(OUTDIR)

λ = collect(300.0:5.0:1000.0)

# Perovskite front face is identical (185 nm) for all three stacks, since the
# layers ahead of it (MgF2/ITO/SnO2/C60) are the same. A coarse in-perovskite
# grid is fine here: A_pvk / T do not depend on the spatial grid.
x_nm = collect(185.0:20.0:985.0)

builders = [("Stack 1 (Si)",    build_stack1),
            ("Stack 2 (HTL,Si)", build_stack2),
            ("Stack 1_Ag",       build_stack1_Ag)]

pA = plot(xlabel = "wavelength (nm)", ylabel = "A_pvk",
          title = "Useful perovskite absorptance", legend = :bottomleft,
          ylims = (0, 1))
pT = plot(xlabel = "wavelength (nm)", ylabel = "T (back side)",
          title = "Transmittance into back medium", legend = :topleft,
          ylims = (0, 1))

results = Dict{String,Any}()
for (name, build) in builders
    stack = build()
    res   = compute_G_matrix(stack, λ, x_nm)
    A_pvk = res.A_per_layer[stack.active_layer, :]
    plot!(pA, λ, A_pvk; lw = 2, label = name)
    plot!(pT, λ, res.T; lw = 2, label = name)
    results[name] = (A_pvk = A_pvk, T = res.T)
    @printf("%-16s  A_pvk peak = %.3f   mean T(600-900nm) = %.3f\n",
            name, maximum(A_pvk),
            sum(res.T[(λ .>= 600) .& (λ .<= 900)]) / count((λ .>= 600) .& (λ .<= 900)))
end

plt = plot(pA, pT; layout = (1, 2), size = (1000, 420))
savefig(plt, joinpath(OUTDIR, "04_compare_stacks.png"))

# --- CSV ---------------------------------------------------------------------
open(joinpath(OUTDIR, "04_compare_stacks.csv"), "w") do io
    println(io, "wavelength_nm," *
        "Apvk_stack1,T_stack1,Apvk_stack2,T_stack2,Apvk_stack1Ag,T_stack1Ag")
    M = hcat(λ,
             results["Stack 1 (Si)"].A_pvk,    results["Stack 1 (Si)"].T,
             results["Stack 2 (HTL,Si)"].A_pvk, results["Stack 2 (HTL,Si)"].T,
             results["Stack 1_Ag"].A_pvk,       results["Stack 1_Ag"].T)
    writedlm(io, M, ',')
end

println("wrote outputs/04_compare_stacks.png and .csv")
