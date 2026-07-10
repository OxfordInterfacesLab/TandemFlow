# examples/01_rta_spectrum.jl
# -----------------------------------------------------------------------------
# Reflectance R(λ), transmittance T(λ) and useful perovskite absorptance
# A_pvk(λ) for Stack 1 (the paper optic: air / MgF2 / ITO / SnO2 / C60 /
# CsFAPbIBr / ITO_inter / Si).
#
# Demonstrates: build a Stack from bundled n,k data, run the TMM over a
# wavelength grid, and read R / T / per-layer absorptance out of one call to
# compute_G_matrix.
#
# Run from the repo root:
#     julia --project=.            # then, at the REPL:
#     julia> include("examples/01_rta_spectrum.jl")
# -----------------------------------------------------------------------------

using PerovskiteTMMlite
ENV["GKSwstype"] = "100"   # headless-safe PNG output (no display needed)
using Plots
using Printf
using DelimitedFiles

const OUTDIR = joinpath(@__DIR__, "outputs")
mkpath(OUTDIR)

# --- wavelength grid (nm). 300–1000 nm is fully covered by all n,k files. ---
λ = collect(300.0:5.0:1000.0)

stack = build_stack1()
pvk   = stack.active_layer          # perovskite layer index (= 6)

# A single spatial node inside the perovskite is enough to get R/T/A; here we
# use a light grid across the whole device so the call is also reusable.
x_nm = collect(0.0:5.0:1000.0)

res = compute_G_matrix(stack, λ, x_nm)   # (; G, A_per_layer, R, T)

R     = res.R
T     = res.T
A_pvk = res.A_per_layer[pvk, :]          # useful absorption in the perovskite
A_tot = vec(sum(res.A_per_layer[2:end-1, :], dims = 1))   # all finite layers

# --- energy-conservation sanity check: R + Σ(A_layers) + T ≈ 1 ---------------
closure = R .+ A_tot .+ T
@printf("Energy closure  R+ΣA+T:  min = %.4f   max = %.4f  (should be ≈ 1)\n",
        minimum(closure), maximum(closure))
@printf("A_pvk peak: %.3f at λ = %d nm\n",
        maximum(A_pvk), Int(round(λ[argmax(A_pvk)])))

# --- plot --------------------------------------------------------------------
plt = plot(λ, R;     label = "R (reflectance)", lw = 2,
           xlabel = "wavelength (nm)", ylabel = "fraction",
           title  = "Stack 1 — R / T / A_pvk", legend = :right, ylims = (0, 1))
plot!(plt, λ, T;      label = "T (into Si)",      lw = 2)
plot!(plt, λ, A_pvk;  label = "A_pvk (useful)",   lw = 2)
plot!(plt, λ, A_tot;  label = "A_total (all layers)", lw = 1, ls = :dash)
savefig(plt, joinpath(OUTDIR, "01_rta_stack1.png"))

# --- CSV out -----------------------------------------------------------------
open(joinpath(OUTDIR, "01_rta_stack1.csv"), "w") do io
    println(io, "wavelength_nm,R,T,A_pvk,A_total")
    writedlm(io, hcat(λ, R, T, A_pvk, A_tot), ',')
end

println("wrote outputs/01_rta_stack1.png and .csv")
