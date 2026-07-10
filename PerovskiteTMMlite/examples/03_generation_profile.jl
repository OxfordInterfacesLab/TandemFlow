# examples/03_generation_profile.jl
# -----------------------------------------------------------------------------
# Spatially- and spectrally-resolved carrier generation inside the perovskite.
#
#   - G(x) at four wavelengths: short λ is absorbed near the front (C60 side),
#     long λ near the bandgap penetrates deeper toward the Si side.
#   - G(x, λ) as a log-scale heatmap over the perovskite depth.
#
# photon_flux defaults to 1 per wavelength, so G here is the optical absorption
# rate a(x,λ) in 1/nm (i.e. generation per unit incident photon flux). Pass a
# real AM1.5G photon flux vector to compute_G_matrix to get carriers/s.
#
#     julia> include("examples/03_generation_profile.jl")
# -----------------------------------------------------------------------------

using PerovskiteTMMlite
ENV["GKSwstype"] = "100"   # headless-safe PNG output (no display needed)
using Plots
using Printf
using DelimitedFiles

const OUTDIR = joinpath(@__DIR__, "outputs")
mkpath(OUTDIR)

stack = build_stack1()
ℓ     = stack.active_layer

# Perovskite front/back positions (nm from the front of the whole stack).
starts    = layer_starts_nm(stack.d_nm)
pvk_front = starts[ℓ]
pvk_back  = pvk_front + stack.d_nm[ℓ]

λ    = collect(300.0:5.0:780.0)                 # stop near the 1.72 eV bandgap
x_nm = collect(pvk_front:2.0:pvk_back)          # depth grid inside perovskite
z_nm = x_nm .- pvk_front                        # 0 at front face

res = compute_G_matrix(stack, λ, x_nm)          # photon_flux = 1 (default)
G   = res.G                                     # [N_x × N_λ], units 1/nm

# --- G(x) at four representative wavelengths ---------------------------------
targets = (400.0, 500.0, 600.0, 700.0)
plt1 = plot(xlabel = "depth into perovskite z (nm)",
            ylabel = "generation a(z) (1/nm)",
            title = "Stack 1 — G(z) vs wavelength", legend = :topright)
for λt in targets
    j = argmin(abs.(λ .- λt))
    plot!(plt1, z_nm, G[:, j]; lw = 2, label = @sprintf("λ = %d nm", Int(λ[j])))
end
savefig(plt1, joinpath(OUTDIR, "03_Gz_stack1.png"))

# --- G(x, λ) heatmap (log10) -------------------------------------------------
Gsafe = clamp.(G, 1e-8, Inf)                    # avoid log10(0)
Z     = permutedims(log10.(Gsafe))              # size (N_λ, N_x): λ on Y, z on X
plt2  = heatmap(z_nm, λ, Z;
                xlabel = "depth into perovskite z (nm)",
                ylabel = "wavelength (nm)",
                title = "Stack 1 — log10 G(z, λ)",
                colorbar_title = "log10 a (1/nm)")
savefig(plt2, joinpath(OUTDIR, "03_Gzlambda_heatmap_stack1.png"))

# --- checkpoint-style sanity: decay from front at 500 nm ---------------------
j500 = argmin(abs.(λ .- 500.0))
@printf("G at 500 nm: front = %.3e,  back = %.3e  (front should exceed back)\n",
        G[1, j500], G[end, j500])

# --- CSV: G(z) at the four wavelengths ---------------------------------------
open(joinpath(OUTDIR, "03_Gz_stack1.csv"), "w") do io
    cols = ["z_nm"]
    data = z_nm
    for λt in targets
        j = argmin(abs.(λ .- λt))
        push!(cols, @sprintf("G_%dnm", Int(λ[j])))
        data = hcat(data, G[:, j])
    end
    println(io, join(cols, ","))
    writedlm(io, data, ',')
end

println("wrote outputs/03_Gz_stack1.png, 03_Gzlambda_heatmap_stack1.png and .csv")
