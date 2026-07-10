# examples/run_all.jl
# Runs every example in order. From the repo root:
#     julia --project=. examples/run_all.jl
# (first run compiles Plots — expect ~30 s before the first figure appears.)

for f in ("01_rta_spectrum.jl",
          "02_layer_absorption.jl",
          "03_generation_profile.jl",
          "04_compare_stacks.jl")
    println("\n=== running $f ===")
    include(joinpath(@__DIR__, f))
end
println("\nAll examples done. See examples/outputs/ for PNGs and CSVs.")
