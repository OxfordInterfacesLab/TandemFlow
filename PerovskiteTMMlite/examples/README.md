# Examples

Four self-contained TMM examples. They use only the bundled `data/nk/*.csv`
files, so no external data is required. Each writes a PNG **and** a CSV into
`examples/outputs/` so results can be inspected with or without a plotting
window.

Run from the repository root:

```julia
julia --project=.
julia> ]                       # enter Pkg mode
(PerovskiteTMMlite) pkg> instantiate   # first time only
(PerovskiteTMMlite) pkg> <backspace>   # back to the REPL
julia> include("examples/01_rta_spectrum.jl")
```

or run everything at once:

```bash
julia --project=. examples/run_all.jl
```

The first plot takes ~30 s while `Plots` precompiles; later runs are fast.

| File | What it shows | Key functions |
|------|---------------|---------------|
| `01_rta_spectrum.jl` | R(λ), T(λ), useful A_pvk(λ) for Stack 1; prints an energy-closure check (R+ΣA+T ≈ 1). | `build_stack1`, `compute_G_matrix` |
| `02_layer_absorption.jl` | Per-layer absorptance — useful (perovskite) vs parasitic (MgF2/ITO/SnO2/C60) vs transmitted into Si. | `compute_G_matrix`, `Stack` fields |
| `03_generation_profile.jl` | G(z) at 400/500/600/700 nm and a log-scale G(z,λ) heatmap inside the perovskite. | `compute_G_matrix`, `layer_starts_nm` |
| `04_compare_stacks.jl` | A_pvk(λ) and back-side T(λ) for all three stacks (Si tandem, +HTL, single-junction Ag). | `build_stack1/2/1_Ag` |

## Expected quick-look numbers (Stack 1)

- Energy closure R + ΣA + T stays within a few ×10⁻³ of 1 across 300–1000 nm.
- A_pvk peaks around 0.8–0.9 in the visible and falls off past ~720 nm
  (CsFAPbIBr, Eg ≈ 1.72 eV).
- In example 3, G(z=front) > G(z=back) at 500 nm (blue light absorbed near the
  C60 side).
- In example 4, the Ag stack has T ≈ 0 (opaque back reflector); the two
  Si-terminated stacks pass a filtered long-λ spectrum on to the silicon.
