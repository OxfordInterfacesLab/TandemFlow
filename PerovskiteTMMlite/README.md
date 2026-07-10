# PerovskiteTMMlite.jl

A small, dependency-light **transfer-matrix-method (TMM)** optics package for
multilayer thin-film solar cells (perovskite and perovskite/silicon tandems),
written in Julia.

It computes, for an arbitrary stack of coherent thin films:

- **R(λ), T(λ)** — reflectance and transmittance into the back medium,
- **per-layer absorptance A_i(λ)** — useful vs parasitic absorption,
- **G(x, λ)** — spatially- and spectrally-resolved carrier-generation rate.

The optical core is a Julia translation of the coherent-TMM routines in
Steven Byrnes' Python [`tmm`](https://github.com/sbyrnes321/tmm) library
(physics/derivations: [arXiv:1603.02720](https://arxiv.org/abs/1603.02720)).
It is the optical front-end used to feed carrier generation into a
drift-diffusion (ChargeTransport.jl) EQE model; that coupling lives in a
separate pipeline and is **not** part of this repository.

## Install

Requires Julia ≥ 1.10.

```bash
git clone <this-repo-url> PerovskiteTMMlite
cd PerovskiteTMMlite
julia --project=.
```
```julia
julia> ]                          # Pkg mode
(PerovskiteTMMlite) pkg> instantiate  # download & build deps (first time only)
```

## Quickstart

```julia
using PerovskiteTMMlite

stack = build_stack1()                 # air/MgF2/ITO/SnO2/C60/CsFAPbIBr/ITO_inter/Si
λ     = collect(300.0:5.0:1000.0)      # nm
x_nm  = collect(0.0:5.0:1000.0)        # depth grid (nm from front face)

res = compute_G_matrix(stack, λ, x_nm) # NamedTuple: (G, A_per_layer, R, T)

res.R                                  # reflectance,  length(λ)
res.T                                  # transmittance into Si, length(λ)
res.A_per_layer[stack.active_layer, :] # perovskite absorptance A_pvk(λ)
res.G                                  # generation matrix, size (length(x_nm), length(λ))
```

A single-wavelength TMM call, if you want the raw field solution:

```julia
n_list = layer_nk_at(stack, 500.0)     # Vector{ComplexF64}, one index per layer
coh    = coh_tmm("s", n_list, stack.d_nm, 0.0, 500.0)
coh["R"], coh["T"]
absorp_in_each_layer(coh)              # per-layer absorptance at 500 nm
```

## Examples

See [`examples/`](examples/) for four runnable scripts (R/T/A spectra,
per-layer absorption, G(x,λ) profiles/heatmap, and a three-stack comparison).
Each writes a PNG and a CSV to `examples/outputs/`.

```bash
julia --project=. examples/run_all.jl
```

## Bundled stacks

| Builder | Stack | Back medium |
|---------|-------|-------------|
| `build_stack1()`    | air / MgF2(100) / ITO(50) / SnO2(15) / C60(20) / CsFAPbIBr(800) / ITO_inter(15) / Si | silicon (∞) |
| `build_stack2()`    | …/ CsFAPbIBr(800) / spiro-TTB(50) / ITO_inter(15) / Si | silicon (∞) |
| `build_stack1_Ag()` | …/ CsFAPbIBr(800) / Ag | silver (∞) |

Perovskite is layer 6 in every stack. `n,k` data (250–1450 nm) for MgF2, ITO,
SnO2, C60, CsFAPbIBr, spiro-TTB, Si and Ag live in `data/nk/`; the AM1.5G
reference spectrum is in `data/AM15G.csv`.

## Public API

```
Optics       coh_tmm, absorp_in_each_layer, position_resolved,
             interface_r, interface_t
Materials    NKData, load_nk_data, interpolate_nk, load_and_interp
Stacks       Stack, build_stack1, build_stack2, build_stack1_Ag,
             layer_nk_at, layer_starts_nm
Generation   compute_G_matrix, tmm_generation_slice, beer_lambert_generation,
             AbsorpAnalyticFn, fill_in!, evaluate
```

Every exported function has a docstring: `julia> ?compute_G_matrix`.

## Data & units

- Wavelengths and layer thicknesses are in **nanometres** throughout the TMM
  (Byrnes' convention: pick one length unit and stay consistent).
- `n,k` CSVs are three columns `wl,n,k` with `wl` in nm.
- `compute_G_matrix` returns `G` in `1/nm` per unit incident photon flux when
  `photon_flux = 1`; pass a photon-flux vector (photons·m⁻²·s⁻¹) to obtain a
  rate. `tmm_generation_slice` returns SI `carriers·m⁻³·s⁻¹` for coupling to a
  device solver.

## Attribution

Optical core ported from Steven Byrnes' `tmm` (MIT). If this optical model
supports published work, please also cite
R. S. Bonilla, *The impact of transparent conducting electrodes on tandem
solar cell efficiency*, **Joule** 9, 102211 (2025).

## License

See [LICENSE](LICENSE).
