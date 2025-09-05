using ChargeTransport

include("../cells/PSC_C60_PVK_NiO.jl")
include("../cells/PSC_C60_PVK_NiO_graded.jl")

using .PSC_3
using .PSC_C60_PVK_NiO_graded

include("../utils/ct_utils.jl")
using .CTUtils

# eps_rs = [23.13, 34.71, 49.26, 59.55]

# for eps_r in eps_rs
#     println("Running for ε_r = $eps_r")
#     PSC_3.p.ε[2] = eps_r
#     IV_pvk = PSC_3.main(n=40, plotting=false, verbose=false)
# end


## graded dielectric
PSC_IV = PSC_C60_PVK_NiO_graded.main(n=40, plotting=false, verbose=false)