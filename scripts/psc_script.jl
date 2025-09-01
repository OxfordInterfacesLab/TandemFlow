using ChargeTransport

include("../cells/PSC_C60_PVK_NiO.jl")

using .PSC_3

include("../utils/ct_utils.jl")
using .CTUtils

IV_pvk = PSC_3.main(n=6, plotting=true, verbose=false)