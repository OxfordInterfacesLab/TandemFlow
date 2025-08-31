include("../cells/Si_TOPCON.jl")
include("../cells/PSC_C60_PVK_PTAA.jl")
include("../cells/PSC_C60_PVK_NiO.jl")

using .Si_TOPCON
using .PSC_2
using .PSC_3

include("../utils/ct_utils.jl")
using .CTUtils

using ChargeTransport

function transmitted_photon_flux(phi0, alpha, d)
    return phi0 .* exp.(-alpha .* d)
end

incidentFlux = 4.3e21 / (m^2 * s) # flux onto perovskite

transmittedFlux = incidentFlux

for ireg in 1:PSC_2.numberOfRegions
    if ireg == 1
        layerHeight = PSC_2.heightLayers[ireg]
    else
        layerHeight = PSC_2.heightLayers[ireg] - PSC_2.heightLayers[ireg - 1]
    end

    global transmittedFlux = transmitted_photon_flux(transmittedFlux, PSC_2.absorption[ireg], layerHeight)
end

# println("Started simulations: PSC")
# IV_pvk = PSC_2.main(n=6, plotting=true, incidentPhotonFlux = incidentFlux)
# println("Finished simulations: PSC\nStarted simulations: TOPCON")
# IV_si = Si_TOPCON.main(n=6, plotting=true, incidentPhotonFlux = transmittedFlux)
# println("Finished simulations: TOPCON")

# char_pvk = get_cell_characteristics(IV_pvk)
# char_si = get_cell_characteristics(IV_si)

# totalPower = char_pvk["Pmax"] + char_si["Pmax"]

# println("----- PVK CELL ------")
# println("Voc: $(char_pvk["Voc"])V")
# println("Jsc: $(char_pvk["Jsc"])A/m^2")
# println("Pmax: $(char_pvk["Pmax"])W/m^2")
# println("---------------------\n\n----- SI CELL ------")
# println("Voc: $(char_si["Voc"])V")
# println("Jsc: $(char_si["Jsc"])A/m^2")
# println("Pmax: $(char_si["Pmax"])W/m^2")
# println("---------------------")

# efficiency = totalPower / ((1000) * (W / (m^2 * s))) * 100

# println("Efficiency: $(efficiency)")

IV_pvk = PSC_3.main(n=6, plotting=true, verbose=true, incidentPhotonFlux=incidentFlux)