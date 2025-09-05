using LessUnitful: @local_unitfactors, @ufac_str, @ph_str
using ChargeTransport

@kwdef struct Params_PSC_C60_PVK_NiO

    #####################################################################
    ############################ parameters ############################

    # physical constants
    constants = ChargeTransport.constants

    # used unit factors
    nm = ufac"nm"
    cm = ufac"cm"
    K = ufac"K"
    m = ufac"m"
    V = ufac"V"
    s = ufac"s"

    eV = constants.q * V

    ########## charge carriers ##########

    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    iphia = 3
    numberOfCarriers = 3 # electrons, holes and anion vacancies

    ########## device geometry ##########

    # region numbers
    regionDonor = 1
    regionIntrinsic = 2
    regionAcceptor = 3
    regions = [regionDonor, regionIntrinsic, regionAcceptor]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregionDonor = 1
    bregionAcceptor = 2
    bregionJ1 = 3
    bregionJ2 = 4
    numberOfBoundaryRegions = 4

    ## length domains
    h_ndoping = 20.0 * nm
    h_intrinsic = 100.0 * nm
    h_pdoping = 13.0 * nm
    h_total = h_ndoping + h_intrinsic + h_pdoping
    heightLayers = [
        h_ndoping,
        h_ndoping + h_intrinsic,
        h_ndoping + h_intrinsic + h_pdoping,
    ]
    ########## physical values ##########

    ## charge numbers
    zn = -1
    zp = 1
    za = 1

    ## temperature
    T = 300.0 * K

    ## band edge energies
    En = [-4.00, -4.0, -2.05] .* eV
    Ep = [-6.00, -5.68, -5.65] .* eV
    Ea = [0.0, 0.0, 0.0] .* eV

    ## effective densities of density of states
    Nn = [1.0e26, 2.75e24, 2.5e26] ./ (m^3)
    Np = [1.0e26, 3.80e24, 2.5e26] ./ (m^3)
    Na = [0.0, 1.21e28, 0.0] ./ (m^3)

    ## mobilities
    μn = [8.0e-6, 1.0e-3, 2.8e-4] .* (m^2) / (V * s)
    μp = [3.5e-7, 1.0e-3, 2.8e-4] .* (m^2) / (V * s)
    μa = [0.0, 1.0e-12, 0.0] .* (m^2) / (V * s)

    ## relative dielectric permittivity
    # 23.13, 34.71, 49.26, 59.55
    ε = [2.5, 23.13, 11.7] .* 1.0

    ## radiative recombination
    r0 = [1.0e-20, 1.0e-20, 1.0e-20] .* cm^3 / s

    ## life times and trap densities
    τn = [1.0e-7, 2.5e-9, 1.0e-6] .* s
    τp = [1.0e-7, 2.5e-9, 1.0e-6] .* s

    Augn = [1.0e-31, 1.0e-31, 1.0e-31] .* (cm^6) / (s)
    Augp = [1.0e-31, 1.0e-31, 1.0e-31] .* (cm^6) / (s)

    ## SRH trap densities
    nTrapDensity = [1.0e20, 4.0e21, 1.0e20] ./ (m^3)
    pTrapDensity = [1.0e20, 4.0e21, 1.0e20] ./ (m^3)

    ## generation
    incidentPhotonFlux = [0.0, 4.3e21, 0.0] ./ (m^2 * s)
    absorption = [0.0, 1.0e7, 0.0] ./ m
    generationPeak = h_ndoping

    generation_uniform = [0.0, 2.64e27, 0.0] ./ (m^3 * s)

    ## doping
    Cn = 1.0e24 / (m^3)
    Cp = 3.0e24 / (m^3)
    Ca = 1.0e24 / (m^3)
end

function Params(p)

    @local_unitfactors m

    params = ChargeTransport.Params(
        p.numberOfRegions,
        p.numberOfBoundaryRegions,
        p.numberOfCarriers
    )

    params.temperature = p.T
    params.chargeNumbers[p.iphin] = p.zn
    params.chargeNumbers[p.iphip] = p.zp
    params.chargeNumbers[p.iphia] = p.za

    params.dielectricConstant = p.ε * p.constants.ε_0

    ## effective DOS, band edge energy and mobilities
    params.densityOfStates[p.iphin, :] = p.Nn
    params.densityOfStates[p.iphip, :] = p.Np
    params.densityOfStates[p.iphia, :] = p.Na

    params.bandEdgeEnergy[p.iphin, :] = p.En
    params.bandEdgeEnergy[p.iphip, :] = p.Ep
    params.bandEdgeEnergy[p.iphia, :] = p.Ea

    params.mobility[p.iphin, :] = p.μn
    params.mobility[p.iphip, :] = p.μp
    params.mobility[p.iphia, :] = p.μa

    for ireg in 1:p.numberOfRegions ## interior region data
        ## recombination parameters
        params.recombinationRadiative[ireg] = p.r0[ireg]
        params.recombinationSRHLifetime[p.iphin, ireg] = p.τn[ireg]
        params.recombinationSRHLifetime[p.iphip, ireg] = p.τp[ireg]
        
        params.recombinationSRHTrapDensity[p.iphin, ireg] = p.nTrapDensity[ireg]
        params.recombinationSRHTrapDensity[p.iphip, ireg] = p.pTrapDensity[ireg]

        params.recombinationAuger[p.iphin, ireg] = p.Augn[ireg]
        params.recombinationAuger[p.iphip, ireg] = p.Augp[ireg]
    end

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[p.iphin, p.bregionJ1] = p.Nn[p.regionIntrinsic]
    params.bDensityOfStates[p.iphip, p.bregionJ1] = p.Np[p.regionIntrinsic]

    params.bDensityOfStates[p.iphin, p.bregionJ2] = p.Nn[p.regionIntrinsic]
    params.bDensityOfStates[p.iphip, p.bregionJ2] = p.Np[p.regionIntrinsic]

    params.bBandEdgeEnergy[p.iphin, p.bregionJ1] = p.En[p.regionIntrinsic]
    params.bBandEdgeEnergy[p.iphip, p.bregionJ1] = p.Ep[p.regionIntrinsic]

    params.bBandEdgeEnergy[p.iphin, p.bregionJ2] = p.En[p.regionIntrinsic]
    params.bBandEdgeEnergy[p.iphip, p.bregionJ2] = p.Ep[p.regionIntrinsic]

    ## no surface recombination
    # params.recombinationSRHvelocity[p.iphin, p.bregionJ1] = 1.0e1 * cm / s
    # params.recombinationSRHvelocity[p.iphip, p.bregionJ1] = 1.0e5 * cm / s

    # params.bRecombinationSRHTrapDensity[p.iphin, p.bregionJ1] = params.recombinationSRHTrapDensity[p.iphin, p.regionIntrinsic]
    # params.bRecombinationSRHTrapDensity[p.iphip, p.bregionJ1] = params.recombinationSRHTrapDensity[p.iphip, p.regionIntrinsic]

    # params.recombinationSRHvelocity[p.iphin, p.bregionJ2] = 1.0e7 * cm / s
    # params.recombinationSRHvelocity[p.iphip, p.bregionJ2] = 1.0e1 * cm / s

    # params.bRecombinationSRHTrapDensity[p.iphin, p.bregionJ2] = params.recombinationSRHTrapDensity[p.iphin, p.regionIntrinsic]
    # params.bRecombinationSRHTrapDensity[p.iphip, p.bregionJ2] = params.recombinationSRHTrapDensity[p.iphip, p.regionIntrinsic]

    ##############################################################

    ## interior doping
    params.doping[p.iphin, p.regionDonor] = p.Cn
    params.doping[p.iphip, p.regionAcceptor] = p.Cp
    params.doping[p.iphin, p.regionIntrinsic] = 1.0e20 / (m^3) # absorber n-doping
    params.doping[p.iphia, p.regionIntrinsic] = p.Ca

    ## no generation

    # parameter which passes the shift information in the Beer-Lambert generation
    # params.generationPeak = p.generationPeak

    ## generation parameters
    # params.generationIncidentPhotonFlux = p.incidentPhotonFlux
    # params.generationAbsorption = p.absorption
    # params.generationUniform = p.generation_uniform

    return params
end