using LessUnitful: @local_unitfactors, @ufac_str, @ph_str
using ChargeTransport

@kwdef struct Params_PSC_C60_PVK_NiO_graded

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
    regionETL1 = 1
    regionETL2 = 2
    regionETL3 = 3
    regionPerovskite = 4
    regionHTL = 5
    regions = [regionETL1, regionETL2, regionETL3, regionPerovskite, regionHTL]
    numberOfRegions = length(regions)

    # boundary region numbers
    bregionETL = 1 # ETL-contact
    bregionHTL = 2 # HTL-contact
    bregionJ1 = 3 # ETL1-ETL2
    bregionJ2 = 4 # ETL2-ETL3
    bregionJ3 = 5 # ETL3-PVK
    bregionJ4 = 6 # PVK-HTL
    numberOfBoundaryRegions = 6

    ## length domains
    h_etl1 = 10.0 * nm
    h_etl2 = 5.0 * nm
    h_etl3 = 5.0 * nm
    h_pvk = 100.0 * nm
    h_htl = 13.0 * nm
    h_total = h_etl1 + h_etl2 + h_etl3 + h_pvk + h_htl
    heightLayers = [
        h_etl1,
        h_etl1 + h_etl2,
        h_etl1 + h_etl2 + h_etl3,
        h_etl1 + h_etl2 + h_etl3 + h_pvk,
        h_etl1 + h_etl2 + h_etl3 + h_pvk + h_htl,
    ]

    ########## physical values ##########

    ## charge numbers
    zn = -1
    zp = 1
    za = 1

    ## temperature
    T = 300.0 * K

    ## band edge energies
    En = [-4.00, -4.00, -4.00, -4.00, -2.05] .* eV
    Ep = [-6.00, -6.00, -6.00, -5.68, -5.65] .* eV
    Ea = [0.0, 0.0, 0.0, 0.0, 0.0] .* eV

    ## effective densities of density of states
    Nn = [1.0e26, 1.0e26, 1.0e26, 2.75e24, 2.5e26] ./ (m^3)
    Np = [1.0e26, 1.0e26, 1.0e26, 3.80e24, 2.5e26] ./ (m^3)
    Na = [0.0, 0.0, 0.0, 1.21e28, 0.0] ./ (m^3)

    ## mobilities
    μn = [8.0e-6, 8.0e-6, 8.0e-6, 1.0e-3, 2.8e-4] .* (m^2) / (V * s)
    μp = [3.5e-7, 3.5e-7, 3.5e-7, 1.0e-3, 2.8e-4] .* (m^2) / (V * s)
    μa = [0.0, 0.0, 0.0, 1.0e-12, 0.0] .* (m^2) / (V * s)

    ## relative dielectric permittivity
    ε = [2.5, 10.0, 20.0, 23.13, 11.7] .* 1.0

    ## radiative recombination
    r0 = [1.0e-20, 1.0e-20, 1.0e-20, 1.0e-20, 1.0e-20] .* cm^3 / s

    ## life times and trap densities
    τn = [1.0e-7, 1.0e-7, 1.0e-7, 2.5e-9, 1.0e-6] .* s
    τp = [1.0e-7, 1.0e-7, 1.0e-7, 2.5e-9, 1.0e-6] .* s

    Augn = [1.0e-31, 1.0e-31, 1.0e-31, 1.0e-31, 1.0e-31] .* (cm^6) / (s)
    Augp = [1.0e-31, 1.0e-31, 1.0e-31, 1.0e-31, 1.0e-31] .* (cm^6) / (s)

    ## SRH trap densities
    nTrapDensity = [1.0e20, 1.0e20, 1.0e20, 4.0e21, 1.0e20] ./ (m^3)
    pTrapDensity = [1.0e20, 1.0e20, 1.0e20, 4.0e21, 1.0e20] ./ (m^3)

    ## doping
    dop_n = [1.0e24, 1.0e24, 1.0e24, 1.0e20, 0.0] ./ (m^3)
    dop_p = [0.0, 0.0, 0.0, 0.0, 3.0e24] ./ (m^3)
    dop_vac = [0.0, 0.0, 0.0, 1.0e24, 0.0] ./ (m^3)
end

function Params(p)

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

        params.doping[p.iphin, ireg] = p.dop_n[ireg]
        params.doping[p.iphip, ireg] = p.dop_p[ireg]
        params.doping[p.iphia, ireg] = p.dop_vac[ireg]
    end

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    for bregion in [p.bregionJ1, p.bregionJ2]
        params.bDensityOfStates[p.iphin, bregion] = p.Nn[p.regionETL2]
        params.bDensityOfStates[p.iphip, bregion] = p.Np[p.regionETL2]

        params.bBandEdgeEnergy[p.iphin, bregion] = p.En[p.regionETL2]
        params.bBandEdgeEnergy[p.iphip, bregion] = p.Ep[p.regionETL2]
    end

    for bregion in [p.bregionJ3, p.bregionJ4]
        params.bDensityOfStates[p.iphin, bregion] = p.Nn[p.regionPerovskite]
        params.bDensityOfStates[p.iphip, bregion] = p.Np[p.regionPerovskite]

        params.bBandEdgeEnergy[p.iphin, bregion] = p.En[p.regionPerovskite]
        params.bBandEdgeEnergy[p.iphip, bregion] = p.Ep[p.regionPerovskite]
    end

    return params
end