#=

Code for simulation of a Perovskite solar cell

TODO: Add the ETL, HTL, and Absorber specs

=#

module PSC_2

using ChargeTransport
using ExtendableGrids
using PyPlot
using CSV
using DataFrames

## import utilities
include("../utils/ct_utils.jl")
using .CTUtils

## modify this dictionary to choose what to plot
toPlot = Dict(
    "grid" => false,
    "generation" => false,
    "dark-sc" => false,
    "light-sc" => false,
    "light-bias" => false,
    "light-oc" => false,
    "iv" => false
)

parameter_file = "../params/Params_PSC_C60_PVK_PTAA.jl"
include(parameter_file)

function main(;
        n = 6, Plotter = PyPlot, plotting = false,
        verbose = false, 
        BeerLambertGeneration = true,
        incidentPhotonFlux = 4.3e21 / (m^2 * s)
    )

    println("--- Define physical parameters and model ---")

    include(parameter_file) # include the parameter file we specified

    ## max contact voltage during I-V scan
    maxVoltage = 1.3 * V

    ## primary data for I-V scan protocol
    scanrate = 0.3 * V / s
    ntsteps = 91 # number of time steps
    tend = maxVoltage / scanrate

    tvalues = range(0, stop = tend, length = ntsteps)

    println("--- Set up grid and regions ---")

    δ = 4 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_n_u = collect(range(0.0, h_ndoping / 2, step = h_ndoping / (0.8 * δ)))
    coord_n_g = geomspace(
        h_ndoping / 2, h_ndoping,
        h_ndoping / (0.7 * δ), h_ndoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        h_ndoping, h_ndoping + h_intrinsic / k,
        h_intrinsic / (5.1 * δ), h_intrinsic / (1.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        h_ndoping + h_intrinsic / k, h_ndoping + h_intrinsic,
        h_intrinsic / (1.1 * δ), h_intrinsic / (5.1 * δ),
        tol = t
    )
    coord_p_g = geomspace(
        h_ndoping + h_intrinsic, h_ndoping + h_intrinsic + h_pdoping / 2,
        h_pdoping / (1.3 * δ), h_pdoping / (0.6 * δ),
        tol = t
    )
    coord_p_u = collect(range(h_ndoping + h_intrinsic + h_pdoping / 2, h_ndoping + h_intrinsic + h_pdoping, step = h_pdoping / (0.8 * δ)))

    coord = glue(coord_n_u, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_p_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0 * μm], [heightLayers[1]], regionDonor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionAcceptor, tol = 1.0e-18)  # p-doped region   = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], bregionDonor, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionAcceptor, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2, tol = 1.0e-18) # second inner interface

    ## plot node grid
    if plotting && toPlot["grid"]
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
        Plotter.show()
    end

    println("--- Define system and fill in information about model ---")

    ## initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers)

    ## choose simulation type
    ## possible choices: Stationary, Transient
    data.modelType = Transient

    ## choose statistics model
    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    ## set relevant recombination mechanisms
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Beer-Lambert or uniform generation
    ## more complex generation models can be implemented via the user-defined generation model
    if !BeerLambertGeneration
        data.generationModel = GenerationUniform
    else
        data.generationModel = GenerationBeerLambert
    end

    ## set interface types
    ## possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionJ2] = InterfaceRecombination
    data.boundaryType[bregionDonor] = OhmicContact

    ## present ionic vacancies in perovskite layer
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])
    
    ## set flux discretization scheme
    ## possible choices: DiffusionEnhanced, DiffusionEnhancedModifiedDrift, ExcessChemicalPotential,
    ## ExcessChemicalPotentialGraded, GeneralizedSG, ScharfetterGummel, ScharfetterGummelGraded
    data.fluxApproximation .= ExcessChemicalPotential

    println("--- Define Params ---")

    ## create params instance - contains all parameters for the cell
    params = Params(numberOfRegions, numberOfRegions + 1, numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp
    params.chargeNumbers[iphia] = za

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg] = ε[ireg] * ε0

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn[ireg]
        params.densityOfStates[iphip, ireg] = Np[ireg]
        params.densityOfStates[iphia, ireg] = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg] = En[ireg]
        params.bandEdgeEnergy[iphip, ireg] = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg] = Ea[ireg]

        params.mobility[iphin, ireg] = μn[ireg]
        params.mobility[iphip, ireg] = μp[ireg]
        params.mobility[iphia, ireg] = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]

        params.recombinationSRHTrapDensity[iphin, ireg] = nTrapDensity[ireg]
        params.recombinationSRHTrapDensity[iphip, ireg] = pTrapDensity[ireg]

        params.recombinationAuger[iphin, ireg] = Augn[ireg]
        params.recombinationAuger[iphip, ireg] = Augp[ireg]

        if BeerLambertGeneration
            params.generationAbsorption[ireg] = absorption[ireg]
        end
    end

    ## set photon flux for each region
    ## assume only the absorber can absorb light, and CTLs are perfectly transparent
    params.generationIncidentPhotonFlux = [0.0, incidentPhotonFlux, 0.0]

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[iphin, bregionJ1] = Nn[regionIntrinsic]
    params.bDensityOfStates[iphip, bregionJ1] = Np[regionIntrinsic]

    params.bDensityOfStates[iphin, bregionJ2] = Nn[regionIntrinsic]
    params.bDensityOfStates[iphip, bregionJ2] = Np[regionIntrinsic]

    params.bBandEdgeEnergy[iphin, bregionJ1] = En[regionIntrinsic]
    params.bBandEdgeEnergy[iphip, bregionJ1] = Ep[regionIntrinsic]

    params.bBandEdgeEnergy[iphin, bregionJ2] = En[regionIntrinsic]
    params.bBandEdgeEnergy[iphip, bregionJ2] = Ep[regionIntrinsic]

    ## surface recombination velocities
    params.recombinationSRHvelocity[iphin, bregionJ1] = 1.0e1 * cm / s
    params.recombinationSRHvelocity[iphip, bregionJ1] = 1.0e5 * cm / s

    params.recombinationSRHvelocity[iphin, bregionJ2] = 1.0e7 * cm / s
    params.recombinationSRHvelocity[iphip, bregionJ2] = 1.0e1 * cm / s

    ## set interface trap densities
    params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    ##############################################################

    ## interior doping
    params.doping[iphin, regionDonor] = Cn # ETL doping
    params.doping[iphip, regionAcceptor] = Cp # HTL doping
    params.doping[iphia, regionIntrinsic] = Ca # initial anion concentration

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)

    println("--- Define control parameters for solver ---")

    ## solver parameters - tweak when convergence is an issue
    control = SolverControl()
    control.verbose = verbose
    control.damp_initial = 0.5
    control.damp_growth = 1.21 # >= 1
    control.maxiters = 1000

    println("--- Solve in equilibrium ---")

    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    ## set axis labels for plots if plotting is set to 'on'
    if plotting
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
        label_density[iphia] = "\$ n_a \$";      label_solution[iphia] = "\$ \\varphi_a\$"
    end

    ## plot carrier densities and energies in short-circuit in the dark
    if plotting && toPlot["dark-sc"]
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Dark Short-Circuit", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Dark Short-Circuit", label_density)
        Plotter.show()
    end

    ## array which defines light intensity at each step as we ramp up illumination
    I = collect(20:-1:0.0)
    LAMBDA = 10 .^ (-I)

    ## set Neumann boundary conditions for anions - otherwise simulation won't know what potential to put anions at
    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJ2] = 0.0

    ## ramp up light intensity
    for istep in 1:(length(I) - 1)
        ## ramp up generation
        ctsys.data.λ2 = LAMBDA[istep + 1]

        println("increase generation with λ2 = $(data.λ2)")

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution
    end

    ## plot carrier densities and energies in short-circuit under illumination
    if plotting && toPlot["light-sc"]
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_density)
        Plotter.show()
    end

    println("--- IV Curve ---")

    ## turn off Neumann boundary conditions for anions
    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 0.0
    ctsys.fvmsys.boundary_values[iphia, bregionJ2] = 0.0

    ## for saving I-V data
    currents = zeros(0) # for current densities (A m^{-2})
    biasValues = zeros(0) # for bias values (V)

    for istep in 2:ntsteps

        t = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep - 1] # Time step size

        ## Apply new voltage (set non-equilibrium values)
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        println("time value: Δt = $(t), bias: Δu = $(Δu)")

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)
        inival = solution

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)

        push!(currents, current)
        push!(biasValues, Δu)
    end 

    ## flip currents - use this based on cell architecture
    currents = -currents

    ## plot IV curve
    if plotting && toPlot["iv"]
        plot_IV(Plotter, biasValues, -currents, "bias \$\\Delta u\$ = $(maxVoltage)")
        show()
    end

    ## return IV data in an IV struct
    return IV(biasValues, currents)
end

end