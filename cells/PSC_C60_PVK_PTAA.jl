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

# import utilities
include("../utils/ct_utils.jl")
using .CTUtils

# Plotting with scienceplots
using PyCall
pyimport("scienceplots")
mpl = PyPlot.matplotlib
plt.style.use(["science","nature"])
PyPlot.rc("figure", figsize=(8, 8))
PyPlot.rc("axes", titlesize=16)
PyPlot.rc("axes", labelsize=14)
PyPlot.rc("xtick", labelsize=12)
PyPlot.rc("ytick", labelsize=12)
PyPlot.rc("legend", fontsize=12)
ENV["MPLBACKEND"] = "qt5agg"

# Modify this dictionary to choose what to plot
toPlot = Dict(
    "grid" => false,
    "generation" => false,
    "dark-sc" => false,
    "light-sc" => false,
    "light-bias" => false,
    "light-oc" => false,
    "iv" => true
)

parameter_file = "../params/Params_PSC_C60_PVK_PTAA.jl"
include(parameter_file)

# you can also use other Plotters, if you add them to the example file
function main(;
        n = 6, Plotter = PyPlot, plotting = false,
        verbose = false, test = false,
        parameter_file = "../params/Params_PSC_C60_PVK_PTAA.jl", # choose the parameter file
        BeerLambertGeneration = true,
        incidentPhotonFlux = 4.3e21 / (m^2 * s)
    )

    println("--- Define physical parameters and model ---")

    include(parameter_file) # include the parameter file we specified

    ## contact voltage
    maxVoltage = 1.3 * V

    ## primary data for I-V scan protocol
    scanrate = 0.3 * V / s
    ntsteps = 91
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

    if plotting && toPlot["grid"]
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
        Plotter.show()
    end

    println("--- Define system and fill in information about model ---")

    ## Initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers)

    ## Possible choices: Stationary, Transient
    data.modelType = Transient

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    if !BeerLambertGeneration
        data.generationModel = GenerationUniform
    else
        data.generationModel = GenerationBeerLambert
    end

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionJ2] = InterfaceRecombination
    data.boundaryType[bregionDonor] = OhmicContact

    ## Present ionic vacancies in perovskite layer
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    data.fluxApproximation .= ExcessChemicalPotential

    println("--- Define Params ---")

    params = Params(grid, numberOfCarriers)

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

        if BeerLambertGeneration
            params.generationAbsorption[ireg] = absorption[ireg]
        end
    end

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

    ## for surface recombination
    params.recombinationSRHvelocity[iphin, bregionJ1] = 1.0e1 * cm / s
    params.recombinationSRHvelocity[iphip, bregionJ1] = 1.0e5 * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    params.recombinationSRHvelocity[iphin, bregionJ2] = 1.0e7 * cm / s
    params.recombinationSRHvelocity[iphip, bregionJ2] = 1.0e1 * cm / s

    params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = params.recombinationSRHTrapDensity[iphin, regionIntrinsic]
    params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = params.recombinationSRHTrapDensity[iphip, regionIntrinsic]

    ##############################################################

    ## interior doping
    params.doping[iphin, regionDonor] = Cn
    params.doping[iphip, regionAcceptor] = Cp
    params.doping[iphia, regionIntrinsic] = Ca

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)

    println("--- Define control parameters for solver ---")

    control = SolverControl()
    control.verbose = verbose
    control.damp_initial = 0.5
    control.damp_growth = 1.21 # >= 1
    control.maxiters = 1000

    println("--- Solve in equilibrium ---")

    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    ### D: ILLUMINATION

    I = collect(20:-1:0.0)
    LAMBDA = 10 .^ (-I)

    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJ2] = 0.0

    for istep in 1:(length(I) - 1)

        ## turn slowly generation on
        ctsys.data.λ2 = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop

    if plotting && toPlot["light-sc"]
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
        label_energy[1, iphia] = "\$E_a-q\\psi\$"; label_energy[2, iphia] = "\$ - q \\varphi_a\$"; label_BEE[iphia] = "\$E_a\$"
        label_density[iphia] = "\$ n_a \$";      label_solution[iphia] = "\$ \\varphi_a\$"

        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_density)
        Plotter.show()
    end

    println("--- IV Curve ---")

    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 0.0
    ctsys.fvmsys.boundary_values[iphia, bregionJ2] = 0.0

    ## for saving I-V data
    currents = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values

    for istep in 2:ntsteps

        t = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep - 1] # Time step size

        ## Apply new voltage (set non-equilibrium values)
        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: Δt = $(t)")
        end

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)
        inival = solution

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)

        push!(currents, current)
        push!(biasValues, Δu)
    end 

    currents = -currents

    if plotting && toPlot["iv"]
        plot_IV(Plotter, biasValues, -currents, "bias \$\\Delta u\$ = $(maxVoltage)")
        show()
    end

    return IV(biasValues, currents)
end

end