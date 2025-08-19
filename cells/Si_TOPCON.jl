using ChargeTransport
using ExtendableGrids
using PyPlot
using CSV
using DataFrames

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
    "light-bias" => true,
    "light-oc" => false,
    "iv" => false
)

function main(;
        n = 6, Plotter = PyPlot, plotting = true,
        verbose = false, test = false,
        #parameter_file = "../parameter_files/Params_PSC_TiO2_MAPI_spiro.jl", # choose the parameter file
        parameter_file = "../params/Params_Si_TOPCON.jl",
    )

    if plotting
        Plotter.close("all")
    end

    println("--- Define physical parameters and model ---")

    include(parameter_file) # include the parameter file we specified

    ## contact voltage
    maxVoltage = 0.95 * V

    ## primary data for I-V scan protocol
    scanrate = 0.3 * V / s
    ntsteps = 151
    tend = maxVoltage / scanrate

    tvalues = range(0, stop = tend, length = ntsteps)

    println("--- Set up grid and regions ---")

    δ = 6 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_em = collect(range(0.0, 3*h_emitter, step = h_emitter / (0.8 * δ))) # p+ emitter

    # n C-Si region
    coord_cz_u = collect(range(3*h_emitter, 0.95 * h_cz, step = h_cz / (0.8 * δ)))
    coord_cz_g = geomspace(
        0.95 * h_cz, h_cz,
        (0.05 * h_cz) / (0.7 * δ), (0.05 * h_cz) / (1.1 * δ),
        tol = t
    )

    # n poly-Si region
    coord_poly_g = geomspace(
        h_cz, h_cz + h_poly / 2,
        h_poly / (1.3 * δ), h_poly / (0.6 * δ),
        tol = t
    )
    coord_poly_u = collect(range(h_cz + h_poly / 2, h_cz + h_poly, step = h_poly / (0.8 * δ)))

    coord = glue(coord_em, coord_cz_u, tol = 10 * t)
    coord = glue(coord, coord_cz_g, tol = 10 * t)
    coord = glue(coord, coord_poly_g, tol = 10 * t)
    coord = glue(coord, coord_poly_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    numberOfNodes = length(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0], [heightLayers[1]], regionCz, tol = 1.0e-18) # n C-Si
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionPoly, tol = 1.0e-18)  # n poly-Si

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], bregionCz, tol = 1.0e-18) # n C-Si outer boundary
    bfacemask!(grid, [h_total], [h_total], bregionPoly, tol = 1.0e-18)  # n poly-Si outer boundary
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18) # junction

    ## Plot node grid
    if toPlot["grid"]
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
        Plotter.show()
    end

    println("--- Define system and fill in information about model ---")

    ## set up generation data
    subg1 = subgrid(grid, [regionCz]); subg2 = subgrid(grid, [regionPoly]);

    generation_file = "simulation_data/scaps/si-topcon-auto.gen"
    generation_rate = generation_from_scaps(generation_file) # function to get generation rate from SCAPS file

    gen1 = generation_rate.(subg1[Coordinates]) # initialize generation in c-Si layer
    # gen2 = zeros(length(subg2[Coordinates]) - 1) # set absorption in poly layer to zero
    gen2 = generation_rate.(subg2[Coordinates]) # initialize generation in poly-Si layer

    generationData = [gen1'; gen2']

    ## Initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers, generationData = generationData)

    data.modelType = Transient # choices: Transient, Stationary
    carrier_stats = Boltzmann # TODO: choices
    data.F = [carrier_stats, carrier_stats]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    data.generationModel = GenerationUserDefined

    data.boundaryType[bregionPoly] = SchottkyContact
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionCz] = SchottkyContact

    data.fluxApproximation .= ExcessChemicalPotential

    println("--- Define Params ---")

    params = Params(grid, numberOfCarriers)
    paramsnodal = ParamsNodal(grid, numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg] = ε[ireg] * ε0

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn[ireg]
        params.densityOfStates[iphip, ireg] = Np[ireg]

        params.bandEdgeEnergy[iphin, ireg] = En[ireg]
        params.bandEdgeEnergy[iphip, ireg] = Ep[ireg]

        params.mobility[iphin, ireg] = μn[ireg]
        params.mobility[iphip, ireg] = μp[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg] = τcz[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τpoly[ireg]

        params.recombinationSRHTrapDensity[iphin, ireg] = 1.0e14 / (m^3)
        params.recombinationSRHTrapDensity[iphip, ireg] = 1.0e14 / (m^3)

        params.recombinationAuger[iphin, ireg] = Augn[ireg]
        params.recombinationAuger[iphip, ireg] = Augp[ireg]
    end

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[iphin, bregionJ1] = Nn[regionPoly]
    params.bDensityOfStates[iphip, bregionJ1] = Np[regionPoly]

    params.bBandEdgeEnergy[iphin, bregionJ1] = En[regionPoly]
    params.bBandEdgeEnergy[iphip, bregionJ1] = Ep[regionPoly]

    # Schottky Barrier
    # TODO: just picked these values - need to justify
    params.SchottkyBarrier[bregionCz] = 1.105 * (eV)
    params.SchottkyBarrier[bregionPoly] = -0.018 * (eV)

    params.bVelocity[iphin, bregionCz] = 1.0e-2 * cm/s
    params.bVelocity[iphip, bregionCz] = 1.0e7 * cm/s

    params.bVelocity[iphin, bregionPoly] = 1.0e7 * cm/s
    params.bVelocity[iphip, bregionPoly] = 1.0e-2 * cm/s

    ## Positive doping corresponds to acceptors
    for icoord in 1:numberOfNodes
        if icoord <= (length(coord_em) + length(coord_cz_u) + length(coord_cz_g) - 2) # n C-Si region
            paramsnodal.doping[icoord] = tanh_diffusion_profile(coord[icoord]; 
                x0 = 0.0,
                x1 = h_emitter,
                ymin = -Ccz,
                ymax = Cem,
                sharpness = 6.0
            )
        else
            paramsnodal.doping[icoord] = -Cpoly
        end
    end

    data.params = params
    data.paramsnodal = paramsnodal

    if plotting
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
    end

    ctsys = System(grid, data, unknown_storage = :sparse)

    println("--- Define control parameters for solver ---")

    control = SolverControl()
    control.verbose = verbose
    control.maxiters = 1000
    control.damp_initial = 0.3
    control.damp_growth = 1.11 # >= 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    solution = equilibrium_solve!(ctsys, control = control)
    inival = solution

    save_cell_profile("simulation_data/chargetransport/si-topcon-schottky-dark-sc.csv", solution, ctsys)

    ipsi = ctsys.fvmsys.physics.data.index_psi
    Vbi = solution[ipsi, end] - solution[ipsi, 1]
    println("\nBuilt-in Voltage: $(Vbi)V\n")

    ### D: ILLUMINATION

    I = collect(20:-0.2:0.0)
    LAMBDA = 10 .^ (-I)

    for istep in 1:(length(I) - 1)

        ## turn slowly generation on
        ctsys.data.λ2 = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop

    # Plot illuminated short-circuit
    if plotting && toPlot["light-sc"]
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_density)
        Plotter.show()
    end

    save_cell_profile("simulation_data/chargetransport/si-topcon-schottky-illuminated-sc.csv", solution, ctsys)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    ## for saving I-V data
    IV = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values
    VocExceeded = [false, false] # first term for if scaps Voc exceeded, second for if current < 0

    if toPlot["light-bias"]
        Plotter.rc("figure", figsize=(18, 8))
        fig = Plotter.figure()
    end
    figure = Plotter.figure()
    for istep in 2:ntsteps
        t = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep - 1] # Time step size

        ## Apply new voltage (set non-equilibrium values)
        set_contact!(ctsys, bregionCz, Δu = Δu)

        if test == false
            println("time value: Δt = $(t)")
        end

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)
        inival = solution

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)


        if plotting && toPlot["light-bias"]
            Plotter.clf()
            
            fig.suptitle("Illuminated, Bias = $(round(Δu, digits=2))", fontsize=16)

            Plotter.subplot(1, 2, 1)
            plot_energies(Plotter, ctsys, solution, "Band Diagram", label_energy, clear=false)
            Plotter.subplot(1, 2, 2)
            plot_densities(Plotter, ctsys, solution, "Carrier Densities", label_density, clear=false)

            fig.tight_layout()
            Plotter.pause(0.1)
            display(gcf())
        end

        if Δu >= 0.7400 && VocExceeded[1] == false
            VocExceeded[1] = true
            save_cell_profile("simulation_data/chargetransport/si-topcon-schottky-illuminated-scaps-oc.csv", solution, ctsys)
        end

        if current < 0.0 && VocExceeded[2] == false
            VocExceeded[2] = true
            # Plot bands and carrier densities at Voc
            save_cell_profile("simulation_data/chargetransport/si-topcon-schottky-illuminated-ct-oc.csv", solution, ctsys)
            println("Graph plotted at V = $(Δu)")

            if plotting && toPlot["light-oc"]
                Plotter.rc("figure", figsize=(18, 8))
                fig = Plotter.figure()

                fig.suptitle("Illuminated, Bias = $(round(Δu, digits=2))", fontsize=16)

                Plotter.subplot(1, 2, 1)
                plot_energies(Plotter, ctsys, solution, "Band Diagram", label_energy, clear=false)
                Plotter.subplot(1, 2, 2)
                plot_densities(Plotter, ctsys, solution, "Carrier Densities", label_density, clear=false)

                fig.tight_layout()
                Plotter.show()
            end
        end

        push!(IV, current)
        push!(biasValues, Δu)
    end # time loop

    # Plot IV curve
    if plotting && toPlot["iv"]
        plot_IV(Plotter, biasValues, -IV, "bias \$\\Delta u\$ = $(maxVoltage)")
        show()
    end
    
    save_iv("simulation_data/chargetransport/si-topcon-schottky-iv.csv", biasValues, IV)

    powerDensity = biasValues .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    Voc = compute_open_circuit_voltage(biasValues, IV)

    IncidentLightPowerDensity = 1000.0 * W / m^2

    efficiency = MaxPD / IncidentLightPowerDensity
    fillfactor = (biasValues[indexPD] * IV[indexPD]) / (IV[1] * Voc)

    println("\nIsc = $(IV[1])")
    println("Voc = $(Voc)")
    println("FF = $(fillfactor)")
    println("PCE = $(efficiency)")
end

# DEBUG
main(n = 12, verbose = false)