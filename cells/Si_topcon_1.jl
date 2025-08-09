ENV["MPLBACKEND"] = "qt5agg"

using ChargeTransport
using ExtendableGrids
using PyPlot
using CSV
using DataFrames

# Plotting with scienceplots
using PyCall
pyimport("scienceplots")
mpl = PyPlot.matplotlib
plt.style.use(["science","nature"])
PyPlot.rc("figure", figsize=(8, 8))
PyPlot.rc("axes", titlesize=16)       # Title of plots
PyPlot.rc("axes", labelsize=14)       # Axis labels
PyPlot.rc("xtick", labelsize=12)      # X-tick labels
PyPlot.rc("ytick", labelsize=12)      # Y-tick labels
PyPlot.rc("legend", fontsize=12)      # Legend text

# you can also use other Plotters, if you add them to the example file
function main(;
        n = 3, Plotter = PyPlot, plotting = true,
        verbose = false,
        parameter_file = "../params/Params_Si_topcon_1.jl"
    )

    if plotting
        Plotter.close("all")
    end

    ### DEBUG
    if test == false
        println("Set up grid and regions")
    end

    include(parameter_file)

    voltageAcceptor = 0.9 * V

    scanrate = 0.3 * V / s
    ntsteps = 91
    vend = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend = vend / scanrate

    tvalues = range(0, stop = tend, length = ntsteps)

    ### DEBUG
    if test == false
        println("--- done ---\n")
        println("\nSet up grid and regions")
    end

    # Constant node density grid
    nodedensity = 100 / (cm)
    
    numberOfNodes = Int(h_total * nodedensity)
    coord = collect(range(0.0, stop = h_total, length = numberOfNodes))
    grid = ExtendableGrids.simplexgrid(coord)

    cellmask!(grid, [0.0 * μm], [heightLayers[1]], regionDonor, tol = 1.0e-18) 
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18)
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionAcceptor, tol = 1.0e-18)

    bfacemask!(grid, [0.0], [0.0], bregionDonor, tol = 1.0e-18)    
    bfacemask!(grid, [h_total], [h_total], bregionAcceptor, tol = 1.0e-18)
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18)
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2, tol = 1.0e-18)

    if test == false
        println("--- done ---\n")
        println("\nDefine System and fill in information about model")
    end

    ## Initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers)

    data.modelType = Transient
    carrier_stats = Boltzmann
    data.F = [carrier_stats, carrier_stats]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    # TODO: Would be better to setup user-defined model - uniform in TLs and then Beer-Lambert in absorber

    data.generationModel = GenerationUniform

    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionJ2] = InterfaceRecombination
    data.boundaryType[bregionDonor] = OhmicContact

    data.fluxApproximation .= ExcessChemicalPotential

    if test == false
        println("--- done ---\n")
        println("\nDefine Params and fill in physical parameters")
    end

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
        params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]

        # ## TODO: Trap densities

        # if ireg == regionAcceptor || ireg == regionDonor
        #     params.recombinationSRHTrapDensity[iphin, ireg] = 0.0 / (m^3)
        #     params.recombinationSRHTrapDensity[iphip, ireg] = 0.0 / (m^3)
        # else
        #     params.recombinationSRHTrapDensity[iphin, ireg] = 2.0e21 / (m^3)
        #     params.recombinationSRHTrapDensity[iphip, ireg] = 2.0e21 / (m^3)
        # end
        
        # params.generationIncidentPhotonFlux[ireg] = incidentPhotonFlux[ireg]
        # params.generationAbsorption[ireg] = absorption[ireg]

        params.generationUniform[ireg] = generation_uniform[ireg]
    end

    # params.generationPeak = generationPeak

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

    # TODO: Recombination parameters for interfaces

    # min_SRV = 1.0e3 * cm/s
    # maj_SRV = 1.0e7 * cm/s

    # ## for surface recombination
    # params.recombinationSRHvelocity[iphin, bregionJ1] = maj_SRV * 5.0
    # params.recombinationSRHvelocity[iphip, bregionJ1] = min_SRV * 5.0

    # params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = 2.0e18 / (m^3)
    # params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = 2.0e18 / (m^3)

    # params.recombinationSRHvelocity[iphin, bregionJ2] = min_SRV
    # params.recombinationSRHvelocity[iphip, bregionJ2] = maj_SRV

    # params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = 1.0e17 / (m^3)
    # params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = 1.0e17 / (m^3)

    ##############################################################

    ## interior doping
    for icoord in 1:numberOfNodes
        paramsnodal.doping[icoord] = 
    end

    # TODO: Determine contact work functions - 4.26eV for Silver on either side, but is it necessary?
    # Contact work functions
    # ctsys.fvmsys.boundary_factors[data.index_psi, bregionDonor] = 1.0e30
    # ctsys.fvmsys.boundary_values[data.index_psi, bregionDonor] = -4.1

    # ctsys.fvmsys.boundary_factors[data.index_psi, bregionAcceptor] = 1.0e30
    # ctsys.fvmsys.boundary_values[data.index_psi, bregionAcceptor] = -5.4

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    control = SolverControl()
    control.verbose = verbose
    control.maxiters = 1000
    control.damp_initial = 0.9
    control.damp_growth = 1.61 # >= 1
    control.max_round = 5

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

    if plotting == true
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        # Plot dark equilibrium
        # Plotter.figure()
        # plot_energies(Plotter, ctsys, solution, "Equilibrium, Dark", label_energy)
        # Plotter.figure()
        # plot_densities(Plotter, ctsys, solution, "Equilibrium, Dark", label_density)
        # Plotter.show()
    end

    ## DARK EQUILIBRIUM
    ipsi = ctsys.fvmsys.physics.data.index_psi
    Vbi = solution[ipsi, end] - solution[ipsi, 1]
    println("Built-in Voltage: $(Vbi)V")

    # Plot equilibrium voltage profile
    # if plotting
    #     Plotter.figure()
    #     plot_solution(Plotter, ctsys, solution, "Equilibrium, Dark", label_solution)

    #     Plotter.show()
    # end

    # save_device_profile_csv("si_1_dark_sc.csv", solution, ctsys)
    # exit()
    
    solution_dark = solution

    ### D: ILLUMINATION

    I = collect(20:-0.5:0.0)
    LAMBDA = 10 .^ (-I)

    # ctsys.fvmsys.boundary_factors[data.index_psi, bregionDonor] = 0
    # ctsys.fvmsys.boundary_values[data.index_psi, bregionDonor] = 0

    # ctsys.fvmsys.boundary_factors[data.index_psi, bregionAcceptor] = 0
    # ctsys.fvmsys.boundary_values[data.index_psi, bregionAcceptor] = 0

    for istep in 1:(length(I) - 1)

        ## turn slowly generation on
        ctsys.data.λ2 = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop

    save_device_profile_csv("si_1_ill_sc.csv", solution, ctsys)
    exit()

    # Plot illuminated short-circuit
    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_density)
        Plotter.show()
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    Isc = -get_current_val(ctsys, solution)

    ## for saving I-V data
    IV = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values
    tested = false

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

        push!(IV, current)
        push!(biasValues, Δu)

        if isapprox(current, 0, atol = 10) && current > 0 && tested == false && plotting
            println("CURRENT: $(current)")
            # determine QFLS near Voc
            subg = subgrid(grid, [regionIntrinsic])
            EFn = view(solution[iphin, :], subg)
            EFp = view(solution[iphip, :], subg)

            QFLS = abs(minimum(EFn) - maximum(EFp))

            # println("EF_n: $(EFn)")
            # println("EF_p: $(EFp)")
            println("QFLS: $(QFLS)")

            # carrier densities at absorber edges 
            e_dens = get_density(solution, regionIntrinsic, ctsys, iphin)
            h_dens = get_density(solution, regionIntrinsic, ctsys, iphip)

            println("\n")
            println("PVK-ETL INTERFACE")
            println("Electron Concentration: $(e_dens[1])")
            println("Hole Concentration: $(h_dens[1])")
            println("\n")
            println("PVK-HTL INTERFACE")
            println("Electron Concentration: $(e_dens[end])")
            println("Hole Concentration: $(h_dens[end])")

            Plotter.figure()
            plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Δu)", label_energy)
            Plotter.figure()
            plot_densities(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Δu)", label_density)
            # Plotter.ylim(1e7, 1e19)            
            show()
            tested = true
        end

    end # time loop

    plot_IV(Plotter, biasValues, IV, "bias \$\\Delta u\$ = $(vend)")
    show()

    IV = -IV

    powerDensity = biasValues .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    Voc = compute_open_circuit_voltage(biasValues, IV)

    IncidentLightPowerDensity = 1000.0 * W / m^2

    efficiency = biasValues[indexPD] * IV[indexPD] / IncidentLightPowerDensity
    fillfactor = (biasValues[indexPD] * IV[indexPD]) / (IV[1] * Voc)

    println("Isc = $(Isc))")
    println("Voc = $(Voc)")
    println("FF = $(fillfactor)")
    println("PCE = $(efficiency)")

    if test == false
        println("*** done\n")
    end

    testval = sum(filter(!isnan, solution)) / length(solution) # when using sparse storage, we get NaN values in solution
    return testval

end

# DEBUG
main()