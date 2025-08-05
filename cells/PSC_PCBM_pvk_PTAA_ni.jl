#=

Code for simulation of a perovskite solar cell with the following layers:
- ETL: PCBM
- Absorber: Cs15FA85Pb(I60Br40)3
- HTL: PTAA

This is an attempt to reproduce the results of the paper:
https://doi.org/10.1038/s41467-023-36141-8

=#

ENV["MPLBACKEND"] = "qt5agg"

using ChargeTransport
using ExtendableGrids
using PyPlot

# Plotting via SciencePlots framework
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

# helper function for determining median energy of a band
function middle_element(A)
    idx = (length(A) + 1) ÷ 2
    return A[idx]
end

function carrier_densities(solution, ctsys, regionIntrinsic, iphin, iphip)
    e_dens = get_density(solution, regionIntrinsic, ctsys, iphin)
    h_dens = get_density(solution, regionIntrinsic, ctsys, iphip)

    println("\nINTERFACE CARRIER DENSITIES")
    println("PVK-HTL INTERFACE")
    println("Electron Concentration: $(e_dens[1])")
    println("Hole Concentration: $(h_dens[1])")
    println("\nPVK-ETL INTERFACE")
    println("Electron Concentration: $(e_dens[end])")
    println("Hole Concentration: $(h_dens[end])\n")
end

function main(;
        n = 6, Plotter = PyPlot, plotting = true,
        verbose = false, test = false,
        parameter_file = "params/Params_PSC_PCBM_pvk_PTAA.jl",
    )

    if plotting
        Plotter.close("all")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    include(parameter_file)

    ## contact voltage
    voltageAcceptor = 1.5 * V

    ## primary data for I-V scan protocol
    scanrate = 0.3 * V / s
    ntsteps = 201
    vend = voltageAcceptor
    tend = vend / scanrate
    tvalues = range(0, stop = tend, length = ntsteps)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## I have used the same adaptive grid as in the PSC examples
    δ = 6 * n       
    t = 0.5 * (cm) / δ 
    k = 1.5  

    coord_p_u = collect(range(0.0, h_pdoping / 2, step = h_pdoping / (0.8 * δ)))
    coord_p_g = geomspace(
        h_pdoping / 2, h_pdoping,
        h_pdoping / (0.7 * δ), h_pdoping / (1.1 * δ),
        tol = t
    )
    coord_i_g1 = geomspace(
        h_pdoping, h_pdoping + h_intrinsic / k,
        h_intrinsic / (5.1 * δ), h_intrinsic / (1.1 * δ),
        tol = t
    )
    coord_i_g2 = geomspace(
        h_pdoping + h_intrinsic / k, h_pdoping + h_intrinsic,
        h_intrinsic / (1.1 * δ), h_intrinsic / (5.1 * δ),
        tol = t
    )
    coord_n_g = geomspace(
        h_pdoping + h_intrinsic, h_pdoping + h_intrinsic + h_ndoping / 2,
        h_ndoping / (1.3 * δ), h_ndoping / (0.6 * δ),
        tol = t
    )
    coord_n_u = collect(range(h_pdoping + h_intrinsic + h_ndoping / 2, h_pdoping + h_intrinsic + h_ndoping, step = h_ndoping / (0.8 * δ)))

    coord = glue(coord_p_u, coord_p_g, tol = 10 * t)
    coord = glue(coord, coord_i_g1, tol = 10 * t)
    coord = glue(coord, coord_i_g2, tol = 10 * t)
    coord = glue(coord, coord_n_g, tol = 10 * t)
    coord = glue(coord, coord_n_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    cellmask!(grid, [0.0 * μm], [heightLayers[1]], regionAcceptor, tol = 1.0e-18)     # n-doped region   = 1
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [heightLayers[2]], [heightLayers[3]], regionDonor, tol = 1.0e-18)  # p-doped region   = 3

    bfacemask!(grid, [0.0], [0.0], bregionAcceptor, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionDonor, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2, tol = 1.0e-18) # second inner interface

    # Plot node grid
    # if plotting
    #     gridplot(grid, Plotter = Plotter, legend = :lt)
    #     Plotter.title("Grid")
    #     #show()
    # end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    data = Data(grid, numberOfCarriers)
    data.modelType = Transient

    carrier_stats = Boltzmann
    data.F = [carrier_stats, carrier_stats, FermiDiracMinusOne]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    data.generationModel = GenerationUniform

    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionJ2] = InterfaceRecombination
    data.boundaryType[bregionDonor] = OhmicContact

    ## The simulation in the paper above uses SCAPS, so to match the results, I've turned off the ionic carrier for now
    # enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    data.fluxApproximation .= ExcessChemicalPotential

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params = Params(grid, numberOfCarriers)

    params.temperature = T
    params.UT = (kB * params.temperature) / q
    params.chargeNumbers[iphin] = zn
    params.chargeNumbers[iphip] = zp
    # params.chargeNumbers[iphia] = za

    for ireg in 1:numberOfRegions

        params.dielectricConstant[ireg] = ε[ireg] * ε0

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg] = Nn[ireg]
        params.densityOfStates[iphip, ireg] = Np[ireg]
        # params.densityOfStates[iphia, ireg] = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg] = En[ireg]
        params.bandEdgeEnergy[iphip, ireg] = Ep[ireg]
        # params.bandEdgeEnergy[iphia, ireg] = Ea[ireg]

        params.mobility[iphin, ireg] = μn[ireg]
        params.mobility[iphip, ireg] = μp[ireg]
        # params.mobility[iphia, ireg] = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg] = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg] = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg] = τp[ireg]

        if ireg == regionAcceptor || ireg == regionDonor
            params.recombinationSRHTrapDensity[iphin, ireg] = 0.0 / (m^3)
            params.recombinationSRHTrapDensity[iphip, ireg] = 0.0 / (m^3)
        else
            params.recombinationSRHTrapDensity[iphin, ireg] = 2.0e21 / (m^3)
            params.recombinationSRHTrapDensity[iphip, ireg] = 2.0e21 / (m^3)
        end
        
        # Beer-Lambert generation for now
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

    min_SRV = 1.0e3 * cm/s
    maj_SRV = 1.0e7 * cm/s

    ## Trap densities and SRH velocities from the following paper: https://doi.org/10.1002/solr.202100219
    params.recombinationSRHvelocity[iphin, bregionJ1] = min_SRV
    params.recombinationSRHvelocity[iphip, bregionJ1] = maj_SRV

    params.bRecombinationSRHTrapDensity[iphin, bregionJ1] = 1.0e17 / (m^3)
    params.bRecombinationSRHTrapDensity[iphip, bregionJ1] = 1.0e17 / (m^3)

    params.recombinationSRHvelocity[iphin, bregionJ2] = maj_SRV * 5.0
    params.recombinationSRHvelocity[iphin, bregionJ2] = min_SRV * 5.0

    params.bRecombinationSRHTrapDensity[iphin, bregionJ2] = 2.0e18 / (m^3)
    params.bRecombinationSRHTrapDensity[iphip, bregionJ2] = 2.0e18 / (m^3)

    ##############################################################

    ## interior doping
    params.doping[iphin, regionDonor] = Cn
    params.doping[iphip, regionAcceptor] = Cp

    # intrinsic carrier concentration in absorber
    params.doping[iphip, regionIntrinsic] = Ci_eff
    # params.doping[iphia, regionIntrinsic] = Ca

    data.params = params
    ctsys = System(grid, data, unknown_storage = :sparse)
    
    # set boundary conditions for potential corresponding to work function of the contacts
    ctsys.fvmsys.boundary_factors[data.index_psi, bregionDonor] = 1.0e30
    ctsys.fvmsys.boundary_values[data.index_psi, bregionDonor] = -4.1

    ctsys.fvmsys.boundary_factors[data.index_psi, bregionAcceptor] = 1.0e30
    ctsys.fvmsys.boundary_values[data.index_psi, bregionAcceptor] = -5.4

    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define control parameters for Solver")
    end
    ################################################################################

    # Solver Control parameters
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

    # COMPARE: Find carrier densities at interfaces and compare to SCAPS results
    println("\nDARK EQUILIBRIUM\n")
    carrier_densities(solution, ctsys, regionIntrinsic, iphin, iphip)

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        ## add labels for anion vacancy
        # label_density[iphia] = "\$ n_a \$";
    end

    ############ DARK EQUILIBRIUM ############

    ipsi = ctsys.fvmsys.physics.data.index_psi
    Vbi = solution[ipsi, end] - solution[ipsi, 1]
    println("Built-in Voltage: $(Vbi)V")

    # Plot equilibrium voltage profile
    # if plotting
    #     Plotter.figure()
    #     plot_solution(Plotter, ctsys, solution, "Equilibrium, Dark", label_solution)

    #     Plotter.show()
    # end

    solution_dark = solution

    # reset Neumann boundary conditions after initial solution is found
    ctsys.fvmsys.boundary_factors[data.index_psi, bregionDonor] = 0
    ctsys.fvmsys.boundary_values[data.index_psi, bregionDonor] = 0

    ctsys.fvmsys.boundary_factors[data.index_psi, bregionAcceptor] = 0
    ctsys.fvmsys.boundary_values[data.index_psi, bregionAcceptor] = 0

    ########## ILLUMINATED EQUILIBRIUM ##########

    I = collect(20:-0.5:0.0)
    LAMBDA = 10 .^ (-I)

    for istep in 1:(length(I) - 1)

        ## ramp up generation
        ctsys.data.λ2 = LAMBDA[istep + 1]

        if test == false
            println("increase generation with λ2 = $(data.λ2)")
        end

        solution = solve(ctsys, inival = inival, control = control)
        inival = solution

    end # generation loop

    # COMPARE: Illuminated equilibrium carrier densities at interfaces
    println("\nILLUMINATED EQUILIBRIUM")
    carrier_densities(solution, ctsys, regionIntrinsic, iphin, iphip)

    # DEBUG: Find peak and minimum of Beer-Lambert generation in the absorber
    # subg = subgrid(grid, [2])
    # println("\nBEER-LAMBERT GENERATION PROFILE INFO")
    # println("PEAK: $(BeerLambert(ctsys, 2, subg[Coordinates])[1])")
    # println("MIN: $(BeerLambert(ctsys, 2, subg[Coordinates])[end])")

    # Plot dark and illuminated short-circuit band energies and carrier densities
    if plotting
        PyPlot.rc("figure", figsize=(10, 10))
        fig = Plotter.figure()

        # 2×2 grid, slot 1
        Plotter.subplot(2, 2, 1)
        plot_energies(Plotter, ctsys, solution_dark, "Equilibrium, Dark", label_energy, clear = false)

        # slot 2
        Plotter.subplot(2, 2, 2)
        plot_densities(Plotter, ctsys, solution_dark, "Equilibrium, Dark", label_density, clear = false)

        # slot 3
        Plotter.subplot(2, 2, 3)
        plot_energies(Plotter, ctsys, solution, "Equilibrium, Illuminated", label_energy, clear = false)

        # slot 4
        Plotter.subplot(2, 2, 4)
        plot_densities(Plotter, ctsys, solution, "Equilibrium, Illuminated", label_density, clear = false)

        # clean up
        fig.tight_layout()
        Plotter.show()

        PyPlot.rc("figure", figsize=(8, 8))
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    # short-circuit current
    Isc = -get_current_val(ctsys, solution)

    # variables for I-V data
    IV = zeros(0)
    biasValues = zeros(0)
    tested = false

    for istep in 2:ntsteps
        t = tvalues[istep]       
        Δu = t * scanrate         
        Δt = t - tvalues[istep - 1] 

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
        
        # DEBUG: Near Voc, determine the quasi-Fermi level splitting (QFLS) and carrier densities at interfaces
        if isapprox(current, 0, atol = 10) && current > 0 && tested == false && plotting
            println("CURRENT: $(current)")
            # determine QFLS near Voc
            subg = subgrid(grid, [regionIntrinsic])
            EFn = view(solution[iphin, :], subg)
            EFp = view(solution[iphip, :], subg)

            QFLS = abs(minimum(EFn) - maximum(EFp))
            println("QFLS: $(QFLS)")

            # carrier densities at absorber edges 
            carrier_densities(solution, ctsys, regionIntrinsic, iphin, iphip)

            Plotter.figure()
            plot_energies(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Δu)", label_energy)
            Plotter.figure()
            plot_densities(Plotter, ctsys, solution, "bias \$\\Delta u\$ = $(Δu)", label_density)        
            show()
            tested = true
        end

    end # time loop

    IV = -IV

    plot_IV(Plotter, biasValues, IV, "bias \$\\Delta u\$ = $(vend)")
    show()

    IV = -IV

    powerDensity = biasValues .* (IV)   
    MaxPD, indexPD = findmax(powerDensity)

    Voc = compute_open_circuit_voltage(biasValues, IV)

    IncidentLightPowerDensity = 1000.0 * W / m^2

    efficiency = biasValues[indexPD] * IV[indexPD] / IncidentLightPowerDensity
    fillfactor = (biasValues[indexPD] * IV[indexPD]) / (IV[1] * Voc)

    println("\nSOLAR CELL PARAMETERS")
    println("Isc = $(Isc))")
    println("Voc = $(Voc)")
    println("FF = $(fillfactor)")
    println("PCE = $(efficiency)")

    if test == false
        println("*** done\n")
    end

end # main
