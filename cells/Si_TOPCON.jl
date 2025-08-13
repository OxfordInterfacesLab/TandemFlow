ENV["MPLBACKEND"] = "qt5agg"

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
PyPlot.rc("axes", titlesize=16)       # Title of plots
PyPlot.rc("axes", labelsize=14)       # Axis labels
PyPlot.rc("xtick", labelsize=12)      # X-tick labels
PyPlot.rc("ytick", labelsize=12)      # Y-tick labels
PyPlot.rc("legend", fontsize=12)      # Legend text

function middle_element(A)
    idx = (length(A) + 1) ÷ 2
    return A[idx]
end

#=
tanh-like drop from `ymax` to `ymin` between x0..x1, then stays at `ymin`.

Arguments:
- x0: start of transition
- x1: end of transition
- ymin, ymax: target bounds
- sharpness: larger => steeper transition
=#
function tanh_plateau(x; x0=0.0, x1=10.0, ymin=0.0, ymax=1.0, sharpness=8.0)
    x1 <= x0 && throw(ArgumentError("require x1 > x0"))
    if x <= x0
        return ymax
    elseif x >= x1
        return ymin
    else
        k   = sharpness / (x1 - x0)
        mid = (x0 + x1)/2
        # tanh gives ~+1 near x0 and ~-1 near x1
        t = (tanh(k*(mid - x)) + 1) / 2   # maps [-1,1] -> [0,1]
        return ymin + (ymax - ymin) * t
    end
end

# save profile as a CSV file
function save_device_profile_csv(filename, solution, ctsys)
    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    x = zeros(0)
    n = zeros(0)
    p = zeros(0)
    Ec = zeros(0)
    Ev = zeros(0)
    EFn = zeros(0)
    EFp = zeros(0)

    for ireg in 1:numberOfRegions
        subg = subgrid(grid, [ireg])

        Ec0 = get_BEE(iphip, ireg, ctsys)
        Ev0 = get_BEE(iphin, ireg, ctsys)
        solpsi = view(solution[data.index_psi, :], subg)
        solp = view(solution[iphip, :], subg)
        soln = view(solution[iphin, :], subg)

        append!(x, subg[Coordinates]')
        append!(n, get_density(solution, ireg, ctsys, iphin))
        append!(p, get_density(solution, ireg, ctsys, iphip))
        append!(Ec, Ec0 ./ q .- solpsi)
        append!(Ev, Ev0 ./ q .- solpsi)
        append!(EFn, -soln)
        append!(EFp, -solp)
    end

    # Build DataFrame
    df = DataFrame(
        x = x,
        n = n,
        p = p,
        Ec = Ec,
        Ev = Ev,
        EFn = EFn,
        EFp = EFp
    )

    CSV.write(filename, df)
end

# you can also use other Plotters, if you add them to the example file
function main(;
        n = 6, Plotter = PyPlot, plotting = true,
        verbose = false, test = false,
        #parameter_file = "../parameter_files/Params_PSC_TiO2_MAPI_spiro.jl", # choose the parameter file
        parameter_file = "../params/Params_Si_TOPCON.jl", # choose the parameter file
    )

    if plotting
        Plotter.close("all")
    end

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    include(parameter_file) # include the parameter file we specified

    ## contact voltage
    voltageAcceptor = 0.75 * V

    ## primary data for I-V scan protocol
    scanrate = 0.3 * V / s
    ntsteps = 91
    vend = voltageAcceptor # bias goes until the given voltage at acceptor boundary
    tend = vend / scanrate

    ## with fixed timestep sizes we can calculate the times a priori
    tvalues = range(0, stop = tend, length = ntsteps)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # TODO: Need density increase around interface

    δ = 6 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_em = collect(range(0.0, 3*h_emitter, step = h_emitter / (0.8 * δ)))
    coord_cz_u = collect(range(3*h_emitter, 0.95 * h_cz, step = h_cz / (0.8 * δ)))
    coord_cz_g = geomspace(
        0.95 * h_cz, h_cz,
        (0.05 * h_cz) / (0.7 * δ), (0.05 * h_cz) / (1.1 * δ),
        tol = t
    )
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
    cellmask!(grid, [0.0], [heightLayers[1]], regionCz, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [heightLayers[1]], [heightLayers[2]], regionPoly, tol = 1.0e-18)  # p-doped region   = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], bregionCz, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [h_total], [h_total], bregionPoly, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1, tol = 1.0e-18) # first  inner interface

    ## Plot node grid
    # if plotting
    #     gridplot(grid, Plotter = Plotter, legend = :lt)
    #     Plotter.title("Grid")
    #     Plotter.show()
    # end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## set up generation data
    subg1 = subgrid(grid, [regionCz]); subg2 = subgrid(grid, [regionPoly]);

    gen1 = zeros(length(subg1[Coordinates])); gen2 = zeros(length(subg2[Coordinates]) - 1)

    for i in 1:length(subg1[Coordinates])
        gen1[i] = generation_from_scaps("simulation_data/scaps/si-topcon-auto.gen", subg1[Coordinates][i])
    end

    generationData = [gen1; gen2]

    ## Initialize Data instance and fill in data
    data = Data(grid, numberOfCarriers, generationData = generationData)

    data.modelType = Transient
    carrier_stats = Boltzmann
    data.F = [carrier_stats, carrier_stats]

    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    data.generationModel = GenerationUserDefined

    data.boundaryType[bregionPoly] = OhmicContact
    data.boundaryType[bregionJ1] = InterfaceRecombination
    data.boundaryType[bregionCz] = OhmicContact

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
    end

    ##############################################################
    ## inner boundary region data (we choose the intrinsic values)
    params.bDensityOfStates[iphin, bregionJ1] = Nn[regionPoly]
    params.bDensityOfStates[iphip, bregionJ1] = Np[regionPoly]

    params.bBandEdgeEnergy[iphin, bregionJ1] = En[regionPoly]
    params.bBandEdgeEnergy[iphip, bregionJ1] = Ep[regionPoly]

    # Schottky Implementation (didn't fix issues)

    # params.SchottkyBarrier[bregionCz] = 1.10 * (eV)
    # params.SchottkyBarrier[bregionPoly] = 0.0 * (eV)

    # params.bVelocity[iphin, bregionCz] = 1.0e-2 * cm/s
    # params.bVelocity[iphip, bregionCz] = 1.0e7 * cm/s

    # params.bVelocity[iphin, bregionPoly] = 1.0e7 * cm/s
    # params.bVelocity[iphip, bregionPoly] = 1.0e-2 * cm/s

    ## Positive doping corresponds to acceptors
    for icoord in 1:numberOfNodes
        if icoord <= (length(coord_em) + length(coord_cz_u) + length(coord_cz_g) - 2) # n C-Si region
            paramsnodal.doping[icoord] = tanh_plateau(coord[icoord]; 
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

    ctsys = System(grid, data, unknown_storage = :sparse)

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
    control.damp_initial = 0.5
    control.damp_growth = 1.51 # >= 1

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

    save_device_profile_csv("simulation_data/chargetransport/si-topcon-dark-sc.csv", solution, ctsys)

    if plotting == true
        ################################################################################
        println("Plot electroneutral potential, band-edge energies and doping")
        ################################################################################
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)

        # Plot dark equilibrium
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Equilibrium, Dark", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Equilibrium, Dark", label_density)
        Plotter.show()
    end

    ## DARK EQUILIBRIUM
    ipsi = ctsys.fvmsys.physics.data.index_psi
    Vbi = solution[ipsi, end] - solution[ipsi, 1]
    println("Built-in Voltage: $(Vbi)V")

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

    # save_device_profile_csv("si_1_ill_sc.csv", solution, ctsys)
    # Plot illuminated short-circuit
    if plotting
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_energy)
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Illuminated Short-Circuit", label_density)
        Plotter.show()
    end

    save_device_profile_csv("simulation_data/chargetransport/si-topcon-illuminated-sc.csv", solution, ctsys)
    
    println("DEBUG/BL-MIN: $(data.generationData[length(gen1)])")
    Plotter.figure()
    Plotter.plot(coord, data.generationData)
    Plotter.show()

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("I-V Measurement Loop")
    end
    ################################################################################

    Isc = get_current_val(ctsys, solution)

    ## for saving I-V data
    IV = zeros(0) # for IV values
    biasValues = zeros(0) # for bias values
    tested = false

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

        push!(IV, current)
        push!(biasValues, Δu)


    end # time loop

    # Plot IV curve
    # plot_IV(Plotter, biasValues, -IV, "bias \$\\Delta u\$ = $(vend)")
    # show()

    IV = IV

    powerDensity = biasValues .* (IV)           # power density function
    MaxPD, indexPD = findmax(powerDensity)

    Voc = compute_open_circuit_voltage(biasValues, IV)

    IncidentLightPowerDensity = 1000.0 * W / m^2

    efficiency = biasValues[indexPD] * IV[indexPD] / IncidentLightPowerDensity
    fillfactor = (biasValues[indexPD] * IV[indexPD]) / (IV[1] * Voc)

    println("\nIsc = $(Isc)")
    println("Voc = $(Voc)")
    println("FF = $(fillfactor)")
    println("PCE = $(efficiency)")
end

# DEBUG
main(n = 12, verbose = true)