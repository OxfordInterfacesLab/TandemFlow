#=

Perovskite solar cell simulation // Graded-dielectric ETL

ETL: C60 
Perovskite: FaCsPbIBr
HTL: NiO

=#

module PSC_C60_PVK_NiO_graded

using ChargeTransport
using ExtendableGrids
using PyPlot
using Interpolations

## import utilities
include("../utils/ct_utils.jl")
using .CTUtils

## modify this dictionary to choose what to plot
toPlot = Dict(
    "grid" => false,
    "generation" => false,
    "dark-sc" => true,
    "light-sc" => false,
    "light-bias" => true,
    "light-oc" => false,
    "iv" => false
)

parameter_file = "../params/Params_PSC_C60_PVK_NiO_graded.jl"
include(parameter_file) # include the parameter file we specified
p = Params_PSC_C60_PVK_NiO_graded() # create an instance of the parameters

function main(;
        n = 6, Plotter = PyPlot, plotting = false,
        verbose = false
    )

    @local_unitfactors μm cm s ns V K ps Hz W m eV

    println("--- Define physical parameters and model ---")

    ## max contact voltage during I-V scan
    maxVoltage = 1.0 * V

    ## primary data for I-V scan protocol
    scanrate = 0.03 * V / s
    ntsteps = 201 # number of time steps
    tend = maxVoltage / scanrate

    tvalues = range(0, stop = tend, length = ntsteps)

    println("--- Set up grid and regions ---")

    δ = 4 * n        # the larger, the finer the mesh
    t = 0.5 * (cm) / δ # tolerance for geomspace and glue (with factor 10)
    k = 1.5        # the closer to 1, the closer to the boundary geomspace

    coord_etl1_u = collect(range(0.0, p.h_etl1 / 2, step = p.h_etl1 / (0.8 * δ))) # uniform ETL1
    coord_etl1_g = geomspace(
        p.h_etl1 / 2, p.heightLayers[1],
        p.h_etl1 / (0.7 * δ), p.h_etl1 / (1.1 * δ),
        tol = t
    ) # geometric ETL1
    coord_etl2_g1 = geomspace(
        p.heightLayers[1], p.heightLayers[1] + p.h_etl2 / 2,
        p.h_etl2 / (1.1 * δ), p.h_etl2 / (0.7 * δ),
        tol = t
    ) # geometric ETL2/1
    coord_etl2_g2 = geomspace(
        p.heightLayers[1] + p.h_etl2 / 2, p.heightLayers[2],
        p.h_etl2 / (0.7 * δ), p.h_etl2 / (1.1 * δ),
        tol = t
    ) # geometric ETL2/2
    coord_etl3_g1 = geomspace(
        p.heightLayers[2], p.heightLayers[2] + p.h_etl3 / 2,
        p.h_etl3 / (1.1 * δ), p.h_etl3 / (0.7 * δ),
        tol = t
    ) # geometric ETL3/1
    coord_etl3_g2 = geomspace(
        p.heightLayers[2] + p.h_etl3 / 2, p.heightLayers[3],
        p.h_etl3 / (0.7 * δ), p.h_etl3 / (1.1 * δ),
        tol = t
    ) # geometric ETL3/2
    coord_pvk_g1 = geomspace(
        p.heightLayers[3], p.heightLayers[3] + p.h_pvk / k,
        p.h_pvk / (5.1 * δ), p.h_pvk / (1.1 * δ),
        tol = t
    ) # geometric Pvk/1
    coord_pvk_g2 = geomspace(
        p.heightLayers[3] + p.h_pvk / k, p.heightLayers[4],
        p.h_pvk / (1.1 * δ), p.h_pvk / (5.1 * δ),
        tol = t
    ) # geometric Pvk/2
    coord_htl_g = geomspace(
        p.heightLayers[4], p.heightLayers[4] + p.h_htl / 2,
        p.h_htl / (1.3 * δ), p.h_htl / (0.6 * δ),
        tol = t
    ) # geometric HTL
    coord_htl_u = collect(range(p.heightLayers[4] + p.h_htl / 2, p.heightLayers[5], step = p.h_htl / (0.8 * δ))) # uniform HTL

    coord = glue(coord_etl1_u, coord_etl1_g, tol = 10 * t)
    coord = glue(coord, coord_etl2_g1, tol = 10 * t)
    coord = glue(coord, coord_etl2_g2, tol = 10 * t)
    coord = glue(coord, coord_etl3_g1, tol = 10 * t)
    coord = glue(coord, coord_etl3_g2, tol = 10 * t)
    coord = glue(coord, coord_pvk_g1, tol = 10 * t)
    coord = glue(coord, coord_pvk_g2, tol = 10 * t)
    coord = glue(coord, coord_htl_g, tol = 10 * t)
    coord = glue(coord, coord_htl_u, tol = 10 * t)
    grid = ExtendableGrids.simplexgrid(coord)

    ## set different regions in grid
    cellmask!(grid, [0.0], [p.heightLayers[1]], p.regionETL1, tol = 1.0e-18)     # ETL1
    cellmask!(grid, [p.heightLayers[1]], [p.heightLayers[2]], p.regionETL2, tol = 1.0e-18)     # ETL2
    cellmask!(grid, [p.heightLayers[2]], [p.heightLayers[3]], p.regionETL3, tol = 1.0e-18) # ETL3
    cellmask!(grid, [p.heightLayers[3]], [p.heightLayers[4]], p.regionPerovskite, tol = 1.0e-18)  # Pvk
    cellmask!(grid, [p.heightLayers[4]], [p.heightLayers[5]], p.regionHTL, tol = 1.0e-18)  # HTL

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0], [0.0], p.bregionETL, tol = 1.0e-18)     # outer left boundary
    bfacemask!(grid, [p.h_total], [p.h_total], p.bregionHTL, tol = 1.0e-18)  # outer right boundary
    bfacemask!(grid, [p.heightLayers[1]], [p.heightLayers[1]], p.bregionJ1, tol = 1.0e-18) # first  inner interface
    bfacemask!(grid, [p.heightLayers[2]], [p.heightLayers[2]], p.bregionJ2, tol = 1.0e-18) # second inner interface
    bfacemask!(grid, [p.heightLayers[3]], [p.heightLayers[3]], p.bregionJ3, tol = 1.0e-18) # third  inner interface
    bfacemask!(grid, [p.heightLayers[4]], [p.heightLayers[4]], p.bregionJ4, tol = 1.0e-18) # fourth inner interface

    ## plot node grid
    if plotting && toPlot["grid"]
        gridplot(grid, Plotter = Plotter, legend = :lt)
        Plotter.title("Grid")
        Plotter.show()
    end

    println("--- Define system and fill in information about model ---")

    ## initialize Data instance and fill in data
    data = Data(grid, p.numberOfCarriers)

    ## choose simulation type
    ## possible choices: Stationary, Transient
    data.modelType = Transient

    ## choose statistics model
    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore
    data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

    ## set relevant recombination mechanisms
    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## set interface types
    ## possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[p.bregionETL] = OhmicContact
    data.boundaryType[p.bregionHTL] = OhmicContact

    ## present ionic vacancies in perovskite layer
    enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionPerovskite])

    ## set flux discretization scheme
    ## possible choices: DiffusionEnhanced, DiffusionEnhancedModifiedDrift, ExcessChemicalPotential,
    ## ExcessChemicalPotentialGraded, GeneralizedSG, ScharfetterGummel, ScharfetterGummelGraded
    data.fluxApproximation .= DiffusionEnhanced

    println("--- Define Params ---")
    
    data.params = Params(p)

    ctsys = System(grid, data, unknown_storage = :sparse)

    println("--- Define control parameters for solver ---")

    ## solver parameters - tweak when convergence is an issue
    control = SolverControl()
    control.verbose = verbose
    control.damp_initial = 0.2
    control.damp_growth = 1.11 # >= 1
    control.maxiters = 1000

    println("--- Solve in equilibrium ---")

    ## determine volume of perovskite layer
    subg = subgrid(grid, [p.regionPerovskite]) # perovskite subgrid

    mOmega = 0.0 # volume of perovskite layer
    for icellVol in subg[CellVolumes]
        mOmega = mOmega + icellVol
    end

    ## bisection search loop algorithm to determine appropriate Ea value
    vacancyDensityRelTol = 1e-2 # relative tolerance for vacancy density
    vacancyDensityMaxIters = 100 # maximum number of iterations for bisection algorithm

    iters = 0
    vacancyDensity = 0.0
    Ea_range = [-4.0, -6.0]
    Ea_mid = (Ea_range[1] + Ea_range[2]) / 2

    solution = nothing
    inival = nothing

    while iters < vacancyDensityMaxIters && (iters == 0 || abs(vacancyDensity - p.dop_vac[p.regionPerovskite]) / vacancyDensity > vacancyDensityRelTol)
        if iters > 0
            if vacancyDensity > p.dop_vac[p.regionPerovskite]
                Ea_range[1] = Ea_mid
            else
                Ea_range[2] = Ea_mid
            end
        end

        Ea_mid = (Ea_range[1] + Ea_range[2]) / 2
        data.params.bandEdgeEnergy[p.iphia, p.regionPerovskite] = Ea_mid * eV

        solution = equilibrium_solve!(ctsys, control = control)
        inival = solution

        intncc = ChargeTransport.integrate(ctsys, storage!, solution)./p.constants.q
        int = intncc[p.iphia, p.regionPerovskite]/mOmega
        vacancyDensity = int
        iters += 1

        println("ITERATION: $(iters) ---- CURRENT RELATIVE DIFFERENCE: $(abs(vacancyDensity - p.dop_vac[p.regionPerovskite]) / vacancyDensity) ---- CURRENT Ea: $(Ea_mid) eV")
    end

    ## set axis labels for plots if plotting is set to 'on'
    if plotting
        label_solution, label_density, label_energy, label_BEE = set_plotting_labels(data)
        label_energy[1, p.iphia] = "\$E_a-q\\psi\$"; label_energy[2, p.iphia] = "\$ - q \\varphi_a\$"; label_BEE[p.iphia] = "\$E_a\$"
        label_density[p.iphia] = "\$ n_a \$";      label_solution[p.iphia] = "\$ \\varphi_a\$"
    end

    ## plot carrier densities and energies in short-circuit in the dark
    if plotting && toPlot["dark-sc"]
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "Dark Short-Circuit", label_energy)
        # savefig("psc-dark-sc-bands.png")
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "Dark Short-Circuit", label_density)
        # savefig("psc-dark-sc-densities.png")
        Plotter.figure()
        plot_Efield(Plotter, ctsys, solution, "Dark Short-Circuit")
        # savefig("psc-dark-sc-Efield.png")
        Plotter.clf()
        # Plotter.show()
    end

    save_cell_profile("sims/PSC_C60_PVK_NiO_graded_sc.csv", solution, ctsys, true, p.iphia)

    println("--- IV Curve ---")

    ## for saving I-V data
    currents = zeros(0) # for current densities (A m^{-2})
    biasValues = zeros(0) # for bias values (V)

    for istep in 2:ntsteps
        t = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep - 1] # Time step size

        ## Apply new voltage (set non-equilibrium values)
        set_contact!(ctsys, p.bregionETL, Δu = Δu)

        solution = solve(ctsys, inival = inival, control = control, tstep = Δt)
        inival = solution

        if verbose
            println("Time: $(round(t, digits=4))s, Voltage: $(round(Δu, digits=4))V")
        end

        ## get I-V data
        current = get_current_val(ctsys, solution, inival, Δt)

        push!(currents, current)
        push!(biasValues, Δu)
    end 

    ## flip currents - use this based on cell architecture
    currents = -currents

    if plotting && toPlot["light-bias"]
        Plotter.figure()
        plot_energies(Plotter, ctsys, solution, "1V Reverse Bias", label_energy)
        # savefig("psc-1Vrb-bands.png")
        Plotter.figure()
        plot_densities(Plotter, ctsys, solution, "1V Reverse Bias", label_density)
        # savefig("psc-1Vrb-densities.png")
        Plotter.figure()
        plot_Efield(Plotter, ctsys, solution, "1V Reverse Bias")
        # savefig("psc-1Vrb-Efield.png")
        Plotter.clf()
        # Plotter.show()
    end

    save_cell_profile("sims/PSC_C60_PVK_NiO_graded_2Vrb.csv", solution, ctsys, true, p.iphia)

    ## plot IV curve
    if plotting && toPlot["iv"]
        plot_IV(Plotter, biasValues, -currents, "bias \$\\Delta u\$ = $(maxVoltage)")
        show()
    end

    ## return IV data in an IV struct
    return IV(biasValues, currents)
end

end