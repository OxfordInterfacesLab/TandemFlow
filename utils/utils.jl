using ChargeTransport
using ExtendableGrids
using DataFrames
using CSV
using PyPlot

# ----------------------------------------------------
#                    DEFINE STRUCTS
# ----------------------------------------------------

"""
A struct for an IV curve
"""
struct IV
    biasValues::Vector{Float64}  # Voltage in V
    current_density::Vector{Float64}  # Current density in mA/cm2
end

"""
A mutable struct for a cell profile
"""
mutable struct CellProfile
    x::Vector{Float64}      # position in μm
    n::Vector{Float64}      # electron density in /cm3
    p::Vector{Float64}      # hole density in /cm3
    Ec::Vector{Float64}     # conduction band edge in eV
    Ev::Vector{Float64}     # valence band edge in eV
    Efn::Vector{Float64}    # electron quasi-Fermi level in eV
    Efp::Vector{Float64}    # hole quasi-Fermi level in eV
    IV::Matrix{Float64}     # IV curve data

    # TODO: Write a better constructor
    has_profile::Bool
    has_IV::Bool

    # constructor for if the cell profile is not given an IV measurement
    function CellProfile(x::Vector{Float64}, n::Vector{Float64}, p::Vector{Float64}, Ec::Vector{Float64}, Ev::Vector{Float64}, Efn::Vector{Float64}, Efp::Vector{Float64})
        IV = zeros(2, 1)

        # translate the bands so that the left edge of electon quasi-Fermi level is at 0 eV
        translationDistance = 0 - Efn[1]
        Ec = Ec .+ translationDistance
        Ev = Ev .+ translationDistance
        Efn = Efn .+ translationDistance
        Efp = Efp .+ translationDistance

        new(x, n, p, Ec, Ev, Efn, Efp, IV)
    end

    # constructor for if the cell profile is only given an IV measurement
    function CellProfile(IV::Matrix{Float64})
        new(zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), IV)
    end
end

"""
A struct for generation profiles
"""
struct GenerationProfile
    x::Vector{Float64}      # position in m
    G::Vector{Float64}    # generation rate in m^-3 s^-1
end 

# ----------------------------------------------------
#                    SIM UTILITIES
# ----------------------------------------------------

"""
tanh function for error function-esque diffusion profile

Profile drops from ymax to ymin between x0 and x1
"""
function tanh_diffusion_profile(x; x0=0.0, x1=10.0, ymin=0.0, ymax=1.0, sharpness=8.0)
    x1 <= x0 && throw(ArgumentError("require x1 > x0"))

    if x <= x0
        return ymax
    elseif x >= x1
        return ymin
    else
        k   = sharpness / (x1 - x0)
        mid = (x0 + x1)/2
        t = (tanh(k*(mid - x)) + 1) / 2
        return ymin + (ymax - ymin) * t
    end

end

"""
Save the cell profile to a CSV file.
"""

#TODO: Add ion functionality
function save_cell_profile(filename, solution, ctsys, ions=false)
    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data
    numberOfRegions = grid[NumCellRegions]

    iphin = data.bulkRecombination.iphin
    iphip = data.bulkRecombination.iphip

    x = zeros(0)
    n = zeros(0)
    p = zeros(0)
    Ec = zeros(0)
    Ev = zeros(0)
    EFn = zeros(0)
    EFp = zeros(0)

    for ireg in 1:numberOfRegions
        subg = subgrid(grid, [ireg])

        Ec0 = get_BEE(iphin, ireg, ctsys) # short-circuit conduction band edge
        Ev0 = get_BEE(iphip, ireg, ctsys) # short-circuit valence band edge
        solpsi = view(solution[data.index_psi, :], subg) # electrical potential
        solp = view(solution[iphip, :], subg) # quasi-Fermi potential for holes
        soln = view(solution[iphin, :], subg) # quasi-Fermi potential for electrons

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

function save_iv(filename::String, biasValues, IV)
    df = DataFrame(
        V = biasValues,
        J = IV
    )

    CSV.write(filename, df)
end

function get_cell_characteristics(IV::IV)
    biasValues = IV.biasValues
    currents = IV.current_density

    powerDensity = biasValues .* (currents)
    MaxPD, indexPD = findmax(powerDensity)

    Voc = compute_open_circuit_voltage(biasValues, currents)

    fillfactor = (biasValues[indexPD] * currents[indexPD]) / (currents[1] * Voc)

    characteristics = Dict(
        "Voc" => Voc,
        "Pmax" => MaxPD,
        "FF" => fillfactor,
    )

    return characteristics
end

# ----------------------------------------------------
#                    SCAPS UTILITIES
# ----------------------------------------------------

"""
This function converts the tables in SCAPS files to DataFrames.
DO NOT CALL THIS FUNCTION DIRECTLY
"""
function table_to_dataframe(lines)
    lines = filter(line -> !isempty(line), lines)
    # split header and data
    header = strip.(split(strip(lines[1]), '\t'))
    data = [split(strip(line), '\t') for line in lines[2:end]]

    data = [map(x -> tryparse(Float64, x) === nothing ? x : parse(Float64, x), row) for row in data] # convert numerical data to floats
    df = DataFrame([getindex.(data, i) for i in 1:length(header)], Symbol.(header)) # create dataframe

    return df
end

"""
This function condenses a cell profile to desired x range.
DO NOT CALL THIS FUNCTION DIRECTLY
"""
function condenseIndices(profile::CellProfile, startIndex, endIndex)
    condensedProfile = CellProfile(
        profile.x[startIndex:endIndex],
        profile.n[startIndex:endIndex],
        profile.p[startIndex:endIndex],
        profile.Ec[startIndex:endIndex],
        profile.Ev[startIndex:endIndex],
        profile.Efn[startIndex:endIndex],
        profile.Efp[startIndex:endIndex],
        profile.IV
    )

    return condensedProfile
end


"""
Convert a SCAPS file to a DataFrame.
"""
function parse_scaps(filename::String)
    !isfile(filename) && throw(ArgumentError("File does not exist: $filename")) # check if file exists

    lines = readlines(filename)

    if endswith(filename, ".eb") # band diagrams, densities
        start_idx = findfirst(x -> occursin(r"(?i)bulk", x), lines)
        end_idx = findfirst(x -> occursin(r"(?i)interface", x), lines)

        if isnothing(start_idx) || isnothing(end_idx) || end_idx <= start_idx
            error("Could not find valid table section.")
        end

        df = table_to_dataframe(lines[(start_idx+1):(end_idx-1)])

        profile = CellProfile(
            df[!, "x(um)"],
            df[!, "n(/cm3)"],
            df[!, "p(/cm3)"],
            df[!, "Ec(eV)"],
            df[!, "Ev(eV)"],
            df[!, "Fn(eV)"],
            df[!, "Fp(eV)"]
        )

        return profile
    elseif endswith(filename, ".gen") # generation
        gen_idx = findfirst(x -> occursin("GEN", x), lines)

        if isnothing(gen_idx) || gen_idx == length(lines)
            error("No GEN line or header found in file.")
        end

        df = table_to_dataframe(lines[(gen_idx+1):end])

        profile = GenerationProfile(
            df[!, "x (um)"] .* 1e-6, # convert from μm to m
            df[!, "Geh (#/cm3.s)"] .* 1e6 # convert from cm^-3 s^-1
        )

        return profile
    elseif endswith(filename, ".iv") # IV curves
        start_idx = findfirst(x -> occursin("IV/START", x), lines)
        end_idx = findfirst(x -> occursin("IV/END", x), lines)

        if isnothing(start_idx) || isnothing(end_idx) || end_idx <= start_idx
            error("Could not find IV/START and IV/END or invalid order.")
        end
        
        df = table_to_dataframe(lines[(start_idx+1):(end_idx-1)])

        IV_matrix = [df[!, "v(V)"];; df[!, "jtot(mA/cm2)"]]
        
        return IV_matrix
    else
        throw(ArgumentError("scaps_to_df() doesn't support this file format: $filename"))
    end

    return profile
end

function set_IV!(profile::CellProfile, IV::Matrix{Float64})
    if size(IV, 1) != 2
        throw(ArgumentError("IV matrix must have 2 rows: [voltage; current density]"))
    end

    if size(IV, 2) < 1
        throw(ArgumentError("IV matrix must have at least one column"))
    end

    profile.IV = IV
end

function set_IV!(profile::CellProfile, filename::String)
    if endswith(filename, ".iv")
        IV = parse_scaps(filename)
    elseif endswith(filename, ".csv")
        IV_df = CSV.read(filename, DataFrame)
        IV = [IV_df[!, "V"];; IV_df[!, "J"]]
    end

    profile.IV = IV
end

"""
Convert a ChargeTransport file to a cell profile - not currently set for IV curves
"""
function ct_to_profile(filename::String)
    df = CSV.read(filename, DataFrame)

    df.x .= df.x .* 1e6  # Convert x from m to μm
    df.n .= df.n .* 1e-6  # Convert n from /m3 to /cm3
    df.p .= df.p .* 1e-6  # Convert p from /m3 to /cm3

    profile = CellProfile(
        df[!, "x"],
        df[!, "n"],
        df[!, "p"],
        df[!, "Ec"],
        df[!, "Ev"],
        df[!, "EFn"],
        df[!, "EFp"]
    )

    return profile
end

#######################################
############# COMPARISONS #############

"""
Plots carrier density values from SCAPS and ChargeTransport
"""
function compare_densities(profile_ct::CellProfile, profile_scaps::CellProfile; xrange=nothing)
    # check if user wants to plot densities within a specific spatial range
    if xrange !== nothing
        xmin, xmax = xrange
        xmin > xmax && throw(ArgumentError("xrange must be a tuple of (xmin, xmax) with xmin < xmax"))

        ct_indices = findfirst(i -> profile_ct.x[i] >= xmin, 1:length(profile_ct.x)), findlast(i -> profile_ct.x[i] <= xmax, 1:length(profile_ct.x))
        scaps_indices = findfirst(i -> profile_scaps.x[i] >= xmin, 1:length(profile_scaps.x)), findlast(i -> profile_scaps.x[i] <= xmax, 1:length(profile_scaps.x))

        if nothing in ct_indices || nothing in scaps_indices
            error("xrange is out of bounds for one of the profiles.")
        end

        if ct_indices[2] == ct_indices[1] || scaps_indices[2] == scaps_indices[1]
            error("xrange is too narrow, no data points to plot.")
        end

        profile_ct = condenseIndices(profile_ct, ct_indices[1], ct_indices[2])
        profile_scaps = condenseIndices(profile_scaps, scaps_indices[1], scaps_indices[2])
    end

    # Plotting
    figure(figsize=(10, 6))
    plot(profile_ct.x, profile_ct.n, label="ChargeTransport (n)", color="blue")
    plot(profile_scaps.x, profile_scaps.n, label="SCAPS (n)", color="blue", linestyle="--")
    plot(profile_ct.x, profile_ct.p, label="ChargeTransport (p)", color="red")
    plot(profile_scaps.x, profile_scaps.p, label="SCAPS (p)", color="red", linestyle="--")

    yscale("log")
    xlabel("x (um)")
    ylabel("n (/cm3)")
    title("Comparison of Carrier Densities")
    legend()

    grid(true)
end

function compare_densities(filename_ct::String, filename_scaps::String; xrange=nothing)
    profile_ct = ct_to_profile(filename_ct)
    profile_scaps = parse_scaps(filename_scaps)

    compare_densities(profile_ct, profile_scaps; xrange=xrange)
end

"""
Plots bands values from SCAPS and ChargeTransport
Solid lines are from ChargeTransport, dashed lines are from SCAPS
"""
function compare_bands(profile_ct::CellProfile, profile_scaps::CellProfile; xrange=nothing)
    if xrange !== nothing
        xmin, xmax = xrange
        xmin > xmax && throw(ArgumentError("xrange must be a tuple of (xmin, xmax) with xmin < xmax"))

        ct_indices = findfirst(i -> profile_ct.x[i] >= xmin, 1:length(profile_ct.x)), findlast(i -> profile_ct.x[i] <= xmax, 1:length(profile_ct.x))
        scaps_indices = findfirst(i -> profile_scaps.x[i] >= xmin, 1:length(profile_scaps.x)), findlast(i -> profile_scaps.x[i] <= xmax, 1:length(profile_scaps.x))

        if nothing in ct_indices || nothing in scaps_indices
            error("xrange is out of bounds for one of the profiles.")
        end

        if ct_indices[2] == ct_indices[1] || scaps_indices[2] == scaps_indices[1]
            error("xrange is too narrow, no data points to plot.")
        end

        profile_ct = condenseIndices(profile_ct, ct_indices[1], ct_indices[2])
        profile_scaps = condenseIndices(profile_scaps, scaps_indices[1], scaps_indices[2])
    end

    figure(figsize=(10, 6))
    plot(profile_ct.x, profile_ct.Ec, color="black")
    plot(profile_ct.x, profile_ct.Ev, color="black")
    plot(profile_ct.x, profile_ct.Efn, color="blue")
    plot(profile_ct.x, profile_ct.Efp, color="red")
    plot(profile_scaps.x, profile_scaps.Ec, color="black", linestyle="--")
    plot(profile_scaps.x, profile_scaps.Ev, color="black", linestyle="--")
    plot(profile_scaps.x, profile_scaps.Efn, color="blue", linestyle="--")
    plot(profile_scaps.x, profile_scaps.Efp, color="red", linestyle="--")

    xlabel("x (um)")
    ylabel("Energy (eV)")
    title("Comparison of Energy Bands")
    grid(true)
end

function compare_bands(filename_ct::String, filename_scaps::String; xrange=nothing)
    profile_ct = ct_to_profile(filename_ct)
    profile_scaps = parse_scaps(filename_scaps)

    compare_bands(profile_ct, profile_scaps; xrange=xrange)
end

"""
Plots IV curves from SCAPS and ChargeTransport
"""
function compare_iv(profile_ct::CellProfile, profile_scaps::DataFrame)
    # Plotting IV curves
    figure(figsize=(10, 6))
    plot(profile_ct.IV[:, 1], profile_ct.IV[:, 2], label="CT", color="blue")
    plot(profile_scaps.IV[:, 1], profile_scaps.IV[:, 2], label="SCAPS", color="red", linestyle="--")

    xlabel("Voltage (V)")
    ylabel("Current Density (mA/cm2)")
    title("IV Curve Comparison")
    legend()
    grid(true)
end

function compare_iv(filename_ct::String, filename_scaps::String)
    profile_ct = ct_to_profile(filename_ct)
    profile_scaps = parse_scaps(filename_scaps)

    compare_iv(profile_ct, profile_scaps)
end

"""
Generates a function to get the generation rate from a SCAPS generation profile
"""
function generation_from_scaps(profile::GenerationProfile)
    return function(x::Float64)
        if x <= profile.x[1]
            return profile.G[1]
        elseif x >= profile.x[end]
            return profile.G[end]
        else
            i = findfirst(i -> profile.x[i] <= x < profile.x[i+1], 1:length(profile.x)-1)
            if isnothing(i)
                error("x is out of bound s or data is not sorted.")
            end
            x1, x2 = profile.x[i], profile.x[i+1]
            G1, G2 = profile.G[i], profile.G[i+1]
            return G1 * exp(log(G2/G1) * (x - x1) / (x2 - x1))
        end
    end
end

function generation_from_scaps(filename::String)
    profile = parse_scaps(filename)

    generation_from_scaps(profile)
end