using DataFrames
using CSV
using PyPlot

struct CellProfile
    x::Vector{Float64}      # position in μm
    n::Vector{Float64}      # electron density in /cm3
    p::Vector{Float64}      # hole density in /cm3
    Ec::Vector{Float64}     # conduction band edge in eV
    Ev::Vector{Float64}     # valence band edge in eV
    Efn::Vector{Float64}    # electron quasi-Fermi level in eV
    Efp::Vector{Float64}    # hole quasi-Fermi level in eV
    IV::Matrix{Float64}     # IV curve data

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

struct GenerationProfile
    x::Vector{Float64}      # position in μm
    G::Vector{Float64}    # generation rate in /cm3.s
end

#######################################
########## FILE CONVERSIONS ###########
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
            df[!, "x (um)"],
            df[!, "Geh (#/cm3.s)"]
        )

        return profile
    elseif endswith(filename, ".iv") # IV curves
        start_idx = findfirst(x -> occursin("IV/START", x), lines)
        end_idx = findfirst(x -> occursin("IV/END", x), lines)

        if isnothing(start_idx) || isnothing(end_idx) || end_idx <= start_idx
            error("Could not find IV/START and IV/END or invalid order.")
        end
        
        df = table_to_dataframe(lines[(start_idx+1):(end_idx-1)])
    else
        throw(ArgumentError("scaps_to_df() doesn't support this file format: $filename"))
    end

    return profile
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
    plot(profile_ct.IV[1], profile_ct.IV[2], label="CT", color="blue")
    plot(profile_scaps.IV[1], profile_scaps.IV[2], label="SCAPS", color="red", linestyle="--")

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

##################################
########### GENERATION ###########

function generation_from_scaps(filename::String)
    profile = parse_scaps(filename)

    return function(x::Float64)
        x = x / (μm)
        if x <= profile.x[1]
            return profile.G[1] * 1e6
        elseif x >= profile.x[end]
            return profile.G[end] * 1e6
        else
            i = findfirst(i -> profile.x[i] <= x < profile.x[i+1], 1:length(profile.x)-1)
            if isnothing(i)
                error("x is out of bounds or data is not sorted.")
            end
            x1, x2 = profile.x[i], profile.x[i+1]
            G1, G2 = profile.G[i], profile.G[i+1]
            return G1 * exp(log(G2/G1) * (x - x1) / (x2 - x1)) * 1e6
        end
    end
end