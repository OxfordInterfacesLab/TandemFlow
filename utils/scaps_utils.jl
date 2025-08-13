using DataFrames
using CSV
using PyPlot

struct CellInfo
    x::Float64
    n::Float64
    p::Float64
    Ec::Float64
    Ev::Float64
    Efn::Float64
    Efp::Float64
end

"""
Convert a SCAPS file to a DataFrame.
"""
function scaps_to_df(filename::String)
    lines = readlines(filename)
    start_idx = findfirst(x -> occursin(r"(?i)bulk", x), lines)
    end_idx = findfirst(x -> occursin(r"(?i)interface", x), lines)

    if isnothing(start_idx) || isnothing(end_idx) || end_idx <= start_idx
        error("Could not find valid table section.")
    end

    table_lines = lines[(start_idx+1):(end_idx-1)]

    # Remove empty lines
    table_lines = filter(x -> !isempty(strip(x)), table_lines)

    # Split header and data using tab delimiter
    header = strip.(split(strip(table_lines[1]), '\t'))

    data = [split(strip(line), '\t') for line in table_lines[2:end]]

    # Convert data to numbers if possible
    data = [map(x -> tryparse(Float64, x) === nothing ? x : parse(Float64, x), row) for row in data]
    df = DataFrame([getindex.(data, i) for i in 1:length(header)], Symbol.(header))

    return df
end



"""
Convert a SCAPS file to a DataFrame.
"""
function ct_to_df(filename::String)
    df = CSV.read(filename, DataFrame)

    df.x .= df.x .* 1e6  # Convert x from cm to um
    df.n .= df.n .* 1e-6  # Convert n from /cm3 to /um3
    df.p .= df.p .* 1e-6  # Convert x from cm to um
    return df
end

"""
Plots carrier density values from SCAPS and ChargeTransport
"""
function compare_densities(df_ct::DataFrame, df_scaps::DataFrame; xrange=nothing)
    df_ct_plot = df_ct
    df_scaps_plot = df_scaps

    if xrange !== nothing
        xmin, xmax = xrange
        df_ct_plot = filter(row -> xmin <= row["x"] <= xmax, df_ct)
        df_scaps_plot = filter(row -> xmin <= row["x(um)"] <= xmax, df_scaps)
    end
    # Plotting
    figure(figsize=(10, 6))
    plot(df_ct_plot[!, "x"], df_ct_plot[!, "n"], label="CT (n)", color="blue")
    plot(df_scaps_plot[!, "x(um)"], df_scaps_plot[!, "n(/cm3)"], label="SCAPS (n)", color="blue", linestyle="--")
    plot(df_ct_plot[!, "x"], df_ct_plot[!, "p"], label="CT (p)", color="red")
    plot(df_scaps_plot[!, "x(um)"], df_scaps_plot[!, "p(/cm3)"], label="SCAPS (p)", color="red", linestyle="--")

    yscale("log")
    xlabel("x(um)")
    ylabel("n(/cm3)")
    title("Comparison of Carrier Densities")
    legend()
    grid(true)
end

"""
Plots carrier density values from SCAPS and ChargeTransport
"""
function compare_bands(df_ct::DataFrame, df_scaps::DataFrame; xrange=nothing)
    # Plotting
    df_ct_plot = df_ct
    df_scaps_plot = df_scaps

    if xrange !== nothing
        xmin, xmax = xrange
        df_ct_plot = filter(row -> xmin <= row["x"] <= xmax, df_ct)
        df_scaps_plot = filter(row -> xmin <= row["x(um)"] <= xmax, df_scaps)
    end

    figure(figsize=(10, 6))
    plot(df_ct_plot[!, "x"], df_ct_plot[!, "Ec"], color="black")
    plot(df_ct_plot[!, "x"], df_ct_plot[!, "Ev"], color="black")
    plot(df_ct_plot[!, "x"], df_ct_plot[!, "EFn"], color="blue")
    plot(df_ct_plot[!, "x"], df_ct_plot[!, "EFp"], color="red")
    plot(df_scaps_plot[!, "x(um)"], df_scaps_plot[!, "Ec(eV)"], color="black", linestyle="--")
    plot(df_scaps_plot[!, "x(um)"], df_scaps_plot[!, "Ev(eV)"], color="black", linestyle="--")
    plot(df_scaps_plot[!, "x(um)"], df_scaps_plot[!, "Fn(eV)"], color="blue", linestyle="--")
    plot(df_scaps_plot[!, "x(um)"], df_scaps_plot[!, "Fp(eV)"], color="red", linestyle="--")

    xlabel("x(um)")
    ylabel("Energy (eV)")
    title("Comparison of Energy Bands")
    grid(true)
end

"""
Create a CellInfo struct at a certain spatial position to learn all key information
"""

# TODO: Combine these functions
function find_info_scaps(df::DataFrame, x::Float64)
    idx = argmin(abs.(df[!, "x(um)"] .- x))
    return CellInfo(
        df[idx, "x(um)"],         # position
        df[idx, "n(/cm3)"],       # electron density
        df[idx, "p(/cm3)"],       # hole density
        df[idx, "Ec(eV)"],            # conduction band edge
        df[idx, "Ev(eV)"],            # valence band edge
        df[idx, "Fn(eV)"],           # electron quasi-Fermi level
        df[idx, "Fp(eV)"]            # hole quasi-Fermi level
    )
end

function find_info_ct(df::DataFrame, x::Float64)
    idx = argmin(abs.(df[!, "x"] .- x))
    return CellInfo(
        df[idx, "x"],         # position
        df[idx, "n"],       # electron density
        df[idx, "p"],       # hole density
        df[idx, "Ec"],            # conduction band edge
        df[idx, "Ev"],            # valence band edge
        df[idx, "EFn"],           # electron quasi-Fermi level
        df[idx, "EFp"]            # hole quasi-Fermi level
    )
end

function scaps_to_df_iv(filename::String)
    lines = readlines(filename)
    # Find IV/START and IV/END
    start_idx = findfirst(x -> occursin("IV/START", x), lines)
    end_idx = findfirst(x -> occursin("IV/END", x), lines)
    if isnothing(start_idx) || isnothing(end_idx) || end_idx <= start_idx
        error("Could not find IV/START and IV/END or invalid order.")
    end
    # The header is the line after IV/START
    header_line = strip(lines[start_idx + 1])
    header = split(header_line)
    # Data lines are between header and IV/END
    data_lines = lines[(start_idx + 2):(end_idx - 1)]
    # Remove empty lines
    data_lines = filter(x -> !isempty(strip(x)), data_lines)
    # Split and parse each data line
    data = [split(strip(line)) for line in data_lines]
    # Convert to Float64 if possible
    data = [map(x -> tryparse(Float64, x) === nothing ? x : parse(Float64, x), row) for row in data]
    df = DataFrame([getindex.(data, i) for i in 1:length(header)], Symbol.(header))
    return df
end

function compare_iv(df_ct::DataFrame, df_scaps::DataFrame)
    # Plotting IV curves
    figure(figsize=(10, 6))
    plot(df_ct[!, "V"], df_ct[!, "J"], label="CT", color="blue")
    plot(df_scaps[!, "v(V)"], df_scaps[!, "jtot(mA/cm2)"], label="SCAPS", color="red", linestyle="--")

    xlabel("Voltage (V)")
    ylabel("Current Density (mA/cm2)")
    title("IV Curve Comparison")
    legend()
    grid(true)
end