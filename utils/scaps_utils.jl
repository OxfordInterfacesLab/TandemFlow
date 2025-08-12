using DataFrames
using CSV
using PyPlot

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
function compare_densities(df_ct::DataFrame, df_scaps::DataFrame)
    # Plotting
    figure(figsize=(10, 6))
    plot(df_ct[!, "x"], df_ct[!, "n"], label="CT (n)", color="blue")
    plot(df_scaps[!, "x(um)"], df_scaps[!, "n(/cm3)"], label="SCAPS (n)", color="blue", linestyle="--")
    plot(df_ct[!, "x"], df_ct[!, "p"], label="CT (p)", color="red")
    plot(df_scaps[!, "x(um)"], df_scaps[!, "p(/cm3)"], label="SCAPS (p)", color="red", linestyle="--")

    yscale("log")
    xlabel("x(um)")
    ylabel("n(/cm3)")
    title("Comparison of Carrier Densities")
    legend()
    grid(true)
end