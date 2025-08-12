using DataFrames

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