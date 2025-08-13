using ChargeTransport
using DataFrames
using CSV
using PyPlot

function scaps_to_df_generation(filename::String)
    lines = readlines(filename)
    # Find the line containing 'GEN'
    gen_idx = findfirst(x -> occursin("GEN", x), lines)
    if isnothing(gen_idx) || gen_idx == length(lines)
        error("No GEN line or header found in file.")
    end
    # The header is the line after 'GEN'
    header_idx = gen_idx + 1
    header = strip.(split(strip(lines[header_idx]), '\t'))
    # Data starts after the header
    data_start = header_idx + 1
    data_lines = filter(x -> !isempty(strip(x)), lines[data_start:end])
    # Split each data line into fields using tab delimiter
    data = [split(strip(line), '\t') for line in data_lines]
    # Convert to numbers if possible
    data = [map(x -> tryparse(Float64, x) === nothing ? x : parse(Float64, x), row) for row in data]
    df = DataFrame([getindex.(data, i) for i in 1:length(header)], Symbol.(header))

    return df
end

function generation_from_scaps(filename::String, xq::Float64)
    df_gen = scaps_to_df_generation(filename)

    x = df_gen[!, "x (um)"]
    Geh = df_gen[!, "Geh (#/cm3.s)"]

    xq = xq / (μm)

    # Find the interval containing xq
    if xq <= x[1]
        return Geh[1] * 1e6
    elseif xq >= x[end]
        return Geh[end] * 1e6
    else
        i = findfirst(i -> x[i] <= xq < x[i+1], 1:length(x)-1)

        if isnothing(i)
            error("xq is out of bounds or data is not sorted.")
        end
        # Exponential interpolation: Geh(x) = Geh1 * exp(log(Geh2/Geh1) * (xq-x1)/(x2-x1))
        x1, x2 = x[i], x[i+1]
        G1, G2 = Geh[i], Geh[i+1]
        return G1 * exp(log(G2/G1) * (xq - x1) / (x2 - x1)) * 1e6
    end
end