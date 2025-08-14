using ChargeTransport
using DataFrames
using CSV
using PyPlot

#######################################################
#######################################################

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

# save profile as a CSV file
function save_device_profile(filename, solution, ctsys)
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

function save_iv(filename::String, biasValues, IV)
    df = DataFrame(
        V = biasValues,
        J = IV
    )

    CSV.write(filename, df)
end