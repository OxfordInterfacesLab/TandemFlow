# src/stacks.jl
#
# Reusable optical-stack constructors. Lifted verbatim from
# test/day7_checkpoint.jl so that day8_eqe_test.jl (and every later script)
# can build a Stack WITHOUT re-running the Day 7 checkpoint as a side effect.
#
# ASSUMES `Stack` (generation.jl) and `load_nk_data`/`NKData` (materials.jl)
# are already in scope. In your driver:
#     include("src/generation.jl")   # Stack, tmm_generation_slice, ...
#     include("src/materials.jl")    # load_nk_data, NKData, interpolate_nk
#     include("src/stacks.jl")       # build_stack1(), ...

const NKDIR = joinpath(@__DIR__, "..", "data", "nk")   # CSV folder
nkfile(name) = joinpath(NKDIR, name * ".csv")

"""
    build_stack1() -> Stack

Stack 1 (Prof. Sebastian's paper optic, Si-terminated, no HTL):
    air / MgF2(100) / ITO(50) / SnO2(15) / C60(20) /
    CsFAPbIBr(800) / ITO_inter(15) / Si(∞)
8 layers, perovskite = layer 6. `air` is the `nothing` vacuum placeholder.
ITO is reused for both the front electrode and the 15 nm interconnect.
"""
function build_stack1()
    mgf2 = load_nk_data(nkfile("MgF2"),      name = "MgF2")
    ito  = load_nk_data(nkfile("ITO"),       name = "ITO")
    sno2 = load_nk_data(nkfile("SnO2"),      name = "SnO2")
    c60  = load_nk_data(nkfile("C60"),       name = "C60")
    pvk  = load_nk_data(nkfile("CsFAPbIBr"), name = "CsFAPbIBr")
    si   = load_nk_data(nkfile("Si"),        name = "Si")

    names = ["air", "MgF2", "ITO", "SnO2", "C60", "CsFAPbIBr", "ITO_inter", "Si"]
    nk    = [nothing, mgf2, ito, sno2, c60, pvk, ito, si]
    d_nm  = [Inf, 100.0, 50.0, 15.0, 20.0, 800.0, 15.0, Inf]
    return Stack(names, nk, d_nm, 6)
end

# Stacks 2 and 3 (spiro-TTB inserted; TMM identical for both) arrive on Day 9/10.
# When you add them, perovskite is still layer 6 in a 9-layer list:
#     air / MgF2 / ITO / SnO2 / C60 / CsFAPbIBr / spiro-TTB(50) / ITO_inter(15) / Si
# function build_stack2() ... return Stack(names, nk, d_nm, 6) end
"""
    build_stack2() -> Stack
 
Stack 2 (realistic device: spiro-TTB HTL inserted, Si-terminated):
    air / MgF2(100) / ITO(50) / SnO2(15) / C60(20) /
    CsFAPbIBr(800) / spiro-TTB(50) / ITO_inter(15) / Si(∞)
9 layers, perovskite STILL = layer 6 (spiro-TTB inserted *after* it at slot 7,
so the perovskite index does not move). ITO_inter -> slot 8, Si -> slot 9.
 
TMM-identical to Stack 3 — Stack 3 differs only in the ChargeTransport back
boundary condition (SchottkyContact, Day 10), not in the optics.
 
nk file: data/nk/SpiroTTB.csv  (loaded under the display name "spiro-TTB").
"""
function build_stack2()
    mgf2  = load_nk_data(nkfile("MgF2"),      name = "MgF2")
    ito   = load_nk_data(nkfile("ITO"),       name = "ITO")
    sno2  = load_nk_data(nkfile("SnO2"),      name = "SnO2")
    c60   = load_nk_data(nkfile("C60"),       name = "C60")
    pvk   = load_nk_data(nkfile("CsFAPbIBr"), name = "CsFAPbIBr")
    spiro = load_nk_data(nkfile("SpiroTTB"),  name = "spiro-TTB")   # <-- NEW
    si    = load_nk_data(nkfile("Si"),        name = "Si")
 
    names = ["air", "MgF2", "ITO", "SnO2", "C60", "CsFAPbIBr",
             "spiro-TTB", "ITO_inter", "Si"]
    nk    = [nothing, mgf2, ito, sno2, c60, pvk, spiro, ito, si]
    d_nm  = [Inf, 100.0, 50.0, 15.0, 20.0, 800.0, 50.0, 15.0, Inf]
    return Stack(names, nk, d_nm, 6)   # perovskite is still layer 6
end

function build_stack1_Ag()
    mgf2 = load_nk_data(nkfile("MgF2"),      name = "MgF2")
    ito  = load_nk_data(nkfile("ITO"),       name = "ITO")
    sno2 = load_nk_data(nkfile("SnO2"),      name = "SnO2")
    c60  = load_nk_data(nkfile("C60"),       name = "C60")
    pvk  = load_nk_data(nkfile("CsFAPbIBr"), name = "CsFAPbIBr")
    ag   = load_nk_data(nkfile("Ag"),        name = "Ag")

    names = ["air", "MgF2", "ITO", "SnO2", "C60", "CsFAPbIBr", "Ag"]
    nk    = [nothing, mgf2, ito, sno2, c60, pvk, ag]
    d_nm  = [Inf, 100.0, 50.0, 15.0, 20.0, 800.0, Inf]
    return Stack(names, nk, d_nm, 6)
end