# Parameters for Si PERC solar cell simulation

iphin = 1 # electron quasi Fermi potential
iphip = 2 # hole quasi Fermi potential
numberOfCarriers = 2 # electrons, holes

########## device geometry ##########

# region numbers
regionCz = 1
regionPoly = 2
regions = [regionCz, regionPoly]
numberOfRegions = length(regions)

# boundary region numbers
bregionCz = 1
bregionPoly = 2
bregionJ1 = 3

## length domains
h_cz = 100.0 * μm
h_poly = 60.0 * nm
h_total = h_cz + h_poly
heightLayers = [
    h_cz,
    h_cz + h_poly
]
########## physical values ##########

## charge numbers
zn = -1
zp = 1

## temperature
T = 300.0 * K

## band edge energies
En = [-4.05, -3.9] .* eV  
Ep = [-5.17, -5.1] .* eV

## effective densities of density of states
Nn = [2.8e25, 1.0e26] ./ (m^3) 
Np = [1.04e25, 1.0e26] ./ (m^3)

## mobilities
μn = [1.4e-1, 1.1e-2] .* (m^2) / (V * s)
μp = [4.5e-2, 4.2e1] .* (m^2) / (V * s)

## relative dielectric permittivity
ε = [11.9, 11.9] .* 1.0  # a-Si and Cz-Si have similar values

## radiative recombination
r0 = [0.0, 0.0] .* cm^3 / s  # Not provided, set to zero

## life times and trap densities
τn = [1.0e-8, 1.0e13] .* s  # a-Si (minority), Cz-Si, a-Si (minority)
τp = [1.0e13, 1.0e-8] .* s  # a-Si (majority), Cz-Si, a-Si (majority)

generation_uniform = [2.50e18, 1.0e15] ./ (m^3 * s)

## doping
Cn = 1.0e22 / (m^3)
Cp = 1.0e23 / (m^3)

# effective intrinsic carrier density
# Ci_eff = 1.0e22 / (m^3)

UT = kB * T / q