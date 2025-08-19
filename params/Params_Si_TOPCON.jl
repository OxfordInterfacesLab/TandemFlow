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
h_emitter = 1.0 * μm
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
r0 = [4.22e-15, 4.22e-15] .* (cm^3) / (s)

## auger recombination
Augn = [6.0e-30, 7.02e-31] .* (cm^6) / (s)
Augp = [9.0e-31, 3.5e-30] .* (cm^6) / (s)

## life times and trap densities
τcz = [1.0, 1.0] .* s  # a-Si (minority), Cz-Si, a-Si (minority)
τpoly = [1.0, 1.0] .* s  # a-Si (majority), Cz-Si, a-Si (majority)

## doping
Cem = 6.0e24 / (m^3)
Ccz = 3.0e21 / (m^3)
Cpoly = 2.0e26 / (m^3)

UT = kB * T / q