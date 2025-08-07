# Parameters for Si PERC solar cell simulation

iphin = 1 # electron quasi Fermi potential
iphip = 2 # hole quasi Fermi potential
numberOfCarriers = 2 # electrons, holes

########## device geometry ##########

# region numbers
regionDonor = 1
regionIntrinsic = 2
regionAcceptor = 3
regions = [regionDonor, regionIntrinsic, regionAcceptor]
numberOfRegions = length(regions)

# boundary region numbers
bregionDonor = 1
bregionAcceptor = 2
bregionJ1 = 3
bregionJ2 = 4

## length domains
h_ndoping = 20.0 * nm
h_absorber = 100.0 * μm
h_pdoping = 10.0 * nm
h_total = h_ndoping + h_absorber + h_pdoping
heightLayers = [
    h_ndoping,
    h_ndoping + h_absorber,
    h_ndoping + h_absorber + h_pdoping,
]
########## physical values ##########

## charge numbers
zn = -1
zp = 1

## temperature
T = 300.0 * K

## band edge energies
En = [-3.9, -4.05, -3.9] .* eV  
Ep = [-5.7, -6.17, -5.7] .* eV

## effective densities of density of states
Nn = [1.0e26, 2.8e26, 1.0e26] ./ (m^3) 
Np = [1.0e26, 1.0e26, 1.0e26] ./ (m^3)

## mobilities
μn = [2.0e-3, 0.04, 2.0e-3] .* (m^2) / (V * s)
μp = [7.5e-3, 0.045, 5.0e-4] .* (m^2) / (V * s)

## relative dielectric permittivity
ε = [11.9, 11.9, 11.9] .* 1.0  # a-Si and Cz-Si have similar values

## radiative recombination
r0 = [0.0, 0.0, 0.0] .* cm^3 / s  # Not provided, set to zero

## life times and trap densities
τn = [1.0e13, 1.0e-4, 1.0e-8] .* s  # a-Si (minority), Cz-Si, a-Si (minority)
τp = [1.0e-8, 1.0e-4, 1.0e13] .* s  # a-Si (majority), Cz-Si, a-Si (majority)

incidentPhotonFlux = [0.0, 4.30e20, 0.0] ./ (m^2 * s)
absorption = [0.0, 1.0e7, 0.0] ./ m
generationPeak = h_ndoping

# generation_uniform = [0.0, 2.64e25, 0.0] ./ (m^3 * s)

## doping
Cn = 1.0e22 / (m^3)
Cp = 1.0e23 / (m^3)

# effective intrinsic carrier density
Ci_eff = 1.0e22 / (m^3)

UT = kB * T / q