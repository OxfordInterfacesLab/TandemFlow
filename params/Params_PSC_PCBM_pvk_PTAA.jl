#=

Parameters from the following papers:
- https://doi.org/10.1038/s41467-023-36141-8
- https://doi.org/10.1002/solr.202100219

=#

#####################################################################
############################ parameters ############################

########## charge carriers ##########

iphin = 1 # electron quasi Fermi potential
iphip = 2 # hole quasi Fermi potential
# iphia = 3
numberOfCarriers = 2 # electrons, holes and anion vacancies

########## device geometry ##########

# region numbers
regionDonor = 3
regionIntrinsic = 2
regionAcceptor = 1
regions = [regionAcceptor, regionIntrinsic, regionDonor]
numberOfRegions = length(regions)

# boundary region numbers
bregionAcceptor = 1
bregionDonor = 2
bregionJ1 = 3
bregionJ2 = 4

## length domains
h_ndoping = 30.0 * nm
h_intrinsic = 400.0 * nm
h_pdoping = 10.0 * nm
h_total = h_ndoping + h_intrinsic + h_pdoping
heightLayers = [
    h_pdoping,
    h_pdoping + h_intrinsic,
    h_pdoping + h_intrinsic + h_ndoping,
]

########## physical values ##########

## charge numbers
zn = -1
zp = 1
# za = 1

## temperature
T = 300.0 * K

## band edge energies
En = [-2.6, -3.9, -4.05] .* eV
Ep = [-5.6, -5.7, -6.05] .* eV
Ea = [0.0, -4.66, 0.0] .* eV

## effective densities of density of states
Nn = [1.0e26, 3.1e24, 1.0e26] ./ (m^3)
Np = [1.0e26, 3.1e24, 1.0e26] ./ (m^3)
Na = [0.0, 1.0e26, 0.0] ./ (m^3)

## mobilities
μn = [1.0e-9, 1.0e-3, 1.0e-6] .* (m^2) / (V * s)
μp = [1.0e-9, 1.0e-3, 1.0e-6] .* (m^2) / (V * s)
μa = [0.0, 1.0e-12, 0.0] .* (m^2) / (V * s)

## relative dielectric permittivity
ε = [3.5, 24.0, 4.99] .* 1.0

## radiative recombination
r0 = [2.30e-18, 2.07e-16, 1.02e-18] .* cm^3 / s

## life times and trap densities
τn = [1e14, 1.0e-6, 0.6e-8] .* s
τp = [0.6e-8, 1.0e-6, 1e14] .* s

## SRH trap energies - not used as I have directly input the trap densities
# EI = [-5.0, -4.6, -4.05] .* eV

## generation
incidentPhotonFlux = [0.0, 4.3e21, 0.0] ./ (m^2 * s)
absorption = [1.0e5, 1.0e7, 0.0] ./ m
generationPeak = h_pdoping

generation_uniform = [0.0, 2.64e25, 0.0] ./ (m^3 * s)

## doping
TLdop = 1.0e11

Cn = TLdop / (m^3)
Cp = TLdop / (m^3)
Ca = 1.0e22 / (m^3)

# effective intrinsic carrier density
Ci_eff = 1.0e16 / (m^3)

UT = kB * T / q
