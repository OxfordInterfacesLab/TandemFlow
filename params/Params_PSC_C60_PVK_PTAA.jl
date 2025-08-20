# TODO: Preamble

#####################################################################
############################ parameters ############################

########## charge carriers ##########

iphin = 1 # electron quasi Fermi potential
iphip = 2 # hole quasi Fermi potential
iphia = 3
numberOfCarriers = 3 # electrons, holes and anion vacancies

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
h_ndoping = 30.0 * nm
h_intrinsic = 400.0 * nm
h_pdoping = 10.0 * nm
h_total = h_ndoping + h_intrinsic + h_pdoping
heightLayers = [
    h_ndoping,
    h_ndoping + h_intrinsic,
    h_ndoping + h_intrinsic + h_pdoping,
]
########## physical values ##########

## charge numbers
zn = -1
zp = 1
za = 1

## temperature
T = 300.0 * K

## band edge energies
En = [-3.9, -3.9, -2.5] .* eV
Ep = [-5.9, -5.53, -5.5] .* eV
Ea = [0.0, -4.66, 0.0] .* eV

## effective densities of density of states
Nn = [1.0e26, 1.0e24, 1.0e26] ./ (m^3)
Np = [1.0e26, 2.2e24, 1.0e26] ./ (m^3)
Na = [0.0, 1.0e27, 0.0] ./ (m^3)


## mobilities
μn = [1.0e-6, 5.0e-4, 1.0e-8] .* (m^2) / (V * s)
μp = [1.0e-6, 2.0e-4, 1.0e-8] .* (m^2) / (V * s)
μa = [0.0, 1.0e-12, 0.0] .* (m^2) / (V * s)


## relative dielectric permittivity
ε = [5.0, 22.0, 3.5] .* 1.0

## radiative recombination
r0 = [6.8e-17, 3.6e-18, 6.3e-17] .* cm^3 / s

## life times and trap densities
τn = [1.0, 4.0e-8, 1.0] .* s
τp = [1.0, 4.0e-8, 1.0] .* s

## SRH trap energies
nTrapDensity = [1.59e9, 4.48e10, 6.33] ./ (m^3)
pTrapDensity = [1.59e9, 4.48e10, 6.33] ./ (m^3)

## generation
incidentPhotonFlux = [0.0, 4.3e21, 0.0] ./ (m^2 * s)
absorption = [0.0, 1.0e7, 0.0] ./ m
generationPeak = h_ndoping

generation_uniform = [0.0, 2.64e27, 0.0] ./ (m^3 * s)

## doping
Cn = 2.09e24 / (m^3)
Cp = 2.09e24 / (m^3)
Ca = 6.0e22 / (m^3)

UT = kB * T / q
