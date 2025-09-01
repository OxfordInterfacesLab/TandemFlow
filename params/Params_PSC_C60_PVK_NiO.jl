using LessUnitful: @local_unitfactors, @ufac_str, @ph_str

constants = ChargeTransport.constants

# used unit factors
nm = ufac"nm"
cm = ufac"cm"
K = ufac"K"
m = ufac"m"
V = ufac"V"
s = ufac"s"

eV = constants.q * V

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
h_ndoping = 20.0 * nm
h_intrinsic = 100.0 * nm
h_pdoping = 13.0 * nm
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
En = [-4.00, -4.0, -2.05] .* eV
Ep = [-6.00, -5.68, -5.65] .* eV
Ea = [0.0, -4.5, 0.0] .* eV

## effective densities of density of states
Nn = [1.0e26, 2.75e24, 2.5e26] ./ (m^3)
Np = [1.0e26, 3.80e24, 2.5e26] ./ (m^3)
Na = [0.0, 4.0e27, 0.0] ./ (m^3)

## mobilities
μn = [8.0e-6, 1.0e-3, 2.8e-4] .* (m^2) / (V * s)
μp = [3.5e-7, 1.0e-3, 2.8e-4] .* (m^2) / (V * s)
μa = [0.0, 1.0e-17, 0.0] .* (m^2) / (V * s)

## relative dielectric permittivity
ε = [2.5, 60.0, 11.7] .* 1.0

## radiative recombination
r0 = [1.0e-20, 1.0e-20, 1.0e-20] .* cm^3 / s

## life times and trap densities
τn = [1.0e-7, 2.5e-9, 1.0e-6] .* s
τp = [1.0e-7, 2.5e-9, 1.0e-6] .* s

Augn = [1.0e-31, 1.0e-31, 1.0e-31] .* (cm^6) / (s)
Augp = [1.0e-31, 1.0e-31, 1.0e-31] .* (cm^6) / (s)

## SRH trap densities
nTrapDensity = [1.0e20, 4.0e21, 1.0e20] ./ (m^3)
pTrapDensity = [1.0e20, 4.0e21, 1.0e20] ./ (m^3)

## generation
incidentPhotonFlux = [0.0, 4.3e21, 0.0] ./ (m^2 * s)
absorption = [0.0, 1.0e7, 0.0] ./ m
generationPeak = h_ndoping

generation_uniform = [0.0, 2.64e27, 0.0] ./ (m^3 * s)

## doping
Cn = 1.0e24 / (m^3)
Cp = 3.0e24 / (m^3)
Ca = 1.6e25 / (m^3)

UT = constants.k_B * T / constants.q
