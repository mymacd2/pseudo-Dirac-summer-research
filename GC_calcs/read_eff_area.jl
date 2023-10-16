using DelimitedFiles, StatsBase, Interpolations

eff_a = readdlm("gc_data/EffA_GC_approx_2.csv", ',')
effarea_non_extrap = Interpolations.interpolate((eff_a[:, 1] ./ 1000,), eff_a[:, 2], Gridded(Linear())) # in m^2
effarea = extrapolate(effarea_non_extrap, 0.0)

# The min and max energy vals (in TeV) we consider based on the limits of the effective area function
emin = 0.7249213596719925
emax = 1099.9079675727165