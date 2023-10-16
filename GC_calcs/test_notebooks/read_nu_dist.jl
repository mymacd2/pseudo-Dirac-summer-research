using DelimitedFiles, StatsBase, Interpolations

# Getting the model-dependent neutrino distribution
νmodel = readdlm("gc_data/GC_model_2.txt", comments=true)
numv = size(νmodel)[1]

# Histogramming and interpolating the distribution
bsize = 0.3
xbins, ybins, zbins = -30:bsize:30, -30:bsize:30, -3:bsize:3
nubins= fit(Histogram, (vec(νmodel[:, 24]), vec(νmodel[:, 25]), vec(νmodel[:, 26])), (xbins, ybins, zbins), closed=:left)

binvals = nubins.weights ./ (numv*bsize^3)
xedges = view(collect(xbins), 1:length(xbins)-1)
yedges = view(collect(ybins), 1:length(ybins)-1)
zedges = view(collect(zbins), 1:length(zbins)-1)

# Interpolating the distribution to get a continuous function, and assuming zero neutrinos outside the sampled window
probdens_non_extrap = Interpolations.interpolate((xedges, yedges, zedges), binvals, Gridded(Linear()))
probdens = extrapolate(probdens_non_extrap, 0.0);