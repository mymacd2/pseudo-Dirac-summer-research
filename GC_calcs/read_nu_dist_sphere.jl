using DelimitedFiles, StatsBase, Interpolations

# Getting the model-dependent neutrino distribution
νmodel = readdlm("gc_data/GC_model_2.txt", comments=true)
numv = size(νmodel)[1]


xs = νmodel[:, 24]
ys = νmodel[:, 25]
zs = νmodel[:, 26]

tuples = [zeros(3) for _ in 1:numv]

for i in 1:numv
        tuples[i] = [xs[i], ys[i], zs[i]]
end

function galcoord(carttuple)
    x, y, z = carttuple
    r = sqrt((8-x)^2 + y^2 + z^2)
    l = atan(y, (8-x))
    b = asin(z / r)
    
    return [r, l, b]

end

for i in 1:numv
    tuples[i] .= galcoord(tuples[i])
end

rs, ls, bs = zeros(numv), zeros(numv), zeros(numv)

for i in 1:numv
    rs[i] = tuples[i][1]
    ls[i] = tuples[i][2]
    bs[i] = tuples[i][3]
end

# We now have what we started with, but in galactic coordinates. Let's bin by r, l, and b:
rbinsize = 0.1
θbinsize = π/100

rbins, lbins, bbins = rbinsize:rbinsize:20+rbinsize, -π:θbinsize:π, -π/2:θbinsize:π/2

nubins_galcoord = fit(Histogram, (rs, ls, bs), (rbins, lbins, bbins), closed=:left)


nuweights = nubins_galcoord.weights
newarray = zeros(size(nuweights))

for i in 1:length(rbins)-1
    for j in 1:length(lbins)-1
        for k in 1:length(bbins)-1
            newarray[i, j, k] = nuweights[i, j, k] * (1/(rbins[i]+0.05)^2)
        end 
    end
end