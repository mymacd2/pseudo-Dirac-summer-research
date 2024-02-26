using DelimitedFiles, StatsBase, Interpolations

function read_models(model1; model2=0, get_es=true)
    # Getting the model-dependent neutrino distribution
    model_1 = readdlm("/Users/millermacdonald/Desktop/Research_shit/GP_nu_sim_data/$model1.txt", comments=true)

    if model2 !== 0
        model_2 = readdlm("/Users/millermacdonald/Desktop/Research_shit/GP_nu_sim_data/$model2.txt", comments=true)
        νmodel = vcat(model_1, model_2)
        numv = size(νmodel)[1]
    else
        νmodel = model_1
        numv = size(νmodel)[1]
    end


    xs = νmodel[:, 4]
    ys = νmodel[:, 5]
    zs = νmodel[:, 6]

    rweights = vec(νmodel[:, 8])

    tuples = [zeros(3) for _ in 1:numv]

    for i in 1:numv
            tuples[i] = [xs[i], ys[i], zs[i]]
    end

    function galcoord(carttuple)
        x, y, z = carttuple
        r = sqrt((8.5-x)^2 + y^2 + z^2)
        l = atan(-y, (8.5-x))
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

    if get_es==true
        logνes = log10.(νmodel[:, 3] ./ 1e12)
        logCRes = log10.(νmodel[:, 7] ./ 1e12)
        return (rs, ls, bs, rweights, logνes, logCRes)
    else
        return (rs, ls, bs, rweights)
    end
end

# nuweights[1:5, :, :] .= 0
#==
newarray = zeros(size(nuweights))

for i in 1:length(rbins)-1
    for j in 1:length(lbins)-1
        for k in 1:length(bbins)-1
            newarray[i, j, k] = nuweights[i, j, k] * (1/(rbins[i]+0.05)^2)
        end 
    end
end
==#