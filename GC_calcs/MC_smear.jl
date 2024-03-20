include("read_energy_res.jl")

using Distributions

##############################################################################################################

# Functions that smear the Monte Carlo in angle

function angsmearMC(ls, bs)
    numν = length(ls)
    ls_sm = zeros(numν)
    bs_sm = zeros(numν)

    for i in 1:numν
        dist = MvNormal([ls[i], bs[i]], 0.122)
        l_sm, b_sm = rand(dist, 1)

        # Changing the coordinates to fit in the \ell, b bounds
        if b_sm > π/2
            b_sm = π - b_sm
            if l_sm > 0
                l_sm -= π
            else
                l_sm += π
            end
        end
        if b_sm < -π/2
            b_sm = -π - b_sm
            if l_sm > 0
                l_sm -= π
            else
                l_sm += π
            end
        end

        if l_sm > π
            l_sm -= 2π
        end
        if l_sm < -π
            l_sm += 2π
        end
            
        ls_sm[i] = l_sm
        bs_sm[i] = b_sm
    end
    return (ls_sm, bs_sm)
end

##############################################################################################################

# Functions that smear the Monte Carlo in (log) energy

# Reshape the original matrix into a 25x8x25x8 array
reshaped_eres= reshape(erestrue, 8, 25, 8, 25)

# Sum along the 2nd and 4th dimensions to collapse the 8x8 blocks
eresbinned = sum(reshaped_eres, dims=(1, 3))[1, :, 1, :]

for i in 1:25
    normconst = sum(eresbinned[:,i])
    eresbinned[:,i] .= eresbinned[:,i]/normconst
end


function esmear(u::Float64, eresmat)

    loges25 = range(log10(5e-1), log10(4e3), 25)

    du = loges25[2] - loges25[1]

    bin_index(z) = ceil((z - loges25[1]) / du)
    i = Int(bin_index(u))
    if i < 1
        i = 1
    end

    e_cdf = ecdf(loges25; weights=Weights(eresmat[:, i]))

    if i < 23
        iend = i + 3
    else
        iend = i
    end

    inv_cdf_non_extrap = Interpolations.interpolate((vcat([0], e_cdf.(loges25[2:iend])),), loges25[1:iend], Gridded(Linear())) # in m^2
    inv_cdf = extrapolate(inv_cdf_non_extrap, Interpolations.Flat())

    return inv_cdf

end

function esmear(i::Int, eresmat)

    loges25 = range(log10(5e-1), log10(4e3), 25)

    e_cdf = ecdf(loges25; weights=Weights(eresmat[:, i]))

    if i < 23
        iend = i + 3
    else
        iend = i
    end

    inv_cdf_non_extrap = Interpolations.interpolate((vcat([0], e_cdf.(loges25[2:iend])),), loges25[1:iend], Gridded(Linear())) # in m^2
    inv_cdf = extrapolate(inv_cdf_non_extrap, Interpolations.Flat())

    return inv_cdf

end


loges25 = range(log10(5e-1), log10(4e3), 25)

esmearvec = [esmear(i, eresbinned) for i in 1:25]


function smearuMC(us, esmearvec)

    loges25 = range(log10(5e-1), log10(4e3), 25)

    du = loges25[2] - loges25[1]

    bin_index(z) = ceil((z - loges25[1]) / du)

    lenus = length(us)

    smeared_us = zeros(lenus)

    for i in 1:lenus
        j = Int(bin_index(us[i]) + 1)
        usm = esmearvec[j](rand())
        smeared_us[i] = usm
    end

    return smeared_us

end