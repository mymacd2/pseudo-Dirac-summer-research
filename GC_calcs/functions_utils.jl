using StaticArrays

function oscprob(et, dm2, leff)

    # Norms squared of the PMNS matrix (exact vals might change but always needs to be unitary)
    u = @SMatrix [0.674743 0.302844 0.0224125;
                  0.0946105 0.360415 0.544974;
                  0.230646 0.33674 0.432613]

    # Conversion factor to go from kpc to 1/eV
    convfactor = 3.086e19 * 5.06773093741 * 1e6

    leff *= convfactor

    # Assuming uniform mass splitting over all three mass states
    osc = (cos((dm2 * leff)/(4*et*1e12)))^2

    prob_ee = osc*((u[1]*u[1]) + (u[4]*u[4]) + (u[7]*u[7]))
    prob_μe = osc*((u[1]*u[2]) + (u[4]*u[5]) + (u[7]*u[8]))

    prob_eτ = osc*((u[3]*u[1]) + (u[6]*u[4]) + (u[9]*u[7]))
    prob_μτ = osc*((u[3]*u[2]) + (u[6]*u[5]) + (u[9]*u[8]))

    prob_eμ = osc*((u[2]*u[1]) + (u[5]*u[4]) + (u[8]*u[7]))
    prob_μμ = osc*((u[2]*u[2]) + (u[5]*u[5]) + (u[8]*u[8]))

    prob_e = 0.333333prob_ee + 0.666666prob_μe
    prob_τ = 0.333333prob_eτ + 0.666666prob_μτ
    prob_μ = 0.333333prob_eμ + 0.666666prob_μμ

    # νμ contribution comes from the 25% chance of a neutral current interaction, which appears as a cascade
    return prob_e + prob_τ + 0.25prob_μ

end


function poissonlog(data, hyp)
    # Avoiding weird /0 scenarios
    if hyp < 1e-20
        hyp = 1e-20
    end
    val = log((hyp^data) * exp(-hyp) / gamma(data+1))
end


function lrt(null, alt)
    altsummand = sum(poissonlog.(null, alt) .- poissonlog.(null, null))
    ts = -2 * altsummand
end

function cartesian(r, l, b)
    x = 8.5 - r*cos(l)*cos(b)
    y = r*sin(l)*cos(b)
    z = r*sin(b)
    return [x, y, z]
end

function cartx(r, l, b)
    x = 8.5 - r*cos(l)*cos(b)
    return x
end

function carty(r, l, b)
    y = r*sin(l)*cos(b)
    return y
end

function cartz(r, l, b)
    z = r*sin(b)
    return z
end

#=
function energy_cut(rs, ls, bs, rws, us, oneweights; umin=log10(umin), umax=log10(umax))
    concatmat = hcat(rs, ls, bs, rws, us, oneweights)
    filtered_matrix_rows = [row for row in eachrow(concatmat) if umin <= row[5] <= umax]
    filtered_mat = hcat(filtered_matrix_rows...)
    fmat = transpose(filtered_mat)
end
=#

function energy_cut(x...; umin=log10(umin), umax=log10(umax))
    concatmat = hcat(x...)
    filtered_matrix_rows = [row for row in eachrow(concatmat) if umin <= row[5] <= umax]
    filtered_mat = hcat(filtered_matrix_rows...)
    return transpose(filtered_mat)
end

function r_cut(x...; rmin=0.3)
    concatmat = hcat(x...)
    filtered_matrix_rows = [row for row in eachrow(concatmat) if rmin <= row[1]]
    filtered_mat = hcat(filtered_matrix_rows...)
    return transpose(filtered_mat)
end

power_law_flux(E, γ; ϕ₀=1, E0=1) = ϕ₀ * (E / E0)^(-γ)