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