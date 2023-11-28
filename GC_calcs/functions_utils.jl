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