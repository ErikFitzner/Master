using LinearAlgebra

function g_r(r::Tuple{Int,Int}; N=1200)
    kx = range(-π, stop=π, length=N)
    ky = range(-π, stop=π, length=N)
    dk = 2π / N
    sumval = 0 + 0im
    for kxi in kx, kyi in ky
        eps = - (cos(kxi) + cos(kyi))  # take J=1, occupied if eps<0 => cos kx + cos ky > 0
        if eps < 0  # occupied region (half filling)
            sumval += exp(im*(kxi*r[1] + kyi*r[2]))
        end
    end
    return sumval * (dk^2) / (4π^2)    # normalization (1/(2π)^2) and dk^2 grid weight
end

for r in ((0,0),(1,0),(1,1),(2,0))
    g = g_r(r; N=1400)
    C = -abs2(g)
    println(C)
end
