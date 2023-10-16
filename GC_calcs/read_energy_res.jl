using DelimitedFiles, StatsBase, Interpolations

# Getting the true energy resolution matrix

true_eres = readdlm("gc_data/eres_matrix.dat")
true_eres[1:25, 175:end] = zeros(25, 26)

true_eres = 10 .^ true_eres

for j in 1:200
    for i in 1:200
        if true_eres[i, j] == 1
            true_eres[i, j] = 0
        end
    end
end

erestrue = transpose(true_eres)

loges = range(log10(5e-1), log10(4e3), 200)
du = loges[2] - loges[1]

for i in 1:200
    normconst = sum(erestrue[:,i]*du)
    erestrue[:,i] .= erestrue[:,i]/normconst
end