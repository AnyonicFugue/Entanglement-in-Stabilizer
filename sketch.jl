# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")


n_stab=4
dof_size=4

stab_generators=rand(Bool,(n_stab,dof_size)) # The first index is the index of the stabilizer. The second index is the index of the d.o.f.
stab_generator_recombined=zeros(Bool,(n_stab,dof_size)) # For each stabilizer the first index is the index of the stabilizer. The second index is the index of the d.o.f.

recomb_mat=zeros(Bool,(n_stab,n_stab)) # A rectangular pivotal matrix. The size is n_stab, the number of stabilizers.
for i in 1:n_stab
    recomb_mat[i,i]=true
    recomb_mat[i,i+1:n_stab]=rand(Bool,n_stab-i)
end

stab_generator_recombined=zeros(Bool,(n_stab,dof_size)) # For each stabilizer the first index is the index of the stabilizer. The second index is the index of the d.o.f.
for i in 1:n_stab
    for j in 1:n_stab
        if(recomb_mat[i,j])
            @views stab_generator_recombined[i,:]=stab_generator_recombined[i,:] .⊻ stab_generators[j,:]
        end
    end
end

println("stab_generators:",stab_generators)
println("recomb_mat:",recomb_mat)
println("stab_generator_recombined:",stab_generator_recombined)