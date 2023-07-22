# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")


OldISG=zeros(Bool,(3,4))
Measurements=zeros(Bool,(3,4))


NewISG_vec=Vector{Vector{Bool}}(undef,0)
println(NewISG_vec)