# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")
include("utils.jl")

stab=zeros(Bool,4)

println(size(stab,1))
