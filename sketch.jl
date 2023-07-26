# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")

import Profile

a=[1,2,3,4,5]

a[1],a[3]=3,5
println(string(64)*string(32))
