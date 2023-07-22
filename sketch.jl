# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")

a=[1,1,1]
OldISG=zeros(Bool,(3,4))
Measurements=zeros(Bool,(3,4))

OldISG[1,:]=[0,1,0,1]
OldISG[2,:]=[1,0,1,1]
OldISG[3,:]=[1,1,0,0]

Measurements[1,:]=[0,1,0,1]
Measurements[2,:]=[1,0,1,0]
Measurements[3,:]=[1,1,0,0]

comm_mat=zeros(Bool,(3,3))


println(calc_comm_matrix!(2,OldISG,Measurements,comm_mat))
println(comm_mat)

