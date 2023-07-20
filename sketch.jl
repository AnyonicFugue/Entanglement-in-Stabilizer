# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")

vec_arr=zeros(Bool,(5,5))

vec_arr[1,:]=[1,0,0,0,0]
vec_arr[2,:]=[1,1,0,0,0]
vec_arr[3,:]=[0,1,1,0,0]
vec_arr[4,:]=[0,0,1,0,0]
vec_arr[5,:]=[0,1,0,0,0]


coeff_arr=zeros(Bool,(5,5)) # This array must be square.

println(gaussian_elimination!(vec_arr,coeff_arr))
println(vec_arr)
println(coeff_arr)



vec_arr=zeros(Bool,(5,3))

vec_arr[1,:]=[1,0,0]
vec_arr[2,:]=[1,1,0]
vec_arr[3,:]=[0,1,0]


coeff_arr=zeros(Bool,(5,5)) # This array must be square.

println(gaussian_elimination!(vec_arr,coeff_arr))
println(vec_arr)
println(coeff_arr)