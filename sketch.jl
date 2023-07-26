# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")


function test(n)
    a=zeros(Bool,(n,n))
    b=ones(Bool,(n,n))
    c=rand(Bool,(n,n))

    # println(gaussian_elimination!(a,b,false))
    println(gaussian_elimination!(c,a,false))
    println(c)
    

end

test(4)