# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")


function test(n)
    a=rand(Bool,(n,n))
    b=rand(Bool,(n,n))
    c=rand(Bool,(n,n))

    @time for i in 1:n
        @views a[i,:]=b[i,:] .⊻ c[i,:]
    end

    @time for i in 1:n
        @views a[i,:]=a[i,:] .⊻ b[i,:]
    end
end

test(4096)