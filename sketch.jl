# Compare the performances of bool and GF(2)
# Conclusion: Use Bool. Bool is faster than GF(2).
include("gaussian_elimination.jl")
include("dynamic_update.jl")
include("utils.jl")

n=10

stab=zeros(Bool,n,2*n)

for i in range(1,n)
    stab[i,2*i-1]=true # Pauli-X on site i

    if(i<n)
        stab[i,2*i+2]=true # Pauli-Z on site i+1
    else
        stab[i,2]=true
    end

    if(i>1)
        stab[i,2*i-2]=true # Pauli-Z on site i-1
    else
        stab[i,2*n]=true
    end
end

for i in range(1,n-1)
    Evaluate_EE(stab,i)
end