#= What should the code do?

Convention: Use Bool type. xor serves as + in F_2.

1. Find the solution space of a linear equation.
    This amounts first to use gaussian elimination to reduce the matrix to an echoleon form (and record the coefficients), then find a maximal linear independent subset of the zero echoleon vectors.
    Proposition. All echoleon vectors must be linear independent combinations of the basis vectors, since the transformation matrix is upper-diagonal.

2. Find a maximally linear independent subset of a set of vectors via gaussian elimination.
    Specifically, use the elimination to reduce to an echoleon form and record the coefficients.

=#
@inline function swap!(v1::Vector{Bool},v2::Vector{Bool})
    # The function swaps the two vectors.
    # The function is inlined to avoid the overhead of function call.
    v1=v1 .⊻ v2
    v2=v1 .⊻ v2
    v1=v1 .⊻ v2
end

function gaussian_elimination!(Vecs::Matrix{Bool},Coefficients::Matrix{Bool})
    # The function overwrites the input vectors with the echoleon vectors and records the coefficients in the second vector.
    # The elimination works in the 2-element finite field. Note that summation is xor.

    # Note that the array Coefficients is assumed to be initialized to all zero.

    l=length(Vecs[1,:]) # The length of the vectors.
    n=length(Vecs[:,1]) # The number of vectors.

    for i in 1:n # Number of vectors.
        Coefficients[i,i]=true # The coefficient of the i-th vector is 1.
    end

    for i in 1:min(n,l)
        # Find the first vector with non-zero i-th component.

        found=false # Whether we have found such a vector. If not found we can directly jump to the next one.

        for j in i:n
            if (Vecs[j,i]==true)&&(found==false)
                swap!(Vecs[i,:],Vecs[j,:])
                swap!(Coefficients[i,:],Coefficients[j,:])
                found=true
            end
        end

        # Eliminate the i-th component of the rest vectors.
        if(found)  
            for j in i+1:n
                if Vecs[j,i]==true
                    Vecs[j,:]=Vecs[j,:] .⊻ Vecs[i,:]
                    Coefficients[j,:]=Coefficients[j,:] .⊻ Coefficients[i,:]
                end
            end
        end
    end
        
    # Compute the rank, i.e. nonzero echoleon vectors.
    rank=0
    for i in 1:n
        if !iszero(Vecs[i,:])
            rank+=1
        end
    end

    return rank
end