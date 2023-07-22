#= What should the code do?

Convention: Use Bool type. xor serves as + in F_2.

1. Find the solution space of a linear equation.
    This amounts first to use gaussian elimination to reduce the matrix to an echoleon form (and record the coefficients), then find a maximal linear independent subset of the zero echoleon vectors.
    Proposition. All echoleon vectors must be linear independent combinations of the basis vectors, since the transformation matrix is upper-diagonal.

2. Find a maximally linear independent subset of a set of vectors via gaussian elimination.
    Specifically, use the elimination to reduce to an echoleon form and record the coefficients.

=#
@inline function swaprows!(v::Matrix{Bool},i::Int64,j::Int64)
    # The function swaps the rows of a matrix, v[i,:] and v[j,:].
    # The function is inlined to avoid the overhead of function call.
    v[i,:]=v[i,:] .⊻ v[j,:]
    v[j,:]=v[i,:] .⊻ v[j,:]
    v[i,:]=v[i,:] .⊻ v[j,:]
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

    cur=1 # The current vector we are working on.
    rank=0

    for i in 1:l
        # Find the first vector with non-zero i-th component.

        found=false # Whether we have found such a vector. If not found we can directly jump to the next one.

        for j in cur:n
            if (Vecs[j,i]==true)&&(found==false)
                if(j!=cur)
                    swaprows!(Vecs,cur,j)
                    swaprows!(Coefficients,cur,j)
                end
                found=true
                cur+=1 # The vector at cur has a nonzero i-th component. So for the next component we start the search from cur+1.
                rank+=1 
                break
            end
        end

        # Break if the matrix is already full-ranked.
        if rank==n
            break
        end

        # Eliminate the i-th component of the rest vectors.
        if(found)  
            for j in cur:n
                if Vecs[j,i]==true
                    Vecs[j,:]=Vecs[j,:] .⊻ Vecs[cur-1,:]
                    Coefficients[j,:]=Coefficients[j,:] .⊻ Coefficients[cur-1,:]
                end
            end
        end
    end

    return rank
end