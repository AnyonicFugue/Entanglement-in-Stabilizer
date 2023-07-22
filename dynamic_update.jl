# Problem: Can we easily write out the stabilizer generators for the ISG? Or must we store the whole group?

# Brute-Force Strategy: For each element of the OldISG, check if it commutes with all elements of Measurements.
# If it does, add it and its products with elements of Measurements to NewISG. Note that two products might produce the same element.

# Elegant Strategy: Keep a list of generators for the ISG. Update it with measurements.
# Idea: Linear Algebra on $Z_2$. Use 0 to represent commute and 1 to represent anti-commute.

function calc_comm_matrix!(LatticeSize::Int64,OldISG::Array{Bool,2},Measurements::Array{Bool,2},CommMatrix::Matrix{Bool})
    # The function calculates the commutation matrix of OldISG with Measurements and stores it in CommMatrix.
    # The function assumes that CommMatrix is initialized to all zero.

    comm_vec=zeros(Bool,LatticeSize)
    tmp_vec=zeros(Bool,LatticeSize)

    
    for i in 1:length(OldISG[:,1])
        for j in 1:length(Measurements[:,1])
            for n in 1:LatticeSize

                tmp_vec[n]=iszero((OldISG[i,2*n-1]+OldISG[i,2*n])*(Measurements[j,2*n-1]+Measurements[j,2*n]))|iszero((2*OldISG[i,2*n-1]+OldISG[i,2*n])-(2*Measurements[j,2*n-1]+Measurements[j,2*n]))
                # (OldISG[2*n-1]+OldISG[2*n])*(Measurements[2*n-1]+Measurements[2*n]) tests whether one of the two operators is identity.
                # (2*OldISG[2*n-1]+OldISG[2*n])-(2*Measurements[2*n-1]+Measurements[2*n]) tests whether the two operators are equal.
                # If one of them holds, then the two operators commute on the n-th site. Mark it as 1.

                comm_vec=.! tmp_vec # The commutation vector. 0 means commute and 1 means anti-commute.
            end
            
            if sum(comm_vec)%2==0
                CommMatrix[i,j]=false  # The two operators commute, set as 0.
            else
                CommMatrix[i,j]=true # The two operators anti-commute, set as 1.
            end
        end
    end

end

function update!(LatticeSize::Int64,OldISG::Array{Int32,2},NewISG::Array{Bool,2},Measurements::Array{Bool,2})
    # The data structure: First index is the index of the stablizer. The second index is the dimension in the vector space.
        # Note that for n^th d.o.f. (starting from 1), 2*n-1 is for Pauli_X, 2*n is for Pauli_Z.


    # OldISG records the generators.
    # The function overwrites the new generators into NewISG, thus it must be a different variable from OldISG.

    # Step 1. Find the surviving element of the OldISG under measurements.
    # This is achieved by finding the zero space of by Gaussian elimination. i.e. (rank+1)-th to n-th vectors constitute a basis of the zero space.

    #  Step 1.1. Find the commutation relations of OldISG with Measurements and record into an array.
    n_old=length(OldISG[:,1])
    n_measure=length(Measurements[:,1])

    commute_matrix=zeros(Bool,(n_old,n_measure)) # The commutation matrix. At the entry (i,j), 0 means the i-th generator commutes with the j-th measurement and 1 for anti-commute.
    calc_comm_matrix!(LatticeSize,OldISG,Measurements,commute_matrix)

    coefficients=zeros(Bool,(n_old,n_old)) # The coefficients of the echoleon vectors. The entry (i,j) is the coefficient of the j-th generator in the i-th echoleon vector.
    rank=0
    rank=gaussian_elimination!(commute_matrix,coeff) # The function overwrites the input vectors with the echoleon vectors and records the coefficients in the second vector.
    # 1-rank are nonzero (i.e. non-commuting) vectors. rank+1 to n are zero vectors, i.e. those commuting with all measurements.




    
end