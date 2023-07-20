# Problem: Can we easily write out the stabilizer generators for the ISG? Or must we store the whole group?

# Brute-Force Strategy: For each element of the OldISG, check if it commutes with all elements of Measurements.
# If it does, add it and its products with elements of Measurements to NewISG. Note that two products might produce the same element.

# Elegant Strategy: Keep a list of generators for the ISG. Update it with measurements.
# Idea: Linear Algebra on $Z_2$. Use 0 to represent commute and 1 to represent anti-commute.

function update!(LatticeSize::Int32,OldISG::Array{Int8,2},NewISG::Array{Int8,2},Measurements::Array{Int8,2})
    # The data structure: First index is the index of the stablizer. The second index is the dimension in the vector space.
        # Note that for n^th d.o.f. (starting from 1), 2*n-1 is for Pauli_X on the d.o.f. while 2*n is for Pauli_Z on it.

    # OldISG records the generators.
    # The function overwrites the new generators into NewISG, thus it must be a different variable from OldISG.

    
end