function update!(OldISG::Vector{Tuple{Vector{Int32},Vector{Int8}}},NewISG::Vector{Tuple{Vector{Int32},Vector{Int8}}},Measurements::Vector{Tuple{Vector{Int32},Vector{Int8}}})
    # Problem: Can we easily write out the stabilizer generators for the ISG? Or must we store the whole group?

    # Brute-Force Strategy: For each element of the OldISG, check if it commutes with all elements of Measurements.
    # If it does, add it and its products with elements of Measurements to NewISG. Note that two products might produce the same element.

    # Elegant Strategy: Keep a list of generators for the ISG. Update it with measurements.
    # Idea: Linear Algebra on $Z_2$. Use 0 to represent commute and 1 to represent anti-commute.
end