function calc_entropy(StabGenerators::Vector{Tuple{Vector{Int32},Vector{Int8}}},region::Vector{Int32})
    # region is a vector of the d.o.f. in the region of interest.
    
    println("calc_entropy")

    # Strategy: Calculate (rank G-rank G_A-rank G_B)/2.
    
    rank_GA=0
    rank_GB=0

    for i in 1:length(StabGenerators)
        if isdisjoint(StabGenerators[i][1],region)
            rank_GB+=1
        end
        if issubset(StabGenerators[i][1],region)
            rank_GA+=1
        end
    end

    println("rank_GA:",rank_GA,", rank_GB:",rank_GB)
    println("entropy:",(length(StabGenerators)-rank_GA-rank_GB)/2)
end