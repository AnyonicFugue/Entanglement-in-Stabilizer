include("calc_plot_and_fit.jl")
include("dynamic_update.jl")


function main()
    lattice_size::Int32=100
    stab_generators=Vector{Tuple{Vector{Int32},Vector{Int8}}}()
    calc_entropy(lattice_size,stab_generators)
end

#=
Step 1: Calculate the entropy for specific codes, without considering the bounary conditions.
=#

