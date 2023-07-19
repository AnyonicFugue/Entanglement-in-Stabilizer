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

function toric_code()
    # Note that the d.o.f. are on the edges, not on the vertices.

    l=16
    lattice_size::Int32=2*l*l

    stab_generators=Vector{Tuple{Vector{Int32},Vector{Int8}}}()

    # Vertex Stabilizers
    for m in 1:l*l # m is the index of the vertex
        y=m%l
        x=Int32((m-y)/l)+1

        x_left=x-1
        if x_left==0
            x_left=l
        end
        m_left=(x_left-1)*l+y

        y_up=y-1
        if y_up==0
            y_up=l
        end
        m_up=(x-1)*l+y_up

        push!(stab_generators,([2*m-1,2*m,2*m_left-1,2*m_up],[3,3,3,3]))
        #=
        Each vertex stabilizer acts on 4 edges: 
        Right edge, 2*m-1; 
        Down edge, 2*m; 
        Left edge (or right edge of vertex on the left). The coordinate of the vertex on the left is (x-1 (l if equals zero),y), 
        Up edge (or down edge of vertex on the top), The coordinate of the vertex on the top is (x,y-1 (l if equals zero)).
        3 denotes Pauli_Z.
        =#
    
    end
    # Plaquette Stabilizers
    for m in 1:l*l # Consider the plaquette at the downright of the vertex in question.
    #=
    Each vertex stabilizer acts on 4 edges: 
    Right edge, 2*m-1; 
    Down edge, 2*m; 
    Right edge of the vertex on the downside. The coordinate of the vertex on the downside is (x,y%l+1).
    Down edge of the vertex on the right. The coordinate of the vertex on the right is (x%l+1,y).
   
    1 denotes Pauli_X.
    =#

        y=m%l
        x=Int32((m-y)/l)+1

        y_down=y%l+1
        x_right=x%l+1
        
        m_down=(x-1)*l+y_down
        m_right=(x_right-1)*l+y

        push!(stab_generators,([2*m-1,2*m,2*m_down-1,2*m_right],[1,1,1,1]))
    end

        
    end


    calc_entropy(lattice_size,stab_generators)
end


# A function to visualize the lattice and the stabilizers.



main()