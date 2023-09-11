include("utils.jl")
include("dynamic_update.jl")

function rough_boundary_X(Width::Int,Depth::Int)
    # The lattice is put on a cylinder, i.e. a periodic boundary condition in the x direction and an open boundary condition in the y direction.
    # The lower boundary is a smooth boundary and the upper boundary is a rough boundary.
        # The lower boundary has y coordinate 1. The upper boundary 

    # The function calculates the boundary stabilizers after measuring all single-site X in the bulk.
    
    # Initially, the state is a cluster state.

    lattice_size=Width*Depth
    Initial_SG=zeros(Bool,lattice_size,2*lattice_size) # One stabilizer each site and we need two elements for each site.
    measurements=zeros(Bool,Width*(Depth-1),2*lattice_size) # The lowest row isn't measured.

    # Initialize the stabilizer group.

    for y in 1:Depth
        for x in 1:Width
            n=Width*(y-1)+x
            Initial_SG[n,2*n-1]=true # Pauli-X
            
            if(y>1)
                n_down=Width*(y-2)+x
                Initial_SG[n,2*n_down]=true # Pauli-Z
            end

            if(y<Depth)
                n_up=Width*y+x
                Initial_SG[n,2*n_up]=true # Pauli-Z

                if(x>1) # The upper boundary is rough, so at the upmost there is not left and right neighbors.
                    n_left=n-1
                    Initial_SG[n,2*n_left]=true # Pauli-Z
                else
                    n_left=Width*(y-1)+Width #leftmost
                    Initial_SG[n,2*n_left]=true # Pauli-Z
                end

                if(x<Width)
                    n_right=Width*(y-1)+x+1
                    Initial_SG[n,2*n_right]=true # Pauli-Z
                else
                    n_right=Width*(y-1)+1 # rightmost
                    Initial_SG[n,2*n_right]=true # Pauli-Z
                end
            end

        end
    end


    # Initialize Measurements
    for y in 2:Depth
        for x in 1:Width
            n=Width*(y-1)+x
            measurements[n-Width,2*n-1]=true # Pauli-X # Set the first site on the second row as 1
        end
    end

    # Update the stabilizer group
    println("Initial_SG")
    println(size(Initial_SG))
    println(ISG_to_string(Initial_SG,1,Width,Depth,true))
    new_SG=update(lattice_size,Initial_SG,measurements)
    println("New_SG")
    println(size(new_SG))
    
    println(ISG_to_string(new_SG,1,Width,Depth,true))

end

rough_boundary_X(5,8)