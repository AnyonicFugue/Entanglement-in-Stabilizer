include("utils.jl")
include("dynamic_update.jl")
import PyPlot

function plot_stabilizer_generators(Width::Int,Depth::Int,Stabilizers::Array{Bool,2},PlotNumber::Int)

    # Use PyPlot library to visualize the stabilizer generators on a square lattice and show a graph for each stabilizer
    # A vertex is green if there is a Pauli-X, blue if a Pauli-Z, and red if both.
    
    lattice_size=Width*Depth
    x=zeros(Int,lattice_size)
    y=zeros(Int,lattice_size)
    color=Array{String}(undef,lattice_size)

    for n in 1:lattice_size
        x[n]=mod(n-1,Width)+1
        y[n]=div(n-1,Width)+1
    end

    
    for j in 1:min(PlotNumber,size(Stabilizers,1))
        # Plot the stabilizers one by one.

        for n in 1:lattice_size
            if(Stabilizers[j,2*n-1] && Stabilizers[j,2*n])
                color[n]="red"
            elseif(Stabilizers[j,2*n-1])
                color[n]="green"
            elseif(Stabilizers[j,2*n])
                color[n]="blue"
            else
                color[n]="grey"
            end
        end
    
        PyPlot.figure(figsize=(Width,Depth))
        PyPlot.scatter(x,y,c=color,s=100)
        PyPlot.xlim(0,Width+1)
        PyPlot.ylim(0,Depth+1)
        PyPlot.xticks([])
        PyPlot.yticks([])
        PyPlot.show()
    end


end

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
    # println(ISG_to_string(Initial_SG,1,Width,Depth,true))
    new_SG=update(lattice_size,Initial_SG,measurements)
    println("New_SG")
    println(size(new_SG))
    
    println(ISG_to_string(new_SG,1,Width,Depth,true))

    plot_stabilizer_generators(Width,Depth,new_SG,4)
end


function smooth_boundary_Y(Width::Int,Depth::Int)
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
            end

            if(x>1) # The upper boundary is smooth, so every point has left and right neighbors.
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


    # Initialize Measurements
    for y in 2:Depth
        for x in 1:Width
            n=Width*(y-1)+x

            # Y-measurements
            measurements[n-Width,2*n-1]=true
            measurements[n-Width,2*n]=true 
        end
    end

    # Update the stabilizer group
    println("Initial_SG")
    println(size(Initial_SG))
    # println(ISG_to_string(Initial_SG,1,Width,Depth,true))
    # plot_stabilizer_generators(Width,Depth,Initial_SG)

    new_SG=update(lattice_size,Initial_SG,measurements)
    println("New_SG")
    println(size(new_SG))
    println(ISG_to_string(new_SG,1,Width,Depth,true))
    plot_stabilizer_generators(Width,Depth,new_SG,4)

end

smooth_boundary_Y(6,8)
smooth_boundary_Y(6,10)
smooth_boundary_Y(6,12)




