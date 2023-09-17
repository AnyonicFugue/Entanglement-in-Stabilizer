include("utils.jl")
include("dynamic_update.jl")
import PyPlot


plot_round=3

function plot_stabilizer_generators(Width::Int,Depth::Int,Stabilizers::Array{Bool,2},PlotNumber::Int,Title::String)

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


    if(false) # if(BoundaryOnly) ; this was initially to show the stabilizers only on the boundary.
        cur=0

        for n in Width+1:lattice_size
            color[n]="grey" # Set all bulk colors to grey
        end

        for j in 1:size(Stabilizers,1)
            # Plot the stabilizers one by one. 
            if(cur==PlotNumber)
                break
            end

            for n in 1:Width # Plot the boundary stabilizers
                if(Stabilizers[j,2*n-1] && Stabilizers[j,2*n])
                    color[n]="red"
                elseif(Stabilizers[j,2*n-1])
                    color[n]="gold" # Pauli-X
                elseif(Stabilizers[j,2*n])
                    color[n]="blue" # Pauli-Z
                else
                    color[n]="grey"
                end
            end

            
            if(sum(Stabilizers[j,2*(Width+1)-1:2*lattice_size])==0) # Acts trivially in the boundary

                PyPlot.figure(figsize=(Width,Depth))
                PyPlot.scatter(x,y,c=color,s=100)
                PyPlot.xlim(0,Width+1)
                PyPlot.ylim(0,Depth+1)
                PyPlot.xticks([])
                PyPlot.yticks([])
                PyPlot.title(Title)
                PyPlot.show()

                cur=cur+1
            end
        end
    else
        for j in 1:min(PlotNumber,size(Stabilizers,1))
            # Plot the stabilizers one by one.
    
            for n in 1:lattice_size
                if(Stabilizers[j,2*n-1] && Stabilizers[j,2*n])
                    color[n]="red"
                elseif(Stabilizers[j,2*n-1])
                    color[n]="gold" # Pauli-X
                elseif(Stabilizers[j,2*n])
                    color[n]="blue" # Pauli-Z
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
            PyPlot.title(Title)
            PyPlot.show()
        end

    end

end


function Calculate_boundary_stabilizers(Width::Int,Depth::Int;BoundaryOnly::Bool=false,SeeFullLattice::Bool=true,Plot::Bool=false,FindLocalGenerators::Bool=true,IsXMeasurement::Bool,IsSmoothBoundary::Bool,PrintISG::Bool=false)
    # The lattice is put on a cylinder, i.e. a periodic boundary condition in the x direction and an open boundary condition in the y direction.
    # The lower boundary is a smooth boundary and the upper boundary is a rough boundary.
        # The lower boundary has y coordinate 1. The upper boundary has y coordinate Depth.

    # If BoundaryOnly, then the returning stabilizer group would be truncated to leave only the boundary part.
    # If SeeFullLattice, then the stabilizer group would be truncated before Find_local_generators in order to show the bulk behavior of the stabilizers.

    lattice_size=Width*Depth
    Initial_SG=zeros(Bool,lattice_size,2*lattice_size) # One stabilizer each site and we need two elements for each site.
    measurements=zeros(Bool,Width*(Depth-1),2*lattice_size) # The lowest row isn't measured.

    # Initialize the stabilizer group.
    if(IsSmoothBoundary)
        initialize_smooth_boundary!(Width,Depth,Initial_SG)
    else
        initialize_rough_boundary!(Width,Depth,Initial_SG)
    end

    # Initialize Measurements
    for y in 2:Depth
        for x in 1:Width
            n=Width*(y-1)+x

            # X-measurements
            measurements[n-Width,2*n-1]=true

            if(!(IsXMeasurement))
                measurements[n-Width,2*n]=true # Y measurements
            end
        end
    end

    # Update the stabilizer group
    # println("Initial_SG")
    println(size(Initial_SG))
    # println(ISG_to_string(Initial_SG,1,Width,Depth,true))
    # plot_stabilizer_generators(Width,Depth,Initial_SG)

    new_SG=update(lattice_size,Initial_SG,measurements)
    # println("New_SG")
    

    if(SeeFullLattice)
        new_SG=new_SG[1:Width,:] # Truncate the stabilizers thus we can see the structure on the full lattice
    end


    if(FindLocalGenerators)
        Find_local_generators!(new_SG)
    end


    if(BoundaryOnly)
        new_SG=new_SG[1:Width,:] # Truncate the stabilizer to leave only the boundary
    end

    if(PrintISG)
        println(ISG_to_string(new_SG,1,Width,Depth))
    end

    title="Smooth Boundary X,Width="*string(Width)*",Depth="*string(Depth)

    if(Plot)
        plot_stabilizer_generators(Width,Depth,new_SG,plot_round,title)
    end

    if(BoundaryOnly)
        return new_SG[:,1:2*Width]
    else
        return new_SG
    end
end


function initialize_smooth_boundary!(Width::Int,Depth::Int,StabArray::Array{Bool,2})
    for y in 1:Depth
        for x in 1:Width
            n=Width*(y-1)+x
            StabArray[n,2*n-1]=true # Pauli-X
            
            if(y>1)
                n_down=Width*(y-2)+x
                StabArray[n,2*n_down]=true # Pauli-Z
            end

            if(y<Depth)
                n_up=Width*y+x
                StabArray[n,2*n_up]=true # Pauli-Z
            end

            if(x>1) # The upper boundary is smooth, so every point has left and right neighbors.
                n_left=n-1
                StabArray[n,2*n_left]=true # Pauli-Z
            else
                n_left=Width*(y-1)+Width #leftmost
                StabArray[n,2*n_left]=true # Pauli-Z
            end

            if(x<Width)
                n_right=Width*(y-1)+x+1
                StabArray[n,2*n_right]=true # Pauli-Z
            else
                n_right=Width*(y-1)+1 # rightmost
                StabArray[n,2*n_right]=true # Pauli-Z
            end

        end
    end
end


function initialize_rough_boundary!(Width::Int,Depth::Int,StabArray::Array{Bool,2})
    for y in 1:Depth
        for x in 1:Width
            n=Width*(y-1)+x
            StabArray[n,2*n-1]=true # Pauli-X
            
            if(y>1)
                n_down=Width*(y-2)+x
                StabArray[n,2*n_down]=true # Pauli-Z
            end

            if(y<Depth)
                n_up=Width*y+x
                StabArray[n,2*n_up]=true # Pauli-Z

                if(x>1) # The upper boundary is rough, so no left or right neighbors.
                    n_left=n-1
                    StabArray[n,2*n_left]=true # Pauli-Z
                else
                    n_left=Width*(y-1)+Width #leftmost
                    StabArray[n,2*n_left]=true # Pauli-Z
                end
    
                if(x<Width)
                    n_right=Width*(y-1)+x+1
                    StabArray[n,2*n_right]=true # Pauli-Z
                else
                    n_right=Width*(y-1)+1 # rightmost
                    StabArray[n,2*n_right]=true # Pauli-Z
                end
            end



        end
    end
end


# smooth_boundary_Y(20,48,false)

function sample_EE(Width::Int,StartDepth::Int,EndDepth::Int,IsX::Bool,IsSmooth::Bool,PrintInfo::Bool=false)
    # The function calculates the entanglement entropy of the boundary stabilizers.
    
    # The entanglement entropy is calculated for the region between StartDepth and EndDepth.
    # The function returns a vector of entanglement entropies for each depth.
    EE_arr=zeros(EndDepth-StartDepth+1)

    Threads.@threads for depth in StartDepth:EndDepth
        # Calculate the entanglement entropy for the region between StartDepth and EndDepth.
        BoundaryISG=Calculate_boundary_stabilizers(Width,depth,BoundaryOnly=true,IsXMeasurement=IsX,IsSmoothBoundary=IsSmooth)

        if(PrintInfo)
            coeff=zeros(Bool,Width,Width)
            print("rank=")
            println(gaussian_elimination!(BoundaryISG,coeff,true))
            println(ISG_to_string(BoundaryISG,1,Width,1))
        end

        EE_arr[depth-StartDepth+1]=Evaluate_EE(BoundaryISG,Int(floor(Width/2)))
    end

    # Plot entanglement entropy vs depth by Pyplot package
    PyPlot.plot(range(start=StartDepth,stop=EndDepth,step=1),EE_arr)
    PyPlot.xlabel("Depth")
    PyPlot.ylabel("Entanglement Entropy")
    PyPlot.title("Entanglement Entropy vs Depth, Width="*string(Width))
    PyPlot.show()


    return EE_arr
end


println(sample_EE(24,16,48,false,true))
