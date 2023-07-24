include("calc_plot_and_fit.jl")
include("dynamic_update.jl")

import CurveFit

function toric_code_static()
    # Note that the d.o.f. are on the edges, not on the vertices.

    l=8
    lattice_size::Int32=2*l*l

    stab_generators=Vector{Tuple{Vector{Int32},Vector{Int8}}}() # Each stabilizer is stored as a tuple. The first component are the d.o.f. it acts on and the second are the Pauli matrices it acts as.

    # Vertex Stabilizers
    for m in 1:l*l # m is the index of the vertex
        y=m%l
        if y==0
            y=l
        end
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
        y=m%l
        if y==0
            y=l
        end
            
        x=Int32((m-y)/l)+1

        y_down=y%l+1
        x_right=x%l+1
        
        m_down=(x-1)*l+y_down
        m_right=(x_right-1)*l+y

        #=
        Each vertex stabilizer acts on 4 edges: 
        Right edge, 2*m-1; 
        Down edge, 2*m; 
        Right edge of the vertex on the downside. The coordinate of the vertex on the downside is (x,y%l+1).
        Down edge of the vertex on the right. The coordinate of the vertex on the right is (x%l+1,y).
    
        1 denotes Pauli_X.
        =#

        push!(stab_generators,([2*m-1,2*m,2*m_down-1,2*m_right],[1,1,1,1]))
    end

    println("toric_code")
    println(length(stab_generators))
    # println(stab_generators)

    # Select rectangular regions with increasing sizes and plot entropy vs region size.

    start=Int32(floor(l/5))

    volume_arr=zeros(Int32,Int32(l/2)-start+1)
    area_arr=zeros(Int32,Int32(l/2)-start+1)
    entropy_arr=zeros(Float32,Int32(l/2)-start+1)

    for s in range(start,Int32(l/2))
        region=Vector{Int32}()
        area_arr[s-start+1]=4*s

        for i in 1:s
            for j in 1:s
                m=(i-1)*l+j # The vertex
                if(i<s) # Not on the rightmost
                    push!(region,2*m-1)
                    volume_arr[s-start+1]+=1
                end

                if(j<s) # Not on the downmost
                    push!(region,2*m)
                    volume_arr[s-start+1]+=1
                end
            end
        end
        entropy_arr[s-start+1]=calc_entropy(stab_generators,region)
    end

    println("area_arr:",area_arr)
    plot_and_fit(entropy_arr,volume_arr,area_arr)
    # plot_and_fit_native(entropy_arr,volume_arr,area_arr)
end

function single_triangle()
    LatticeSize=3

    ISG=zeros(Bool,(3,6))
    measurement_rounds=zeros(Bool,(3,1,6)) # We perform 3 rounds of measurement. Each round has 1 measurement, which is the product of XX, YY, ZZ on 3 edges respectively.
    measurement_rounds[1,1,:]=[1,0,1,0,0,0] #X1X2
    measurement_rounds[2,1,:]=[0,0,0,1,0,1] #Z2Z3
    measurement_rounds[3,1,:]=[1,1,0,0,1,1] #Y1Y3

    for i in 1:3
        ISG=update(LatticeSize,ISG,measurement_rounds[i,:,:])
        println("ISG:",ISG)
        println("ISG_to_string:",ISG_to_string(ISG))
    end
end

function ISG_to_string(ISG::Array{Bool,2})
    # Convert the ISG to a string.
    # The first index is the index of the stabilizer. The second index is the index of the d.o.f.
    # The d.o.f. are ordered as X1,Z1,X2,Z2,...,Xn,Zn.

    n=length(ISG[:,1])
    lattice_size=Int(length(ISG[1,:])/2)
    ISG_string=Vector{String}(undef,0)
    
    for i in 1:n
        str=String[]
       
        for k in 1:lattice_size
            if (ISG[i,2*k-1])&&(!ISG[i,2*k])
                push!(str,"X")
                push!(str,string(k))
            elseif (!ISG[i,2*k-1])&&(ISG[i,2*k])
                push!(str,"Z")
                push!(str,string(k))
            elseif (ISG[i,2*k-1])&&(ISG[i,2*k])
                push!(str,"Y")
                push!(str,string(k))
            end  
        end
       
        push!(ISG_string,join(str))
    end
    return ISG_string
end

function single_hexagon_3rounds(cycle)
    # Perform 3 rounds of measurement on a hexagon.
    lattice_size=6
    round=3

    ISG=zeros(Bool,(round,2*lattice_size))
    measurement_rounds=zeros(Bool,(round,2,2*lattice_size)) # We perform 3 rounds of measurement. Each round has 1 measurement, which is the product of XX, YY, ZZ on 3 edges respectively.
    # Round 1: Measure Z1Z6 and Z3Z4
    measurement_rounds[1,1,2]=1
    measurement_rounds[1,1,12]=1
    measurement_rounds[1,2,6]=1
    measurement_rounds[1,2,8]=1

    # Round 2: Measure X1X2 and X4X5
    measurement_rounds[2,1,1]=1
    measurement_rounds[2,1,3]=1
    measurement_rounds[2,2,7]=1
    measurement_rounds[2,2,9]=1

    # Round 3: Measure Y2Y3 and Y5Y6
    measurement_rounds[3,1,3]=1
    measurement_rounds[3,1,4]=1
    measurement_rounds[3,1,5]=1
    measurement_rounds[3,1,6]=1
    
    measurement_rounds[3,2,9]=1
    measurement_rounds[3,2,10]=1
    measurement_rounds[3,2,11]=1
    measurement_rounds[3,2,12]=1

    for i in 1:cycle*round
        ISG=update(lattice_size,ISG,measurement_rounds[(i-1)%round+1,:,:])
        println("ISG:",ISG)
        println("ISG_to_string:",ISG_to_string(ISG))
    end
end

function single_hexagon_2rounds(cycle)
    # Perform 2 rounds of measurement on a hexagon.
    lattice_size=6
    round=2

    ISG=zeros(Bool,(round,2*lattice_size))
    measurement_rounds=zeros(Bool,(round,3,2*lattice_size)) # We perform 2 rounds of measurement. Each round has 1 measurement, which is the product of XX, YY, ZZ on 3 edges respectively.
 
    # Round 1: Measure X1X2, Z3Z4, Y5Y6
    measurement_rounds[1,1,1]=1
    measurement_rounds[1,1,3]=1

    measurement_rounds[1,2,6]=1
    measurement_rounds[1,2,8]=1

    measurement_rounds[1,3,9]=1
    measurement_rounds[1,3,10]=1
    measurement_rounds[1,3,11]=1
    measurement_rounds[1,3,12]=1

    # Round 2: Measure Z1Z6, X4X5, Y2Y3
    measurement_rounds[2,1,2]=1
    measurement_rounds[2,1,12]=1

    measurement_rounds[2,2,7]=1
    measurement_rounds[2,2,9]=1

    measurement_rounds[2,3,3]=1
    measurement_rounds[2,3,4]=1
    measurement_rounds[2,3,5]=1
    measurement_rounds[2,3,6]=1

    for i in 1:cycle*round
        ISG=update(lattice_size,ISG,measurement_rounds[(i-1)%round+1,:,:])
        println("ISG:",ISG)
        println("ISG_to_string:",ISG_to_string(ISG))
    end
end

function test_parallel_comm_matrix(LatticeSize::Int64,ISGSize::Int64)
    # Assume measurements are of the same order of ISG.

    empty=zeros(Bool,(2,2))
    calc_comm_matrix!(1,empty,empty,empty,false)
    calc_comm_matrix!(1,empty,empty,empty,true) # Precompile the function.

    OldISG=rand(Bool,(ISGSize,2*LatticeSize))
    Measurements=rand(Bool,(ISGSize,2*LatticeSize))
    CommMatrix=zeros(Bool,(ISGSize,ISGSize))

    println("Serial Version")
    @time calc_comm_matrix!(LatticeSize,OldISG,Measurements,CommMatrix,false)

    OldISG=rand(Bool,(ISGSize,2*LatticeSize))
    Measurements=rand(Bool,(ISGSize,2*LatticeSize))
    CommMatrix=zeros(Bool,(ISGSize,ISGSize))

    println("Parallel Version")
    @time calc_comm_matrix!(LatticeSize,OldISG,Measurements,CommMatrix,true)

end
# single_hexagon_3rounds(3)
# single_triangle()
#single_hexagon_2rounds(3)

function snake_growing_lattice(Cycle::Int,LatticeSize::Int)
    round=8
end

stab_generators=zeros(Bool,(1,1))
sampling_square(4,3,1,stab_generators)