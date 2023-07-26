include("calc_plot_and_fit.jl")
include("dynamic_update.jl")

import CurveFit
import Profile
import StatProfilerHTML

function toric_code_static(l::Int64,Parallel::Bool)
    # Note that the d.o.f. are on the edges, not on the vertices.
    # The stabilizers are ordered as X1,Z1,X2,Z2,...,Xn,Zn.

    lattice_size::Int64=2*l*l

    # For a lattice of length l, the size is 2*l*l, since each vertex corresponds to two edges.
    # There is one stabilizer on each vertex and each plaquette, so the total number of stabilizers is 2*l*l.

    stab_generators=zeros(Bool,(2*l*l,2*lattice_size)) # For each stabilizer the first index is the index of the stabilizer. The second index is the index of the d.o.f.

    # Vertex Stabilizers
    for m in 1:l*l # m is the index of the vertex
        y=m%l
        if y==0
            y=l
        end
        x=Int64((m-y)/l)+1

        # Relation between the coordinate and the index of the vertex: m=(x-1)*l+y
        # The edge on the right of the vertex m is 2*m-1. The edge on the down of the vertex is 2*m.

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

        stab_generators[m,2*(2*m-1)]=true
        stab_generators[m,2*(2*m)]=true
        stab_generators[m,2*(2*m_left-1)]=true
        stab_generators[m,2*(2*m_up)]=true
        #=
        Each vertex stabilizer acts on 4 edges as Pauli_Z:
        Right edge, 2*m-1; 
        Down edge, 2*m; 
        Left edge (or right edge of vertex on the left). The coordinate of the vertex on the left is (x-1 (l if equals zero),y), 
        Up edge (or down edge of vertex on the top), The coordinate of the vertex on the top is (x,y-1 (l if equals zero)).
        2*n denotes Pauli_Z on the n-th d.o.f.
        =#
    end

    # Plaquette Stabilizers
    for m in 1:l*l # Consider the plaquette at the downright of the vertex in question.
        y=m%l
        if y==0
            y=l
        end
            
        x=Int64((m-y)/l)+1

        y_down=y%l+1
        x_right=x%l+1
        
        m_down=(x-1)*l+y_down
        m_right=(x_right-1)*l+y

        stab_generators[l*l+m,2*(2*m-1)-1]=true
        stab_generators[l*l+m,2*(2*m)-1]=true
        stab_generators[l*l+m,2*(2*m_down-1)-1]=true
        stab_generators[l*l+m,2*(2*m_right)-1]=true
        #=
        Each vertex stabilizer acts on 4 edges: 
        Right edge, 2*m-1; 
        Down edge, 2*m; 
        Right edge of the vertex on the downside. The coordinate of the vertex on the downside is (x,y%l+1).
        Down edge of the vertex on the right. The coordinate of the vertex on the right is (x%l+1,y).
    
        2*n-1 denotes Pauli_X on the n-th d.o.f.
        =#
    end

    println("toric_code")
    println(length(stab_generators[:,1]))
    # println(stab_generators)

    # Select rectangular regions with increasing sizes and plot entropy vs region size.

    start=Int64(floor(l/5))

    volume_arr=zeros(Int64,Int64(l/2)-start+1)
    area_arr=zeros(Int64,Int64(l/2)-start+1)
    entropy_arr=zeros(Float32,Int64(l/2)-start+1)

    region=zeros(Int64,(Int64(l/2)-start+1,2*Int64(l/2)^2))
    cur=zeros(Int64,Int64(l/2)-start+1)

    for s in range(max(start,2),Int64(l/2))
        cur[s-start+1]=1
        area_arr[s-start+1]=4*s

        for i in 1:s
            for j in 1:s
                m=(i-1)*l+j # The vertex
                if(i<s) # Not on the rightmost
                    region[s-start+1,cur[s-start+1]]=2*m-1
                    cur[s-start+1]+=1
                    volume_arr[s-start+1]+=1
                end

                if(j<s) # Not on the downmost
                    region[s-start+1,cur[s-start+1]]=2*m
                    volume_arr[s-start+1]+=1
                    cur[s-start+1]+=1
                end
            end
        end

        # println(region[s-start+1,:]) # Only for debugging
        # entropy_arr[s-start+1]=calc_entropy(stab_generators,region[s-start+1,1:cur-1],Parallel)
    end

    for s in range(max(start,2),Int64(l/2)) # The line is separated for testing the performance of the function.
        entropy_arr[s-start+1]=calc_entropy(stab_generators,view(region,s-start+1,:),Parallel,cur[s-start+1]-1)
    end

    println("area_arr:",area_arr)
    println("entropy_arr:",entropy_arr)
    ## plot_and_fit(entropy_arr,volume_arr,area_arr) 
        # Disabled when timing.
    fit(entropy_arr,volume_arr,area_arr)
    # plot_and_fit_native(entropy_arr,volume_arr,area_arr)
end

function toric_code_static_multisampling(l::Int64,Parallel::Bool)
    # Note that the d.o.f. are on the edges, not on the vertices.
    # The stabilizers are ordered as X1,Z1,X2,Z2,...,Xn,Zn.

    lattice_size::Int64=2*l*l

    # For a lattice of length l, the size is 2*l*l, since each vertex corresponds to two edges.
    # There is one stabilizer on each vertex and each plaquette, so the total number of stabilizers is 2*l*l.

    stab_generators=zeros(Bool,(2*l*l,2*lattice_size)) # For each stabilizer the first index is the index of the stabilizer. The second index is the index of the d.o.f.

    # Vertex Stabilizers
    for m in 1:l*l # m is the index of the vertex
        y=m%l
        if y==0
            y=l
        end
        x=Int64((m-y)/l)+1

        # Relation between the coordinate and the index of the vertex: m=(x-1)*l+y
        # The edge on the right of the vertex m is 2*m-1. The edge on the down of the vertex is 2*m.

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

        stab_generators[m,2*(2*m-1)]=true
        stab_generators[m,2*(2*m)]=true
        stab_generators[m,2*(2*m_left-1)]=true
        stab_generators[m,2*(2*m_up)]=true
        #=
        Each vertex stabilizer acts on 4 edges as Pauli_Z:
        Right edge, 2*m-1; 
        Down edge, 2*m; 
        Left edge (or right edge of vertex on the left). The coordinate of the vertex on the left is (x-1 (l if equals zero),y), 
        Up edge (or down edge of vertex on the top), The coordinate of the vertex on the top is (x,y-1 (l if equals zero)).
        2*n denotes Pauli_Z on the n-th d.o.f.
        =#
    end

    # Plaquette Stabilizers
    for m in 1:l*l # Consider the plaquette at the downright of the vertex in question.
        y=m%l
        if y==0
            y=l
        end
            
        x=Int64((m-y)/l)+1

        y_down=y%l+1
        x_right=x%l+1
        
        m_down=(x-1)*l+y_down
        m_right=(x_right-1)*l+y

        stab_generators[l*l+m,2*(2*m-1)-1]=true
        stab_generators[l*l+m,2*(2*m)-1]=true
        stab_generators[l*l+m,2*(2*m_down-1)-1]=true
        stab_generators[l*l+m,2*(2*m_right)-1]=true
        #=
        Each vertex stabilizer acts on 4 edges: 
        Right edge, 2*m-1; 
        Down edge, 2*m; 
        Right edge of the vertex on the downside. The coordinate of the vertex on the downside is (x,y%l+1).
        Down edge of the vertex on the right. The coordinate of the vertex on the right is (x%l+1,y).
    
        2*n-1 denotes Pauli_X on the n-th d.o.f.
        =#
    end

    println("toric_code")
    println(length(stab_generators[:,1]))
    # println(stab_generators)

    # Select rectangular regions with increasing sizes and plot entropy vs region size.

    start=Int64(floor(l/5))

    volume_arr=zeros(Int64,Int64(l/2)-start+1)
    area_arr=zeros(Int64,Int64(l/2)-start+1)
    entropy_arr=zeros(Float32,Int64(l/2)-start+1)


    for s in range(max(start,2),Int64(l/2))
        
        area_arr[s-start+1]=4*s
        volume_arr[s-start+1]=2*s*s
        
        entropy_arr[s-start+1]=sample_squares(l,s,8,stab_generators)
    end

    # println("area_arr:",area_arr)
    ## plot_and_fit(entropy_arr,volume_arr,area_arr) 
        # Disabled when timing.
    fit(entropy_arr,volume_arr,area_arr)
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

function single_hexagon_3rounds(cycle,Comm_Parallel::Bool,Elimination_Parallel::Bool)
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
        ISG=update(lattice_size,ISG,measurement_rounds[(i-1)%round+1,:,:],Comm_Parallel,Elimination_Parallel)
        println("ISG:",ISG)
        println("ISG_to_string:",ISG_to_string(ISG))
    end
end

function single_hexagon_2rounds(cycle,Comm_Parallel::Bool,Elimination_Parallel::Bool)
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
        ISG=update(lattice_size,ISG,measurement_rounds[(i-1)%round+1,:,:],Elimination_Parallel,Comm_Parallel)
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

function test_comm_matrix(LatticeSize::Int64,ISGSize::Int64,cycle=3)

    for i in 1:cycle
        ISG=rand(Bool,(ISGSize,2*LatticeSize))
        Measurements=rand(Bool,(ISGSize,2*LatticeSize))
        CommMatrix=zeros(Bool,(ISGSize,ISGSize))

        calc_comm_matrix!(LatticeSize,ISG,Measurements,CommMatrix,false)
        
        println("ISG:",ISG_to_string(ISG))
        println("Measurements:",ISG_to_string(Measurements))
        println("CommMatrix:",CommMatrix)
    end
end

function snake_growing_lattice(Cycle::Int,LatticeSize::Int)
    round=8
end

function test_gaussian_elimination(n)
    small=zeros(Bool,(1,1))

    vecs=rand(Bool,(n,n))
    coeff=zeros(Bool,(n,n))

    gaussian_elimination!(small,small)
    #StatProfilerHTML.@profilehtml gaussian_elimination!(vecs,coeff)
    @time gaussian_elimination!(vecs,coeff)
    # Profile.@profile gaussian_elimination!(vecs,coeff)
    # Profile.print()
end

function test_dynamic_update_performance(LatticeSize::Int64,ISGSize::Int64,cycle=3)
    # Assume measurements are of the same order of ISG.

    ISG=rand(Bool,(ISGSize,2*LatticeSize))
    Measurements=rand(Bool,(ISGSize,2*LatticeSize))

    @time update(LatticeSize,ISG,Measurements,false)
end


single_hexagon_3rounds(2,false,false)
println("False False")
single_hexagon_3rounds(2,true,false)
println("True False")
single_hexagon_3rounds(2,false,true)
println("False True")
single_hexagon_3rounds(2,true,true)
println("True True")

# test_dynamic_update_performance(2,2)
# test_dynamic_update_performance(4096,1024)