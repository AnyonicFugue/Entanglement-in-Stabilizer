include("calc_plot_and_fit.jl")
include("dynamic_update.jl")
include("utils.jl")

import CurveFit
# import Profile
# import StatProfilerHTML


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

function toric_code_static_recombination(l::Int64,Parallel::Bool)
    # Note that the d.o.f. are on the edges, not on the vertices.
    # The stabilizers are ordered as X1,Z1,X2,Z2,...,Xn,Zn.

    lattice_size::Int64=2*l*l

    # For a lattice of length l, the size is 2*l*l, since each vertex corresponds to two edges.
    # There is one stabilizer on each vertex and each plaquette, so the total number of stabilizers is 2*l*l, i.e. lattice_size.

    stab_generators=zeros(Bool,(lattice_size,2*lattice_size)) # For each stabilizer the first index is the index of the stabilizer. The second index is the index of the d.o.f.

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

    println("toric code recombination")
    println(length(stab_generators[:,1]))
    # println(stab_generators)

    # Generate Pivotal Matrices
    recomb_mat=zeros(Bool,(lattice_size,lattice_size)) # A rectangular pivotal matrix. The size is lattice_size, the number of stabilizers.
    for i in 1:lattice_size
        recomb_mat[i,i]=true
        recomb_mat[i,i+1:lattice_size]=rand(Bool,lattice_size-i)
    end

    stab_generator_recombined=zeros(Bool,(lattice_size,2*lattice_size)) # For each stabilizer the first index is the index of the stabilizer. The second index is the index of the d.o.f.
    for i in 1:lattice_size
        for j in 1:lattice_size
            if(recomb_mat[i,j])
                @views stab_generator_recombined[i,:]=stab_generator_recombined[i,:] .⊻ stab_generators[j,:]
            end
        end
    end

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
        entropy_arr[s-start+1]=calc_entropy(stab_generator_recombined,view(region,s-start+1,:),Parallel,cur[s-start+1]-1)
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

function straight_growing_ISG(Cycle::Int,LatticeSideLength::Int,Parallel::Bool)
    round=6

    # The vertices are 3-colored. We start from the ISG of single edges on the left.
    # The ISG evolves like this: We first measure a edge connecting the existing one, then measure the vertex between the new edge and the old edge it connects.
        # Note that the vertex anti-commute with a single edge and commute with 2 edges, thus it preserves the large string but eliminate individual edges.
        # Measure Z for vertices and X for edges.

    n_vertex=LatticeSideLength*LatticeSideLength
    LatticeSize=2*n_vertex # The number of edges is twice the number of vertices.
    n_period=Int(floor(LatticeSideLength/3))

    ISG=zeros(Bool,(LatticeSideLength,2*LatticeSize)) 
        # Each edge has a d.o.f., thus we need two indices, one for X and one for Z.
        # We start by single edges on the left.
    Measurement=zeros(Bool,(round,n_period*LatticeSideLength,2*LatticeSize)) 
        # The first index is the round. The second index is the measurements in the round. The third index records the actions of the measurements.
        # In column each round, we measure one in each period.
    
    # Initialize ISG
    for i in 1:LatticeSideLength
        # The vertex has precisely number i!
        n=2*i-1 # The edge on the right of the vertex. 
        ISG[i,2*n-1]=true # Measure X on the edge.
    end

    if(debug)
        println("ISG",ISG_to_string(ISG))
    end

    # Set measurements on each round
    for i in 1:round
        # Measure edges in odd rounds.
        cur=1

        if(i%2==1)
            k=Int((i+1)/2) # k=1 means measure B edges, k=2 means C, k=3 means A.

            for j in 1:n_period # Which period is the measurement in.
                if (j==1 && k==3) # The leftmost A-edge should not be measured.
                    continue
                end

                for y in 1:LatticeSideLength 
                    x=3*(j-1)+(k%3+1) # e.g. When j=1,k=1, x=2, which is the B edge (right edge of the x=2 vertex) in the first period. 
                    m=(x-1)*LatticeSideLength+y # The vertex

                    if(x<LatticeSideLength) # Not on the rightmost
                        Measurement[i,cur,2*(2*m-1)-1]=true # Measure X on the edge. Right edge is 2*m-1 and X on the edge is 2*(2*m-1)-1.
                        cur+=1
                    end
                end
            end
        end

        # Measure vertices in even rounds.
        if(i%2==0)
            k=Int(i/2) # k=1 means measure B vertices (the vertex on the left of B edges), k=2 means C, k=3 means A.

            for j in 1:n_period
                if (j==1 && k==3) # The leftmost A-vertex should not be measured.
                    continue
                end

                for y in 1:LatticeSideLength 
                    x=3*(j-1)+(k%3+1) # e.g. When j=1,k=1, x=2, which is the B vertex in the first period. 
                    m=(x-1)*LatticeSideLength+y # The vertex

                    left_m=(x-2)*LatticeSideLength+y # The vertex on the left of the vertex m.
                    up_m=(x-1)*LatticeSideLength+y-1 # The vertex on the top of the vertex m.

                    # The left edge, i.e. the right edge of the left vertex, is always measured, since we've already excluded the leftmost vertex.
                    Measurement[i,cur,2*(2*left_m-1)]=true # Left edge is 2*left_m-1 and Z on the edge is 2*(2*left_m-1).

                    if(x<LatticeSideLength) # Not on the rightmost
                        Measurement[i,cur,2*(2*m-1)]=true # Measure Z on the right edge of the vertex. Right edge is 2*m-1 and Z on the edge is 2*(2*m).
                    end

                    if(y<LatticeSideLength) # Not on the downmost
                        Measurement[i,cur,2*(2*m)]=true # Measure Z on the down edge of the vertex. Down edge is 2*m and Z on the edge is 2*(2*m-1).
                    end

                    if(y>1) # Not on the upmost
                        Measurement[i,cur,2*(2*up_m)]=true # Measure Z on the up edge of the vertex. Up edge (down edge of the upper vertex) is 2*up_m and Z on the edge is 2*(2*up_m).
                    end

                    cur+=1
                end
            end
        end
    end

    if(debug)
        for i in 1:round
            println("Measurement in round ",i,":",ISG_to_string(Measurement[i,:,:]))
        end
        println()
    end

    measurement_vector=Vector{Matrix{Bool}}(undef,0)
    for i in 1:round
        push!(measurement_vector,Measurement[i,:,:])
    end


    # Start updating
    for i in 1:Cycle

        println("Cycle ",i)

        for j in 1:round
            @views ISG=update(LatticeSize,ISG,measurement_vector[j],Parallel,Parallel)
            # println("ISG:",ISG)
            println("ISG_to_string in round ",j,": ",ISG_to_string(ISG))
        end
        println()
    end
end

function snake_growing_ISG(Cycle::Int,LatticeSideLength::Int)
    round=8

    # The lattice is 4-colored. We start from the ISG of single edges on the left.
    # The ISG evolves like this: We first measure a edge connecting the existing one, then measure the vertex between the new edge and the old edge it connects.
        # Note that the vertex anti-commute with a single edge and commute with 2 edges, thus it preserves the large string but eliminate individual edges.

    n_vertex=LatticeSideLength*LatticeSideLength
    LatticeSize=2*n_vertex # The number of edges is twice the number of vertices.

    ISG=zeros(Bool,(round,2*LatticeSize)) # Each edge has a d.o.f., thus we need two indices, one for X and one for Z.
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

debug=false


