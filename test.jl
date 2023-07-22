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
    measurement_rounds[1,1,:]=[true,false,true,false,false,false] #X1X2
    measurement_rounds[2,1,:]=[false,false,false,true,false,true] #Z2Z3
    measurement_rounds[3,1,:]=[true,true,false,false,true,true] #Y1Y3

    for i in 1:3
        ISG=update(LatticeSize,ISG,measurement_rounds[1,:,:])
        println("ISG:",ISG)
    end

    
    
end

function floquet_code()
    
end

single_triangle()