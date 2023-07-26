include("calc_plot_and_fit.jl")
include("dynamic_update.jl")
include("utils.jl")

import Profile

function straight_growing_ISG(Cycle::Int,LatticeSideLength::Int,SampleInterval::Int64,Parallel::Bool,debug::Bool=false) # Each cycle the string grows length 3. When it reaches the right boundary it would be killed by the rightmost vertex.
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

    measurement_vector=Vector{Matrix{Bool}}(undef,0) # Use a vector to store measurements, in order to avoid dynamic typing for SubArray.
    for i in 1:round
        push!(measurement_vector,Measurement[i,:,:])
    end


    start=Int64(floor(LatticeSideLength/5))

    volume_arr=zeros(Int64,Int64(LatticeSideLength/2)-start+1)
    area_arr=zeros(Int64,Int64(LatticeSideLength/2)-start+1)
    entropy_arr=zeros(Float32,Int64(LatticeSideLength/2)-start+1)

    r2_arealaw_arr=zeros(Float32,Cycle)
    r2_volumelaw_arr=zeros(Float32,Cycle)

    println("Size of ISG:",size(ISG))

    # Start updating
    for i in 1:Cycle
        if(debug)
            println("Cycle ",i)
        end
        
        for j in 1:round
            @views ISG=update(LatticeSize,ISG,measurement_vector[j],Parallel,Parallel)
            #= println("ISG:",ISG)
            if(debug)
                println("ISG_to_string in round ",j,": ",ISG_to_string(ISG))
                println("Round ",j," finished.")
            end
            =#

        end

        println("Size of ISG:",size(ISG))

        for l in range(max(start,2),Int64(LatticeSideLength/2)) # l is the side length of the sampling region.
    
            area_arr[l-start+1]=4*l
            volume_arr[l-start+1]=2*l*l
            
            entropy_arr[l-start+1]=sample_squares(LatticeSideLength,l,SampleInterval,ISG,Parallel)
        end

        # plot_and_fit(entropy_arr,volume_arr,area_arr)
        r2_arealaw_arr[i],r2_volumelaw_arr[i] = fit(entropy_arr,volume_arr,area_arr)

        # println(ISG_to_string(ISG,Int64(LatticeSideLength/4)))
    end

    println("Area law r2: ",r2_arealaw_arr)
    println("Volume law r2: ",r2_volumelaw_arr)
end


function main()
    
    straight_growing_ISG(10,24,16,true,true)
    
    #=
    profile_io=open("profile.txt","w")
    Profile.@profile straight_growing_ISG(10,16,16,true,true)
    println("Please enter the min count in profiling:")
    min_count=parse(Int64,readline())

    Profile.print(profile_io,mincount=min_count)
    close(profile_io)
    =#
end

main()