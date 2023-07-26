

function ISG_to_string(ISG::Array{Bool,2},MinLength::Int)
    # Convert the ISG to a string.
    # The first index is the index of the stabilizer. The second index is the index of the d.o.f.
    # The d.o.f. are ordered as X1,Z1,X2,Z2,...,Xn,Zn.

    n=length(ISG[:,1])
    lattice_size=Int(length(ISG[1,:])/2)
    ISG_string=Vector{String}(undef,0)
    
    for i in 1:n
        str=String[]
        current_len=0
        for k in 1:lattice_size
            if (ISG[i,2*k-1])&&(!ISG[i,2*k])
                push!(str,"X")
                push!(str,string(k))
                current_len+=1
            elseif (!ISG[i,2*k-1])&&(ISG[i,2*k])
                push!(str,"Z")
                push!(str,string(k))
                # current_len+=1
            elseif (ISG[i,2*k-1])&&(ISG[i,2*k])
                push!(str,"Y")
                push!(str,string(k))
                current_len+=1
            end  
        end
        
        if(current_len>=MinLength)
            push!(ISG_string,join(str)*",l="*string(current_len))
        end
    end
    return ISG_string
end

function sample_squares(LatticeSideLength::Int64,RegionSideLength::Int64,SampleInterval::Int64,StabGenerators::Array{Bool,2},Parallel::Bool)

    # region_arr=Vector{Vector{Int64}}(undef,0) # The array to store the region of different start positions.
    n_sample_1D=Int64(floor((LatticeSideLength-RegionSideLength)/SampleInterval)) # The number of samples in a single dimension.
        # e.g. LatticeSideLength=10, RegionSideLength=3, SampleInterval=2, then n_sample=4.
    n_sample=n_sample_1D*n_sample_1D # The number of samples in total.

    region_arr=zeros(Int64,(n_sample,2*RegionSideLength*RegionSideLength)) # The array to store the region of different start positions.
        # The first index is the number of the sample.
        # The number of vertices in the region is RegionSideLength*RegionSideLength, so the number of edges is 2*RegionSideLength*RegionSideLength.
    cur=zeros(Int64,n_sample) # The number of vertices in the region of different start positions.
    entropy_arr=zeros(Float32,n_sample) # The entropy of different start positions.

    n=0 # The cursor of the present start_pos.

    for x_start_pos in range(start=1,step=SampleInterval,length=n_sample_1D) # Generate the regions to sample over
        for y_start_pos in range(start=1,step=SampleInterval,length=n_sample_1D)
            # Take a rectangular region starting from start_pos, of length RegionSideLength.
            n+=1
            cur[n]=1

            for i in range(start=x_start_pos,step=1,length=RegionSideLength)
                for j in range(start=y_start_pos,step=1,length=RegionSideLength)
                    m=(i-1)*LatticeSideLength+j # The vertex
                    if(i<x_start_pos+RegionSideLength-1) # Not on the rightmost
                        region_arr[n,cur[n]]=2*m-1
                        cur[n]+=1
                    end

                    if(j<y_start_pos+RegionSideLength-1) # Not on the downmost
                        region_arr[n,cur[n]]=2*m
                        cur[n]+=1
                    end
                    
                end
            end
        end
    end

    if(Parallel)
        Threads.@threads for i in 1:n_sample
            @views entropy_arr[i]=calc_entropy(StabGenerators,region_arr[i,1:cur[i]-1],false,cur[i]-1)
            # @views entropy_arr[i]=calc_entropy_old(StabGenerators,region_arr[i,1:cur[i]-1],false,cur[i]-1)
        end
    else
        for i in 1:n_sample
            @views entropy_arr[i]=calc_entropy(StabGenerators,region_arr[i,1:cur[i]-1],false,cur[i]-1)
            # @views entropy_arr[i]=calc_entropy_old(StabGenerators,region_arr[i,1:cur[i]-1],false,cur[i]-1)
        end
    end

    return sum(entropy_arr)/n_sample

end