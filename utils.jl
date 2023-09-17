

function ISG_to_string(ISG::Array{Bool,2},MinLength::Int,Width::Int,Depth::Int,Display_Coordinate::Bool=true)
    # Convert the ISG to a string.
    # The first index is the index of the stabilizer. The second index is the index of the d.o.f.
    # The d.o.f. are ordered as X1,Z1,X2,Z2,...,Xn,Zn.

    n=length(ISG[:,1])
    lattice_size=Int(length(ISG[1,:])/2)
    ISG_string=Vector{String}(undef,0)

    if(Display_Coordinate)
        for i in 1:n
            str=String[]
            current_len=0
            for k in 1:lattice_size

                # k = Width*(y-1) + x
                # Calculate x and y from k
                x=(k-1)%Width+1
                y=Int((k-x)/Width)+1

                if (ISG[i,2*k-1])&&(!ISG[i,2*k])
                    push!(str,"X")

                    push!(str,"(")
                    push!(str,string(x))
                    push!(str,",")
                    push!(str,string(y))
                    push!(str,")")
                    current_len+=1
                elseif (!ISG[i,2*k-1])&&(ISG[i,2*k])
                    push!(str,"Z")

                    push!(str,"(")
                    push!(str,string(x))
                    push!(str,",")
                    push!(str,string(y))
                    push!(str,")")
                    current_len+=1
                elseif (ISG[i,2*k-1])&&(ISG[i,2*k])
                    push!(str,"Y")

                    push!(str,"(")
                    push!(str,string(x))
                    push!(str,",")
                    push!(str,string(y))
                    push!(str,")")
                    current_len+=1
                end  
            end
            
            if(current_len>=MinLength)
                push!(ISG_string,join(str)*",l="*string(current_len))
            end
        end
    else
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
                    current_len+=1
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

function N_nontrivial_sites(g::Array{Bool})
    # This function calculates the number of sites where action of g is nontrivial.

    n=0
    l=Int(size(g,1)/2)

    for i in 1:l
        if(g[2*i-1]||g[2*i])
            n+=1
        end
    end

    return n
end

function Find_local_generators!(ISG::Array{Bool,2})
    # ISG stores the generators for the stabilizer group, which could be viewed as a basis for the vector space of the 2-element field. The stabilizer group is therefore a vector space spanned by these generators (i.e. basis elements).
    # The entries of the basis elements are either 0 or 1, with the addition rule 1+1=0.
    # The function performs elementary transformations to the basis vectors (i.e. add one vector to another) to find a basis such that the total number of entries with value 1 is minimized.

    n=size(ISG,1) # The number of stabilizers.

    # First, calculate the total number of 1 in all basis elements.
    total_1=0
    for i in 1:n
        total_1+=sum(ISG[i,:])
    end

    # Then, perform elementary transformations to the basis vectors to minimize the total number of 1.

    flipped=true


    while(flipped)

        flipped=false

        for i in 1:n
            for j in 1:n
                if(i!=j)
                    after_generator=zeros(Bool,1,size(ISG,2))
                    after_generator[1,:]=ISG[j,:] .⊻ ISG[i,:] # The basis element after adding the i-th basis element to the j-th basis element.
                    after_weight=sum(after_generator)

                    if(sum(ISG[j,:])>after_weight)
                        # If the total number of 1 in the j-th basis element is larger than that in the j-th basis element after adding the i-th basis element, then perform the transformation.
                        ISG[j,:] = after_generator # Add the i-th basis element to the j-th basis element.
                        flipped=true
                    elseif(sum(ISG[j,:])==after_weight)
                        # If the total number of 1 in the j-th basis element is equal to that in the j-th basis element after adding the i-th basis element, then compare the number of nontrivial sites.
                        if(N_nontrivial_sites(ISG[j,:])>N_nontrivial_sites(after_generator[1,:]))
                            @views ISG[j,:] = after_generator[1,:] # Add the i-th basis element to the j-th basis element.
                            flipped=true
                        end
                    end

                end
            end
        end
    end

end

function Evaluate_EE(ISG::Array{Bool,2},NSite::Int,Parallel::Bool=false)
    # The function uses Gaussian elimination to calculate the rank of the subgroup of ISG that acts trivially out of sites {1,...,n}.

    double_Width=size(ISG,2) # 2*Width
    Width=Int(double_Width/2) # The width of the lattice.
    n_stab=size(ISG,1)

    # EE = N_A-|G_A|=N_B-|G_B|. Calculate N_B-|G_B| for simplicity since NSite is smaller than half.

    ISG_A=ISG[:,1:2*NSite] # The stabilizer group inside the region.
    
    coefficients=zeros(Bool,n_stab,n_stab) # The coefficient matrix to be used in gaussian_elimination.
    rank=gaussian_elimination!(ISG_A,coefficients,Parallel) # The rank equals the total rank minus the dimension of the zero space of truncated partial stabilizers in A (the number of stabilizers has support only in the region B).
    # rank = |G|-|G_B|=Width-|G_B|
    EE_A=rank-NSite # EE = N_A-|G_A|=N_B-|G_B|

    return EE_A


end