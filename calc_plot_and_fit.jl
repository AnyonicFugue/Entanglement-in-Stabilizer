import PyPlot as plt
import CurveFit as cf

function calc_entropy_old(StabGenerators::Array{Bool,2},Region,Parallel::Bool,RegionSize::Int64)
    # region is a vector of the d.o.f. in the region of interest.

    # Strategy: Calculate (rank G-rank G_A-rank G_B)/2.

    # 'Density' means the distance between different samplings.


    @views n_stab=length(StabGenerators[:,1])
    stab_weights_total=zeros(Int64,n_stab) # The weight of each stabilizer. Note that X and Z are both weight 1, while Y is weight 2.
    stab_weights_region=zeros(Int64,n_stab) # The weight of each stabilizer in the region.
    stab_in_A_or_B=zeros(Bool,n_stab) # Whether the stabilizer is in A or B.

    if(Parallel)
        Threads.@threads for i in 1:n_stab
            @views stab_weights_total[i]=sum(StabGenerators[i,:])

            for j in 1:RegionSize
                # Calculate the weight in the region.
                stab_weights_region[i]+=StabGenerators[i,2*Region[j]-1]
                stab_weights_region[i]+=StabGenerators[i,2*Region[j]]
            end
            
            if (stab_weights_region[i]==0) | (stab_weights_region[i]==stab_weights_total[i]) # The stabilizer is in A or B.
                stab_in_A_or_B[i]=true
            end
        end
    
    else
        for i in 1:n_stab
            @views stab_weights_total[i]=sum(StabGenerators[i,:])

            for j in 1:RegionSize
                # Calculate the weight in the region.
                stab_weights_region[i]+=StabGenerators[i,2*Region[j]-1]
                stab_weights_region[i]+=StabGenerators[i,2*Region[j]]
            end
            
            if (stab_weights_region[i]==0) | (stab_weights_region[i]==stab_weights_total[i]) # The stabilizer is in A or B.
                stab_in_A_or_B[i]=true
            end
        end
    end

    entropy=(n_stab-sum(stab_in_A_or_B))/2
    return entropy
end

function calc_entropy(StabGenerators::Array{Bool,2},Region,Parallel::Bool,RegionSize::Int64)
    # region is a vector of the d.o.f. in the region of interest.

    # Strategy: Calculate (rank G-rank G_A-rank G_B)/2.

    # rank G_A and rank G_B are calculated by the rank of the submatrix of G_A and G_B corresponding to the region.
        # Specifically, rank G_A = # columns - rank P_A(G)

    rank_GA=0
    rank_GB=0

    @views n_stab=length(StabGenerators[:,1])
    @views n_dof=length(StabGenerators[1,:]) # The number of d.o.f. in the region.
    partial_stab_A=zeros(Bool,(n_stab,2*RegionSize)) # The partial matrices in the region A.
    partial_stab_B=zeros(Bool,(n_stab,2*(n_dof-RegionSize))) # The partial matrices in the region B.

    coefficients=zeros(Bool,(n_stab,n_stab)) # The coefficients to call the gaussian_elimination.

    cur_A=zeros(Int64,n_stab) # The number of vertices in the region of different start positions.
    cur_B=zeros(Int64,n_stab) # The number of vertices in the region of different start positions.

    # Construct the partial matrices.
    if(Parallel)
        Threads.@threads for i in 1:n_stab
            cur_A[i]=1
            cur_B[i]=1

            for j in 1:Int64(n_dof/2) # The number of vertices.
                if j in Region
                    partial_stab_A[i,cur_A[i]]=StabGenerators[i,2*j-1]
                    cur_A[i]+=1
                    partial_stab_A[i,cur_A[i]]=StabGenerators[i,2*j]
                    cur_A[i]+=1
                else
                    partial_stab_B[i,cur_B[i]]=StabGenerators[i,2*j-1]
                    cur_B[i]+=1
                    partial_stab_B[i,cur_B[i]]=StabGenerators[i,2*j]
                    cur_B[i]+=1
                end
            end
        end
    else
        for i in 1:n_stab
            cur_A[i]=1
            cur_B[i]=1

            for j in 1:Int64(n_dof/2) # The number of vertices.
                if j in Region
                    partial_stab_A[i,cur_A[i]]=StabGenerators[i,2*j-1]
                    cur_A[i]+=1
                    partial_stab_A[i,cur_A[i]]=StabGenerators[i,2*j]
                    cur_A[i]+=1
                else
                    partial_stab_B[i,cur_B[i]]=StabGenerators[i,2*j-1]
                    cur_B[i]+=1
                    partial_stab_B[i,cur_B[i]]=StabGenerators[i,2*j]
                    cur_B[i]+=1
                end
            end
        end
    end

   
    rank_GB=n_stab-gaussian_elimination!(partial_stab_A,coefficients,Parallel) # Mathematically this is the dimension of the zero space of A-submatrices. Each basis element of the zero space corresponds to a generator of G_B.
    rank_GA=n_stab-gaussian_elimination!(partial_stab_B,coefficients,Parallel)

    entropy=(n_stab-rank_GA-rank_GB)/2
    # println("rank_GA:",rank_GA,", rank_GB:",rank_GB)
    # println("entropy:",entropy)

    return entropy
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

function plot_and_fit(EntropyArr::Vector{Float32},VolumeArr::Vector{Int64},AreaArr::Vector{Int64})
    # Plot entropy vs area and volume.


    bv,av=cf.linear_fit(VolumeArr,EntropyArr)
    ba,aa=cf.linear_fit(AreaArr,EntropyArr)

    println("Area-law fitting: y=",aa,"x+",ba)
    plt.plot(AreaArr,EntropyArr,AreaArr,aa*AreaArr.+ba)
    
    plt.title("Entropy vs Area")
    plt.grid()
    plt.show()
    
    println("Volume-law fitting: y=",av,"x+",bv)
    plt.plot(VolumeArr,EntropyArr,VolumeArr,av*VolumeArr.+bv)
    plt.title("Entropy vs Volume")
    plt.grid()
    plt.show()
end

function fit(EntropyArr::Vector{Float32},VolumeArr::Vector{Int64},AreaArr::Vector{Int64})
    # Plot entropy vs area and volume.


    bv,av=cf.linear_fit(VolumeArr,EntropyArr)
    ba,aa=cf.linear_fit(AreaArr,EntropyArr)

    println("Area-law fitting: y=",aa,"x+",ba)

    println("Volume-law fitting: y=",av,"x+",bv)

end