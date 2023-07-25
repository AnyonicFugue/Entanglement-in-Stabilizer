import PyPlot as plt
import CurveFit as cf

# import Plots


function calc_entropy(StabGenerators::Array{Bool,2},Region,Parallel::Bool,RegionSize::Int64)
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
    ## println("rank_GA:",rank_GA,", rank_GB:",rank_GB)
    ## println("entropy:",entropy)
    return entropy
end

function sample_squares(LatticeSideLength::Int64,RegionSideLength::Int64,SampleInterval::Int64,StabGenerators::Array{Bool,2})

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


    Threads.@threads for i in 1:n_sample
        @views entropy_arr[i]=calc_entropy(StabGenerators,region_arr[i,1:cur[i]-1],false,cur[i]-1)
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

function plot_and_fit_native(EntropyArr::Vector{Float32},VolumeArr::Vector{Int64},AreaArr::Vector{Int64})
    # Plot entropy vs area and volume.


    bv,av=cf.linear_fit(VolumeArr,EntropyArr)
    ba,aa=cf.linear_fit(AreaArr,EntropyArr)

    println("Area-law fitting: y=",aa,"x+",ba)
    Plots.plot(AreaArr,EntropyArr,label="Entropy vs Area")
    Plots.plot(AreaArr,aa*AreaArr.+ba,label="Fitted")
    ## plt.plot(AreaArr,EntropyArr,AreaArr,aa*AreaArr.+ba)
    
    #Plots.title("Entropy vs Area")
    #Plots.grid()
    ## plt.show() Special solution for github codespace
    # readline()
    
    println("Volume-law fitting: y=",av,"x+",bv)
    Plots.plot(VolumeArr,EntropyArr,label="Entropy vs Volume")
    Plots.plot(VolumeArr,av*VolumeArr.+bv,label="Fitted")
    # lots.title("Entropy vs Volume")
    # Plots.grid()
    # readline()
    ## plt.show()
end


function fit(EntropyArr::Vector{Float32},VolumeArr::Vector{Int64},AreaArr::Vector{Int64})
    # Plot entropy vs area and volume.


    bv,av=cf.linear_fit(VolumeArr,EntropyArr)
    ba,aa=cf.linear_fit(AreaArr,EntropyArr)

    println("Area-law fitting: y=",aa,"x+",ba)

    println("Volume-law fitting: y=",av,"x+",bv)

end