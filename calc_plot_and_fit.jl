import PyPlot as plt
import CurveFit as cf
# import Plots

function calc_entropy(StabGenerators::Vector{Tuple{Vector{Int32},Vector{Int8}}},region::Vector{Int32})
    # region is a vector of the d.o.f. in the region of interest.
    
    ## println("calc_entropy")

    # Strategy: Calculate (rank G-rank G_A-rank G_B)/2.
    
    rank_GA=0
    rank_GB=0

    for i in 1:length(StabGenerators)
        if isdisjoint(StabGenerators[i][1],region)
            rank_GB+=1
        end
        if issubset(StabGenerators[i][1],region)
            rank_GA+=1
        end
    end

    entropy=(length(StabGenerators)-rank_GA-rank_GB)/2
    ## println("rank_GA:",rank_GA,", rank_GB:",rank_GB)
    ## println("entropy:",entropy)
    return entropy
end

function plot_and_fit(EntropyArr::Vector{Float32},VolumeArr::Vector{Int32},AreaArr::Vector{Int32})
    # Plot entropy vs area and volume.


    bv,av=cf.linear_fit(VolumeArr,EntropyArr)
    ba,aa=cf.linear_fit(AreaArr,EntropyArr)

    println("Area-law fitting: y=",aa,"x+",ba)
    plt.plot(AreaArr,EntropyArr,AreaArr,aa*AreaArr.+ba)
    
    plt.title("Entropy vs Area")
    plt.grid()
    ## plt.show() Special solution for github codespace
    readline()
    
    println("Volume-law fitting: y=",av,"x+",bv)
    plt.plot(VolumeArr,EntropyArr,VolumeArr,av*VolumeArr.+bv)
    plt.title("Entropy vs Volume")
    plt.grid()
    readline()
    ## plt.show()
end

function plot_and_fit_native(EntropyArr::Vector{Float32},VolumeArr::Vector{Int32},AreaArr::Vector{Int32})
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


function fit(EntropyArr::Vector{Float32},VolumeArr::Vector{Int32},AreaArr::Vector{Int32})
    # Plot entropy vs area and volume.


    bv,av=cf.linear_fit(VolumeArr,EntropyArr)
    ba,aa=cf.linear_fit(AreaArr,EntropyArr)

    println("Area-law fitting: y=",aa,"x+",ba)

    println("Volume-law fitting: y=",av,"x+",bv)

end