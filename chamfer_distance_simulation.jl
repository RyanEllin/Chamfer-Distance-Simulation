using StatsBase
using Plots
using Distributions

#Used to simulate Chamfer-Squared Distance for points uniformly on hypersphere. 
#Can be found exactly with order-statistics, but they're computationally ugly
#Computationally efficient method of estimation using chi-squared distribution:
#Assume that the first point in question is at (1,0,0,...,0) by spherical symmetry
#Chi-Squared to simulate the component not in the e_1 direction
#Component in e_1 direction seperately sampled from normal  
#Then normalize, etc

#size_s2 refers to the number of sample points being compared at once. Badly named
#beta, even worse, is the chi-squared distribution

function single_distance_squared(r, size_s2, nrm, beta)
    dist_list = []
    for j in 1:size_s2
        x_hat = rand(nrm)
        y_hat = rand(beta)
        lgth = sqrt(x_hat^2+y_hat)
        x = x_hat/ lgth
        y = y_hat / (lgth^2)
        dist = (x-1)^2+y 
        append!(dist_list, dist)
    end
    return r^2*minimum(dist_list)
end



function simulator(dim, r, size_s2, trials)
    sum=0
    nrm=Normal()
    beta=Chisq(dim)
    for i in 1:trials
        sum += single_distance_squared(r, size_s2, nrm, beta)
    end
    return sum/trials
end


function main(dim, r, size_low, size_high, trials)
    ys=[]
    for k in size_low:size_high
        append!(ys,simulator(dim,r, k , trials) )
        if mod(k,50) == 0
            print("$k Complete\n")
        end
    end
    xs = [k for k in size_low:size_high]
    pic = plot(xs,ys, size=(1200, 800))
    title!("Chamfer Squared Distance vs Number of Sample Points in $dim Dimensions", titlefontsize=10)
    display(pic)
    readline()
end


main(200, 1, 50, 5000,50)