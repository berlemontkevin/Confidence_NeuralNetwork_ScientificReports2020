# Computing some distributions for discrete values
using Distributions

function make_distributions(c_min,c_max,number,N)

values = linspace(c_min,c_max,number)

    return sample(values,N)

end



function invertLinear(t,xmin,xmax,ymin,ymax)
    # special case, if Pol1
    ymin == ymax && return xmin + t*(xmax-xmin)
    # otherwise, Pol2
    A = (ymax-ymin)/(xmax-xmin)/2;
    B = ymin;
    C = - t*(ymin*(xmax-xmin)+(xmax-xmin)*(ymax-ymin)/2);
    x = xmin + (-B + sqrt(B^2-4*A*C))/(2*A)
    return x
end


function getRandom(func, xmin, xmax, Nvalues, Ndiv)
# Sampling from a specific distribution
    x = linspace(xmin,xmax,Ndiv)
    P = func.(x)
    # remove interval where the function is negative
    P = [v<0 ? 0 : v for v in P]
    # weight is given by integral
    Pa = [(P[i]+P[i+1])/2 for i in 1:Ndiv-1]
    Pa = Pa/sum(Pa)
    # generate set on random bin indexes
    inds = StatsBase.sample(1:Ndiv-1, Weights(Pa), Nvalues)
    # convert the set of indexes to random variables inside a bin
    return [invertLinear(rand(),x[i],x[i+1],P[i],P[i+1]) for i in inds]
end

function myhist(data, min, max, nbins)
# Compute the normalized histogram
  N = length(data)             # How many elements in the input vector 'data' ?
  delta = (max-min)/nbins      # Bin size is inferred here from the maximal, minimal, and bin number
  out = zeros(nbins)           # Let's initialize the output data structures for the bin count
  bin = zeros(nbins)           # and for the bin centres...

  start = min                  # Left edge
  for k=1:nbins
    stop   = start + delta   # Right edge
    out[k] = length(findall((data .>= start) .& (data .< stop))) # Count how many elements are between left and right
    bin[k] = start + delta/2. # Centre of the bin
    start  = stop            # New left edge
   end
   return out, bin
  end


