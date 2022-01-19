#######################################
# Analysiss by the model with confidence as difference in FR
########################################

# necessary packages

using DataFrames,CSV

rand(123456789)
#plotlyjs()
print("Go")
using Distributions
#using LaTeXStrings              # nice LaTeX strings
#using FileIO

using Bootstrap
print("Go")
using StatsBase

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
print("Go")
function getRandom(func, xmin, xmax, Nvalues, Ndiv)
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


#%% Load the data
include(".\\decision_making.jl")
include(".\\various_functions.jl")

function sim_mf_decisionmaking_pattern_fastv2(coh_temp::Array{Float64},Ir::Float64,If::Float64,
  param::Dict{String,Float64},fixed_param::Dict{String,Float64},simu_param::Dict{String,Float64})



  # FI curve parameters - adapted from Wanf fit
  a=fixed_param["a"]
  b=fixed_param["b"]
  d=fixed_param["d"]
  falset= param["false_time"]
  thr = param["Threshold"]


  # time variables of the  simulation

  dt=simu_param["dt"]
  time_wind=simu_param["time_wind"] # mean window of the variables
  slide_wind = simu_param["slide_wind"]
  false_time=param["false_time"]
  falsetIf=param["false_timeIf"]


  #---- Intialise and vectorise variables to be used in loops below ------

  # Final Variables
  res=zeros(length(coh_temp))
  reaction_time=zeros(length(coh_temp))
  diffS=zeros(length(coh_temp))
  S1=zeros(length(coh_temp))
  S2=zeros(length(coh_temp))
  R1=zeros(length(coh_temp))
  R2=zeros(length(coh_temp))

  ###### Simu variables aprameters
  noise_amp=param["noise_amp"]

  mu0=  param[  "mu0"]

  I0E1p=  param[  "I0E1"]
  I0E2p =  param[  "I0E2"]
  JN11 =  param[ "JN11"]
  JN22  = param[ "JN22"]
  JN12 =  param[ "JN12"]
  JN21 =  param[ "JN21"]
  JAext  = param[ "JAext"]
  Tnmda = param[  "Tnmda"]
  gamma = param[  "gamma"]
  Tstim = param[  "Tstim"]
  Tampa = param[  "Tampa"]
  coh=  param[ "coh"]

  false_timeIf=  param[  "false_timeIf"]

  #####

  # Temporary variables

    Isyn1=0.0
  Isyn2=0.0
  Trest=param["T_rest"]
  s1=0.1
  s2=0.1
  phi1=0.0
  phi2=0.0
  I_eta1=noise_amp * randn()
  I_eta2=noise_amp * randn()

  nu1=Float64[2.0]
  nu2=Float64[2.0]
  k=0 # variable for mod 4 in the mean
  th=false
  result=true

  for (j,coh) in enumerate(coh_temp) # loop on trials number

    th=false
    result=true
    t=0

    nu1=Float64[nu1[end]]
    nu2=Float64[nu2[end]]
    k=0
    while th==false # while we do not cross the threshold
      t=t+dt
      k=k+1

      if t<Trest && result==true
        I0E1=I0E1p -Ir*exp(-t/falset)
        I0E2=I0E2p -Ir*exp(-t/falset)

        I_stim1=0
        I_stim2=0

      elseif t<Trest
        I0E1=I0E1p -Ir*exp(-t/falset) -If*exp(-t/falsetIf)
        I0E2=I0E2p -Ir*exp(-t/falset) -If*exp(-t/falsetIf)
        I_stim1=0
        I_stim2=0

      elseif t<Trest+dt/2 && t>Trest-dt/2
        I0E1=I0E1p
        I0E2=I0E2p
        I_stim1 = JAext*mu0*(1+coh/100)
        I_stim2 = JAext*mu0*(1-coh/100)


      else
        I0E1=I0E1p
        I0E2=I0E2p
        I_stim1 = JAext*mu0*(1+coh/100)
        I_stim2 = JAext*mu0*(1-coh/100)
      end


      Isyn1=JN11*s1-JN12*s2+I_stim1 + I_eta1
      Isyn2=JN22*s2-JN21*s1+I_stim2 + I_eta2

      phi1=(a*Isyn1-b)/(1-exp(-d*(a*Isyn1-b)))
      phi2=(a*Isyn2-b)/(1-exp(-d*(a*Isyn2-b)))

      #---- Dynamical equations -------------------------------------------
      # Mean NMDA-mediated synaptic dynamics updating
      s1 = s1 + dt*(-s1/Tnmda + (1-s1)*gamma*nu1[end]/1000);


      s2 = s2 + dt*(-s2/Tnmda + (1-s2)*gamma*nu2[end]/1000);


      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1 = I_eta1 + (dt/Tampa)*(I0E1-I_eta1) + sqrt(dt/Tampa)*noise_amp*randn() ;
      I_eta2 = I_eta2 + (dt/Tampa)*(I0E2-I_eta2) + sqrt(dt/Tampa)*noise_amp*randn() ;

      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1 < 0
        push!(nu1,0)
        phi1 = 0;
      else
        push!(nu1,phi1)
      end;
      if phi2 < 0
        push!(nu2,0)
        phi2 = 0;
      else
        push!(nu2,phi2)
      end;

      if mod(k,4)==0

        # Test the crossing of a threshold
        @views m1= mean(nu1[k-3:k])
        @views m2 = mean(nu2[k-3:k])
        if  t>Trest&& m1>thr
          th=true
          if coh>=0
          @inbounds  res[j]=1
            result=true
          else
            @inbounds  res[j]=-1
            result=false
          end

        elseif t>Trest&& m2>thr
          th = true
          if coh>=0
            @inbounds  res[j]=-1
            result=false
          else
            @inbounds  res[j]=1
            result=true
          end

        end
      end
    end
  @inbounds  reaction_time[j]=t-Trest
  @inbounds    diffS[j]=s2-s1
  @inbounds    R1[j]=nu1[end]
  @inbounds  R2[j]=nu2[end]
  @inbounds    S1[j]=s1
  @inbounds  S2[j]=s2
    th=false
  end

  result=DataFrame()
  result[:RTs]=reaction_time
  result[:NotError]=res
  result[:Distri]=coh_temp
  result[:DiffS]=diffS
  result[:S1]=S1
  result[:S2]=S2
  result[:R1]=R1
  result[:R2]=R2
  result[:LastTrial]=vcat([0],res[1:end-1])
  return result
end


function sim_mf_decisionmaking_pattern_save(coh_temp::Array{Float64},Ir::Float64,If::Float64,
  param::Dict{String,Float64},fixed_param::Dict{String,Float64},simu_param::Dict{String,Float64})



  # FI curve parameters - adapted from Wang fit
  a=fixed_param["a"]
  b=fixed_param["b"]
  d=fixed_param["d"]
  falset= param["false_time"]
  thr = param["Threshold"]


  # time variables of the  simulation

  dt=simu_param["dt"]
  time_wind=simu_param["time_wind"] # mean window of the variables
  slide_wind = simu_param["slide_wind"]
  false_time=param["false_time"]
  falsetIf=param["false_timeIf"]


  #---- Intialise and vectorise variables to be used in loops below ------

  # Final Variables
  res=zeros(length(coh_temp))
  reaction_time=zeros(length(coh_temp))
  diffS=zeros(length(coh_temp))
  S1=zeros(length(coh_temp))
  S2=zeros(length(coh_temp))
  R1=zeros(length(coh_temp))
  R2=zeros(length(coh_temp))

  FR1=[]
  FR2 = []

  ###### Simu variables aprameters
  noise_amp=param["noise_amp"]

  mu0=  param[  "mu0"]

  I0E1p=  param[  "I0E1"]
  I0E2p =  param[  "I0E2"]
  JN11 =  param[ "JN11"]
  JN22  = param[ "JN22"]
  JN12 =  param[ "JN12"]
  JN21 =  param[ "JN21"]
  JAext  = param[ "JAext"]
  Tnmda = param[  "Tnmda"]
  gamma = param[  "gamma"]
  Tstim = param[  "Tstim"]
  Tampa = param[  "Tampa"]
  coh=  param[ "coh"]

  false_timeIf=  param[  "false_timeIf"]

  #####

  # Temporary variables

    Isyn1=0.0
  Isyn2=0.0
  Trest=param["T_rest"]
  s1=0.1
  s2=0.1
  phi1=0.0
  phi2=0.0
  I_eta1=noise_amp * randn()
  I_eta2=noise_amp * randn()

  nu1=Float64[2.0]
  nu2=Float64[2.0]
  k=0 # variable for mod 4 in the mean
  th=false
  result=true
  tempFR1=Float64[]
  tempFR2=Float64[]
  for (j,coh) in enumerate(coh_temp) # loop on trials number

    th=false
    result=true
    t=0

    nu1=Float64[nu1[end]]
    nu2=Float64[nu2[end]]
    k=0
    if j!=1
      push!(FR1,tempFR1)
      push!(FR2,tempFR2)
    end
    tempFR1=Float64[]
    tempFR2=Float64[]
    while th==false # while we do not cross the threshold
      t=t+dt
      k=k+1

      if t<Trest && result==true
        I0E1=I0E1p -Ir*exp(-t/falset)
        I0E2=I0E2p -Ir*exp(-t/falset)

        I_stim1=0
        I_stim2=0

      elseif t<Trest
        I0E1=I0E1p -Ir*exp(-t/falset) -If*exp(-t/falsetIf)
        I0E2=I0E2p -Ir*exp(-t/falset) -If*exp(-t/falsetIf)
        I_stim1=0
        I_stim2=0

      elseif t<Trest+dt/2 && t>Trest-dt/2
        I0E1=I0E1p
        I0E2=I0E2p
        I_stim1 = JAext*mu0*(1+coh/100)
        I_stim2 = JAext*mu0*(1-coh/100)


      else
        I0E1=I0E1p
        I0E2=I0E2p
        I_stim1 = JAext*mu0*(1+coh/100)
        I_stim2 = JAext*mu0*(1-coh/100)
      end


      Isyn1=JN11*s1-JN12*s2+I_stim1 + I_eta1
      Isyn2=JN22*s2-JN21*s1+I_stim2 + I_eta2

      phi1=(a*Isyn1-b)/(1-exp(-d*(a*Isyn1-b)))
      phi2=(a*Isyn2-b)/(1-exp(-d*(a*Isyn2-b)))

      #---- Dynamical equations -------------------------------------------
      # Mean NMDA-mediated synaptic dynamics updating
      s1 = s1 + dt*(-s1/Tnmda + (1-s1)*gamma*nu1[end]/1000);


      s2 = s2 + dt*(-s2/Tnmda + (1-s2)*gamma*nu2[end]/1000);


      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1 = I_eta1 + (dt/Tampa)*(I0E1-I_eta1) + sqrt(dt/Tampa)*noise_amp*randn() ;
      I_eta2 = I_eta2 + (dt/Tampa)*(I0E2-I_eta2) + sqrt(dt/Tampa)*noise_amp*randn() ;

      push!(tempFR1,s1)
      push!(tempFR2,s2)
      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1 < 0
        push!(nu1,0)
        phi1 = 0;
      else
        push!(nu1,phi1)
      end;
      if phi2 < 0
        push!(nu2,0)
        phi2 = 0;
      else
        push!(nu2,phi2)
      end;

      if mod(k,4)==0

        # Test the crossing of a threshold
        @views m1= mean(nu1[k-3:k])
        @views m2 = mean(nu2[k-3:k])
        if  t>Trest&& m1>thr
          th=true
          if coh>=0
          @inbounds  res[j]=1
            result=true
          else
            @inbounds  res[j]=-1
            result=false
          end

        elseif t>Trest&& m2>thr
          th = true
          if coh>=0
            @inbounds  res[j]=-1
            result=false
          else
            @inbounds  res[j]=1
            result=true
          end

        end
      end
    end
  @inbounds  reaction_time[j]=t-Trest
  @inbounds    diffS[j]=s2-s1
  @inbounds    R1[j]=nu1[end]
  @inbounds  R2[j]=nu2[end]
  @inbounds    S1[j]=s1
  @inbounds  S2[j]=s2
    th=false
  end

  result=DataFrame()
  result[:RTs]=reaction_time
  result[:NotError]=res
  result[:Distri]=coh_temp
  result[:DiffS]=diffS
  result[:S1]=S1
  result[:S2]=S2
  result[:R1]=R1
  result[:R2]=R2
  result[:LastTrial]=vcat([0],res[1:end-1])
  return result,FR1,FR2
end

rand(123456789)
#%% Param of the subjects

print("Go")
weightsL = [0.12, 0.18, 0.2, 0.2, 0.18,0.12]
df2 = CSV.read(".\\Manip3.csv")



param = Dict(
    "mu0"=>30.0,
    "noise_amp"=>0.02,
    "I0E1"=>0.3255,
    "I0E2"=>0.3255,
    "JN11"=>0.2609,
    "JN22"=>0.2609,
    "JN12"=>0.0497,
    "JN21"=>0.0497,
    "JAext" => 0.00052,
    "Tnmda" => 100.0,
    "Tampa" => 2.0,
    "gamma" => 0.641,
    "Tstim"=>5000.0,
    "T_rest"=>700.0,
    "coh"=>5.0,
    "false_time"=>150.0,
    "false_timeIf"=>500.0,
    "Threshold"=>12.0
    )
    fixed_param=Dict("a"=>270.0,"b"=>108.0,"d"=>0.1540)
  simu_param=Dict("dt"=>0.5,"time_wind"=>4.0,"slide_wind"=>4.0) #4=2/dt
Name=[1,3,5,6,7,9]
Threshold=[13.08,13.70,14.95,12.56,12.55,14.65]
c02=[3.2,2.68,1.32,2.35,1.91,3.9]
c08=[14.95,11.81,7.1,8.15,13.49,14.1]
c16=[34.5,18.04,17.5,30.01,50.21,28.1]
Variable_mu = [0.422602, 0.390527,
 0.510844, 0.529926, 0.436503,  0.433333]
Variable_gamma = [3599.17,  657.64,
 323.661, 699.619, 1147.85,  1268.4]
Variable_lambda =[6.98586,  10.6735,  5.06607,
3.11227, 5.63951,  6.51043]

vmufinal=[0.17476,
0.10728,
0.193,
0.268,
0.203,
0.148]


total = 0
for n=1:6

t=length(df[:AngleValue][df[:Name].==Name[n]])
total=t+total
end

# Compute our histogram completely normalized
function myhist(data, min, max, nbins)
  N = length(data)             # How many elements in the input vector 'data' ?
  delta = (max-min)/nbins      # Bin size is inferred here from the maximal, minimal, and bin number
  out = zeros(nbins)           # Let's initialize the output data structures for the bin count
  bin = zeros(nbins)           # and for the bin centres...

  start = min                  # Left edge
  for d=1:nbins
    stop   = start + delta   # Right edge
    out[d] = length(find((data .>= start) .& (data .< stop))) # Count how many elements are between left and right
    bin[d] = start + delta/2. # Centre of the bin
    start  = stop            # New left edge
   end
   return out, bin
  end

dftot=DataFrame()

NameS=zeros(total)
Pconf=zeros(total)
RT=zeros(total)
Acc=zeros(total)
Distri=zeros(total)
ltemp=1

n=1
for n=1:6
  vgamma2=Variable_gamma[n]
  vmu2=Variable_mu[n]
  vlambda=Variable_lambda[n]
cohvalue= [-c16[n], -c08[n], -c02[n], c02[n], c08[n], c16[n]]
If=0.0
  distri=df[:AngleValue][df[:Name].==Name[n]]

  distri[distri.==1.6]=c16[n]
  distri[distri.==-1.6]=-c16[n]
  distri[distri.==-0.8]=-c08[n]
  distri[distri.==0.8]=c08[n]
  distri[distri.==-0.2]=-c02[n]
  distri[distri.==0.2]=c02[n]
distri2=zeros(length(distri))
for i=1:length(distri)
distri2[i]=distri[i]
end
distri=distri2

  ltempdis=length(distri)

Ir=0.033
param["Threshold"]=Threshold[n]
param["mu0"]=35.0
resultM3= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If, param,fixed_param,simu_param)


Rateconf=abs.(resultM3[:R1]-resultM3[:R2])
h=StatsBase.fit(Histogram,Rateconf,nbins=200,closed=:right)

cum_function=Float64[0] # Cumulative function of DeltaS
for w in h.weights
        push!(cum_function,w+cum_function[end])
end

out, bin = myhist(df[:Resp2][df[:Name].==Name[n]], 0, 9, 10)  # Compute the histogram, between -10 and 10, with 100 bins

normalised = out ./ (length(df[:Resp2][df[:Name].==Name[n]]) )    # Normalised
cprime=cum_function./cum_function[end]
param_step=Float64[]
for i in 1:10
        push!(param_step,(searchsortedlast(cprime,sum(normalised[1:i]))))
end


rlim2=zeros(length(cprime))
for x=0:length(cprime)-1
    rlim2[x+1]=x*(h.edges[1][2]-h.edges[1][1])+h.edges[1][1]

end

rlim=rlim2[Int64.(param_step)]

function conf_non_parametrique(x,param)
        #  Directly compute the hist. with Laughlin entropy (non parametric)
        # param : subdivision of DeltaS space
        rep = Int[]
        for xtemp in x
                push!(rep,searchsortedlast(param,xtemp))
        end
        return rep
end

dfName1=conf_non_parametrique(Rateconf,rlim)
resultM3[:Conf]=dfName1
print(h.edges)
NameS[Int(ltemp):Int(ltempdis+ltemp-1)]=n.*ones(ltempdis)
Pconf[Int(ltemp):Int(ltempdis+ltemp-1)]=dfName1
RT[Int(ltemp):Int(ltempdis+ltemp-1)]=resultM3[:RTs]
Acc[Int(ltemp):Int(ltempdis+ltemp-1)] = resultM3[:NotError]
Distri[Int(ltemp):Int(ltempdis+ltemp-1)] = resultM3[:Distri]

ltemp=ltemp+ltempdis

end


dftot[:Acc]=Acc
dftot[:Name]=NameS

dftot[:Pconf]=Pconf
dftot[:Distri]=Distri
dftot[:RTs]=RT

CSV.write(".\\dataNetwork.csv",dftot)
#<<<<<<<<<<<<<<
