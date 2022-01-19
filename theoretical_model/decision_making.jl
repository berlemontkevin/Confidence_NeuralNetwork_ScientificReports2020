# Implementation of a decision making model with feedback of responses or confidence

# Some part of this implementation are extended from the original model.
# The idea is to model decision making with two competitve population
# in Mean field framework
using Distributions
using DataFrames, CSV



function sim_mf_decisionmaking_pattern_confidence(coh_temp,Ir,IF,If_function,conf_function,
  param,fixed_param,simu_param)
  # Description : Return the trajectory of the dynamical system in the case of
  # decision making with two competitive populations
  # and for a serie of stimulus.
  # Here the reset is done to coincide with the confidence task of Jerome


  # Args:
  # List_rest : list of the time of reset between stimuli
  # List_coh : coherence of the several stimuli
  # List_time : duration of the stimuli
  #time : dictionnary with all the time specificity of our model
  #simu_param : Dictionary of the specifity of simulation
  # param : general param of the model
  # Dic_spec : specific param to the simulation
  # If_function : Mapping between DeltaS,Confidence and If
  # Ir : value of reset current
  # If : max value of inhibitory reward current

  # Output :
  # Reaction time with half direction selectivity
  # Return 1 or 2 r1 r2
  # Save param et time
  # return timing spike (for example fixed at the value in Wang)
  # return the true list of decision, and maybe rest distribution ????
  # Wrong or Right list
  # S1 S2 R1 R2
  # Confidence list : Just note that the first term do not have any sense : Just init
  # of the list in order to have the same length


  # FI curve parameters - adapted from Wanf fit
  a=fixed_param["a"]
  b=fixed_param["b"]
  d=fixed_param["d"]
  th = 20


  # time variables of the  simulation

  dt=simu_param["dt"]
  time_wind=simu_param["time_wind"] # mean window of the variables

  slide_wind = simu_param["slide_wind"]

  false_time=param["false_time"]
    false_timeIf=param["false_timeIf"]

  deltaSstim=[]

  res=zeros(length(coh_temp))
  reaction_time=zeros(length(coh_temp))
  diffS=zeros(length(coh_temp))
  S1=zeros(length(coh_temp))
  S2=zeros(length(coh_temp))
 R1=zeros(length(coh_temp))
  R2=zeros(length(coh_temp))

 conflist=zeros(length(coh_temp))

  #---- Intialise and vectorise variables to be used in loops below ------

  # Final Variables



  Istim1 = []
  # Temporary variables

  Isyn1=[0.0]
  Isyn2=[0.0]
  Trest=param["T_rest"]
  s1=[0.1]
  s2=[0.1]
  phi1=[0.0]
  phi2=[0.0]
  I_eta1=[param["noise_amp"] * randn()]
  I_eta2=[param["noise_amp"] * randn()]
  nu1=[2.0]
  nu2=[2.0]
  I=[0.0]

  k=0 # variable for mod 4 in the mean
  th=false
  result=true
  for j=1:length(coh_temp) # loop on trials number
    coh=coh_temp[j]
    t=0

    # reinit
    Isyn1=[Isyn1[end]]
    Isyn2=[Isyn2[end]]
    s1=[s1[end]]
    s2=[s2[end]]
    phi1=[phi1[end]]
    phi2=[phi2[end]]
    I_eta1=[I_eta1[end]]
    I_eta2=[I_eta2[end]]
    nu1=[nu1[end]]
    nu2=[nu2[end]]
    k=0
    confidence = conf_function(nu1[end],nu2[end])
    conflist[j]=confidence
    while th==false # while we do not cross the threshold
      t=t+dt
      k=k+1

      if t<Trest
        I0E1=0.3255 -Ir*exp(-t/false_time) - If_function(IF,0,9,confidence)*exp(-t/false_timeIf) #100ms time duration
        I0E2=0.3255 -Ir*exp(-t/false_time)- If_function(IF,0,9,confidence)*exp(-t/false_timeIf) #100ms time duration
        #push!(I,I0E1)
        I_stim1=0
        I_stim2=0
        cht=0
      elseif t<Trest+dt/2 && t>Trest-dt/2
        I0E1=0.3255#dic_param["I0E1"]
        #push!(I,I0E1) #assuming input on 1 and 2 are the same
        I0E2=0.3255#dic_param["I0E2"]
        I_stim1 = param["JAext"]*param["mu0"]*(1+coh/100)
        I_stim2 = param["JAext"]*param["mu0"]*(1-coh/100)
        cht=coh
        #push!(deltaSstim,abs(s1[end] - s2[end]))

      else
        I0E1=0.3255#dic_param["I0E1"]
        #push!(I,I0E1) #assuming input on 1 and 2 are the same
        I0E2=0.3255#dic_param["I0E2"]
        I_stim1 = param["JAext"]*param["mu0"]*(1+coh/100)
        I_stim2 = param["JAext"]*param["mu0"]*(1-coh/100)
        cht=coh

      end


      push!(Isyn1,param["JN11"]*s1[end]-param["JN12"]*s2[end]+I_stim1 + I_eta1[end])
      push!(Isyn2,param["JN22"]*s2[end]-param["JN21"]*s1[end]+I_stim2 + I_eta2[end])

      push!(phi1,(a*Isyn1[end]-b)/(1-exp(-d*(a*Isyn1[end]-b))))
      push!(phi2,(a*Isyn2[end]-b)/(1-exp(-d*(a*Isyn2[end]-b))))
      #---- Dynamical equations -------------------------------------------

      # Mean NMDA-mediated synaptic dynamics updating
      s1t = s1[end] + dt*(-s1[end]/param["Tnmda"] + (1-s1[end])*param["gamma"]*nu1[end]/1000);
      push!(s1,s1t)
      s2t = s2[end] + dt*(-s2[end]/param["Tnmda"] + (1-s2[end])*param["gamma"]*nu2[end]/1000);
      push!(s2,s2t)
      # Ornstein-Uhlenbeck generation of noise in pop1 and 2
      I_eta1t = I_eta1[end] + (dt/param["Tampa"])*(I0E1-I_eta1[end]) + sqrt(dt/param["Tampa"])*param["noise_amp"]*randn() ;
      I_eta2t = I_eta2[end] + (dt/param["Tampa"])*(I0E2-I_eta2[end]) + sqrt(dt/param["Tampa"])*param["noise_amp"]*randn() ;
      push!(I_eta1,I_eta1t)
      push!(I_eta2,I_eta2t)
      # To ensure firing rates are always positive. Large noise amplitude
      # may result in unwanted negative values
      if phi1[end] < 0
        push!(nu1,0)
        phi1[end] = 0;
      else
        push!(nu1,phi1[end])
      end;
      if phi2[end] < 0
        push!(nu2,0)
        phi2[end] = 0;
      else
        push!(nu2,phi2[end])
      end;

      if mod(k,4)==0

        # Test the crossing of a threshold
        if mean(nu1[k-3:k])>20 && t>100
          th=true
          if coh>=0
            res[j]=1
            result=true
          else
           res[j]=-1
            result=false
          end

        elseif mean(nu2[k-3:k])[end]>20 && t>100
          th = true
          if coh>=0
           res[j]=-1
            result=false
          else
            res[j]=1
            result=true
          end

        end
      end
    end

   reaction_time[j]=t-Trest
    diffS[j]=s2[end]-s1[end]
    R1[j]=nu1[end]
    R2[j]=nu2[end]
 S1[j]=s1[end]
    S2[j]=s2[end]
    th=false
  end


  result=DataFrame()
  result[:RTs]=reaction_time
  result[:NotError]=res
  result[:Distri]=coh_temp
  result[:DiffS]=diffS
  result[:R1]=R1
  result[:R2]=R2
 result[:S1]=S1
  result[:S2]=S2
 result[:confidence]=conflist
  result[:LastTrial]=vcat([0],res[1:end-1])

  return result
end



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

  k=0

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

      if t<Trest

        I0E1=I0E1p -Ir*exp(-t/falset)
        I0E2=I0E2p -Ir*exp(-t/falset)

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
