# This script does the figures conrresponding to the fit of the decision-making model.

rand(123456789)
function list_elements(l,n=4)
      a = Float64[]
        i=1
        while (length(a) <n)
           if any(a.==abs(l[i]))

            else
                push!(a,abs(l[i]))
            end
            i=i+1
        end
        return sort(a)
    end

include(".\\decision_making.jl")
include(".\\various_functions.jl")
using CSV,Query,DataFrames
using Distributions, StatsBase


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
    "false_time"=>150.0,#200
    "false_timeIf"=>500.0,
    "Threshold"=>12.0
    )
    df= CSV.read(".\\data\\Manip1.csv")

fixed_param=Dict("a"=>270.0,"b"=>108.0,"d"=>0.1540)
simu_param=Dict("dt"=>0.5,"time_wind"=>4.0,"slide_wind"=>4.0) #4=2/dt
Name=[1,3,5,6,7,9]
Threshold=[10.78,12.69,14.07,12.80,10.05,12.89]
Drift=[35.0,35.0,35.0,35.0,35.0,35.0]

c02=[1.22,1.57,0.93,0.87,0.76,1.10]
c05=[3.88,5.39,4.94,2.82,3.88,4.2]

c08=[7.86,8.25,6.65,2.94,5.85,9.1]
c16=[13.5,12.41,13.95,6.82,13.25,13.9]
Ir=0.033
If=0.0
Variable_mu = [0.402, 0.369,
 0.476, 0.422, 0.425,  0.415]
Variable_gamma = [159.301,  372,
 202, 293, 335,  308]
Variable_lambda =[14.37,  15.23,  6.68,
9.92, 6.97,  9.52]
vmufinal=[0.21893340000000003,
0.188122,
0.2672158,
0.23036479999999998,
0.1653072,
0.16317199999999998
]
using Bootstrap
using Plots
plotlyjs()


function theta_correspondance(angle,name)
    # it returns depending on the block 1 or 3
    # careful angle in radian
    if name == "Sujet 1"

        return (1.46*angle-0.01*angle*angle,2.5*angle +0.092*angle*angle)

    elseif name == "Sujet 2"

        return (1.92*angle-0.068*angle*angle,2.78*angle-0.098*angle*angle)

    elseif name =="Sujet 3"

        return (1.38*angle,1.67*angle)

    elseif name == "Sujet 4"
        return (0.71*angle-0.004*angle*angle,0.50*angle+0.15*angle*angle)

    elseif name == "Sujet 5"

        return (1.02*angle+0.029*angle*angle,0.48*angle+0.34*angle*angle)

    elseif name == "Sujet 6"

        return (1.45*angle,2.71*angle)

    else
        return (0,0)
    end
end

pyplot()

upscale = 2 #8x upscaling in resolution
fntsm = Plots.font("Times New Roman", 15.0*upscale)
fntlg = Plots.font("Times New Roman", 60.0*upscale)
default(titlefont=fntlg,guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
fig=plot(layout=(4,6),legend=false,size=(800*2*upscale,upscale*600*2))
fig2=plot(layout=6,legend=false,size=(800*upscale,upscale*600))
ylabel!(fig[1],"Response times (s)")
ylabel!(fig[7],"Accuracy")
ylabel!(fig[13],"Response times (s)")
ylabel!(fig[19],"Accuracy")
for n=1:6
rtsuj=Float64[]
rtsujup=Float64[]

rtsujdown=Float64[]

accsujup=Float64[]

accsujdown=Float64[]

#suj 2 fait un truc trop chelou en RTs
dfsuj1rt=df[:Rt1][df[:Name].==Name[n]]
dfsuj1av=df[:AngleValue][df[:Name].==Name[n]]
dfsuj1acc=df[:Acc][df[:Name].==Name[n]]
accsuj=Float64[]

rt=df[:Rt1][df[:Name].==Name[n]]

    rtf = zeros(length(rt))
    conff=zeros(length(rt))
    accf = zeros(length(rt))

    for i = 1:length(rt)
    rtf[i] = rt[i]
    conff[i] = abs(dfsuj1av[i])
    accf[i]=dfsuj1acc[i]
    end
    meanconf = zeros(10)
    meanconf2 = zeros(10)

    n_boot = 10000

     ## Basic bootstrap


     ## basic CI
     for a in list_elements(dfsuj1av)


         bs1 = bootstrap(filter(!isnan,rtf[conff.==a]), mean, BasicSampling(n_boot))
         cil = 0.99;
     bci1 = Bootstrap.ci(bs1, BasicConfInt(cil));
    push!(rtsuj,original(bs1)[1])#0.5+0.5.*mean(rtf[conff.==i])
    push!(rtsujdown,bci1[1][1]-bci1[1][2])
    push!(rtsujup,bci1[1][3]-bci1[1][1])


    bs1 = bootstrap(filter(!isnan,accf[conff.==a]), mean, BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.ci(bs1, BasicConfInt(cil));
    push!(accsuj,0.5+0.5.*original(bs1)[1])#0.5+0.5.*mean(rtf[conff.==i])
    push!(accsujdown,0.5*bci1[1][1]-0.5*bci1[1][2])
    push!(accsujup,0.5*bci1[1][3]-0.5*bci1[1][1])



    end



#cohvalue= [-c16[n], -c08[n],-c05[n], -c02[n], c02[n],c05[n], c08[n], c16[n]]

param["Threshold"]=Threshold[n]
param["mu0"]=Drift[n]
#weights = [0.05, 0.1, 0.15, 0.2, 0.2, 0.15,0.1,0.05]
  #distri = sample(cohvalue,Weights(weights),10000)

cohvalue = zeros(20)
for i=1:10
    (c,c2) = theta_correspondance(0.1*i*2.0*2*pi,"Sujet $n")
    cohvalue[10-i+1] = -c
    cohvalue[10+i] = c
end
weights= 1.0./20.0.*ones(20)
distri = sample(cohvalue,Weights(weights),200000)

res= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If, param,fixed_param,simu_param)

rt=Float64[]
acc=Float64[]
for a in list_elements(res[:Distri],10)
  push!(rt,mean(res[:RTs][abs.(res[:Distri]).==a])./1000)
  push!(acc,mean(res[:NotError][abs.(res[:Distri]).==a]))

end

cohvalue= [-c16[n], -c08[n],-c05[n], -c02[n], c02[n],c05[n], c08[n], c16[n]]

param["Threshold"]=Threshold[n]
param["mu0"]=Drift[n]
weights = [0.05, 0.1, 0.15, 0.2, 0.2, 0.15,0.1,0.05]
distri = sample(cohvalue,Weights(weights),10000)

res= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If, param,fixed_param,simu_param)

plot!(fig[n],0.1:0.2:2.0,rt.+ (mean(dfsuj1rt) .- mean(res[:RTs]./1000)),linewidth=4*upscale,linestyle=:solid,color=:darkblue,alpha=0.9)#marker=:diamond,markersize=6*upscale

scatter!(fig[n],[0.2,0.5,0.8,1.6],rtsuj
,yerr=(rtsujdown,rtsujup),markersize=8*upscale,marker=:diamond,color=:red)

ylims!(fig[n],(0.37,0.72))
xlims!(fig[n],(0,1.7))


plot!(fig[n+6],0.1:0.2:2.0,0.5.+0.5.*acc,linewidth=4*upscale,color=:darkblue,linestyle=:solid,alpha=0.9)#marker=:diamond,markersize=6*upscale,

scatter!(fig[n+6],[0.2,0.5,0.8,1.6],accsuj,yerr=(accsujdown,accsujup),markersize=8*upscale,marker=:diamond,color=:red)
title!(fig[n],"Participant $n")
#title!(fig[n+6],"Participant $n")


ylims!(fig[n+6],(0.5,1.02))
xlims!(fig[n+6],(0,1.7))


plot!(fig2[n],0.1:0.2:2.0,0.5.+0.5.*acc,linewidth=4*upscale,color=:darkblue,linestyle=:solid,alpha=0.9)#marker=:diamond,markersize=6*upscale,

scatter!(fig2[n],[0.2,0.5,0.8,1.6],accsuj,yerr=(accsujdown,accsujup),markersize=8*upscale,marker=:diamond,color=:red)

ylims!(fig2,(0.5,1.02))
xlims!(fig2,(0,1.7))

print(Variable_mu[n] - mean(res[:RTs])./1000)
print("\n")
end



fixed_param=Dict("a"=>270.0,"b"=>108.0,"d"=>0.1540)
simu_param=Dict("dt"=>0.5,"time_wind"=>4.0,"slide_wind"=>4.0) #4=2/dt
Name=[1,3,5,6,7,9]
Threshold=[13.08,13.70,14.95,12.96,12.55,14.65]
Drift = [35.0,35.0,35.0,35.0,35.0,35.0]
Threshold=[13.08,13.70,14.95,12.96,12.55,14.65]
c02=[3.2,2.68,1.32,2.35,1.91,3.9]
c08=[14.95,11.81,7.1,8.15,13.49,14.1]
c16=[34.5,18.04,17.5,30.01,50.21,28.1]
Variable_mu = [0.422602, 0.390527,
 0.510844, 0.529926, 0.436503,  0.433333]
Variable_gamma = [3599.17,  657.64,
 323.661, 699.619, 1147.85,  1268.4]
Variable_lambda =[6.98586,  10.6735,  5.06607,
3.11227, 5.63951,  6.51043]

function list_elements(l)
      a = Float64[]
        i=1
        while (length(a) <3)
           if any(a.==abs(l[i]))

            else
                push!(a,abs(l[i]))
            end
            i=i+1
        end
        return sort(a)
    end
vmufinal=[0.17476,
0.10728,
0.193,
0.268,
0.203,
0.148]

using Bootstrap
using Plots
df= CSV.read(".\\data\\Manip3.csv")


for n=1:6
rtsuj=Float64[]
rtsujup=Float64[]

rtsujdown=Float64[]

accsujup=Float64[]

accsujdown=Float64[]

dfsuj1rt=df[:Rt1][df[:Name].==Name[n]]
dfsuj1av=df[:AngleValue][df[:Name].==Name[n]]
dfsuj1acc=df[:Acc][df[:Name].==Name[n]]
accsuj=Float64[]

rt=df[:Rt1][df[:Name].==Name[n]]

    rtf = zeros(length(rt))
    conff=zeros(length(rt))
    accf = zeros(length(rt))

    for i = 1:length(rt)
    rtf[i] = rt[i]
    conff[i] = abs(dfsuj1av[i])
    accf[i]=dfsuj1acc[i]
    end
    meanconf = zeros(10)
    meanconf2 = zeros(10)

    n_boot = 10000

     ## Basic bootstrap


     ## basic CI
     for a in list_elements(dfsuj1av,3)


         bs1 = bootstrap(filter(!isnan,rtf[conff.==a]), mean, BasicSampling(n_boot))
         cil = 0.99;
     bci1 = Bootstrap.ci(bs1, BasicConfInt(cil));
    push!(rtsuj,original(bs1)[1])
    push!(rtsujdown,bci1[1][1]-bci1[1][2])
    push!(rtsujup,bci1[1][3]-bci1[1][1])


    bs1 = bootstrap(filter(!isnan,accf[conff.==a]), mean, BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.ci(bs1, BasicConfInt(cil));
    push!(accsuj,0.5+0.5.*original(bs1)[1])
    push!(accsujdown,0.5*bci1[1][1]-0.5*bci1[1][2])
    push!(accsujup,0.5*bci1[1][3]-0.5*bci1[1][1])



    end



    param["Threshold"]=Threshold[n]
    param["mu0"]=Drift[n]
    #weights = [0.05, 0.1, 0.15, 0.2, 0.2, 0.15,0.1,0.05]
      #distri = sample(cohvalue,Weights(weights),10000)

    cohvalue = zeros(20)
    for i=1:10
        (c,c2) = theta_correspondance(0.1*i*2.0*2*pi,"Sujet $n")
        cohvalue[10-i+1] = -c2
        cohvalue[10+i] = c2
    end
    weights= 1.0./20.0.*ones(20)
    distri = sample(cohvalue,Weights(weights),200000)

    res= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If, param,fixed_param,simu_param)

    rt=Float64[]
    acc=Float64[]
    for a in list_elements(res[:Distri],10)
      push!(rt,mean(res[:RTs][abs.(res[:Distri]).==a])./1000)
      push!(acc,mean(res[:NotError][abs.(res[:Distri]).==a]))

    end

cohvalue= [-c16[n], -c08[n], -c02[n], c02[n], c08[n], c16[n]]

param["Threshold"]=Threshold[n]
param["mu0"]=Drift[n]
  weights = [0.12, 0.18, 0.2, 0.2, 0.18,0.12]

  distri = sample(cohvalue,Weights(weights),10000)
res= sim_mf_decisionmaking_pattern_fastv2(distri,Ir,If, param,fixed_param,simu_param)


plot!(fig[n+12],0.1:0.2:2.0,rt.+ (mean(dfsuj1rt) .- mean(res[:RTs]./1000)),linewidth=4*upscale,linestyle=:solid,color=:darkblue,alpha=0.9)#marker=:diamond,markersize=6*upscale

scatter!(fig[n+12],[0.2,0.8,1.6],rtsuj
,yerr=(rtsujdown,rtsujup),markersize=8*upscale,marker=:diamond,color=:red)

ylims!(fig[n+12],(0.41,0.98))
xlims!(fig[n+12],(0,1.7))
#title!(fig[n+12],"Participant $n")
plot!(fig[n+18],0.1:0.2:2.0,0.5.+0.5.*acc,linewidth=4*upscale,color=:darkblue,linestyle=:solid,alpha=0.9)#marker=:diamond,markersize=6*upscale,

scatter!(fig[n+18],[0.2,0.8,1.6],accsuj,yerr=(accsujdown,accsujup),markersize=8*upscale,marker=:diamond,color=:red)

ylims!(fig[n+18],(0.5,1.02))
xlims!(fig[n+18],(0,1.7))
#title!(fig[n+18],"Participant $n")


plot!(fig2[n],0.1:0.2:2.0,0.5.+0.5.*acc,linewidth=4*upscale,color=:darkblue,linestyle=:solid,alpha=0.9)#marker=:diamond,markersize=6*upscale,

scatter!(fig2[n],[0.2,0.8,1.6],accsuj,yerr=(accsujdown,accsujup),markersize=8*upscale,marker=:diamond,color=:red)

ylims!(fig2,(0.5,1.02))
xlims!(fig2,(0,1.7))

print(Variable_mu[n] - mean(res[:RTs])./1000)
print("\n")
end


fig

fig2
