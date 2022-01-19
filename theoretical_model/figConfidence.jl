####
# RT and Acc fig
######

rand(123456789)
Name=[1,3,5,6,7,9]
using Plots
gr()
using Statistics
using DataFrames
using CSV,Query

df = CSV.read(".\\data\\Manip3.csv")

pyplot()
gr()


fig=plot(layout=6,size=(800,600))


meanconf = zeros(10)
meanconf2 = zeros(10)
errconfdown=zeros(10)
errconfup=zeros(10)
Name=[1,3,5,6,7,9]
using Bootstrap
n=0
i=0
for n=1:6
    rt=df[:Rt1][df[:Name].==Name[n]].*1000
    conf = df[:Resp2][df[:Name].==Name[n]]

    rtf = zeros(length(rt))
    conff=zeros(length(rt))

    for i = 1:length(rt)
        rtf[i] = rt[i]
        conff[i] = conf[i]
    end
    meanconf = zeros(10)
    meanconf2 = zeros(10)
    for i=0:9
        n_boot = 10000

        ## Basic bootstrap


        ## basic CI
        if n==4 && i==0
        elseif n==1 && i==9
        elseif n==1 && i==0
        elseif n==2 && i==0
        elseif n==3 && i==0
        elseif n==4 && i==1
        else
            bs1 = bootstrap(filter(isfinite,rtf[conff.==i]), mean, BasicSampling(n_boot))
            cil = 0.99;
            bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
            meanconf2[i+1] =original(bs1)[1]
            errconfdown[i+1]=bci1[1][1]-bci1[1][2]
            errconfup[i+1]=bci1[1][3]-bci1[1][1]

        end
    end

    datar= CSV.read(".\\data\\RT$n.csv")

    colorlist=[:blue,:red,:darkorange,:green,:violet,:black]

    dataRT=DataFrame()
    dataRT[:C0]=collect(Missings.replace(datar[:C0], 1))
    dataRT[:C1]=collect(Missings.replace(datar[:C1], 1))
    dataRT[:C2]=collect(Missings.replace(datar[:C2], 1))
    dataRT[:C3]=collect(Missings.replace(datar[:C3], 1))
    dataRT[:C4]=collect(Missings.replace(datar[:C4], 1))
    dataRT[:C5]=collect(Missings.replace(datar[:C5], 1))
    dataRT[:C6]=collect(Missings.replace(datar[:C6], 1))
    dataRT[:C7]=collect(Missings.replace(datar[:C7], 1))
    dataRT[:C8]=collect(Missings.replace(datar[:C8], 1))
    dataRT[:C9]=collect(Missings.replace(datar[:C9], 1))

    meanlist=[NaN,
    NaN,
    mean(dataRT[:C2][isfinite.(dataRT[:C2])]),
    mean(dataRT[:C3][isfinite.(dataRT[:C3])]),
    mean(dataRT[:C4][isfinite.(dataRT[:C4])]),
    mean(dataRT[:C5][isfinite.(dataRT[:C5])]),
    mean(dataRT[:C6][isfinite.(dataRT[:C6])])
    ,mean(dataRT[:C7][isfinite.(dataRT[:C7])]),
    mean(dataRT[:C8][isfinite.(dataRT[:C8])]),
    NaN]

    varlist=[var(dataRT[:C0][isfinite.(dataRT[:C0])]),
    var(dataRT[:C1][isfinite.(dataRT[:C1])]),
    var(dataRT[:C2][isfinite.(dataRT[:C2])])
    ,var(dataRT[:C3][isfinite.(dataRT[:C3])]),
    var(dataRT[:C4][isfinite.(dataRT[:C4])]),
    var(dataRT[:C5][isfinite.(dataRT[:C5])]),
    var(dataRT[:C6][isfinite.(dataRT[:C6])])
    ,var(dataRT[:C7][isfinite.(dataRT[:C7])]),
    var(dataRT[:C8][isfinite.(dataRT[:C8])]),
    var(dataRT[:C9][isfinite.(dataRT[:C9])])]

    ribbonup=zeros(10)
    ribbondown=zeros(10)
    n_boot = 1000
    print(n)
    if n==4
    elseif n==3
    elseif n==1

    else

        bs1 = bootstrap( mean,filter(isfinite,dataRT[:C0][1:20]), BasicSampling(n_boot))
        cil = 0.99;
        bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
        meanlist[1]=bci1[1][1]
        ribbonup[1]=bci1[1][3]-bci1[1][1]
        ribbondown[1]=bci1[1][1]-bci1[1][2]
    end
    ############################################""
    if n==4

    elseif n==1

    else
        bs1 = bootstrap(mean,filter(isfinite,dataRT[:C1][1:20]),  BasicSampling(n_boot))
        cil = 0.99;
        bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
        meanlist[2]=bci1[1][1]
        ribbonup[2]=bci1[1][3]-bci1[1][1]
        ribbondown[2]=bci1[1][1]-bci1[1][2]
    end
    ############################################""

    bs1 = bootstrap(mean,filter(isfinite,dataRT[:C2][1:20]),  BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[3]=bci1[1][1]
    ribbonup[3]=bci1[1][3]-bci1[1][1]
    ribbondown[3]=bci1[1][1]-bci1[1][2]

    ############################################""

    bs1 = bootstrap( mean,filter(isfinite,dataRT[:C3][1:20]), BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[4]=bci1[1][1]
    ribbonup[4]=bci1[1][3]-bci1[1][1]
    ribbondown[4]=bci1[1][1]-bci1[1][2]
    ############################################""

    bs1 = bootstrap( mean,filter(isfinite,dataRT[:C4][1:20]), BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[5]=bci1[1][1]
    ribbonup[5]=bci1[1][3]-bci1[1][1]
    ribbondown[5]=bci1[1][1]-bci1[1][2]
    ############################################""

    bs1 = bootstrap(mean,filter(isfinite,dataRT[:C5][1:20]),  BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[6]=bci1[1][1]
    ribbonup[6]=bci1[1][3]-bci1[1][1]
    ribbondown[6]=bci1[1][1]-bci1[1][2]
    ############################################""

    bs1 = bootstrap( mean,filter(isfinite,dataRT[:C6][1:20]), BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[7]=bci1[1][1]
    ribbonup[7]=bci1[1][3]-bci1[1][1]
    ribbondown[7]=bci1[1][1]-bci1[1][2]
    ############################################""

    bs1 = bootstrap(mean,filter(isfinite,dataRT[:C7][1:20]),  BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[8]=bci1[1][1]
    ribbonup[8]=bci1[1][3]-bci1[1][1]
    ribbondown[8]=bci1[1][1]-bci1[1][2]
    ############################################""

    bs1 = bootstrap(mean,filter(isfinite,dataRT[:C8][1:20]),  BasicSampling(n_boot))
    cil = 0.99;
    bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
    meanlist[9]=bci1[1][1]
    ribbonup[9]=bci1[1][3]-bci1[1][1]
    ribbondown[9]=bci1[1][1]-bci1[1][2]

    ############################################""
    if n==1
    else
        bs1 = bootstrap(mean,filter(isfinite,dataRT[:C9][1:20]),  BasicSampling(n_boot))
        cil = 0.99;
        bci1 = Bootstrap.confint(bs1, BasicConfInt(cil));
        meanlist[10]=bci1[1][1]
        ribbonup[10]=bci1[1][3]-bci1[1][1]
        ribbondown[10]=bci1[1][1]-bci1[1][2]
    end
    ############################################""

    #ribbon=(ribbondown,ribbonup)
    plot!(fig[n],[0,1,2,3,4,5,6,7,8,9],meanlist,ribbon=(ribbondown,ribbonup),
    linewidth=6,color=colorlist[n])
    print(meanlist)
    scatter!(fig[n],[0,1,2,3,4,5,6,7,8,9],meanconf2,yerr=(errconfdown,errconfup),color=colorlist[n]
    ,markersize=8)
    ylims!((300,1250))
end


plot!(legend=false)
fig




savefig(fig,".\\confidencepersubjectRT.svg")
