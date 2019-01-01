import numpy as np
from matplotlib.pyplot import *
from statsmodels.tsa.filters.hp_filter import hpfilter
import json


# returns index starting from 1
def gridlookup2(x0,xgrid):

    nx = xgrid.shape[0]

    ix = 0
    for jx in range(nx):
        if x0<=xgrid[jx]:
            break
        ix = ix+1

    ix = min(max(1,ix),nx-1)

    return ix


def calcDH(knotsm,mpmat0,pmat0,izvec,Kvec,statflag):

    tott = izvec.shape[0]

    mp = Kvec[0]
    mpvec = np.zeros(tott)
    pvec = np.zeros(tott)

    for tt in range(tott):

        if statflag==1:
            mnow = Kvec[tt]
        else:
            mnow = mp

        # NOTE: im and iz are indices starting from 1, not 0 as in Fortran
        iz = izvec[tt]
        im = gridlookup2(np.log(mnow),np.log(knotsm))
        # print(np.array([im,iz]))
        wm = np.log(knotsm[im]/mnow)/np.log(knotsm[im]/knotsm[im-1])
        mp = np.exp(wm*np.log(mpmat0[im-1,iz-1]) + (1-wm)*np.log(mpmat0[im,iz-1]))
        cnow = 1/np.exp(wm*np.log(pmat0[im-1,iz-1]) + (1-wm)*np.log(pmat0[im,iz-1]))

        mpvec[tt] = mp
        pvec[tt] = 1/cnow

    return mpvec,pvec


def calcaccurat(json,simT,drop,filename):


    knotsm = np.array(json["input"]["knotsm"])
    izvec = np.array(json["input"]["izvec"])
    mpmat0 = np.array(json["output"]["mpmat0"])
    pmat0 = np.array(json["output"]["pmat0"])
    Kvec = np.array(json["output"]["Kvec"])
    Kpvec = np.array(json["output"]["Kpvec"])
    Cvec = np.array(json["output"]["Cvec"])

    nm = np.array(json["input"]["knotsm"]).shape[0]
    nz = np.array(json["input"]["Gz"]).shape[0]
    mpmat = np.zeros((nm,nz))
    pmat = np.zeros((nm,nz))

    for iz in range(nz):
        for im in range(nm):

            index = nm*iz+im
            # print(index)
            mpmat[im,iz] = mpmat0[index]
            pmat[im,iz] = pmat0[index]

    # dynamic forecasts
    mpvec0,pvec0 = calcDH(knotsm, mpmat, pmat, izvec, Kvec, 0)
    # static forecasts
    mpvec1,pvec1 = calcDH(knotsm, mpmat, pmat, izvec, Kvec, 1)

    izvec = izvec[drop:simT+drop]
    Kvec = Kvec[drop:simT+drop]
    Kpvec = Kpvec[drop:simT+drop]
    Cvec = Cvec[drop:simT+drop]
    mpvec0 = mpvec0[drop:simT+drop]
    pvec0 = pvec0[drop:simT+drop]
    mpvec1 = mpvec1[drop:simT+drop]
    pvec1 = pvec1[drop:simT+drop]

    DHmax = np.zeros((nz,2))
    DHmean = np.zeros((nz,2))
    RMSE = np.zeros((nz,2))
    Rsquared = np.zeros((nz,2))

    for iz in range(nz):

        DHmax[iz,0] = 100*np.max(np.abs(np.log(mpvec0[izvec==iz+1]/Kpvec[izvec==iz+1])))
        DHmax[iz,1] = 100*np.max(np.abs(np.log(pvec0[izvec==iz+1]*Cvec[izvec==iz+1])))
        DHmean[iz,0] = 100*np.sum(np.abs(np.log(mpvec0[izvec==iz+1]/Kpvec[izvec==iz+1])))/simT
        DHmean[iz,1] = 100*np.sum(np.abs(np.log(pvec0[izvec==iz+1]*Cvec[izvec==iz+1])))/simT

        y0 = np.log(Kpvec[izvec==iz+1])
        y0 = np.reshape(y0,(y0.shape[0],1))
        # print(y0.shape)
        y1 = -np.log(Cvec[izvec==iz+1])
        y1 = np.reshape(y1,(y1.shape[0],1))
        # how to obtain row vectors?
        y = np.hstack([y0,y1])
        # below does not work
        # y = np.hstack([np.log(Kpvec[izvec==iz+1]),-np.log(Cvec[izvec==iz+1].T)])
        # print(y.shape)
        ymean = np.mean(y,axis=0)
        # print(ymean.shape)
        ydev = y-ymean
        # ydev0 = y[:,0]-ymean[0]
        # ydev1 = y[:,1]-ymean[1]
        # print(ydev.shape)
        e0 = np.log(mpvec1[izvec==iz+1]/Kpvec[izvec==iz+1])
        e1 = np.log(pvec1[izvec==iz+1]*Cvec[izvec==iz+1])

        RMSE[iz,0] = 100*np.sqrt(np.mean(np.power(e0,2)))
        RMSE[iz,1] = 100*np.sqrt(np.mean(np.power(e1,2)))
        Rsquared[iz,0] = 1-np.dot(e0.T,e0)/np.dot(ydev[:,0].T,ydev[:,0])
        Rsquared[iz,1] = 1-np.dot(e1.T,e1)/np.dot(ydev[:,1].T,ydev[:,1])

    data = np.hstack([DHmax, DHmean, RMSE, Rsquared])
    print(data)
    np.savetxt(filename + ".csv", data, delimiter=",", fmt='%.4f')

    return mpvec0, pvec0, mpvec1, pvec1


def calcstochss(json,irfdrop,filename):

    Kirvec = np.array(json["irf"]["Kvec"])
    Zirvec = np.array(json["irf"]["Zvec"])
    Yirvec = np.array(json["irf"]["Yvec"])
    Iirvec = np.array(json["irf"]["Ivec"])
    Nirvec = np.array(json["irf"]["Nvec"])
    Cirvec = np.array(json["irf"]["Cvec"])

    ikirvec = np.zeros((Kirvec.shape[0],7))
    ikirvec[:,0] = np.array(json["irf"]["ikmean"])
    ikirvec[:,1] = np.array(json["irf"]["ikstddev"])
    ikirvec[:,2] = np.array(json["irf"]["ikinaction"])
    ikirvec[:,3] = np.array(json["irf"]["ikspikepos"])
    ikirvec[:,4] = np.array(json["irf"]["ikspikeneg"])
    ikirvec[:,5] = np.array(json["irf"]["ikpos"])
    ikirvec[:,6] = np.array(json["irf"]["ikneg"])

    mnow = Kirvec[irfdrop-1]
    ynow = Yirvec[irfdrop-1]
    nnow = Nirvec[irfdrop-1]
    cnow = Cirvec[irfdrop-1]

    data1 = np.array([mnow/ynow, nnow, 1/cnow, 0, 0, 0, 0])
    data2 = ikirvec[irfdrop-1,:]
    data = np.vstack([data1,data2])
    print(data)
    np.savetxt(filename + ".csv", data, delimiter=",", fmt='%.4f')

    return Yirvec,Iirvec,Nirvec,Kirvec,Cirvec,Zirvec


def calcaggstat(json,simT,drop,filename):

    Kvec = np.array(json["output"]["Kvec"])
    Zvec = np.array(json["output"]["Zvec"])
    Yvec = np.array(json["output"]["Yvec"])
    Ivec = np.array(json["output"]["Ivec"])
    Nvec = np.array(json["output"]["Nvec"])
    Cvec = np.array(json["output"]["Cvec"])
    Kpvec = np.array(json["output"]["Kpvec"])
    Xvec = np.array(json["output"]["Xvec"])

    Kvec = Kvec[drop:simT + drop]
    Zvec = Zvec[drop:simT + drop]
    Yvec = Yvec[drop:simT + drop]
    Ivec = Ivec[drop:simT + drop]
    Nvec = Nvec[drop:simT + drop]
    Cvec = Cvec[drop:simT + drop]
    Kpvec = Kpvec[drop:simT + drop]

    # BUG in the use of np.hstack and hpfilter???
    y = np.log(np.hstack((Yvec, Cvec, Ivec, Nvec, Kvec, Zvec)))
    # print(y.shape)
    yd, yf = hpfilter(y, 100)
    yd = np.reshape(yd, (simT, 6), 'F')
    # yf = np.reshape(yf, (simT, 6), 'F')
    yd = yd[8:simT - 8, :]
    std0 = np.std(yd, axis=0)
    print('Standard deviation')
    # print('Y C I N K Z')
    data1 = np.array([std0[0]*100, std0[1]/std0[0], std0[2]/std0[0], std0[3]/std0[0],
                    std0[4]/std0[0], std0[5]/std0[0]])
    print(data1)

    corr0 = np.corrcoef(yd, rowvar=False)
    print('Output correlation')
    data2 = corr0[0, 0:6]
    print(data2)

    print('Aggregate investment rate')
    xmom = np.zeros(4)
    xmom[0] = np.sum(Xvec[drop:simT + drop]) / simT
    xmom[1] = np.sum(np.power(Xvec[drop:simT + drop] - xmom[0], 2)) / simT
    xmom[2] = np.sum(np.power(Xvec[drop:simT + drop] - xmom[0], 3)) / simT
    xmom[3] = np.sum(np.power(Xvec[drop:simT + drop] - xmom[0], 4)) / simT
    # persistence
    rho1 = np.dot(Xvec[drop:simT + drop - 1] - xmom[0], Xvec[drop + 1:simT + drop] - xmom[0])
    rho1 = rho1 / np.sum(np.power(Xvec[drop:simT + drop - 1] - xmom[0], 2))
    # standard deviation
    sig1 = np.power(xmom[1], 0.5)
    # skewness
    g1 = xmom[2] / np.power(xmom[1], 1.5)
    # excess kurtosis
    g2 = xmom[3] / np.power(xmom[1], 2) - 3.0
    data3 = np.array([rho1, sig1, g1, g2, 0, 0])
    print(data3)

    data = np.vstack([data1,data2,data3])
    np.savetxt(filename + ".csv", data, delimiter=",", fmt='%.4f')

    return Yvec,Ivec,Nvec,Kvec,Cvec,Kpvec,Zvec


########################################
# main part of the script
########################################
# for KT 2008
jsonKS = json.load(open('./results_extend_KS.json','r'))
jsonXpa = json.load(open('./results_extend_Xpa.json','r'))
fname1 = "simall_extend"
fname2 = "irfall_extend"
fname3 = "DHall_extend"
fname4 = "eptime_extend"
# for KT 2003
# jsonKS = json.load(open('./results_traditional_KS.json','r'))
# jsonXpa = json.load(open('./results_traditional_Xpa.json','r'))
# fname1 = "simall_trad"
# fname2 = "irfall_trad"
# fname3 = "DHall_trad"
# fname4 = "eptime_trad"

# the below is common between KS and Xpa
drop = jsonXpa["input"]["drop"]
simT = jsonXpa["input"]["simT"]
irfdrop = jsonXpa["input"]["irfdrop"]
nm = np.array(jsonXpa["input"]["knotsm"]).shape[0]
nz = np.array(jsonXpa["input"]["Gz"]).shape[0]
GAMY = jsonXpa["input"]["param"][0]
DELTA = jsonXpa["input"]["param"][2]
# print(DELTA)


## unconditional moments
YvecXpa,IvecXpa,NvecXpa,KvecXpa,CvecXpa,KpvecXpa,Zvec = calcaggstat(jsonXpa,simT,drop,fname1+'Xpa')
YvecKS,IvecKS,NvecKS,KvecKS,CvecKS,KpvecKS,Zvec = calcaggstat(jsonKS,simT,drop,fname1+'KS')

figure()
st = 1001
ed = 1050
time = np.linspace(st,ed,50)

subplot(231)
plot(time,np.log(YvecXpa[st-1:ed]),'b-o',markerfacecolor='none',label='Xpa')
plot(time,np.log(YvecKS[st-1:ed]),'k-x',label='KS')
title('Output')
ylabel('Log')
xlim(st,ed)
ylim(-0.8,-0.1)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))
legend(loc="upper right", edgecolor="black")

subplot(232)
plot(time,np.log(IvecXpa[st-1:ed]),'b-o',markerfacecolor='none')
plot(time,np.log(IvecKS[st-1:ed]),'k-x')
title('Investment')
xlim(st,ed)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))

subplot(233)
plot(time,np.log(NvecXpa[st-1:ed]),'b-o',markerfacecolor='none')
plot(time,np.log(NvecKS[st-1:ed]),'k-x')
title('Labor')
xlim(st,ed)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))

subplot(234)
plot(time,np.log(KvecXpa[st-1:ed]),'b-o',markerfacecolor='none')
plot(time,np.log(KvecKS[st-1:ed]),'k-x')
title('Capital')
xlabel('Year')
ylabel('Log')
xlim(st,ed)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))

subplot(235)
plot(time,np.log(CvecXpa[st-1:ed]),'b-o',markerfacecolor='none')
plot(time,np.log(CvecKS[st-1:ed]),'k-x')
title('Consumption')
xlabel('Year')
xlim(st,ed)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))

subplot(236)
plot(time,np.log(Zvec[st-1:ed]),'r-')
title('TFP')
xlabel('Year')
xlim(st,ed)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))

subplots_adjust(wspace=0.4, hspace=0.4)
# subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
# show()
savefig(fname1+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0)


## conditional moments (IRFs)
YirvecXpa,IirvecXpa,NirvecXpa,KirvecXpa,CirvecXpa,Zirvec = calcstochss(jsonXpa,irfdrop,fname2+'Xpa')
YirvecKS,IirvecKS,NirvecKS,KirvecKS,CirvecKS,Zirvec = calcstochss(jsonKS,irfdrop,fname2+'KS')

figure()
st = irfdrop+1
ed = irfdrop+22
time = np.linspace(st,ed,22)

subplot(231)
plot(time,100*np.log(YirvecXpa[st-1:ed]/YirvecXpa[st-2]),'b-o',markerfacecolor='none',label='Xpa')
plot(time,100*np.log(YirvecKS[st-1:ed]/YirvecKS[st-2]),'k-x',label='KS')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('Output')
ylabel('Percent')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1
legend(loc="upper right", edgecolor="black")

subplot(232)
plot(time,100*np.log(IirvecXpa[st-1:ed]/IirvecXpa[st-2]),'b-o',markerfacecolor='none')
plot(time,100*np.log(IirvecKS[st-1:ed]/IirvecKS[st-2]),'k-x')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('Investment')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1

subplot(232)
plot(time,100*np.log(IirvecXpa[st-1:ed]/IirvecXpa[st-2]),'b-o',markerfacecolor='none')
plot(time,100*np.log(IirvecKS[st-1:ed]/IirvecKS[st-2]),'k-x')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('Investment')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1

subplot(233)
plot(time,100*np.log(NirvecXpa[st-1:ed]/NirvecXpa[st-2]),'b-o',markerfacecolor='none')
plot(time,100*np.log(NirvecKS[st-1:ed]/NirvecKS[st-2]),'k-x')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('Labor')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1

subplot(234)
plot(time,100*np.log(KirvecXpa[st-1:ed]/KirvecXpa[st-2]),'b-o',markerfacecolor='none')
plot(time,100*np.log(KirvecKS[st-1:ed]/KirvecKS[st-2]),'k-x')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('Capital')
xlabel('Year')
ylabel('Percent')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1

subplot(235)
plot(time,100*np.log(CirvecXpa[st-1:ed]/CirvecXpa[st-2]),'b-o',markerfacecolor='none')
plot(time,100*np.log(CirvecKS[st-1:ed]/CirvecKS[st-2]),'k-x')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('Consumption')
xlabel('Year')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1

subplot(236)
plot(time,100*np.log(Zirvec[st-1:ed]/Zirvec[st-2]),'r-')
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0)
title('TFP')
xlabel('Year')
xlim(st,ed)
xticks([(st+1)+0,(st+1)+5,(st+1)+10,(st+1)+15,(st+1)+20],[0,5,10,15,20]) # st is set to -1

subplots_adjust(wspace=0.4, hspace=0.4)
# subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
# show()
savefig(fname2+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0)


## accuracy statistics

mpvec0Xpa,pvec0Xpa,mpvec1Xpa,pvec1Xpa = calcaccurat(jsonXpa,simT,drop,fname3+'Xpa')
mpvec0KS,pvec0KS,mpvec1KS,pvec1KS = calcaccurat(jsonKS,simT,drop,fname3+'KS')

figure()
st = 1001
ed = 1050
time = np.linspace(st,ed,50)

subplot(221)
plot(time,np.log(mpvec0KS[st-1:ed]),'m-^',markerfacecolor='none',label='Dynamic')
plot(time,np.log(mpvec1KS[st-1:ed]),'c-o',markerfacecolor='none',label='Static')
plot(time,np.log(KpvecKS[st-1:ed]),'k-x',label='Actual')
title('KS:Capital')
ylabel('Log')
xlim(st,ed)
ylim(-0.2,0.4)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))
legend(loc="lower right", edgecolor="black")

subplot(222)
plot(time,-np.log(pvec0KS[st-1:ed]),'m-^',markerfacecolor='none',label='Dynamic')
plot(time,-np.log(pvec1KS[st-1:ed]),'c-o',markerfacecolor='none',label='Static')
plot(time,np.log(CvecKS[st-1:ed]),'k-x',label='Actual')
title('KS:Consumption')
ylabel('Log')
xlim(st,ed)
ylim(-0.95,-0.75)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))
# legend(loc="lower right", edgecolor="black")

subplot(223)
plot(time,np.log(mpvec0Xpa[st-1:ed]),'m-^',markerfacecolor='none',label='Dynamic')
plot(time,np.log(mpvec1Xpa[st-1:ed]),'c-o',markerfacecolor='none',label='Static')
plot(time,np.log(KpvecXpa[st-1:ed]),'k-x',label='Actual')
title('Xpa:Capital')
ylabel('Log')
xlim(st,ed)
ylim(-0.2,0.4)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))
# legend(loc="lower right", edgecolor="black")

subplot(224)
plot(time,-np.log(pvec0Xpa[st-1:ed]),'m-^',markerfacecolor='none',label='Dynamic')
plot(time,-np.log(pvec1Xpa[st-1:ed]),'c-o',markerfacecolor='none',label='Static')
plot(time,np.log(CvecXpa[st-1:ed]),'k-x',label='Actual')
title('Xpa:Consumption')
ylabel('Log')
xlim(st,ed)
ylim(-0.95,-0.75)
xticks(np.array([st+1,st+10,st+20,st+30,st+40,st+50]),np.array([1,10,20,30,40,50]))
# legend(loc="lower right", edgecolor="black")

subplots_adjust(wspace=0.4, hspace=0.4)
# subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
show()
# savefig(fname3+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0)
savefig('DHall_extend.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0)


# elapsed time
eptimeinKS = np.array(jsonKS["output"]["eptimein"])
eptimeoutKS = np.array(jsonKS["output"]["eptimeout"])
numloopKS = eptimeoutKS.shape[0]
eptimeKS = np.sum(eptimeinKS)+np.sum(eptimeoutKS)
eptimeinXpa = np.array(jsonXpa["output"]["eptimein"])
eptimeoutXpa = np.array(jsonXpa["output"]["eptimeout"])
numloopXpa = eptimeoutXpa.shape[0]
eptimeXpa = np.sum(eptimeinXpa)+np.sum(eptimeoutXpa)

data1 = np.array([np.mean(eptimeinKS),np.mean(eptimeoutKS),numloopKS,eptimeKS])
data2 = np.array([np.mean(eptimeinXpa),np.mean(eptimeoutXpa),numloopXpa,eptimeXpa])
data3 = data2/data1
data = np.vstack([data1,data2,data3])
print(data)
np.savetxt(fname4 + ".csv", data, delimiter=",", fmt='%.4f')
