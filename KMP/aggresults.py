import numpy as np;
from matplotlib.pyplot import *
import json;


# returns index starting from 1
def gridlookup2(x0,xgrid):

    nx = xgrid.shape[0];

    ix = 0;
    for jx in range(nx):
        if x0<=xgrid[jx]:
            break;
        ix = ix+1;

    ix = min(max(1,ix),nx-1);

    return ix;


def calcDH(knotsm,mpmat0,izvec,Kvec,statflag):

    tott = izvec.shape[0];

    mp = Kvec[0];
    mpvec = np.zeros(tott);
    # pvec = np.zeros(tott);

    for tt in range(tott):

        if statflag==1:
            mnow = Kvec[tt];
        else:
            mnow = mp;

        # NOTE: im and iz are indices starting from 1, not 0 as in Fortran
        iz = izvec[tt];
        im = gridlookup2(np.log(mnow),np.log(knotsm));
        # print(np.array([im,iz]));
        wm = np.log(knotsm[im]/mnow)/np.log(knotsm[im]/knotsm[im-1]);
        mp = np.exp(wm*np.log(mpmat0[im-1,iz-1]) + (1-wm)*np.log(mpmat0[im,iz-1]));
        # cnow = 1/np.exp(wm*np.log(pmat0[im-1,iz-1]) + (1-wm)*np.log(pmat0[im,iz-1]));

        mpvec[tt] = mp
        # pvec[tt] = 1/cnow;

    return mpvec;


def calcaccurat(json,simT,drop,filename):


    knotsm = np.array(json["input"]["knotsm"]);
    izvec = np.array(json["input"]["izvec"]);
    mpmat0 = np.array(json["output"]["mpmat0"]);
    # pmat0 = np.array(json["output"]["pmat0"]);
    Kvec = np.array(json["output"]["Kvec"]);
    Kpvec = np.array(json["output"]["Kpvec"]);
    # Cvec = np.array(json["output"]["Cvec"]);

    nm = np.array(json["input"]["knotsm"]).shape[0];
    nz = np.array(json["input"]["Gz"]).shape[0];
    mpmat = np.zeros((nm,nz));
    pmat = np.zeros((nm,nz));

    for iz in range(nz):
        for im in range(nm):

            index = nm*iz+im;
            # print(index);
            mpmat[im,iz] = mpmat0[index];
            # pmat[im,iz] = pmat0[index];

    # dynamic forecasts
    mpvec0 = calcDH(knotsm, mpmat, izvec, Kvec, 0);
    # static forecasts
    mpvec1 = calcDH(knotsm, mpmat, izvec, Kvec, 1);

    izvec = izvec[drop:simT+drop];
    Kvec = Kvec[drop:simT+drop];
    Kpvec = Kpvec[drop:simT+drop];
    # Cvec = Cvec[drop:simT+drop];
    mpvec0 = mpvec0[drop:simT+drop];
    # pvec0 = pvec0[drop:simT+drop];
    mpvec1 = mpvec1[drop:simT+drop];
    # pvec1 = pvec1[drop:simT+drop];

    DHmax = np.zeros((nz,1));
    DHmean = np.zeros((nz,1));
    RMSE = np.zeros((nz,1));
    Rsquared = np.zeros((nz,1));

    for iz in range(nz):

        DHmax[iz,0] = 100*np.max(np.abs(np.log(mpvec0[izvec==iz+1]/Kpvec[izvec==iz+1])));
        # DHmax[iz,1] = 100*np.max(np.abs(np.log(pvec0[izvec==iz+1]*Cvec[izvec==iz+1])));
        DHmean[iz,0] = 100*np.sum(np.abs(np.log(mpvec0[izvec==iz+1]/Kpvec[izvec==iz+1])))/simT;
        # DHmean[iz,1] = 100*np.sum(np.abs(np.log(pvec0[izvec==iz+1]*Cvec[izvec==iz+1])))/simT;

        y0 = np.log(Kpvec[izvec==iz+1]);
        y0 = np.reshape(y0,(y0.shape[0],1));
        # print(y0.shape);
        # y1 = -np.log(Cvec[izvec==iz+1]);
        # y1 = np.reshape(y1,(y1.shape[0],1));
        # how to obtain row vectors?
        # y = np.hstack([y0,y1]);
        # below does not work
        # y = np.hstack([np.log(Kpvec[izvec==iz+1]),-np.log(Cvec[izvec==iz+1].T)]);
        # print(y.shape);
        y = y0;
        ymean = np.mean(y,axis=0);
        # print(ymean.shape);
        ydev = y-ymean;
        # ydev0 = y[:,0]-ymean[0];
        # ydev1 = y[:,1]-ymean[1];
        # print(ydev.shape);
        e0 = np.log(mpvec1[izvec==iz+1]/Kpvec[izvec==iz+1]);
        # e1 = np.log(pvec1[izvec==iz+1]*Cvec[izvec==iz+1]);

        RMSE[iz,0] = 100*np.sqrt(np.mean(np.power(e0,2)));
        # RMSE[iz,1] = 100*np.sqrt(np.mean(np.power(e1,2)));
        Rsquared[iz,0] = 1-np.dot(e0.T,e0)/np.dot(ydev[:,0].T,ydev[:,0]);
        # Rsquared[iz,1] = 1-np.dot(e1.T,e1)/np.dot(ydev[:,1].T,ydev[:,1]);

    data = np.hstack([DHmax, DHmean, RMSE, Rsquared]);
    print(data);
    np.savetxt(filename + ".csv", data, delimiter=",", fmt='%.4f')

    return mpvec0, mpvec1;


def calcstochss(json,irfdrop,filename):

    Kirvec = np.array(json["irf"]["Kvec"]);
    Zirvec = np.array(json["irf"]["Zvec"]);
    Yirvec = np.array(json["irf"]["Yvec"]);
    Iirvec = np.array(json["irf"]["Ivec"]);
    Nirvec = np.array(json["irf"]["Nvec"]);
    Cirvec = np.array(json["irf"]["Cvec"]);

    shareWirvec = np.zeros((Kirvec.shape[0],8));
    shareWirvec[:,0] = np.array(json["irf"]["shareW1"]);
    shareWirvec[:,1] = np.array(json["irf"]["shareW2"]);
    shareWirvec[:,2] = np.array(json["irf"]["shareW3"]);
    shareWirvec[:,3] = np.array(json["irf"]["shareW4"]);
    shareWirvec[:,4] = np.array(json["irf"]["shareW5"]);
    shareWirvec[:,5] = np.array(json["irf"]["shareW9095"]);
    shareWirvec[:,6] = np.array(json["irf"]["shareW9599"]);
    shareWirvec[:,7] = np.array(json["irf"]["shareWT1"]);
    gini = np.array(json["irf"]["gini"]);

    mnow = Kirvec[irfdrop-1];
    ynow = Yirvec[irfdrop-1];
    # nnow = Nirvec[irfdrop-1];
    # cnow = Cirvec[irfdrop-1];

    data1 = shareWirvec[irfdrop-1,:];
    data2 = np.array([gini[irfdrop-1], mnow/ynow, 100*np.log(Cirvec[irfdrop+1]/Cirvec[irfdrop])]);
    data = np.hstack([data1,data2]);
    print(data);
    np.savetxt(filename + ".csv", data, delimiter=",", fmt='%.4f')

    return Yirvec,Iirvec,Nirvec,Kirvec,Cirvec,Zirvec;


def calcaggstat(json,simT,drop,filename):

    Kvec = np.array(json["output"]["Kvec"]);
    Zvec = np.array(json["output"]["Zvec"]);
    Yvec = np.array(json["output"]["Yvec"]);
    Ivec = np.array(json["output"]["Ivec"]);
    Nvec = np.array(json["output"]["Nvec"]);
    Cvec = np.array(json["output"]["Cvec"]);
    Kpvec = np.array(json["output"]["Kpvec"]);
    # Xvec = np.array(json["output"]["Xvec"]);

    Kvec = Kvec[drop:simT + drop];
    Zvec = Zvec[drop:simT + drop];
    Yvec = Yvec[drop:simT + drop];
    Ivec = Ivec[drop:simT + drop];
    Nvec = Nvec[drop:simT + drop];
    Cvec = Cvec[drop:simT + drop];
    Kpvec = Kpvec[drop:simT + drop];

    KYmean = np.mean(Kvec/Yvec);
    corr1 = np.corrcoef(Yvec,Cvec);
    std1 = np.std(Ivec);
    ac1 = np.corrcoef(Yvec[0:simT-3],Yvec[3:simT]);

    data = np.array([KYmean,corr1[0,1],std1,ac1[0,1]]);
    print(data);
    np.savetxt(filename + ".csv", np.reshape(data,(1,data.shape[0])), delimiter=",", fmt='%.4f')

    return Yvec,Ivec,Nvec,Kvec,Cvec,Kpvec,Zvec;


########################################
# main part of the script
########################################
# for KMP calibration
jsonKS = json.load(open('./results_calibKMP_KS.json','r'));
jsonXpa = json.load(open('./results_calibKMP_Xpa.json','r'));
fname1 = "simall_KMP";
fname2 = "irfall_KMP";
fname3 = "DHall_KMP";
fname4 = "eptime_KMP";
# for KS calibration
# jsonKS = json.load(open('./results_calibKS_KS.json','r'));
# jsonXpa = json.load(open('./results_calibKS_Xpa.json','r'));
# fname1 = "simall_KS";
# fname2 = "irfall_KS";
# fname3 = "DHall_KS";
# fname4 = "eptime_KS";
# for hety calibration
# jsonKS = json.load(open('./results_calibhety_KS.json','r'));
# jsonXpa = json.load(open('./results_calibhety_Xpa.json','r'));
# fname1 = "simall_hety";
# fname2 = "irfall_hety";
# fname3 = "DHall_hety";
# fname4 = "eptime_hety";

drop = jsonXpa["input"]["drop"];
simT = jsonXpa["input"]["simT"];
irfdrop = jsonXpa["input"]["irfdrop"];
nm = np.array(jsonXpa["input"]["knotsm"]).shape[0];
nz = np.array(jsonXpa["input"]["Gz"]).shape[0];


## unconditional moments
YvecXpa,IvecXpa,NvecXpa,KvecXpa,CvecXpa,KpvecXpa,Zvec = calcaggstat(jsonXpa,simT,drop,fname1+'Xpa');
YvecKS,IvecKS,NvecKS,KvecKS,CvecKS,KpvecKS,Zvec = calcaggstat(jsonKS,simT,drop,fname1+'KS');

figure();
st = 1001;
ed = 1200;
time = np.linspace(st,ed,200);

subplot(221);
plot(time,np.log(KvecXpa[st-1:ed]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,np.log(KvecKS[st-1:ed]),'k-x',label='KS');
title('Capital');
ylabel('Log');
xlim(st,ed);
xticks(np.array([st+1,st+40,st+80,st+120,st+160,st+200]),np.array([1,40,80,120,160,200]));
# legend(loc="lower right", edgecolor="black");

subplot(222);
plot(time,np.log(CvecXpa[st-1:ed]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,np.log(CvecKS[st-1:ed]),'k-x',label='KS');
title('Consumption');
ylabel('Log');
xlim(st,ed);
xticks(np.array([st+1,st+40,st+80,st+120,st+160,st+200]),np.array([1,40,80,120,160,200]));
legend(loc="lower right", edgecolor="black");

subplot(223);
plot(time,np.log(YvecXpa[st-1:ed]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,np.log(YvecKS[st-1:ed]),'k-x',label='KS');
title('Output');
ylabel('Log');
xlim(st,ed);
xticks(np.array([st+1,st+40,st+80,st+120,st+160,st+200]),np.array([1,40,80,120,160,200]));
# legend(loc="lower right", edgecolor="black");

subplot(224);
plot(time,np.log(IvecXpa[st-1:ed]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,np.log(IvecKS[st-1:ed]),'k-x',label='KS');
title('Investment');
ylabel('Log');
xlim(st,ed);
xticks(np.array([st+1,st+40,st+80,st+120,st+160,st+200]),np.array([1,40,80,120,160,200]));
# legend(loc="lower right", edgecolor="black");

subplots_adjust(wspace=0.4, hspace=0.4);
subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
# show();
savefig(fname1+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0);


## conditional moments (IRFs)
YirvecXpa,IirvecXpa,NirvecXpa,KirvecXpa,CirvecXpa,Zirvec = calcstochss(jsonXpa,irfdrop,fname2+'Xpa');
YirvecKS,IirvecKS,NirvecKS,KirvecKS,CirvecKS,Zirvec = calcstochss(jsonKS,irfdrop,fname2+'KS');

figure();
st = irfdrop+1;
ed = irfdrop+7;
time = np.linspace(st,ed,7);

subplot(221);
plot(time,100*np.log(Zirvec[st-1:ed]/Zirvec[st-2]),'r-');
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0);
title('TFP');
ylabel('Percent');
xlim(st,ed);
xticks([(st+1)+0,(st+1)+1,(st+1)+2,(st+1)+3,(st+1)+4,(st+1)+5],[0,1,2,3,4,5]); # st is set to -1
# legend(loc="upper right", edgecolor="black");

subplot(222);
plot(time,100*np.log(CirvecXpa[st-1:ed]/CirvecXpa[st-2]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,100*np.log(CirvecKS[st-1:ed]/CirvecKS[st-2]),'k-x',label='KS');
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0);
title('Consumption');
ylabel('Percent');
xlim(st,ed);
xticks([(st+1)+0,(st+1)+1,(st+1)+2,(st+1)+3,(st+1)+4,(st+1)+5],[0,1,2,3,4,5]); # st is set to -1
legend(loc="lower right", edgecolor="black");

subplot(223);
plot(time,100*np.log(YirvecXpa[st-1:ed]/YirvecXpa[st-2]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,100*np.log(YirvecKS[st-1:ed]/YirvecKS[st-2]),'k-x',label='KS');
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0);
title('Output');
ylabel('Percent');
xlim(st,ed);
xticks([(st+1)+0,(st+1)+1,(st+1)+2,(st+1)+3,(st+1)+4,(st+1)+5],[0,1,2,3,4,5]); # st is set to -1
# legend(loc="lower right", edgecolor="black");

subplot(224);
plot(time,100*np.log(IirvecXpa[st-1:ed]/IirvecXpa[st-2]),'b-o',markerfacecolor='none',label='Xpa');
plot(time,100*np.log(IirvecKS[st-1:ed]/IirvecKS[st-2]),'k-x',label='KS');
plot(np.array([st,ed]),np.array([0,0]),'k-',linewidth=1.0);
title('Investment');
ylabel('Percent');
xlim(st,ed);
xticks([(st+1)+0,(st+1)+1,(st+1)+2,(st+1)+3,(st+1)+4,(st+1)+5],[0,1,2,3,4,5]); # st is set to -1
# legend(loc="lower right", edgecolor="black");

subplots_adjust(wspace=0.4, hspace=0.4);
subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
# show();
savefig(fname2+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0);


## accuracy statistics

mpvec0Xpa,mpvec1Xpa = calcaccurat(jsonXpa,simT,drop,fname3+'Xpa');
mpvec0KS,mpvec1KS = calcaccurat(jsonKS,simT,drop,fname3+'KS');

figure();
st = 1001;
ed = 1500;
time = np.linspace(st,ed,500);

subplot(211);
plot(time,np.log(mpvec0KS[st-1:ed]),'m-.',markerfacecolor='none',label='Dynamic');
plot(time,np.log(mpvec1KS[st-1:ed]),'c--',markerfacecolor='none',label='Static');
plot(time,np.log(KpvecKS[st-1:ed]),'k-',label='Actual');
title('KS:Capital');
ylabel('Log');
xlim(st,ed);
legend(loc="lower right", edgecolor="black");

subplot(212);
plot(time,np.log(mpvec0Xpa[st-1:ed]),'m-.',markerfacecolor='none',label='Dynamic');
plot(time,np.log(mpvec1Xpa[st-1:ed]),'c--',markerfacecolor='none',label='Static');
plot(time,np.log(KpvecXpa[st-1:ed]),'k-',label='Actual');
title('Xpa:Capital');
ylabel('Log');
xlim(st,ed);
# legend(loc="lower right", edgecolor="black");

subplots_adjust(wspace=0.4, hspace=0.4);
subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
# show();
savefig(fname3+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0);


# elapsed time
eptimeinKS = np.array(jsonKS["output"]["eptimein"]);
eptimeoutKS = np.array(jsonKS["output"]["eptimeout"]);
numloopKS = eptimeoutKS.shape[0];
eptimeKS = np.sum(eptimeinKS)+np.sum(eptimeoutKS);
eptimeinXpa = np.array(jsonXpa["output"]["eptimein"]);
eptimeoutXpa = np.array(jsonXpa["output"]["eptimeout"]);
numloopXpa = eptimeoutXpa.shape[0];
eptimeXpa = np.sum(eptimeinXpa)+np.sum(eptimeoutXpa);

data1 = np.array([np.mean(eptimeinKS),np.mean(eptimeoutKS),numloopKS,eptimeKS]);
data2 = np.array([np.mean(eptimeinXpa),np.mean(eptimeoutXpa),numloopXpa,eptimeXpa]);
data3 = data2/data1;
data = np.vstack([data1,data2,data3]);
print(data);
np.savetxt(fname4 + ".csv", data, delimiter=",", fmt='%.4f')
