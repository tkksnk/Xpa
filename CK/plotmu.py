import numpy as np
from matplotlib.pyplot import *
import json


json = json.load(open('./results_CKnx5na51nk501.json','r'))
# fname = 'density_trad'
# json = json.load(open('./results_extend_KS.json','r'))
# fname = 'density_extend'
nx = np.array(json["input"]["xgrid"]).shape[0]
nk = np.array(json["input"]["kgrid"]).shape[0]
# GAMY = json["input"]["param"][0]
# DELTA = json["input"]["param"][2]
mumat = np.array(json["ss"]["mumat"])
# mpmat = np.array(json["ss"]["mpmat"])
# ymat = np.array(json["ss"]["ymat"])
# nmat = np.array(json["ss"]["nmat"])
knotsb = np.array(json["input"]["kgrid"])

# mnow = 0.0
# mp = 0.0
# ynow = 0.0
# nnow = 0.0
#
mumat1 = np.zeros((nk,nx))
# mpmat1 = np.zeros((nb,ne))
# ymat1 = np.zeros((nb,ne))
# nmat1 = np.zeros((nb,ne))
# ikmat1 = np.zeros((nb,ne))
#
for ix in range(nx):

    for ik in range(nk):

        mumat1[ik,ix] = mumat[nk*ix+ik]
        # mpmat1[ib,ie] = mpmat[nb*ie+ib]
        # ymat1[ib,ie]  = ymat[nb*ie+ib]
        # nmat1[ib,ie]  = nmat[nb*ie+ib]
        # ikmat1[ib,ie] = GAMY*mpmat1[ib,ie]/knotsb[ib]-(1-DELTA)

    # mnow = mnow + np.dot(mumat1[:,ie].T,knotsb)
    # mp   = mp   + np.dot(mumat1[:,ie].T,mpmat1[:,ie])
    # ynow = ynow + np.dot(mumat1[:,ie].T,ymat1[:,ie])
    # nnow = nnow + np.dot(mumat1[:,ie].T,nmat1[:,ie])
#
# inow = GAMY*mp-(1-DELTA)*mnow
# cnow = ynow-inow
#
# data1 = np.array([mnow/ynow, nnow, 1/cnow, 0, 0, 0, 0])
#
# ikvec = np.zeros(7)
# ikvec[0] = np.sum(np.sum(mumat1*ikmat1))

# for ie in range(ne):
#
#     for ib in range(nb):
#
#         ikvec[1] = ikvec[1] + mumat1[ib,ie]*np.power(ikmat1[ib,ie]-ikvec[0],2)
#         if np.abs(ikmat1[ib,ie])<0.01:
#             ikvec[2] = ikvec[2] + mumat1[ib,ie]
#         if ikmat1[ib,ie]>0.20:
#             ikvec[3] = ikvec[3] + mumat1[ib,ie]
#         if ikmat1[ib,ie]<-0.20:
#             ikvec[4] = ikvec[4] + mumat1[ib,ie]
#         if ikmat1[ib,ie]>=0.01:
#             ikvec[5] = ikvec[5] + mumat1[ib,ie]
#         if ikmat1[ib,ie]<=-0.01:
#             ikvec[6] = ikvec[6] + mumat1[ib,ie]
#
# data2 = ikvec
# data = np.vstack((data1,data2))
# print(data)
# np.savetxt(fname + ".csv", data, delimiter=",", fmt='%.4f')

print(np.sum(mumat1,0))
print(mumat1[:,2])
nb = nk

figure()
# if (nx==1):
plot(knotsb,np.sum(mumat1,1))
# xlim(knotsb[0],knotsb[nb-1])
# ylim(0,np.max(mumat1[:,0])*1.1)
xlabel('Capital')
ylabel('Density')
# else:
#     subplot(231)
#     plot(knotsb,mumat1[:,0])
#     # xlim(knotsb[0],knotsb[nb-1])
#     # xlim(knotsb[0],10.0)
#     # ylim(0,np.max(mumat1[:,0])*1.1)
#     # xlabel('Capital')
#     ylabel('Density')
#     title(r'$\epsilon_1$')
#
#     subplot(232)
#     plot(knotsb,mumat1[:,1])
#     # xlim(knotsb[0],knotsb[nb-1])
#     # xlim(knotsb[0],10.0)
#     # ylim(0,np.max(mumat1[:,1])*1.1)
#     # xlabel('Capital')
#     # ylabel('Density')
#     title(r'$\epsilon_2$')
#
#     subplot(233)
#     plot(knotsb,mumat1[:,2])
#     # xlim(knotsb[0],knotsb[nb-1])
#     # xlim(knotsb[0],10.0)
#     # ylim(0,np.max(mumat1[:,2])*1.1)
#     # xlabel('Capital')
#     # ylabel('Density')
#     title(r'$\epsilon_3$')
#
#     subplot(234)
#     plot(knotsb,mumat1[:,3])
#     # xlim(knotsb[0],knotsb[nb-1])
#     # ylim(0,np.max(mumat1[:,3])*1.1)
#     xlabel('Capital')
#     ylabel('Density')
#     title(r'$\epsilon_4$')
#
#     subplot(235)
#     plot(knotsb,mumat1[:,4])
#     # xlim(knotsb[0],knotsb[nb-1])
#     # ylim(0,np.max(mumat1[:,4])*1.1)
#     xlabel('Capital')
#     # ylabel('Density')
#     title(r'$\epsilon_5$')
#
# subplots_adjust(wspace=0.6, hspace=0.4)
# subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
# show()
# savefig(fname+'.eps', format='eps', dpi=600, bbox_inches='tight', pad_inches=0)
show()