import numpy as np
from matplotlib.pyplot import *
import json


# json = json.load(open('./results_calibKMP_KS.json','r'))
# fname = 'density_KMP'
json = json.load(open('./results_calibKS_KS.json','r'))
fname = 'density_KS'
# json = json.load(open('./results_calibhety_KS.json','r'))
# fname = 'density_hety'
Ge = np.array(json["input"]["Ge"])
knotsb = np.array(json["input"]["knotsb"])
mux = np.array(json["input"]["mux"])
ne = Ge.shape[0]
nb = knotsb.shape[0]
nx = mux.shape[0]
mu0 = np.array(json["ss"]["muss"])
mpmat = np.array(json["ss"]["mpmat"])
shareWvec = np.array(json["ss"]["shareWvec"])
gini = np.array(json["ss"]["gini"])

mnow = 0.0
mp = 0.0

mumat1 = np.zeros((nb,nx))
mpmat1 = np.zeros((nb,nx))

for ix in range(nx):

    for ib in range(nb):

        mumat1[ib,ix] = mu0[nb*ix+ib]
        mpmat1[ib,ix] = mpmat[nb*ix+ib]

    mnow = mnow + np.dot(mumat1[:,ix].T,knotsb)
    mp   = mp   + np.dot(mumat1[:,ix].T,mpmat1[:,ix])

if (nx>2):

    muu = np.sum(mumat1[:,0:nx-2:2],1)
    mue = np.sum(mumat1[:,1:nx-1:2],1)

else:

    muu = mumat1[:,0]
    mue = mumat1[:,1]

lnow = np.sum(muu)*Ge[0] + np.sum(mue)*Ge[1]
# print(lnow)

data1 = shareWvec
data2 = np.array([gini,np.power(mnow/lnow,1-0.36)])
data = np.hstack((data1,data2))
print(data)
np.savetxt(fname + ".csv", data, delimiter=",", fmt='%.4f')
