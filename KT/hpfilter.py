import scipy as sp
from matplotlib.pyplot import *
from scipy import linalg as la
from scipy import sparse
# from numpy.linalg import inv;
# from scipy.sparse import lil_matrix, csc_matrix;
# from scipy.sparse.linalg import inv;
# import time;


def hp_filter(y, w):
    # make sure the inputs are the right shape
    m, n = y.shape
    if m < n:
        y = y.T
        m = n

    a = sp.array([w, -4 * w, ((6 * w + 1) / 2.)])
    d = sp.tile(a, (m, 1))

    d[0, 1] = -2. * w
    d[m - 2, 1] = -2. * w
    d[0, 2] = (1 + w) / 2.
    d[m - 1, 2] = (1 + w) / 2.
    d[1, 2] = (5 * w + 1) / 2.
    d[m - 2, 2] = (5 * w + 1) / 2.

    B = sparse.spdiags(d.T, [-2, -1, 0], m, m)
    B = B + B.T

    # report the filtered series, s
    s = sp.dot(la.inv(B.todense()), y)
    return s


def testhp():
    # read in and assign variables
    data = sp.loadtxt('investment.txt')
    data = sp.log(data)
    y = sp.atleast_2d(data)
    s = hp_filter(y, 1600)
    devs = y.T - s

    # plot the data, the filtered series, and abs. deviation
    subplot(211)
    plot(data, 'k')
    plot(s, 'r-')
    title("Data and HP filter")
    ylabel("log of investment")

    subplot(212)
    plot(devs, 'k')
    title("Stationary series")
    ylabel("log of investment")
    xlabel("Quarters from 1947Q1-2011Q4")
    show()

# yd = yd[8:simT-8,0];

# lam = 100;
#
# d = y.shape[0];
# start = time.time();
# # a = np.zeros((d,d));
# a = lil_matrix((d,d));
# for i in range(2,d-2):
#     a[i,i] = 6*lam+1;
#     a[i,i+1] = -4*lam;
#     a[i,i+2] = lam;
#     a[i,i-2] = lam;
#     a[i,i-1] = -4*lam;
#
# a[1,1] = 1+5*lam;
# a[1,2] = -4*lam;
# a[1,3] = lam;
# a[1,0] = -2*lam;
# a[0,0] = 1+lam;
# a[0,1] = -2*lam;
# a[0,2] = lam;
#
# a[d-2,d-2] = 1+5*lam;
# a[d-2,d-3] = -4*lam;
# a[d-2,d-4] = lam;
# a[d-2,d-1] = -2*lam;
# a[d-1,d-1] = 1+lam;
# a[d-1,d-2] = -2*lam;
# a[d-1,d-3] = lam;
#
# # inva = inv(a);
# inva = inv(csc_matrix(a));
# # yt = np.dot(inva,y);
# # yf = y-yt;
# elapsed_time = time.time() - start;
# print ("elapsed_time: {0}".format(elapsed_time) + " [sec]");
# yt = np.linalg.inv(a)*y;
# yf = y-yt;
# so = sp.atleast_2d(so);
# sd = hpfilter(so,100);
# sd = so-sf;
