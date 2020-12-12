import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

"""
Follows
"Dynamic Critical Behavior of the Chayes–Machta Algorithm for the Random-Cluster Model, I. Two Dimensions"
chapter 4, i.e. Sokal's windowing method.
"""


def autocorrelation_function_fft(ts_):
    ts = np.asarray(ts_, dtype=np.float64)
    ts -= ts.mean()
    ts_len = len(ts)
    cov = fftconvolve(ts, ts[::-1], mode='full')
    idx = int((len(cov) - 1) / 2)
    ints = np.arange(0, idx + 1)
    cov = cov[idx:]
    cov /= (ts_len - ints)
    rho = cov / cov[0]
    return rho, cov[0]


def tau_int_cal_rho(rho, c=6.):
    rho[0] *= 0.5
    taus = np.cumsum(rho)
    taus *= c
    ms = np.arange(0, len(rho))
    idx = np.argmax(ms >= taus)
    tau = taus[idx] / c
    m = ms[idx]
    # m=0
    # while m<c*taus[m]:
    #    m+=1
    #tau = taus[m]
    err_tau = np.sqrt(2. * (2. * m + 1) / len(rho)) * tau
    return tau, err_tau
