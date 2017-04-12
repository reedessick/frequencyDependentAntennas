description = "a module that contains functions for computing the frequency-dependent antenna patterns of an IFO"
author = "Reed Essick"

#-------------------------------------------------

import numpy as np

#-------------------------------------------------

def antenna_response( theta, phi, psi, q, T=1., freqs=0. ):
    '''
    the input variables mean things
    theta, phi, psi, and freqs should all be the same length np.ndarray objects if they are not floats

    return F+, Fx
    '''
    ### compute trigonometric functions only once!
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)

    cosPhi = np.cos(phi)
    sinPhi = np.sin(phi)

    cos2psi = np.cos(2*psi)
    sin2psi = np.cos(2*psi)

    ### define relevant coordinate unit vectors for convenience
    nx = sinTheta*cosPhi
    ny = sinTheta*sinPhi
#    nz = cosTheta
    
    ### compute terms common to each antenna pattern
    freqsT = freqs*T
    twoFreqsT = 2*freqsT
    exp_twoFreqsT = np.exp(-twoFreqsT)

    dV_prefact = exp_twoFreqsT/(twoFreqsT*(1-q*exp_twoFreqsT))

    a = 1-nx
    b = 1+nx
    dV_nx = (1-np.exp(-a*freqsT))/a - exp_twoFreqsT*(1-np.exp(b*freqsT))/b

    a = 1-ny
    b = 1+ny
    dV_ny = (1-np.exp(-a*freqsT))/a - exp_twoFreqsT*(1-np.exp(b*freqsT))/b

    ### assemble the derivatives
    dV_dhp = dV_prefact*(dV_nx*(nx**2 - sinPhi**2) - dV_ny*(ny**2 - cosPhi**2)
    dV_dhx = -dV_prefact*(dV_nx + dV_ny)*cosTheta*2*sinPhi*cosPhi

    return dV_dhp*cos2psi - dV_dhc*sin2psi, dV_dhp*sin2psi + dV_dhc*cos2psi
