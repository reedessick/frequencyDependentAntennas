description = """\
a module that contains functions for computing the frequency-dependent antenna patterns of an IFO.

Expressions are taken LIGO-T060237 in a coordinate-dependent form. 
These were then modified to a coordinate independent form by the module's author.

This included validation against Eqn. B7 from Anderson, et all PhysRevD 63(04) 2003
"""
author = "Reed Essick"

#-------------------------------------------------

import numpy as np

#-------------------------------------------------

#--- helper functions

def __dV__( q, exp_twoFreqsT ):
    '''
    helper function that determines the prefactor for the antenna patterns.
    this may typically be included in the PSD rather than the antenna patterns, so we may want to set this to 1
    '''
    return exp_twoFreqsT/(1-q*exp_twoFreqsT)

def __D__( freqsT, twoFreqsT, exp_twoFreqsT, N ):
    '''
    helper function that returns the part of the frequency dependence that depends on the arm's directions
    '''
    if isinstance(freqsT, (int, float)):
        if freqsT==0:
            return 1

    elif np.all(freqsT==0):
        return 1

    a = 1-N
    b = 1+N

    return ((1-np.exp(-a*twoFreqsT))/a - exp_twoFreqsT*(1-np.exp(b*freqsT))/b)/twoFreqsT

#--- actual responses

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
    
    ### compute terms common to each antenna pattern
    freqsT = 2*np.pi*freqs*T
    twoFreqsT = 2*freqsT
    exp_twoFreqsT = np.exp(-twoFreqsT)

    dV = __dV__(q, exp_twoFreqsT)
    dV_xx = __D__(freqsT, twoFreqsT, exp_twoFreqsT, nx)
    dV_yy = __D__(freqsT, twoFreqsT, exp_twoFreqsT, ny)

    ### assemble the derivatives
    dV_dhp = dV*(dV_xx*(nx**2 - sinPhi**2) - dV_yy*(ny**2 - cosPhi**2))
    dV_dhx = -dV*(dV_xx + dV_yy)*cosTheta*2*sinPhi*cosPhi

    return dV_dhp*cos2psi - dV_dhc*sin2psi, dV_dhp*sin2psi + dV_dhc*cos2psi

def coordIndep_antenna_response( theta, phi, psi, ex, ey, q, T=1., freqs=0. ):
    '''
    the input variables mean things
    theta, phi, psi, and freqs should all be the same length np.ndarray objects if they are not floats

    ex and ey are the unit vectors describing the directions of the detector arms in an Earth-fixed frame (should be np.ndarray objects)

    return F+, Fx
    '''
    ### compute the trigonometric functions only once
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)

    cosPhi = np.cos(phi)
    sinPhi = np.sin(phi)

    cosPsi = np.cos(psi)
    sinPsi = np.sin(psi)

    ### compute unit vectors in Earth-Fixed frame from the wave-frame
    ex_wave = np.array([sinPhi*cosPsi - sinPsi*cosPhi*cosTheta, -cosPhi*cosPsi - sinPsi*sinPhi*cosTheta, sinPhi*sinTheta])
    ey_wave = np.array([-sinPhi*sinPsi - cosPsi*cosPhi*cosTheta, cosPhi*sinPsi - cosPsi*sinPhi*cosTheta, sinTheta*cosPsi])

    ### compute cartesian vector for the line-of-sight to the source
    n = np.array([sinTheta*cosPhi, sinTheta*sinPhi, cosTheta])

    ### compute detector matrix
    freqsT = 2*np.pi*freqs*T
    twoFreqsT = 2*freqsT
    exp_twoFreqsT = np.exp(-twoFreqsT)

    dV = __dV__(q, exp_twoFreqsT)
    dV_yy = dV * __D__(freqsT, twoFreqsT, exp_twoFreqsT, np.sum(n*ey))
    dV_xx = dV * __D__(freqsT, twoFreqsT, exp_twoFreqsT, np.sum(n*ex))

    ### assemble these parts into antenna responses
    Fp = 0.
    Fx = 0.
    for i in xrange(3): ### FIXME? may be able to do this more efficiently vi np.array manipulations?
        for j in xrange(3):
            detector = dV_xx * ex[i]*ex[j] - dV_yy * ey[i]*ey[j] ### compute matrix element for detector matrix
            Fp += detector * (ex_wave[i]*ex_wave[j] - ey_wave[i]*ey_wave[j]) ### multiply by matrix element from wave polarizations
            Fx += detector * (ex_wave[i]*ey_wave[j] + ey_wave[i]*ex_wave[j])

    return Fp, Fx
