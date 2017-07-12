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

def __D__( freqsT, N ):
    '''
    helper function that returns the part of the frequency dependence that depends on the arm's directions
    '''
    if isinstance(freqsT, (int, float, complex)):
        if freqsT==0:
            return np.ones_like(N, dtype=complex)

    elif np.all(freqsT==0):
        return np.ones((len(freqsT), len(N)), dtype=complex)

    #--- matt's condensed version
    n = np.outer(np.ones_like(freqsT), N)
    freqsT = np.outer(freqsT, np.ones_like(N))
    phi = freqsT/1j

    ans = np.empty_like(freqsT)

    truth = n==1
    ans[truth] = 0.5*(1+np.exp(freqsT[truth])*np.sinc(phi[truth]/np.pi))

    truth = n==-1
    ans[truth] = np.exp(-2*freqsT[truth])*0.5*(1 + np.exp(freqsT[truth])*np.sinc(phi[truth]/np.pi))

    truth = freqsT==0
    ans[truth] = 1

    truth = (np.abs(n)!=1)*(freqsT!=0)
    ans[truth] = np.exp(-freqsT[truth])/(1-n[truth]**2)*(np.sinc(phi[truth]/np.pi) + (n[truth]/freqsT[truth])*(np.cos(phi[truth]) - np.exp(freqsT[truth]*n[truth])))

    return ans

#def __D__( freqsT, N ):
#    '''
#    helper function that returns the part of the frequency dependence that depends on the arm's directions
#
#    DEPRECATED: known to be buggy because it doesn not presever |D(.,n)|==|D(.,-n)|
#    '''
#    if isinstance(freqsT, (int, float, complex)):
#        if freqsT==0:
#            return np.ones_like(N, dtype=complex)
#
#    elif np.all(freqsT==0):
#        return np.ones((len(freqsT), len(N)), dtype=complex)
#
#    a = 1-np.outer(np.ones_like(freqsT), N)
#    b = 2-a
#
#    freqsT = np.outer(freqsT, np.ones_like(N))
#    twoFreqsT = 2*freqsT
#    exp_twoFreqsT = np.exp(-twoFreqsT)
#
#    ### assemble things carefully so we don't have to worry about nans later
#    ans = np.zeros_like(a)
#
#    truth = twoFreqsT==0 ### this kills all terms
#    ans[truth] = 1
#
#    truth = np.logical_not(truth) ### flip it so it's everything else, use this to further filter array assignment
#
#    # where a is non-zero
#    this_truth = truth*(a!=0)
#    ans[this_truth] += (1 - np.exp( -a[this_truth]*freqsT[this_truth] ) ) / ( a[this_truth]*twoFreqsT[this_truth] )
#
#    # wehre a is zero
#    this_truth = truth*(a==0)
#    ans[this_truth] += 0.5
#
#    # where b is non-zero
#    this_truth = truth*(b!=0)
#    ans[this_truth] -= exp_twoFreqsT[this_truth] * (1 - np.exp( b[this_truth]*freqsT[this_truth] ) ) / ( b[this_truth]*twoFreqsT[this_truth] )
#
#    # where b is zero
#    this_truth = truth*(b==0)
#    ans[this_truth] -= 0.5*exp_twoFreqsT[this_truth]
#
#    return ans

#--- actual responses

def antenna_response( theta, phi, psi, ex, ey, T=1., freqs=0. ):
    '''
    the input variables mean things
    theta, phi, psi, and freqs should all be the same length np.ndarray objects if they are not floats

    ex and ey are the unit vectors describing the directions of the detector arms in an Earth-fixed frame (should be np.ndarray objects)

    return F+, Fx
    '''
    if isinstance(theta, (int, float)):
        theta = [theta]
        phi = [phi]
        psi = [psi]

    ### compute the trigonometric functions only once
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    cosPhi = np.cos(phi)
    sinPhi = np.sin(phi)
    cosPsi = np.cos(psi)
    sinPsi = np.sin(psi)

    ### compute unit vectors in Earth-Fixed frame from the wave-frame
    ex_wave = np.array([sinPhi*cosPsi - sinPsi*cosPhi*cosTheta, -cosPhi*cosPsi - sinPsi*sinPhi*cosTheta, sinPsi*sinTheta])
    ey_wave = np.array([-sinPhi*sinPsi - cosPsi*cosPhi*cosTheta, cosPhi*sinPsi - cosPsi*sinPhi*cosTheta, sinTheta*cosPsi])

    ### compute cartesian vector for the line-of-sight to the source
    n = np.array([sinTheta*cosPhi, sinTheta*sinPhi, cosTheta])

    ### compute detector matrix
    freqsT = 2j*np.pi*freqs*T ### this convention should match what is in LAL : x(f) = \int dt e^{-2\pi i f t} x(t)

    # factor of 1/2 is for normalization
    dV_xx = 0.5 * __D__(freqsT, np.sum(np.outer(ex, np.ones_like(theta))*n, axis=0))
    dV_yy = 0.5 * __D__(freqsT, np.sum(np.outer(ey, np.ones_like(theta))*n, axis=0))

    ### assemble these parts into antenna responses
    Fp = 0.
    Fx = 0.
    for i in xrange(3): ### FIXME? may be able to do this more efficiently vi np.array manipulations?
        exi = ex[i]
        eyi = ey[i]
        ex_wavei = ex_wave[i]
        ey_wavei = ey_wave[i]
        for j in xrange(3):
            ex_wavej = ex_wave[j]
            ey_wavej = ey_wave[j]
            Dij = dV_xx * exi*ex[j] - dV_yy * eyi*ey[j] ### compute matrix element for detector matrix

            Fp += Dij * (ex_wavei*ex_wavej - ey_wavei*ey_wavej) ### multiply by matrix element from wave polarizations
            Fx += Dij * (ex_wavei*ey_wavej + ey_wavei*ex_wavej)

    return Fp, Fx

#--- from bayesburst (known to be correct, supported for cross-validation purposes)

def const_antenna_response(theta, phi, psi, nx, ny):
    """
    computes the antenna patterns for detector arms oriented along nx and ny (cartesian vectors). 

    Antenna patterns are computed accoring to Eqn. B7 from Anderson, et all PhysRevD 63(04) 2003
    """
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)

    X = (sin_phi*cos_psi - sin_psi*cos_phi*cos_theta, -cos_phi*cos_psi - sin_psi*sin_phi*cos_theta, sin_psi*sin_theta)
    Y = (-sin_phi*sin_psi - cos_psi*cos_phi*cos_theta, cos_phi*sin_psi - cos_psi*sin_phi*cos_theta, sin_theta*cos_psi)

    ### iterate over x,y,z to compute F+ and Fx
    Fp = 0.
    Fx = 0.
    for i in xrange(3):
        nx_i = nx[i]
        ny_i = ny[i]
        Xi = X[i]
        Yi = Y[i]
        for j in xrange(3):
            Xj = X[j]
            Yj = Y[j]
            Dij = 0.5*(nx_i*nx[j] - ny_i*ny[j])
            Fp += (Xi*Xj - Yi*Yj)*Dij
            Fx += (Xi*Yj + Yi*Xj)*Dij

    return Fp, Fx
