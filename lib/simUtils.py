description = "a module that houses simulation and utility routines and classes"
author = "Reed Essick"

#-------------------------------------------------

from lal.lal import C_SI as c

import numpy as np

import freqDepAntennas as ant
import psds

#-------------------------------------------------

class PSD(object):
    """
    a represenation of a PSD that is callable
    """
    def __init__(self, freqs, vals):
        self.freqs = freqs
        self.vals = vals

    def __call__(self, freqs):
        return np.interp(freqs, self.freqs, self.vals)

#------------------------
### define common PSDs
aLIGO_design = PSD(
    psds.aLIGO_design['freqs'], 
    psds.aLIGO_design['vals'],
)
aVirgo_design = PSD(
    psds.aVirgo_design['freqs'], 
    psds.aVirgo_design['vals'],
)

#-------------------------------------------------

class Detector(object):
    """
    a representation of a detector
    knows how to project strains into it's output
    """

    def __init__(self, name, ex, ey, r, L, PSD):
        self.name = name
        self.ex = ex
        self.ey = ey
        self.r = r
        self.L = L
        self.T = L/c
        self.PSD = PSD

    def project(self, freqs, hpf, hxf, theta, phi, psi, zeroFreq=False):
        """
        project strains into this detector
        if zeroFreq, we use const_antenna_response instead of antenna_response
        """
        ### angular dependence from antenna patterns
        if zeroFreq:
            Fp, Fx = ant.const_antenna_response(theta, phi, psi, self.ex, self.ey)
        else:
            Fp, Fx = ant.antenna_response(theta, phi, psi, self.ex, self.ey, T=self.T, freqs=freqs)
            Fp = Fp[:,0]
            Fx = Fx[:,0]
        ### overall phase delay from extra time-of-flight
        ### r is measured in seconds
        ### FIXME: there could be a sign error here depending on our convention for FFT's...
        sinTheta = np.sin(theta)
        n = -np.array([np.cos(phi)*sinTheta, np.sin(phi)*sinTheta, np.cos(theta)])
        dt = np.sum(self.r*n)
        phs = 2j*np.pi*freqs*dt

        return (Fp*hpf + Fx*hxf)*np.exp(phs)

    def drawNoise(self, freqs):
        """
        simulate Gaussian noise from this detector's PSD
        """
        amp = np.random.randn(*freqs.shape)*self.PSD(freqs)**0.5
        phs = np.random.rand(*freqs.shape)
        return amp*np.cos(phs) + 1j*amp*np.sin(phs)

#------------------------
### define common detectors
LLO = Detector(
    name = "L",
    ex = np.array((-0.9546, -0.1416, -0.2622)),
    ey = np.array((+0.2977, -0.4879, -0.8205)),
    r = np.array((-0.074276, -5.496284, +3.224257))*1e6/c,
    L = 4e3,
    PSD = aLIGO_design,
)
LHO = Detector(
    name = "H",
    ex = np.array((-0.2239, +0.7998, +0.5569)),
    ey = np.array((-0.9140, +0.0261, -0.4049)),
    r = np.array((-2.161415, -3.834695, +4.600350))*1e6/c,
    L = 4e3,
    PSD = aLIGO_design,
)
Virgo = Detector(
    name = "V",
    ex = np.array((-0.7005, +0.2085, +0.6826)),
    ey = np.array((-0.0538, -0.9691, +0.2408)),
    r = np.array((+4.546374, +0.842990, +4.378577))*1e6/c,
    L = 3e3,
    PSD = aVirgo_design,
)

#-------------------------------------------------

def h2pol( h, iota, distance=1. ):
    '''
    map a waveform into polarization components and normalize by distance
    '''
    cos_iota = np.cos(iota)
    return h*0.5*(1+cos_iota**2)/distance, 1j*h*cos_iota/distance

#-------------------------------------------------

def snr( freqs, detector, data, hpf, hxf, theta, phi, psi, zeroFreq=False ):
    """
    computes the SNR of this template against this data
    does NOT maximize over the phase at coalescence
    """
    template = detector.project(freqs, hpf, hxf, theta, phi, psi, zeroFreq=zeroFreq)
    PSD = detector.PSD(freqs)
    deltaF = freqs[1]-freqs[0]
    return np.sum(deltaF*np.conjugate(data)*template/PSD).real / np.sum(deltaF*np.conjugate(template)*template/PSD).real**0.5

def cumsum_snr(freqs, detector, data, hpf, hxf, theta, phi, psi, zeroFreq=False ):
    """
    returns the cumulative sum of the snr as a function of frequency
    does NOT maximize over the phase at coalescence
    """
    template = detector.project(freqs, hpf, hxf, theta, phi, psi, zeroFreq=zeroFreq)
    PSD = detector.PSD(freqs)
    deltaF = freqs[1]-freqs[0]
    return np.cumsum(deltaF*np.conjugate(data)*template/PSD).real / np.sum(deltaF*np.conjugate(template)*template/PSD).real**0.5

#------------------------

def lnLikelihood( (theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=False, **kwargs ):
    """
    log(likelihood) of this template against this data
    """
    hpf, hxf = h2pol(h, iota)
    return np.sum([0.5*snr(freqs, detector, datum, hpf, hxf, theta, phi, psi, zeroFreq=zeroFreq)**2 for detector, datum in zip(detectors, data)])

def likelihood( (theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=False, **kwargs ):
    return np.exp(lnLikelihood((theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=zeroFreq))

#---

def lnPrior( (theta, phi, psi, iota, distance), minDistance=0, maxDistance=1000, **kwargs ):
    """
    log(prior) of extrinsic parameters
    """
    if (theta<0) or (theta>np.pi):
        return -np.infty
    if (phi<0) or (phi>2*np.pi):
        return -np.infty
    if (psi<0) or (psi>2*np.pi):
        return -np.infty
    if (iota<0) or (iota>np.pi):
        return -np.infty
    if (distance<minDistance) or (distance>maxDistance):
        return -np.infty

    return np.log(np.sin(theta)) + np.log(np.sin(iota)) + 2*np.log(distance)

def prior( (theta, phi, psi, iota, distance), minDistance=0, maxDistance=1000, **kwargs ):
    return np.exp(lnPrior((theta, phi, psi, iota, distance), minDistance=minDistance, maxDistance=maxDistance))

#---

def lnPosterior( (theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=False, minDistance=0, maxDistance=1000, **kwargs ):
    """
    log(posterior) of this template and extrinsic params against this data
    NOTE: this is not strictly normalized, so it isn't exactly the posterior
    """
    return lnLikelihood((theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=zeroFreq) + lnPrior((theta, phi, psi, iota, distance), minDistance=minDistance, maxDistance=maxDistance)

def posterior( (theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=False, minDistance=0, maxDistance=1000, **kwargs ):
    return np.exp(lnPosterior((theta, phi, psi, iota, distance), freqs, data, h, detectors, zeroFreq=zeroFreq, minDistance=minDistance, maxDistance=maxDistance))

#---

def load_ensemble( path ):
    """
    a single place where we standardize how we load in data from ensemble.out files produced by localize
    """
    samples = np.genfromtxt(path, skiprows=8, names=True)
    samples['phi'][samples['phi']>np.pi] -= 2*np.pi ### wrap for plotting in mollweide projection

    return samples
