description = "a module that houses simulation and utility routines and classes"
author = "Reed Essick"

#-------------------------------------------------

from lal.lal import C_SI as c
from lal import lalsimulation as lalsim

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

    def project(self, freqs, hpf, hxf, theta, phi, psi):
        """
        project strains into this detector
        """
        Fp, Fx = ant.antenna_resopnse( theta, phi, psi, detector.ex, detector.ey, T=detector.T, freqs=freqs)
        return Fp*hpF + Fx*hxf

    def drawNoise(self, freqs):
        """
        simulate Gaussian noise from this detector's PSD
        """
        amp = np.random.randn(freqs.shape)*self.PSD(freqs)**0.5
        phs = np.random.rand(freqs.shape)
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

def snr( freqs, detector, data, hpf, hxf, theta, phi, psi ):
    """
    computes the SNR of this template against this data
    """
    template = detector.project(freqs, hpf, hxf, theta, phi, psi)
    PSD = detector.PSD(freqs)
    deltaF = freqs[1]-freqs[0]
    return np.sum(deltaF*np.conjugate(data)*template/PSD).real / np.sum(deltaF*np.conjugate(template)*template/PSD)**0.5

def cumsum_snr(freqs, detector, data, hpf, hxf, theta, phi, psi ):
    """
    returns the cumulative sum of the snr as a function of frequency
    """
    return np.cumsum(deltaF*np.conjugate(data)*template/PSD).real / np.sum(deltaF*np.conjugate(template)*template/PSD)**0.5

def lnLikelihood( freqs, detector, data, hpf, hxf, theta, phi, psi ):
    """
    log(likelihood) of this template against this data
    """
    return 0.5*snr(freqs, PSD, data, hpf, hxf, theta, phi, psi)**2
