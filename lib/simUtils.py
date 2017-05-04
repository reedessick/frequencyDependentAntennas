description = "a module that houses simulation and utility routines and classes"
author = "Reed Essick"

#-------------------------------------------------

from lal.lal import C_SI as c

import numpy as np
from scipy.special import erfinv
import healpy as hp

import freqDepAntennas as ant
import psds

import pickle

#-------------------------------------------------

twopi = 2*np.pi
twoIpi = 2j*np.pi

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
### defiine common PSDs

known_psds = dict((name, PSD(psd['freqs'], psd['vals'])) for name, psd in \
    [
        ('aLIGO', psds.aLIGO),
        ('aLIGO_O1', psds.aLIGO_O1),
        ('aLIGO_O2', psds.aLIGO_O2),
        ('aLIGO_O3', psds.aLIGO_O3),
        ('aLIGO_design', psds.aLIGO_design), 
        ('aPlus', psds.aPlus),
        ('aPlus_sqzonly', psds.aPlus_sqzonly),
        ('aVirgo', psds.aVirgo),
        ('aVirgo_sqz', psds.aVirgo_sqz),
        ('aVirgo_wb', psds.aVirgo_wb),
        ('CE', psds.CE),
        ('CE_wb', psds.CE_wb),
        ('ET', psds.ET),
        ('Voyager', psds.Voyager),
    ]
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
        phs = twoIpi*freqs*self.dt(theta, phi)

        return (Fp*hpf + Fx*hxf)*np.exp(phs)

    def dt(self, theta, phi):
        """
        time delay relative to geocenter
        """
        sinTheta = np.sin(theta)
        n = -np.array([np.cos(phi)*sinTheta, np.sin(phi)*sinTheta, np.cos(theta)])
        return np.sum(self.r*n)

    def drawNoise(self, freqs):
        """
        simulate Gaussian noise from this detector's PSD
        """
        amp = np.random.randn(*freqs.shape)*self.PSD(freqs)**0.5
        phs = np.random.rand(*freqs.shape)
        return amp*np.cos(phs) + 1j*amp*np.sin(phs)

#------------------------
### define common detectors

known_detectors = []

for name in ['aLIGO', 'aLIGO_O1', 'aLIGO_O2', 'aLIGO_O3', 'aLIGO_design', 'aPlus', 'aPlus_sqzonly', 'CE', 'CE_wb', 'Voyager']:
    known_detectors += [
        Detector(
            name = "L-"+name,
            ex = np.array((-0.9546, -0.1416, -0.2622)),
            ey = np.array((+0.2977, -0.4879, -0.8205)),
            r = np.array((-0.074276, -5.496284, +3.224257))*1e6/c,
            L = 4e3,
            PSD = known_psds[name],
        ),
        Detector(
            name = "Llong-"+name,
            ex = np.array((-0.9546, -0.1416, -0.2622)),
            ey = np.array((+0.2977, -0.4879, -0.8205)),
            r = np.array((-0.074276, -5.496284, +3.224257))*1e6/c,
            L = 4e4,
            PSD = known_psds[name],
        ),
        Detector(
            name = "H-"+name,
            ex = np.array((-0.2239, +0.7998, +0.5569)),
            ey = np.array((-0.9140, +0.0261, -0.4049)),
            r = np.array((-2.161415, -3.834695, +4.600350))*1e6/c,
            L = 4e3,
            PSD = known_psds[name],
        ),
        Detector(
            name = "Hlong-"+name,
            ex = np.array((-0.2239, +0.7998, +0.5569)),
            ey = np.array((-0.9140, +0.0261, -0.4049)),
            r = np.array((-2.161415, -3.834695, +4.600350))*1e6/c,
            L = 4e4,
            PSD = known_psds[name],
        ),
    ]

for name in ['aVirgo', 'aVirgo_sqz', 'aVirgo_wb', 'ET']:
    known_detectors += [
        Detector(
            name = "V-"+name,
            ex = np.array((-0.7005, +0.2085, +0.6826)),
            ey = np.array((-0.0538, -0.9691, +0.2408)),
            r = np.array((+4.546374, +0.842990, +4.378577))*1e6/c,
            L = 3e3,
            PSD = known_psds[name],
        ),
    ]

known_detectors = dict((detector.name, detector) for detector in known_detectors)

#-------------------------------------------------

def h2hAtT(freqs, h, t0):
    """
    add in phasing for time-at-coalescence
    """
    return h*np.exp(twoIpi*freqs*t0)

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

def symmetries( theta, phi, psi, iota, distance, t0 ):
    """
    enforce symmetries in parameters
    """
    ### enforce symmetries
    theta = theta%twopi
    if theta > np.pi:
        theta = twopi-theta 
        phi += np.pi
    phi = phi%twopi
    psi = psi%twopi
    iota = iota%twopi
    if iota > np.pi:
        iota = twopi-iota

    return theta, phi, psi, iota, distance, t0

def array_symmetries( theta, phi, psi, iota, distance, t0 ):
    '''
    a version of symmetries that works with np.ndarray objects
    '''
    theta = theta%twopi
    truth = theta>np.pi
    theta[truth] = twopi-theta[truth]
    phi[truth] += np.pi

    phi = phi%twopi

    psi = psi%twopi

    iota = iota%twopi
    truth = iota>np.pi
    iota[truth] = twopi-iota[truth]

    return theta, phi, psi, iota, distance, t0    

#------------------------

def IFOs2pole( ifo1, ifo2 ):
    '''
    return the line-of-sight polar direction in Earth-fixed coordinates.
    return theta, phi
    '''
    pole = known_detectors[ifo1].r - known_detectors[ifo2].r
    pole /= np.sum(pole**2)**0.5
    return np.arccos(pole[2]), np.arctan2(pole[1], pole[0]) ### store this in polar coordinates

def lineOfSight2ThetaPhi( theta_los, phi_los, pole=None ):
    """
    convert a direction specified in the line of sight coordinates into Earth-Fixed coordinates
    in Earth-fixed coordinates, the line-of-sight's polar direction is specified by pole=(thetaN, phiN)
    theta_los, phi_los are defined in the line-of-sight frame

    We define the transformation from the Earth-Fixed frame to the line-of-sight frame as 
        - a rotation about the z-axis by -phiN (to bring pole into the x-z plane)
        - a rotation about the y-axis by thetaN (to rotate z into N)
    this choice determines the phi_los=0 plane uniquely
    """
    if pole==None:
        return theta_los, phi_los

    thetaN, phiN = pole

    cosThetaN = np.cos(thetaN)
    sinThetaN = np.sin(thetaN)
    cosPhiN = np.cos(phiN)
    sinPhiN = np.sin(phiN)

    cosThetaLOS = np.cos(theta_los)
    sinThetaLOS = np.sin(theta_los)
    cosPhiLOS = np.cos(phi_los)
    sinPhiLOS = np.sin(phi_los)

    xLOS = sinThetaLOS*cosPhiLOS
    yLOS = sinThetaLOS*sinPhiLOS

    x = cosPhiN*cosThetaN*xLOS - sinPhiN*yLOS + cosPhiN*sinThetaN*cosThetaLOS
    y = sinPhiN*cosThetaN*xLOS + cosPhiN*yLOS + sinPhiN*sinThetaN*cosThetaLOS
    z = -sinThetaN*xLOS + cosThetaN*cosThetaLOS

    return np.arccos(z), np.arctan2(y, x)

def ThetaPhi2LineOfSight( theta, phi, pole=None ):
    """
    should perform the inverse of lineOfSight2ThetaPhi
    """
    if pole==None:
        return theta, phi

    thetaN, phiN = pole

    cosThetaN = np.cos(thetaN)
    sinThetaN = np.sin(thetaN)
    cosPhiN = np.cos(phiN)
    sinPhiN = np.sin(phiN)

    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    cosPhi = np.cos(phi)
    sinPhi = np.sin(phi)

    x = sinTheta*cosPhi
    y = sinTheta*sinPhi

    xLOS = cosPhiN*cosThetaN*x + sinPhiN*cosThetaN*y - sinThetaN*cosTheta
    yLOS = -sinPhiN*x + cosPhiN*y
    zLOS = sinThetaN*cosPhiN*x + sinThetaN*sinPhiN*y + cosThetaN*cosTheta

    return np.arccos(zLOS), np.arctan2(yLOS, xLOS)

#------------------------

def lnLikelihood( (theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=False, pole=None, verbose=False, **kwargs ):
    """
    log(likelihood) of this template against this data
    """
    ### enforce symmetries
    theta, phi, psi, iota, distance, t0 = symmetries(theta, phi, psi, iota, distance, t0)

    if pole!=None:
        ### theta, phi are supplied in the line-of-sight frame. We need to rotate back into Earth-fixed coordinates
        theta, phi = lineOfSight2ThetaPhi(theta, phi, pole)

    hpf, hxf = h2pol(h2hAtT(freqs, h, t0), iota, distance=distance)
    snrs = np.array([snr(freqs, detector, datum, hpf, hxf, theta, phi, psi, zeroFreq=zeroFreq) for detector, datum in zip(detectors, data)])
    if verbose: ### print the SNRs
        for snr2, detector in zip(snrs, detectors):
            print detector.name, snr2
        print ""
    return np.sum(0.5*snrs**2)

def likelihood( (theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=False, pole=None, **kwargs ):
    return np.exp(lnLikelihood((theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=zeroFreq, pole=pole))

#---

def lnPrior( (theta, phi, psi, iota, distance, t0), minDistance=0, maxDistance=1000, minT0=-1., maxT0=1., **kwargs ):
    """
    log(prior) of extrinsic parameters

    we allow some flexibility in angles (extend beyond the strict periodic boundaries), but we do enforce boundaries to keep them from running to infinity
    """
    if (theta<-np.pi) or (theta>twopi):
        return -np.infty
    if (phi<-twopi) or (phi>twopi):
        return -np.infty
    if (psi<-twopi) or (psi>twopi):
        return -np.infty
    if (iota<-np.pi) or (iota>twopi):
        return -np.infty
    if (distance<minDistance) or (distance>maxDistance):
        return -np.infty
    if (t0<minT0) or (t0>maxT0):
        return -np.infty

    theta, phi, psi, iota, distance, t0 = symmetries(theta, phi, psi, iota, distance, t0)
    return np.log(np.sin(theta)) + np.log(np.sin(iota)) + 2*np.log(distance)

def prior( (theta, phi, psi, iota, distance, t0), minDistance=0, maxDistance=1000, minT0=-1., maxT0=1., **kwargs ):
    return np.exp(lnPrior((theta, phi, psi, iota, distance), minDistance=minDistance, maxDistance=maxDistance, minT0=minT0, maxT0=maxT0))

#---

def lnPosterior( (theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=False, minDistance=0, maxDistance=1000, minT0=-1., maxT0=1., pole=None, **kwargs ):
    """
    log(posterior) of this template and extrinsic params against this data
    NOTE: this is not strictly normalized, so it isn't exactly the posterior
    """
    return lnLikelihood((theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=zeroFreq, pole=pole) + lnPrior((theta, phi, psi, iota, distance, t0), minDistance=minDistance, maxDistance=maxDistance, minT0=minT0, maxT0=maxT0)

def posterior( (theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=False, minDistance=0, maxDistance=1000, minT0=-1., maxT0=1., pole=None, **kwargs ):
    return np.exp(lnPosterior((theta, phi, psi, iota, distance, t0), freqs, data, h, detectors, zeroFreq=zeroFreq, minDistance=minDistance, maxDistance=maxDistance, minT0=minT0, maxT0=maxT0, pole=pole))

#---

def load_ensemble( path, header=False, min_sample=0, spacing=1 ):
    """
    a single place where we standardize how we load in data from ensemble.out files produced by localize
    """
    skiprows = 8 ### the number of rows in the header

    samples = np.genfromtxt(path, skiprows=skiprows, names=True)
    samples['theta'], samples['phi'], samples['psi'], samples['iota'], samples['distanceMpc'], samples['timeAtCoalescence'] = \
        array_symmetries(samples['theta'], samples['phi'], samples['psi'], samples['iota'], samples['distanceMpc'], samples['timeAtCoalescence'])

    ### discard everything befor min_sample
    inds = np.arange(len(samples))
    keys = 'k lnprob theta phi psi iota distanceMpc timeAtCoalescence'.split()
    ans = dict((key, np.array([])) for key in keys)
    for truth in [samples['k']==walker for walker in [int(_) for _ in sorted(set(samples['k']))]]:
        truth[np.arange(len(truth))<min_sample] = False
        these_inds = inds[truth][::spacing]
        for key in keys:
            ans[key] = np.concatenate((ans[key], samples[key][these_inds]))

    if header:
        file_obj = open(path, 'r')
        header = [file_obj.readline().strip() for _ in xrange(skiprows)]
#        return samples, header
        return ans, header

    else:
#        return samples
        return ans

def resume_ensemble( path ):
    """
    pick out the last points from the existing samples and return them as a formatted sampler
    """
    samples = load_ensemble(path)

    N = np.arange(len(samples['k'])) ### total length of the current array

    walkers = sorted(set(samples['k'])) ### figure out how many walkers there are
    Nwalkers = len(walkers)

    p0 = np.empty((Nwalkers, 6), dtype=float) ### figure out where chains ended
    for walker in walkers:
        ind = N[samples['k']==walker][-1] ### take the most recent sample for this walker
        p0[walker] = np.array([
            samples['theta'][ind], 
            samples['phi'][ind], 
            samples['psi'][ind], 
            samples['iota'][ind], 
            samples['distanceMpc'][ind], 
            samples['timeAtCoalescence'][ind],
        ])

    return p0, Nwalkers

#---

def dump_setup( path, args, kwargs):
    """
    write data needed for setting up a sampler to disk
    """
    file_obj = open(path,'w')
    pickle.dump(args, file_obj)
    pickle.dump(kwargs, file_obj)
    file_obj.close()

def load_setup( path ):
    """
    load data needed for setting up a sampler from disk
    """
    file_obj = open(path,'r')
    args = pickle.load(file_obj)
    kwargs = pickle.load(file_obj)
    file_obj.close()
    return args, kwargs

#-------------------------------------------------

def lineOfSightGrid( theta, phi, Ntheta, Nphi, sigmaTheta, pole=None ):
    '''
    assumes theta, phi are in Earth fixed coordinates.
    generates a grid in line of sight coordinates around the corresponding theta_los, phi_los
    returns that grid in Earth fixed coords
    return thetaGRID, phiGRID
    '''
    theta_los, phi_los = ThetaPhi2LineOfSight(theta, phi, pole=pole)

    # phi is easy, just wrap it around the azimuth
    phiGRID = np.linspace(0, 2*np.pi, Nphi+1)[:-1]

    # theta is a bit ad hoc, but should be a gaussian centered on the injected location
    if Ntheta%2:
        Ntheta += 1
    thetaGRID = erfinv(np.linspace(0,1,Ntheta/2+1)[:-1])
    thetaGRID = theta_los + sigmaTheta*np.concatenate((-thetaGRID[::-1], thetaGRID[1:]))
    thetaGRID = thetaGRID[(thetaGRID>=0)*(thetaGRID<=np.pi)]

    ### rotate back to Earth-fixed
    thetaGRID, phiGRID = np.meshgrid(thetaGRID, phiGRID)
    thetaGRID, phiGRID = lineOfSight2ThetaPhi(thetaGRID.flatten(), phiGRID.flatten(), pole=pole)

    phiGRID[phiGRID<0] += 2*np.pi ### we want these to be positive

    return thetaGRID, phiGRID

def psiGrid( Npsi ):
    '''
    returns a grid over psi
    '''
    return np.linspace(0, np.pi, Npsi+1)[:-1]

def timeAtCoalescenceGrid( t0, theta_inj, phi_inj, detector, theta, phi, Ntac, minT0=-1., maxT0=1., **kwargs ):
    '''
    returns grid over timeAtCoalescence
    we shift the grid around t0 to get the grid to be centered on the time-of-arival at detector when the signal comes from (theta,phi) instead of (theta_inj,phi_inj)
    '''
    ###     time at detector from injection       shift back to geocenter      sample around that
    return (t0+detector.dt(theta_inj, phi_inj)) - detector.dt(theta, phi) + np.linspace(minT0, maxT0, Ntac)

def iotaDistanceGrid( iota, distance, Niota, Ndistance, minDistance=1, maxDistance=1000, padding=2., **kwargs ):
    '''
    returns a grid over iota and distance

    return iotaGRID, distanceGRID
    '''
    cosIotaGRID = np.outer(np.linspace(-1, 1, Niota), np.ones(Ndistance)) ### spacing of inclinations

    ### maximum allowed "scaling constant" for placing distance grid
    ### take into account the detectability (malmquist prior -> gets rid of a lot of otherwise very densely sampled parameter-space
    cosIota2 = np.cos(iota)**2 
    dM3 = ( (padding*distance)**2 / ((1+cosIota2)**2+cosIota2) )**(3./2)

    ### minimum allowed "scaling constant" for distance grid
    dm3 = (minDistance/5**0.5)**3

    do = np.outer(np.ones(Niota), np.linspace(dm3, dM3, Ndistance)**(1./3)) ### constants for scaling relation

    cosIota2GRID = cosIotaGRID**2
    distanceGRID = do*((1+cosIota2GRID)**2 + cosIota2GRID)**0.5

    distanceGRID = distanceGRID.flatten()
    truth = (distanceGRID>=minDistance)*(distanceGRID<=maxDistance) ### exclude points outside of prior bounds

    return np.arccos(cosIotaGRID).flatten()[truth], distanceGRID[truth]
