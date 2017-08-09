import numpy as np

class Disk(object):
    '''Disk class
    
    Contains methods to generate radial flux and surface brightness profiles
    for a simple power-law axisymmetric optically-thin disk. 
    
    Parameters
    ----------
    lstar : float
        Stellar luminosity in L_sun.
    dist : float
        Distance to star in pc.
    alpha : float
        Disk surface density power-law index.
    z : float
        Number of zodis in disk, simple multiplier.
    wav : float
        Wavelength of interest in um.
    rin : float
        Disk inner radius in AU.
    rout : float
        Disk outer radius in AU.
    r0 : float
        Disk reference radius in AU.
    nr : int
        Number of radial locations in disk.
        
    Attributes
    ----------
    r : array
        Calculates value at centre of each radial bin.
    rarcs : array
        r in arcseconds.
    tbb : array
        Black body temperature of each bin centre r.
    sigzodi : float
        1 zodi surface density by defintion after integrating Kelsall model
        g(xi) over column to get multiplier of 0.629991 for 1.13e-7 volume 
        density.
    bigsig : array
        Disk flux density. AU^2/AU^2 in each of the bins.
    aream : array
        Area in AU^2 of each bin.
    sig : array
        Total dust area in each of the r bins in AU^2.
    '''        
    
    
    
    def __init__(self,lstar=None,dist=None,alpha=0.34,z=1.0,
               wav=11.0,rin=0.034422617777777775,rout=10.0,r0=1.0,nr=500):
        '''Initialise, default alpha=0.34, z=1.0, wav=11.0, 
        rin=0.034422617777777775, rout=10.0, r0=1.0, nr=500
        '''
        
        if lstar == None:
            raise Exception('Must give lstar')
            
        if dist == None:
            raise Exception('Must give dist')
        
        
        self.lstar = float(lstar)
        self.dist = float(dist)
        self.alpha = float(alpha)
        self.z = float(z)
        self.wav = float(wav)
        
        self.rin = float(rin)*np.sqrt(lstar)
        self.rout = float(rout)*np.sqrt(lstar)
        self.r0 = float(r0)*np.sqrt(lstar)

        nrb = nr + 1
        rb = np.linspace(self.rin,self.rout,nrb)
        drm = np.diff(rb)       
        self.r = 0.5*(rb[1:nrb]+rb[0:-1])
        self.rarcs = self.r/self.dist
        
        self.tbb = 278.3*(self.lstar**0.25)/np.sqrt(self.r)
        self.sigzodi = 7.11889e-8
        self.bigsig = self.sigzodi*self.z*(self.r/self.r0)**(-self.alpha)
        self.aream = 2*np.pi*self.r*drm
        self.sig = self.aream*self.bigsig
    
    
    
    def bnu_Jy_sr(self,wav,t):
        '''A function that returns the black body spectral radiance in Jy/sr 
        for a given temperature and at one or more given wavelengths (in um).
        
        Parameters
        ----------
        wav : float
            Wavelength in microns.
        t : array
            Temperature array in Kelvin.
            
        Returns
        -------
        ret : array
            Bnu in Jy/sr for given temperature and wavelength.
        '''
        
        k1 = 3.9728949e19     
        k2 = 14387.69         
        wav = float(wav)
        t = np.array(t,dtype=float)
        fact1 = k1/(wav**3)
        fact2 = k2/(t*wav)
        ret = np.array(fact1/(np.exp(fact2)-1.0))
        return ret
    
    
    
    def fnu_disk(self):
        '''Returns the flux profile for a simple power-law axisymmetric 
        optically-thin disk.
        '''
        
        fnudisk = (2.95e-10)*self.bnu_Jy_sr(self.wav,self.tbb)*(0.25*self.sig/np.pi)/self.dist**2        
        return fnudisk
    
    
    def snu_disk(self):
        '''Returns the surface brightness profile for a simple power-law
        optically thin disk.
        '''
        
        fnudisk = self.fnu_disk()
        snudisk = fnudisk/(self.aream/self.dist**2)
        return snudisk
