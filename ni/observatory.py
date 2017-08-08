import numpy as np

class Observatory(object):
    '''Observatory class
    
    This class is used to contain the attributes of the telescope 
    interferometer being used and methods used to calculate the 
    transmission of an azimuthally integrated radial disk profile through
    the interferometer fringe pattern.
    
    The position angle is the disk position angle relative to the projected 
    fringes, so a disk with PA=0deg will be parallel to the transmission
    pattern and if highly inclined may not be detected.
    
    Parameters
    ----------
    wav : float, optional
        Wavelength to compute transmission at in microns. Default = 11.0
    observatory : string, optional
        'lbti' or 'keck'. Choosing observatory will intialise the appropriate
        observatory mirror spacing and lattitude.
        
    Attributes
    ----------
    mirror_spacing : float
        Mirror spacing of interferometer
    lat : float
        Lattitude of the given observatory

    Methods
    -------
    transmission(phi=None,n_theta=100,inc=None,omegalbti=None)
        Returns the transmission array of an azimuthally integrated 
        radial disk profile through the LBTI fringe pattern.
    t_null(phiproj)
        Returns the transmission of a projected point on the disk through
        the interferometer fringe pattern.
    fringe_spacing()
        Calculates the distance to the first transmission peak using the 
        initialised attributes.
    altitude(has,dec)
        Returns the altitude of the target star over a range of hour 
        angles.
    parallactic_angle(has,dec):
        Returns the position angles of the target star in degrees.        
    '''
    
    def __init__(self,wav=11.0,observatory="lbti"):
        '''Initialise, default wav=11.0, observatory="lbti"'''
        
        self.wav = wav
        
        if observatory is "lbti":
            self.mirror_spacing = 14.4
            self.lat = 32.7016
        elif observatory is "keck":
            self.mirror_spacing = 85.0
            self.lat = 19.286
        else:
            raise Exception('Observatory not recognised')
            
        self.phi_null = self.fringe_spacing()
    
        
    def altitude(self,has,dec):
        """Returns the altitude of the target star over a range of hour 
        angles in degrees.
        
        Parameters
        ----------
        has : array
            Input hour angles in degrees as an array. 
        dec : float
            Star declination in degrees.
            
        Returns
        -------
        h : array
            Altitude of star in degrees.
        """
        delta = np.deg2rad(dec)
        has = np.deg2rad(has)
        latrad = np.deg2rad(self.lat)        
        h = np.arcsin(np.sin(latrad)*np.sin(delta) + np.cos(latrad) * 
                      np.cos(delta)*np.cos(latrad))
        h = np.rad2deg(h)
        return h
        
    
    def parallactic_angle(self,has,dec,has_input='deg'):
        '''Returns the position angles of the target star in degrees.
        default hasinput='deg'
        
        Parameters
        ----------
        has : array
            Input hour angles in degrees as an array.
        dec : float
            Star declination in degrees.
        has_input : string, optional
            User defines whether the hour angles (has) input is in degrees 'deg',
            or in hour angles 'ha'.
            
        Returns
        -------
        va : array
            Position angles of star in degrees.
        '''
        
        if has_input == 'ha':
            # Convert hour angles to degrees unless given in degrees.
            has *= 360/24.0

        delta = np.deg2rad(dec)
        h = self.altitude(has,dec)
        has = np.deg2rad(has)
        latrad = np.deg2rad(self.lat)
        cosva = (np.sin(latrad)*np.cos(delta) - np.cos(latrad) * 
                 np.sin(delta)*np.cos(has))/np.sqrt(1-np.sin(h)**2)
        sinva = np.sin(has) * np.cos(latrad)/np.sqrt(1-np.sin(h)**2)
        va = np.rad2deg(np.arctan2(sinva,cosva))
        return va
    
    
    def fringe_spacing(self):
        '''Calculates the distance to the first transmission peak'''
        
        phi_null = (648000.0/np.pi)*(self.wav)*(1e-6)/(2*self.mirror_spacing)
        self.phi_null = phi_null
        return phi_null

    
    def t_null(self,phiproj):
        '''Returns the transmission of a projected point on the disk through
        the interferometer fringe pattern. Equation (8).
        
        This function is called within the transmission method with the 
        calculated projected radial points to return their transmission as an 
        array.
        
        By calling this function with just the phi array, it will return the 
        sin^2 transmission pattern of the interferometer.
        
        
        Parameters
        ----------
        phiproj : array
            Input array of projected radial points.
            
        Returns
        -------
        T_inter : array
            Transmission at radial points on a disk through the fringe pattern.
        '''
        
        T_inter = np.sin((np.pi*phiproj)/(2*self.phi_null))**2
        return T_inter
    
    
    def transmission(self,phi=None,n_theta=100,inc=None,omegalbti=None):
        '''Default, n_theta=100.
        
        Returns the transmission array of an azimuthally integrated 
        radial disk profile through the LBTI fringe pattern.
        
        Parameters
        ----------
        phi : array
            'numpy' array of radii over distance to star. Normally defined
            phi = Disk().rarcs
        n_theta : int
            Number of points between 0 and 2*pi.
        inc : array
            'numpy' array of length 1 or greater for the inclination of the 
            disk. Must have same length as omegalbti.
        omegalbti : array
            'numpy' array of length 1 or greater for the position angle of 
            disk relative to fringe pattern. Must have same length as inc.
            
        Returns
        -------
        Tarray : array
            Transmission array of an azimuthally integrated radial disk profile
            through the LBTI fringe pattern. Dimensions [n_phi,n_orientations].
        '''
        
        if phi is None:
            raise Exception('phi: Must enter a numpy array of length 1 or greater')
            
        if n_theta is None:
            raise Exception('n_theta: Must enter an integer value for the number of \
                            equally seperated points between 0 and 2*pi')
            
        if inc is None:
            raise Exception('inc: Must enter a numpy array of length 1 or greater')
            
        if omegalbti is None:
            raise Exception('omegalbti: Must enter a numpy array of length 1 or greater')   
        
        
        self.phi = phi
        self.n_theta = n_theta
        #self.theta = np.linspace(0,2*np.pi,self.theta)
        self.theta = (2*np.pi)/(self.n_theta) * np.linspace(0,self.n_theta-1,self.n_theta)
        self.inc = inc
        self.omegalbti = omegalbti
        
        # Tile the phi array with dimensions n_phi to a square matrix with 
        # dimensions [n_phi,n_orientations]
        phi_proj_part1_tile1 = np.tile(phi,(len(self.inc),1))
        
        # Tile again to create phi matrix with dimensions 
        # [n_phi,n_orientations,n_theta]
        phi_proj_part1_tile2 = np.tile(phi_proj_part1_tile1,
                                       (self.theta.size,1,1)).transpose() 
        
        # Bracket part of eq. (6) 
        phi_proj_part2 = np.add(np.outer(np.sin(self.omegalbti),np.cos(self.theta)),np.outer(np.multiply(np.cos(self.omegalbti),np.cos(self.inc)),np.sin(self.theta)))
        
        # Tile equation 6 part to match phi matrix dimensions, 
        # [n_phi,n_orientations,n_theta]
        phi_proj_part2_tile = np.tile(phi_proj_part2,(len(self.phi),1,1))
        
        # phi matrix times the bracket part of eq. (6) 
        phi_proj =  np.multiply(phi_proj_part1_tile2,phi_proj_part2_tile)
        
        # Insert projected radial points phi_proj into eq. (8) in t_null 
        # method
        Tnull = self.t_null(phi_proj)
        
        # Avergae each orientation over the theta values (0 to 2*pi) to give
        # transmission as a function of radial points phi for each orientation
        # Matrx dimensions [n_r,n_orientation].
        Tarray = np.mean(Tnull,axis=2)
        
        return Tarray

