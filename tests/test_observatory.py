import pytest
import numpy as np

from .context import ni
    
def lbti_mirror_spacing_test():
    observatory = "lbti"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.mirror_spacing == 14.4)
    
    
def keck_mirror_spacing_test():
    observatory = "keck"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.mirror_spacing == 85.0)


def lbti_lattitude_test():
    observatory = "lbti"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.lat == 32.7016)

    
def keck_lattitude_test():
    observatory = "keck"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.lat == 19.286)

    
def Tarray_single_orientation_shape_test():
    obs = ni.observatory.Observatory()
    Tarray = obs.transmission(phi=np.array([1,2,3,4,5]),n_theta=100,
                              inc=np.array([0]),omegalbti=np.array([0]))
    radii_length = len(Tarray[:,0])
    orientation_length = len(Tarray[0,:])
    assert(radii_length == len(obs.phi))
    assert(orientation_length == len(obs.inc))

    
def Tarray_multiple_orientation_shape_test():
    obs = ni.observatory.Observatory()
    Tarray = obs.transmission(phi=np.array([1,2,3,4,5]),n_theta=100,
                              inc=np.array([1,2,3,4]),omegalbti=np.array([1,2,3,4]))
    radii_length = len(Tarray[:,0])
    orientation_length = len(Tarray[0,:])
    assert(radii_length == len(obs.phi))
    assert(orientation_length == len(obs.inc))
    
