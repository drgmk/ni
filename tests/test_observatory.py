import pytest
import numpy as np

from .context import ni
    
def test_lbti_mirror_spacing():
    observatory = "lbti"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.mirror_spacing == 14.4)
    
    
def test_keck_mirror_spacing():
    observatory = "keck"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.mirror_spacing == 85.0)


def test_lbti_lattitude():
    observatory = "lbti"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.lat == 32.7016)

    
def test_keck_lattitude():
    observatory = "keck"
    obs = ni.observatory.Observatory(observatory=observatory)
    assert(obs.lat == 19.286)

    
def test_Tarray_single_orientation_shape():
    obs = ni.observatory.Observatory()
    Tarray = obs.transmission(phi=np.array([1,2,3,4,5]),n_theta=100,
                              inc=np.array([0]),omegalbti=np.array([0]))
    radii_length = len(Tarray[:,0])
    orientation_length = len(Tarray[0,:])
    assert(radii_length == 5)
    assert(orientation_length == 1)

    
def test_Tarray_multiple_orientation_shape():
    obs = ni.observatory.Observatory()
    Tarray = obs.transmission(phi=np.array([1,2,3,4,5]),n_theta=100,
                              inc=np.array([1,2,3,4]),omegalbti=np.array([1,2,3,4]))
    radii_length = len(Tarray[:,0])
    orientation_length = len(Tarray[0,:])
    assert(radii_length == 5)
    assert(orientation_length == 4)
    
