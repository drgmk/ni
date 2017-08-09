import pytest
import numpy as np

from .context import ni  

def test_number_of_radii():
    nr = 100
    dsk = ni.disk.Disk(lstar=1.0,dist=10.0,nr=nr)
    fnudisk = dsk.fnu_disk()
    snudisk = dsk.snu_disk()
    assert(len(dsk.r)==nr)
    assert(len(dsk.rarcs)==nr)
    assert(len(dsk.tbb)==nr)
    assert(len(dsk.bigsig)==nr)
    assert(len(dsk.aream)==nr)
    assert(len(dsk.sig)==nr)
    assert(len(fnudisk)==nr)
    assert(len(snudisk)==nr)
    
def test_attributes_are_floats():
    dsk = ni.disk.Disk(lstar=1,dist=10,alpha=1,z=1,wav=11,nr=100)
    assert(dsk.lstar == float(dsk.lstar))
    assert(dsk.dist == float(dsk.dist))
    assert(dsk.alpha == float(dsk.alpha))
    assert(dsk.z == float(dsk.z))
    assert(dsk.wav == float(dsk.wav))
    assert(dsk.rin == float(dsk.rin))
    assert(dsk.rout == float(dsk.rout))
    assert(dsk.r0 == float(dsk.r0))
    
    
def test_Bnu_Jy_sr_outputs_float():
    dsk = ni.disk.Disk(lstar=1,dist=10)
    ret = dsk.bnu_Jy_sr(dsk.wav,dsk.tbb)
    assert(ret.dtype == float)
    
    
def test_fnu_disk_outputs_float():
    dsk = ni.disk.Disk(lstar=1,dist=10)
    fnudisk = dsk.fnu_disk()
    assert(fnudisk.dtype == float)
    
    
def test_snu_disk_outputs_float():
    dsk = ni.disk.Disk(lstar=1,dist=10)
    snudisk = dsk.snu_disk()
    assert(snudisk.dtype == float)
    

    
    
    
