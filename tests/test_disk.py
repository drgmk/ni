import pytest
import numpy as np

from .context import ni

def number_of_bins_test():
    nr = 100
    nrb = nr + 1
    dsk = ni.disk.Disk(lstar=1.0,tstar=5800.0,dist=10.0,nr=nr)
    assert(len(dsk.nrb)==nrb)
    assert(len(dsk.rb)==nrb)
    assert(len(dsk.radii)==nrb)
    

def number_of_radii_test():
    nr = 100
    dsk = ni.disk.Disk(lstar=1.0,tstar=5800.0,dist=10.0,nr=nr)
    fnudisk = dsk.fnu_disk()
    snudisk = dsk.snu_disk()
    assert(len(dsk.nr)==nr)
    assert(len(dsk.rm)==nr)
    assert(len(dsk.drm)==nr)
    assert(len(dsk.rarcs)==nr)
    assert(len(dsk.tbb)==nr)
    assert(len(dsk.bigsig)==nr)
    assert(len(dsk.aream)==nr)
    assert(len(dsk.sig)==nr)
    assert(len(fnudisk)==nr)
    assert(len(snudisk)==nr)
    
def attributes_are_floats_test():
    dsk = ni.disk.Disk(lstar=1,tstar=5800,dist=10,alpha=1,z=1,wav=11,nr=100)
    assert(np.dtype(dsk.lstar)==float)
    assert(np.dtype(dsk.tstar)==float)
    assert(np.dtype(dsk.dist)==float)
    assert(np.dtype(dsk.alpha)==float)
    assert(np.dtype(dsk.z)==float)
    assert(np.dtype(dsk.wav)==float)
    assert(np.dtype(dsk.rin)==float)
    assert(np.dtype(dsk.rout)==float)
    assert(np.dtype(dsk.r0)==float)
    
    
def Bnu_Jy_sr_outputs_float_test():
    dsk = ni.disk.Disk(lstar=1.0,tstar=5800.0,dist=10.0)
    ret = dsk.bnu_Jy_sr(dsk.wav,dsk.tbb)
    assert(np.dtype(ret)==float)
    
    
def fnu_disk_outputs_float_test():
    dsk = ni.disk.Disk(lstar=1.0,tstar=5800.0,dist=10.0)
    fnudisk = dsk.fnu_disk()
    assert(np.dtype(fnudisk)==float)
    
    
def snu_disk_outputs_float_test():
    dsk = ni.disk.Disk(lstar=1.0,tstar=5800.0,dist=10.0)
    snudisk = dsk.snu_disk()
    assert(np.dtype(snudisk)==float)
    

    
    
    
