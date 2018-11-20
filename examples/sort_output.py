"""
Example attempts to recreate Fig. 10 histogram using the orientations function 
and averaging over hour angles.
"""


import numpy as np
from orientation import orientations
import ni
import matplotlib.pyplot as plt

dsk = ni.disk.Disk(lstar=1.0,dist=10.0,nr=100)
obs = ni.observatory.Observatory()

omega = np.linspace(0,np.pi,50)
inc = np.arccos(np.linspace(0,1.0,50))
#HA = np.linspace(-30,30,10)
HA = np.array([0])
HA= obs.parallactic_angle(HA,60)

(incall,omegalbtiall) = orientations(omega,inc,HA)

fnudisk = dsk.fnu_disk()
fnudisk_tile = np.tile(fnudisk,(len(inc)*len(omega),1)).transpose()
tstar=5800.0
fnustar = 1.77 * dsk.bnu_Jy_sr(obs.wav,tstar) * dsk.lstar * tstar**(-4)/(dsk.dist**2)

Tarray = obs.transmission(phi=dsk.rarcs,inc=incall,omegalbti=omegalbtiall)

Tarray_HAmean = np.mean(Tarray.reshape(len(dsk.rarcs),len(omega)*len(inc),len(HA)),axis=2)

FT = Tarray_HAmean * fnudisk_tile / fnustar
FTmean = np.mean(FT,axis=0)

z = (1e-4)/FTmean


weights = np.ones_like(z)/float(len(z))
plt.hist(z,weights=weights,bins=300)
plt.xlim([0,5000])