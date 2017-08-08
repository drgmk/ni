import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


from ni.observatory import Observatory
from ni.disk import Disk


class NI_Examples(object):
    
    
    def null_transmission(self):
        print('Figure 3: Projected transmission')
        obs = Observatory()
        PHInull = obs.phi_null
        PHI = np.linspace(0,10*PHInull,1000)
        Tnull = obs.t_null(PHI)
        plt.figure()
        plt.plot(PHI/PHInull,Tnull,'k-')
        plt.title('Fringe Pattern')
        plt.xlabel(r'$\frac{\phi}{\phi_{null}}$')
        plt.ylabel('Transmission')
        
        
    def face_on_transmission(self):
        print('Figure 3: Transmission for face-on')
        obs = Observatory()
        PHInull = obs.phi_null
        Tarray = obs.transmission(phi=np.linspace(0,10*PHInull,1000),
                                  n_theta=100,
                                  inc=np.array([0]),omegalbti=np.array([0]))
        PHI = obs.phi
        plt.figure()
        plt.plot(PHI/PHInull,Tarray,'b-')
        plt.title('Face-on Disk Transmission')
        plt.xlabel(r'$\frac{\phi}{\phi_{null}}$')
        plt.ylabel('Transmission')    
        
    
    def random_transmission(self):
        print('Figure 3: Transmission for random orientations')
        obs = Observatory()
        PHInull = obs.phi_null
        Tarray = obs.transmission(phi=np.linspace(0,10*PHInull,1000),
                                  n_theta=100,
                                  inc=np.arccos(np.random.uniform(0,1.0,1000)),
                                  omegalbti=np.random.uniform(0,np.pi,1000))
        PHI = obs.phi
        plt.figure()
        plt.plot(PHI/PHInull,np.mean(Tarray,axis=1),'r-')
        plt.title('Random Orientation Disk Transmission')
        plt.xlabel(r'$\frac{\phi}{\phi_{null}}$')
        plt.ylabel('Transmission')      
        
        
    def total_fnu_disk(self):
        print('Figure 4: Total flux density per unit radius')
        dsk = Disk(lstar=1.0,dist=10.0,nr=2000)
        fnudisk = dsk.fnu_disk()
        rm = dsk.rm
        plt.figure()
        plt.plot(rm,fnudisk/((1e-3)*np.mean(np.diff(dsk.rarcs))),'-')
        plt.title('Face-on Total Flux')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Flux in radial bins')
        plt.xlim([0,5])
        
        
    def transmitted_fnu_disk(self):
        print('Figure 4: Transmitted flux density per unit radius')
        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=2000)
        Tarray = obs.transmission(phi=dsk.rarcs,
                                  n_theta=100,
                                  inc=np.array([0]),omegalbti=np.array([0]))

        fnudisk = dsk.fnu_disk()
        rm = dsk.rm
        rarcs = dsk.rarcs
        plt.figure()
        plt.plot(rm,Tarray[:,0]*fnudisk/((1e-3)*np.mean(np.diff(rarcs))),'-')
        plt.title('Face-on Transmitted Flux')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Flux in radial bins')
        plt.xlim([0,5])
        
        
    def total_snu_disk(self):
        print('Figure 1(c): Total surface brightness per unit radius')
        dsk = Disk(lstar=1.0,dist=10.0,nr=1000)
        snudisk=dsk.snu_disk()
        rm = dsk.rm
        plt.figure()
        plt.loglog(rm,snudisk,'-')
        plt.title('Model Surface Brightness')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Surface Brightness')
        
        
    def transmitted_snu_disk(self):
        print('Figure 1(c): Transmitted surface brightness per unit radius')
        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=2000)
        Tarray = obs.transmission(phi=dsk.rarcs,
                                  n_theta=100,
                                  inc=np.array([0]),omegalbti=np.array([0]))

        snudisk=dsk.snu_disk()
        rm = dsk.rm 
        plt.figure()
        plt.loglog(rm,(np.mean(Tarray,axis=1)*snudisk),'-')
        plt.title('Model Surface Brightness')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Surface Brightness')
    
    
    def flux_orientations(self):
        print('Figure 5: Transmitted fluxes for random orientations')
        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=1000)
        fnudisk = dsk.fnu_disk()
        FTarray = []
        
        for i in range(0,10000):
            Tarray = obs.transmission(phi=dsk.rarcs,n_theta=100,
                                      inc=np.arccos(np.random.uniform(0,1.0,1)),
                                      omegalbti=np.random.uniform(0,np.pi,1))
            FT = fnudisk * Tarray[:,0]
            FTsum = np.sum(FT)
            FTarray.append(FTsum)
            
        plt.figure()
        FTarray = np.array(FTarray)
        weights = np.ones_like(FTarray)/float(len(FTarray))
        plt.hist(FTarray/(1e-3),bins=21,weights=weights)
        plt.xlabel('Transmitted flux (mJy)')
        plt.ylabel('Fraction of orientations')
            
        
    def HA_transmission(self):
        print('Figure 7: Instantaneous transmission through range of hour angles')
        obs = Observatory()
        inclination = np.deg2rad([0,30,45,60,75,90])
        has = obs.parallactic_angle(has=np.linspace(-2.0,2.0,1000),dec=60,has_input='ha')
        plt.figure()
        for i in inclination:
            
            dsk = Disk(lstar=1.0,dist=10.0,nr=500)
            Tarray = obs.transmission(phi=dsk.rarcs,n_theta=100,
                                      inc=np.ones(1000)*i,
                                      omegalbti=np.deg2rad(has))

            fnudisk = dsk.fnu_disk()
            FTarray = []
            for j in range(0,len(obs.inc)):
                
                FT = fnudisk * Tarray[:,j]
                FTsum = sum(FT)
                FTarray.append(FTsum)
                
            if i == 0:
                FTmax = max(FTarray)
            
            FTarray = np.array(FTarray)
            ax = plt.subplot()
            ax.plot(np.linspace(-2.0,2.0,1000),FTarray/FTmax,'.')
        
        plt.legend(['0','30','45','60','75','90'])
        plt.xlabel('HA (h)')
        plt.ylabel('Relative transmission')
        
                
    def inc_PA(self):
        print('Figure 6(a): Instantaneous sensitivty relative to face on disk')
        obs = Observatory()
        inc = np.linspace(0,np.pi/2,100)
        omegalbti=np.deg2rad(np.linspace(-90,90,100))
        [INC,OMEGALBTI] = np.meshgrid(inc,omegalbti)
        FTarray = np.empty([len(inc),len(omegalbti)])
        dsk = Disk(lstar=1.0,dist=10.0,nr=100)
        fnudisk = dsk.fnu_disk()
        for a in range(0,len(inc)):
                
            T = obs.transmission(phi=dsk.rarcs,
                                      n_theta=100,inc=np.array([inc[a]]),
                                      omegalbti=omegalbti)
            FT = np.tile(fnudisk,(len(inc),1)).transpose() * T
            FTsum = sum(FT)
            FTarray[:,a] = FTsum
                
        plt.figure()
        levels = np.linspace(0, 1, 101)
        plt.contourf(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTarray/np.amax(FTarray),interpolation='bilinear', origin='lower',cmap=cm.hot,levels=levels)
        plt.colorbar()
        plt.xlabel('Position angle (deg)')
        plt.ylabel('Inclination (deg)')
        
        z = np.array([0.3,0.60,0.80,0.95,0.99])
        contour = plt.contour(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTarray/np.amax(FTarray),levels=z,colors='k')
        plt.clabel(contour)

        
    def inc_PA_average(self):
        print('Figure 6(b): 4 hour average relative sensitivity')
        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0)
        fnudisk = dsk.fnu_disk()
        HA = np.linspace(-30,30,10)
        PA = obs.parallactic_angle(HA,dec=60)
        omega = np.deg2rad(np.linspace(-90,90,30))
        inc = np.linspace(0,np.pi/2,30)
        FTtotalarray = np.empty([len(inc),len(omega),len(HA)])
        for h in range(0,len(HA)):
            FTarray = np.empty([len(inc),len(omega)])
            for a in range(0,len(inc)):
                
                T = obs.transmission(phi=dsk.rarcs,
                                          n_theta=100,inc=np.array([inc[a]]),
                                          omegalbti=omega + np.deg2rad(PA[h]))
                FT = np.tile(fnudisk,(len(inc),1)).transpose() * T
                FTsum = sum(FT)
                FTarray[:,a] = FTsum
            
            FTtotalarray[:,:,h] = FTarray
            
   
        FTtotalmean = np.mean(FTtotalarray,axis=2)
        
        [INC,OMEGALBTI] = np.meshgrid(inc,omega)
        plt.figure()
        levels = np.linspace(0, 1, 101)
        plt.contourf(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTtotalmean/np.amax(FTtotalmean),interpolation='bilinear', origin='lower',cmap=cm.hot,levels=levels)
        plt.colorbar()
            
        z = np.array([0.60,0.80,0.95,0.99])
        contour = plt.contour(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTtotalmean/np.amax(FTtotalmean),levels=z,colors='k')
        plt.clabel(contour)
        plt.xlabel('Position angle (deg)')
        plt.ylabel('Inclination (deg)')


    def Rin(self):
        print('Figure 9(a): Total and transmitted flux as function of rin')
        obs = Observatory()

        plt.figure()
        totarray = []
        transarray = []
        Rin = np.logspace(-2,1,100)
        for r in range(0,len(Rin)):
            
            dsk = Disk(lstar=1.0,dist=10.0,rin=Rin[r],nr=500)
            Tarray = obs.transmission(phi=dsk.rarcs,n_theta=100,
                                      inc=np.array([0]),omegalbti=np.array([0]))
            fnudisk = dsk.fnu_disk()
            transarray.append(sum(fnudisk*Tarray[:,0]))
            totarray.append(sum(fnudisk))
            
        plt.loglog(Rin,totarray,'-')
        plt.loglog(Rin,transarray,'-')
        plt.axis([0.01,10,10e-8,10e-4])      
        plt.xlabel('$R_{in} (AU)$')
        plt.ylabel('Total or transmitted disk flux (Jy)')
        
        
    def Rout(self):
        print('Figure 9(b): Total and transmitted flux as function of rout')
        obs = Observatory()
        plt.figure()
        totarray = []
        transarray = []
        Rout = np.logspace(-2.1,2,100)
        for r in range(0,len(Rout)):
            dsk = Disk(lstar=1.0,dist=10.0,rout=Rout[r],nr=5000)
            Tarray = obs.transmission(phi=dsk.rarcs,n_theta=100,
                                      inc=np.array([0]),omegalbti=np.array([0]))
            fnudisk = dsk.fnu_disk()
            transarray.append(sum(fnudisk*Tarray[:,0]))
            totarray.append(sum(fnudisk))
        plt.loglog(Rout,totarray,'-')
        plt.loglog(Rout,transarray,'-')
        plt.axis([0.05,10,10e-8,10e-4])
        plt.xlabel('$R_{out} (AU)$')
        plt.ylabel('Total or transmitted disk flux (Jy)')


    def zodi_orientations(self):
        print('Figure 10: Distribution of zodi levels for 5000 random orientations')
        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=1000)
        fnudisk = dsk.fnu_disk()                
        FTarray = []
        tstar = 5800.0
        for i in range(0,5000):
            Tarray = obs.transmission(phi=dsk.rarcs,n_theta=100,
                                      inc=np.arccos(np.random.uniform(0,1.0,1)),
                                      omegalbti=np.random.uniform(0,np.pi,1))

            fnustar = 1.77*dsk.bnu_Jy_sr(obs.wav,tstar)*dsk.lstar*tstar**(-4)/(dsk.dist**2)
            FT = (fnudisk * Tarray[:,0])
            FTsum = np.sum(FT)/fnustar
            FTarray.append(FTsum)
        FTarray = np.array(FTarray)
        z = (1e-4)/FTarray
        plt.figure()
        weights = np.ones_like(z)/float(len(z))
        plt.hist(z,weights=weights)
        plt.ylabel('Fraction of orientations')
        plt.xlabel('Zodi')
        
et = NI_Examples()

et.null_transmission()
et.face_on_transmission()
et.random_transmission()
et.total_fnu_disk()
et.transmitted_fnu_disk()
et.total_snu_disk()
et.transmitted_snu_disk()
et.flux_orientations()
et.HA_transmission()
et.inc_PA()
et.inc_PA_average()
et.Rin()
et.Rout()
et.zodi_orientations()
