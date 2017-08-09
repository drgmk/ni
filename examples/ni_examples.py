import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


from ni.observatory import Observatory
from ni.disk import Disk


class NI_Examples(object):
    '''Series of example plots taken from:
    Kennedy, Grant M., et al. "EXO-zodi modeling for the large binocular 
    telescope interferometer." The Astrophysical Journal Supplement Series 
    216.2 (2015): 23.
    '''
    
    def null_transmission(self):
        '''Figure 3: Projected null transmission'''
        
        obs = Observatory()
        PHInull = obs.phi_null
        PHI = np.linspace(0,10*PHInull,1000)
        Tnull = obs.t_null(PHI)
        fig = plt.figure()
        plt.plot(PHI/PHInull,Tnull,'k-')
        plt.title('Fringe Pattern')
        plt.xlabel(r'$\frac{\phi}{\phi_{null}}$')
        plt.ylabel('Transmission')
        fig.savefig('fig3_null.png')
        
        
    def face_on_transmission(self):
        '''Figure 3: Transmission for face-on disk'''
        
        obs = Observatory()
        PHInull = obs.phi_null
        Tarray = obs.transmission(phi=np.linspace(0,10*PHInull,1000),
                                  n_theta=100,
                                  inc=np.array([0]),omegalbti=np.array([0]))
        PHI = obs.phi
        fig = plt.figure()
        plt.plot(PHI/PHInull,Tarray,'b-')
        plt.title('Face-on Disk Transmission')
        plt.xlabel(r'$\frac{\phi}{\phi_{null}}$')
        plt.ylabel('Transmission')
        fig.savefig('fig3_faceon.png')
        
    
    def random_transmission(self):
        '''Figure 3: Average transmission for 1000 random orientations'''
        
        obs = Observatory()
        PHInull = obs.phi_null
        Tarray = obs.transmission(phi=np.linspace(0,10*PHInull,1000),
                                  n_theta=100,
                                  inc=np.arccos(np.random.uniform(0,1.0,1000)),
                                  omegalbti=np.random.uniform(0,np.pi,1000))
        PHI = obs.phi
        fig = plt.figure()
        plt.plot(PHI/PHInull,np.mean(Tarray,axis=1),'r-')
        plt.title('Random Orientation Disk Transmission')
        plt.xlabel(r'$\frac{\phi}{\phi_{null}}$')
        plt.ylabel('Transmission')      
        fig.savefig('fig3_rand.png')
        
        
    def total_fnu_disk(self):
        '''Figure 4: Total flux density per unit radius'''
        
        dsk = Disk(lstar=1.0,dist=10.0,nr=2000)
        fnudisk = dsk.fnu_disk()
        r = dsk.r
        fig = plt.figure()
        plt.plot(r,fnudisk/((1e-3)*np.mean(np.diff(dsk.rarcs))),'-')
        plt.title('Face-on Total Flux')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Flux in radial bins')
        plt.xlim([0,5])
        fig.savefig('fig4_total_flux.png')
        
        
    def transmitted_fnu_disk(self):
        '''Figure 4: Transmitted flux denisty per unit radius'''
        
        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=2000)
        Tarray = obs.transmission(phi=dsk.rarcs,
                                  n_theta=100,
                                  inc=np.array([0]),omegalbti=np.array([0]))

        fnudisk = dsk.fnu_disk()
        r = dsk.r
        rarcs = dsk.rarcs
        fig = plt.figure()
        plt.plot(r,Tarray[:,0]*fnudisk/((1e-3)*np.mean(np.diff(rarcs))),'-')
        plt.title('Face-on Transmitted Flux')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Flux in radial bins')
        plt.xlim([0,5])
        fig.savefig('fig4_trans_flux.png')
        
        
    def total_snu_disk(self):
        '''Figure 1(c): Total surface brightness per unit radius'''

        dsk = Disk(lstar=1.0,dist=10.0,nr=1000)
        snudisk=dsk.snu_disk()
        r = dsk.r
        fig = plt.figure()
        plt.loglog(r,snudisk,'-')
        plt.title('Model Surface Brightness')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Surface Brightness')
        fig.savefig('fig1_total_SB.png')
        
        
    def transmitted_snu_disk(self):
        '''Figure 1(c): Transmitted surface brightness per unit radius'''

        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=2000)
        Tarray = obs.transmission(phi=dsk.rarcs,
                                  n_theta=100,
                                  inc=np.array([0]),omegalbti=np.array([0]))

        snudisk=dsk.snu_disk()
        r = dsk.r
        fig = plt.figure()
        plt.loglog(r,(np.mean(Tarray,axis=1)*snudisk),'-')
        plt.title('Model Surface Brightness')
        plt.xlabel('Radius (AU)')
        plt.ylabel('Surface Brightness')
        fig.savefig('fig1_trans_SB.png')
    
    
    def flux_orientations(self):
        '''Figure 5: Transmitted fluxes for random orientations'''

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
            
        fig = plt.figure()
        FTarray = np.array(FTarray)
        weights = np.ones_like(FTarray)/float(len(FTarray))
        plt.hist(FTarray/(1e-3),bins=21,weights=weights)
        plt.xlabel('Transmitted flux (mJy)')
        plt.ylabel('Fraction of orientations')
        fig.savefig('fig5.png')
            
        
    def HA_transmission(self):
        '''Figure 7: Instantaneous transmission through range of hour angles'''

        obs = Observatory()
        inclination = np.deg2rad([0,30,45,60,75,90])
        has = obs.parallactic_angle(has=np.linspace(-2.0,2.0,1000),dec=60,has_input='ha')
        fig = plt.figure()
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
            ax.plot(np.linspace(-2.0,2.0,1000),FTarray/FTmax,'-')
        
        plt.legend(['0','30','45','60','75','90'])
        plt.xlabel('HA (h)')
        plt.ylabel('Relative transmission')
        fig.savefig('fig7.png')
        
                
    def inc_PA(self):
        '''Figure 6(a): Instantaneous sensitivity relative to face on disk'''

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
                
        fig = plt.figure()
        levels = np.linspace(0, 1, 101)
        plt.contourf(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTarray/np.amax(FTarray),interpolation='bilinear', origin='lower',cmap=cm.hot,levels=levels)
        plt.colorbar()
        plt.xlabel('Position angle (deg)')
        plt.ylabel('Inclination (deg)')
        
        z = np.array([0.3,0.60,0.80,0.95,0.99])
        contour = plt.contour(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTarray/np.amax(FTarray),levels=z,colors='k')
        plt.clabel(contour)
        fig.savefig('fig6a.png')
        
    def inc_PA_average(self):
        '''Figure 6(b): 4 hour average relative sensitivity'''

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
        fig = plt.figure()
        levels = np.linspace(0, 1, 101)
        plt.contourf(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTtotalmean/np.amax(FTtotalmean),interpolation='bilinear', origin='lower',cmap=cm.hot,levels=levels)
        plt.colorbar()
            
        z = np.array([0.60,0.80,0.95,0.99])
        contour = plt.contour(np.rad2deg(OMEGALBTI),np.rad2deg(INC),FTtotalmean/np.amax(FTtotalmean),levels=z,colors='k')
        plt.clabel(contour)
        plt.xlabel('Position angle (deg)')
        plt.ylabel('Inclination (deg)')
        fig.savefig('fig6b.png')

    def Rin(self):
        '''Figure 9(a): Total and transmitted flux as function of rin'''

        obs = Observatory()

        fig = plt.figure()
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
        fig.savefig('fig9a.png')
        
        
    def Rout(self):
        '''Figure 9(b): Total and transmitted flux as function of rout'''

        obs = Observatory()
        fig = plt.figure()
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
        fig.savefig('fig9b.png')


    def zodi_orientations(self):
        '''Figure 10: Distribution of zodi levels for 5000 random orientations'''

        obs = Observatory()
        dsk = Disk(lstar=1.0,dist=10.0,nr=1000)
        fnudisk = dsk.fnu_disk()                
        FTarray = []
        tstar = 5800.0
        HA = np.linspace(-30,30,100)
        PA = np.deg2rad(obs.parallactic_angle(HA,60))
        
        for i in range(0,200):
            incinput = np.arccos(np.random.uniform(0,1.0,1))
            omegalbtiinput = np.random.uniform(0,np.pi,1)
            Tarray = obs.transmission(phi=dsk.rarcs,n_theta=100,
                                      inc=incinput,
                                      omegalbti=omegalbtiinput + PA)
            Tarray = np.mean(Tarray,axis=1)
            print(i)
            fnustar = 1.77*dsk.bnu_Jy_sr(obs.wav,tstar)*dsk.lstar*tstar**(-4)/(dsk.dist**2)
            FT = (fnudisk * Tarray)/fnustar
            FTsum = np.sum(FT)
            FTarray.append(FTsum)
        FTarray = np.array(FTarray)
        z = (1e-4)/FTarray
        fig = plt.figure()
        weights = np.ones_like(z)/float(len(z))
        plt.hist(z,weights=weights)
        plt.ylabel('Fraction of orientations')
        plt.xlabel('Zodi')
        fig.savefig('fig10.png')
        
        
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
