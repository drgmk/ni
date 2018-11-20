import numpy as np

def orientations(omega,inc,HA):
    '''Function that generates two long arrays to be inputted into the 
    observation transmission method such that the output orientations accounts for 
    every hour angle across each omega across every inclination.
    '''
    omegalbtilen = len(omega)*len(HA)
    omegalbti = []
    for o in range(0,len(omega)):
        for p in range(0,len(HA)):
            omegalbti.append(omega[o] + HA[p])

    incall = np.repeat(inc,omegalbtilen,axis=0)            
    omegalbtiall = np.tile(np.array(omegalbti),len(inc))

    return (incall,omegalbtiall)

