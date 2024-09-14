import numpy as np
import math
import surfaces

class Telescope:
    def __init__(self, telname='magellan'):
        self.Dtel        = 650.00
        self.obscuration = 0.07    # Secondary obscuration
        self.Atel        = np.pi * (self.Dtel/2.0)**2 * (1.0-self.obscuration)
        self.nsurf       = 3
        self.surfaces    = []
        # These surfaces follow measurements reported for Clay by Osip/Palunas in March 2018
        self.surfaces.append(surfaces.OpticalSurface('M1','BareAl','clay_m1.txt','LCO','model'))
        self.surfaces.append(surfaces.OpticalSurface('M2','BareAl','clay_m2.txt','LCO','model'))
        self.surfaces.append(surfaces.OpticalSurface('M3','ZeCoat','clay_m3.txt','LCO','model'))

    def throughput(self, input_wave):
        composite_throughput = np.ones(len(input_wave))
        for i in range(self.nsurf):
            self.surfaces[i].loadThroughputTab()
            composite_throughput *= self.surfaces[i].throughput(input_wave)
        return composite_throughput
    
