import numpy as np
import os

class OpticalSurface:
    def __init__(self,name,glass,coating_file,vendor,status):
        self.name             = name
        self.glass            = glass
        self.coating_file     = coating_file
        self.vendor           = vendor
        self.status           = status  # Values: placeholder,model,heritage,witness,asbuilt
        self.wv_tab           = np.array(0)
        self.transmission_tab = np.array(0)
        self.reflectance_tab  = np.array(0)

    def loadThroughputTab(self):
        if (os.path.isfile(self.coating_file)):
            coating_data = np.genfromtxt(self.coating_file,usecols=[0,1],
                                         dtype=None,names=['wv','throughput'])
        elif (os.path.isfile(os.environ['COATINGS_PATH']+self.coating_file)):
            full_path = os.environ['COATINGS_PATH']+self.coating_file
            coating_data = np.genfromtxt(full_path,usecols=[0,1],
                                         dtype=None,names=['wv','throughput'])
        self.wv_tab  = coating_data['wv']
        self.transmission_tab = coating_data['throughput']
        if (self.vendor == "ECI"): 
            self.transmission_tab = (100.-self.transmission_tab)/100.
        if (self.name == "Collimator"): 
            self.transmission_tab = self.transmission_tab/100.
        self.reflectance_tab  = 1-self.transmission_tab

    def throughput(self,input_wave):
        # Should perform checks to make sure that coatings have been loaded, etc
        # Should perform checks to make sure input and tabulated waves are in nm
        interp_output = np.interp(input_wave,self.wv_tab,self.transmission_tab)
        return interp_output

class OpticalElement:
    def __init__(self,name,nsurf,glass,coating_file,vendor,status):    
        self.name         = name
        self.nsurf        = nsurf
        self.glass        = glass
        self.coating_file = coating_file
        self.vendor       = vendor
        self.surfaces     = []
        for i in range(nsurf):
            print("   ", name, glass,coating_file,vendor,status)
            self.surfaces.append(OpticalSurface(name,glass,coating_file,vendor,status))
            self.surfaces[i].loadThroughputTab()

    def throughput(self, input_wave):
        composite_throughput = np.ones(len(input_wave))
        for i in range(self.nsurf):
            composite_throughput *= self.surfaces[i].throughput(input_wave)
        return composite_throughput

class DiffractionGrating:
    def __init__(self,name,ruling,cwl,alpha,blaze_file,vendor,status):
        print("   ", name,ruling,cwl,alpha,blaze_file,vendor,status)
        self.name       = name
        self.ruling     = ruling
        self.cwl        = cwl
        self.alpha      = alpha      # The AOI of incoming beam
        self.blaze_file = blaze_file
        self.vendor     = vendor
        self.status     = status
        self.blaze      = OpticalSurface(name,'VPH',blaze_file,vendor,status)
        self.blaze.loadThroughputTab()

