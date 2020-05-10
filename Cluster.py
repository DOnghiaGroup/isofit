'''
Cameren Swiggum (swiggum2@wisc.edu)

Class representing a stellar cluster
Meant for Gaia DR2 data 
'''

import numpy as np
import pandas as pd

class cluster:

    def __init__(self, data : pd.DataFrame, name : str):
        self.data = data
        self.name = name
        self.apply_cuts = False

    def __str__(self):
        return '{} cluster with {} stars'.format(self.name,len(self.photometry))

    @property
    def photometry(self) -> pd.DataFrame :
        if self.apply_cuts == False:
            g = self.data['phot_g_mean_mag'].values
            b = self.data['phot_bp_mean_mag'].values
            r = self.data['phot_rp_mean_mag'].values
            parallax = self.data['parallax'].values
            label = self.data['label'].values
            
            M_g = g + 5 - (5*np.log10(1000/parallax))
            b_r = b - r

            cluster_dict = {"abs_mag" : M_g, "color" : b_r, 'label' : label}
            return pd.DataFrame(cluster_dict)

        else:
            return pd.DataFrame(self.data)

    @photometry.setter
    def photometry(self, data : pd.DataFrame):
        self.data = data 
    
    def photometric_cuts(self, cut = None):
        df = self.photometry
        self.apply_cuts = True
        if cut == None:
            cut = (((df.color > 0.3 ) & (df.color < 1.8) & (df.abs_mag > df.color*2.46+2.76)) |
                        ((df.color > 1.8 ) & (df.abs_mag > df.color*2.8 +2.16)) |
                        ((df.color > -0.5) & (df.color < 4  ) & (df.abs_mag < df.color*2.14-0.57)) |
                        ((df.color > 1.2 ) & (df.color < 3  ) & (df.abs_mag < df.color*1.11+0.66)))
        
        self.photometry = df[np.invert(cut)]



    
         




