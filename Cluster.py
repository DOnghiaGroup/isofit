'''
Cameren Swiggum (swiggum2@wisc.edu)

Class representing a stellar cluster
Meant for Gaia DR2 data 
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import scipy.interpolate as interp

class cluster:

    def __init__(self, data : pd.DataFrame, name : str):
        self.data = data
        self.name = name
        self.apply_cuts = False

    def __str__(self):
        return '{} cluster with {} stars'.format(self.name,len(self.photometry))

    @property
    def photometry(self) -> pd.DataFrame :
        '''
        Reformats dataframe to contain just relevant info
        '''
        if self.apply_cuts == False:
            g = self.data['phot_g_mean_mag'].values
            b = self.data['phot_bp_mean_mag'].values
            r = self.data['phot_rp_mean_mag'].values
            parallax = self.data['parallax'].values
            label = self.data['label'].values
            try:
                metallicity = self.data['m_h'].values
                extinction = self.data['a_g_val'].values
                color_excess = self.data['e_bp_min_rp_val'].values
            except:
                print('Metallicty, extinction, and color excess values not included')
            
            M_g = g + 5 - (5*np.log10(1000/parallax))
            b_r = b - r

            cluster_dict = {"abs_mag" : M_g, "color" : b_r, 'metallicity': metallicity, 'extinction' : extinction, 'color_excess' : color_excess, 'label' : label}
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

    def plot_cluster_cmd(self, show = True):
        color = self.photometry.color.values
        abs_mag = self.photometry.abs_mag.values
        xy = np.vstack([color,abs_mag])
        number_density = gaussian_kde(xy)(xy)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(color,abs_mag,c=number_density,s=1)
        ax.set_title(self.name + ' CMD')
        ax.set_xlabel(r'$G_{Bp} - G_{Rp}$')
        ax.set_ylabel(r'$M_g$')
        ax.set_xlim(min(color) - .15, max(color) + .15)
        ax.set_ylim(min(abs_mag) - .15, max(abs_mag) + .15)
        plt.gca().invert_yaxis()
        if show == True:
            plt.show()
        return fig, ax 

