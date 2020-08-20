'''
Cluster.py

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

    def __init__(self, data : pd.DataFrame, name : str, error_type : str):
        """__init__.

        Class to contain attributes and functions for Gaia DR2 data of 
        a passed in stellar population. Loaded in from pandas.DataFrame
        of relevant data.

        Parameters
        ----------
        data : pd.DataFrame
            Stellar cluster data.
        name : str
            Name of stellar cluster.
        error_type : str 
            Uncertainty propagation to use : `astrometric` or `astrometric_photometric`
        """

        self.data = data
        self.name = name
        self.apply_cuts = False # Keeps track if cluster has photometric cuts applied
        self.age = None # Assigned age from best fit isochrone
        self.error_type = error_type # Fitting method used to data

    def __str__(self):
        """__str__.
        """
        return '{} cluster with {} stars'.format(self.name,len(self.photometry))


    @property
    def photometry(self) -> pd.DataFrame :
        """photometry.

        Reformats dataframe to contain just relevant info

        Parameters
        ----------
        self : Cluster

        Returns
        -------
        pd.DataFrame
        """
        if self.apply_cuts == False:
            g = self.data['phot_g_mean_mag'].values
            b = self.data['phot_bp_mean_mag'].values
            r = self.data['phot_rp_mean_mag'].values
            parallax = self.data['parallax'].values
            parallax_error = self.data['parallax_error'].values
            label = self.data['label'].values
            try:
                metallicity = self.data['m_h'].values
                extinction = self.data['a_g_val'].values
                color_excess = self.data['e_bp_min_rp_val'].values
            except:
                print('Metallicty, extinction, and color excess values not included')
            
            M_g = g + 5 - (5*np.log10(1000/parallax))
            g_flux = self.data['phot_g_mean_flux'].values
            g_flux_err = self.data['phot_g_mean_flux_error'].values
           
            if self.error_type == 'astrometric':
                error_M_g = error_abs_mag_parallax(M_g, parallax, parallax_error)
            elif self.error_type == 'astrometric_photometric':
                error_M_g = error_abs_mag(g_flux, g_flux_err, parallax, parallax_error)

#            error_M_g = 0.434*(parallax_error/parallax)
            b_r = b - r

            cluster_dict = {"abs_mag" : M_g, 
                    "color" : b_r, 
                    'metallicity': metallicity, 
                    'extinction' : extinction, 
                    'parallax' : parallax, 
                    'parallax_error' : parallax_error, 
                    'color_excess' : color_excess,
                    'abs_mag_err' : error_M_g, 
                    'label' : label}

            return pd.DataFrame(cluster_dict)

        else:
            return pd.DataFrame(self.data)

    @photometry.setter
    def photometry(self, data : pd.DataFrame):
        """photometry.

        Parameters
        ----------
        data : pd.DataFrame
            Data to be set for photometry.
        """
        self.data = data 
    
    def photometric_cuts(self, cut = None):
        """photometric_cuts.
        Apply photometric cuts

        Parameters
        ----------
        cut : 
            Personalized cut in same format as default one defiend
        """

        df = self.photometry
        self.apply_cuts = True
        if cut == None:
            cut = (((df.color > -.02 ) & (df.color < 1.8) & (df.abs_mag > df.color*2.46+2.76)) |
                        ((df.color > 1.8 ) & (df.abs_mag > df.color*2.8 +2.16)) |
                        ((df.color > -0.5) & (df.color < 4  ) & (df.abs_mag < df.color*2.14-0.57)) |
                        ((df.color > 1.2 ) & (df.color < 3  ) & (df.abs_mag < df.color*1.11+0.66)))
        
        self.photometry = df[np.invert(cut)]

    def plot_cluster_cmd(self, show = True):
        '''
        Plots the cluster
        @param show -> bool: If function should display plot or not
        @return fig, ax : figure and axis of plot
        '''
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
        #ax.set_xlim(min(color) - .15, max(color) + .15)
        #ax.set_ylim(min(abs_mag) - .15, max(abs_mag) + .15)
        plt.gca().invert_yaxis()
        if show == True:
            plt.show()
        return fig, ax 

def error_abs_mag_parallax(Mg, parallax, parallax_error):
    """error_abs_mag_parallax.
    Propagated error to be used in likelihood fit.

    Parameters
    ----------
    Mg : np.ndarray
        Absolute magnitude in the Gaia G band
    parallax : np.ndarray
        Parallax values
    parallax_error : np.ndarray
        Parallax errors

    Returns
    -------
    error_abs_mag : np.ndarray
        Propagated error
    """

    error_abs_mag = np.abs(Mg)*np.sqrt((parallax_error/parallax)**2) 
    return error_abs_mag 



# TODO : This is wrong?
def error_abs_mag(gflux, gflux_error, parallax, parallax_error):

    """error_abs_mag.

    Propagates errors in magnitudes (flux) and parallax 
    to get estimate on absolute G-band mag. Possibly not
    the best estimation since magnitudes are not symmetric
    in similarity to flux-space.

    Parameters
    ----------
    gflux : np.ndarray
        flux from g band
    gflux_error : np.ndarray
        flux error from g band
    parallax : np.ndarray 
        parallax
    parallax_error : np.ndarray
        parralax error

    Returns
    -------
    sigma_abs_mag : np.ndarray
        Propagated error

    References: (1) https://arxiv.org/pdf/1804.09368.pdf
                (2) http://gaia.ari.uni-heidelberg.de/gaia-workshop-2018/files/Gaia_DR2_photometry.pdf

    """
    sigma_g_zp = 0.0018 # Uncertainty in G mag zeropoint (1)
    sigma_g = np.sqrt((1.086*gflux_error/gflux)**2 + sigma_g_zp**2)

    # Hopefully this error propagation is correct
    # Obtained by propagating errors for M_g = G + 5 - 5log10(1000/parallax)
    sigma_abs_mag = np.sqrt(parallax_error**2 + ((5*sigma_g)/parallax/np.log(10)))
    #plt.figure()
    #plt.hist(sigma_abs_mag, edgecolor = 'black', bins = int(np.sqrt(len(sigma_abs_mag))))
    #plt.show()
    return sigma_abs_mag
