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
    
    '''
    Class to contain attributes and functions for Gaia DR2 data of 
    a passed in stellar population. Loaded in from pandas.DataFrame
    of relevant data.
    
    Attributes:
        data -> pd.DataFrame: DataFrame of stellar population. Needed Keys
            found in README.md
        name -> str: Given name of population
        apply_cuts -> bool: Keeps track if cluster has photometric cuts applied
        photometry -> pd.DataFrame: Reformatted container of data
    '''

    def __init__(self, data : pd.DataFrame, name : str, error_type):
        """__init__.

        Parameters
        ----------
        data : pd.DataFrame
            Stellar cluster data.
        name : str
            Name of stellar cluster.
        error_type : TODO 
            error_type
        """

        self.data = data
        self.name = name
        self.apply_cuts = False
        self.age = None
        self.error_type = error_type

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
        ax.set_xlim(min(color) - .15, max(color) + .15)
        ax.set_ylim(min(abs_mag) - .15, max(abs_mag) + .15)
        plt.gca().invert_yaxis()
        if show == True:
            plt.show()
        return fig, ax 

    def plot_metallicity_distribution(self):
        '''
        Shows cluster's metallicity distribution. Meant to be used
        before downloading isochrone grid so that a good estimation
        of M_H is entered.
        '''

        metallicity = self.photometry.metallicity.values
        metallicity = metallicity[~np.isnan(metallicity)]
        metallicity = metallicity[np.where(metallicity != -9999.)]
        N = len(metallicity)

        plt.figure()
        plt.title(self.name + ' N = ' + str(N))
        plt.xlabel('Metallicty [dex]')
        plt.ylabel('Count')
        plt.hist(metallicity, bins = int(np.sqrt(N)), edgecolor = 'black')
        plt.show()
    
    def plot_extinction_distribution(self):
        '''
        Shows cluster's distribution in Av to select
        a good approximation for downloading the
        isochrone grid
        '''
        
        Av = self.photometry.extinction.values
        Av = Av[~np.isnan(Av)]
        Av = Av[np.where(Av != -9999.)]
        N = len(Av)

        plt.figure()
        plt.title(self.name + ' N = ' + str(N)) 
        plt.xlabel('Extinction [mag]')
        plt.ylabel('Count')
        plt.hist(Av, bins = int(np.sqrt(N)), edgecolor = 'black')
        plt.show()

    
                  
    def plot_abs_mag_distribution(self):
        '''
        Plots the distribution of abs_mags
        '''
        abs_mag = self.photometry.abs_mag.values
        plt.figure()
        plt.xlabel(r'$M_G$')
        plt.ylabel('Count')
        plt.hist(abs_mag, bins = int(np.sqrt(len(abs_mag))), edgecolor = 'black')
        plt.show()

    def output_cluster(self, directory):
        self.photometry.to_csv(directory + self.name)


'''
Helper functions
'''

def error_abs_mag_parallax(Mg, parallax, parallax_error):
    '''
    Returns error on abs_mag solely based on parallax errors

    @param Mg : Absolute magnitudes in the G band
    @param parallax : Parallax values
    @param parallax_error : Errors in parallax
    @return error_abs_mag : Errors in absolute magntidue of G band
    '''
    error_abs_mag = np.abs(Mg)*np.sqrt((parallax_error/parallax)**2) 
    return error_abs_mag 



# TODO : This is wrong!
def error_abs_mag(gflux, gflux_error, parallax, parallax_error):
    '''
    
    Propagates errors in magnitudes (flux) and parallax 
    to get estimate on absolute G-band mag. Possibly not
    the best estimation since magnitudes are not symmetric
    in similarity to flux-space.

    References: (1) https://arxiv.org/pdf/1804.09368.pdf
                (2) http://gaia.ari.uni-heidelberg.de/gaia-workshop-2018/files/Gaia_DR2_photometry.pdf
    '''
    
    sigma_g_zp = 0.0018 # Uncertainty in G mag zeropoint (1)
    sigma_g = np.sqrt((1.086*gflux_error/gflux)**2 + sigma_g_zp**2)

    # Hopefully this error propagation is correct
    # Obtained by propagating errors for M_g = G + 5 - 5log10(1000/parallax)
    sigma_abs_mag = np.sqrt(parallax_error**2 + ((5*sigma_g)/parallax/np.log(10)))
    #plt.figure()
    #plt.hist(sigma_abs_mag, edgecolor = 'black', bins = int(np.sqrt(len(sigma_abs_mag))))
    #plt.show()
    return sigma_abs_mag
