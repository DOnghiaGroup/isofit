'''
Cameren Swiggum (swiggum2@wisc.edu)

Class to represent data from Parsec/Colibri 
Isochroens. Fount at http://stev.oapd.inaf.it/cgi-bin/cmd
'''
import numpy as np
import pandas as pd

class isochrone:

    def __init__(self, color : np.ndarray, abs_mag : np.ndarray, age : float, metallicity : float):
        self.color = color 
        self.abs_mag = color 
        self.age = age 
        self.metallicity = metallicity
        self.best_fit = False # Will switch to True if isochrone represents the cluster data
        self.error = None 
        
    @property
    def error(self):
        return self.__error
    
    @error.setter
    def error(self, value):
        if self.__error != None:
            self.__error = value
        else:
            raise TypeError("Isochrone error attribute already had a value")
    
    def __str__(self):
        out_str = 'Isochrone with log age {}, metallicity {}, and is'.format(self.age, self.metallicity)

        if self.best_fit:
            out_str += ' THE BEST FIT!'
        elif not self.best_fit:
            out_str += ' not the best fit.'

        return out_str
    

