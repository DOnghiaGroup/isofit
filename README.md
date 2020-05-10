# isofit.py
A simple single-age stellar population fitting procedure that uses supplied Gaia DR2 photometry/astrometry and fits a grid of Parsec-Colibri isochrones to 
the data. Example provided in example.py. Fit so far is only performed by finding the isochrone that minimizes square error. Plans are to incorporate photometric/astrometric error for a chi-square minimization. Grid of isochrones can be downloaded at http://stev.oapd.inaf.it/cgi-bin/cmd_3.3. Not tested thoroughly. 

#Cluster.py
A class to hold given data. Takes in a pandas.DataFrame of Gaia DR2 data with required columns with names: parallax, phot_g_mean_mag, phot_bp_mean_mag, 
phot_rp_mean_mag, label. Since this module currently can't automate the process of obtaining an isochrone grid, there are plotting procedures in this class that 
can help a user determine the extinction (A_v) and metallicty (M_H) values to enter at http://stev.oapd.inaf.it/cgi-bin/cmd_3.3 to obtain the appropriate isochrones.

 
