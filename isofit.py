import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

# Functions

# splits isochrone data files by its metallicities
def isochrone_grid(isochrone_data : pd.DataFrame,groupby_name):
    isochrones_grouped = isochrone_data.groupby([groupby_name])
    # List of dataframes
    isochrone_list = [isochrone for name, isochrone in isochrones_grouped]
    return isochrone_list

# tabulates output isochrone data into useful columns based on selection

def tabulate_isochrone_data(data,columns_select,names):
    data = data.T
    data_to_tabulate = data[[columns_select]]
    data_dict = {}
    for column, name in zip(data_to_tabulate,names):
        data_dict[name] = column
    return pd.DataFrame(data_dict)

def tabulate_my_data(df):
    abs_mag = df['phot_g_mean_mag'] + 5 - (5*np.log10(1000./df['parallax']))
    select = (((df.bp_rp > 0.3 ) & (df.bp_rp < 1.8) & (abs_mag > df.bp_rp*2.46+2.76)) |
                        ((df.bp_rp > 1.8 ) &                    (abs_mag > df.bp_rp*2.8 +2.16)) |
                        ((df.bp_rp > -0.5) & (df.bp_rp < 4  ) & (abs_mag < df.bp_rp*2.14-0.57)) |
                        ((df.bp_rp > 1.2 ) & (df.bp_rp < 3  ) & (abs_mag < df.bp_rp*1.11+0.66)))
    df_bad = df[select]

    df = df[np.invert(select)]

    g = df.phot_g_mean_mag.values
    g_flux_err = df['phot_g_mean_flux'].values
    parallax = df.parallax.values
    parallax_error = df['parallax_error'].values
    Mg = g + 5 - (5*np.log10(1000/parallax))
    b = df['phot_bp_mean_mag'].values
    r = df['phot_rp_mean_mag'].values
    b_r = b - r
    label = df['label']

    data_dict={'g':g,'b':b,'r':r,'Mg':Mg,'b_r':b_r, 'parallax_error':parallax_error,'g_flux_err':g_flux_err, 'label' : label}
    data_df = pd.DataFrame(data_dict)
    return data_df, df_bad

# plots grid of isochrones on top of my data
def plot_analysis(df, isochrones, df_bad=None):

    g = df.g.values
    b_r = df.b_r.values
    Mg = df.Mg.values

    # Plotting my data
    plt.figure()
    xy = np.vstack([b_r,Mg])
    numb_density = gaussian_kde(xy)(xy)
    plt.scatter(b_r,Mg,c=numb_density,s=3)
    plt.gca().invert_yaxis()

    parallax_bad = df_bad['parallax'].values
    b_bad = df_bad['phot_bp_mean_mag'].values
    r_bad = df_bad['phot_rp_mean_mag'].values
    g_bad = df_bad['phot_g_mean_mag'].values
    Mg_bad = g_bad + 5 - (5*np.log10(1000/parallax_bad))
    plt.scatter(b_bad - r_bad, Mg_bad, s=1, alpha=.4, c = 'gray')

    # Plotting isochrones
    for isochrone in isochrones:
        # selection = .loc[(isochrone['G_BP'] - isochrone['G_RP'] > np.min(b_r)) & (isochrone['G_BP'] - isochrone['G_RP'] < np.max(b_r))].values
        g = isochrone['Gmag']
        b = isochrone['G_BP']
        r = isochrone['G_RP']
        plt.plot(b-r,g,linewidth=.2)
    plt.title(r'$\lambda_{ori}$' + ' with Isochrones')
    plt.xlabel(r'$G_b - G_r$')
    plt.ylabel(r'$M_g$')
    plt.show()

def calc_square_error(data,isochrones):
    '''
    Interpolates isochrones to grid of data and computes sum of square in the residuals
    '''
    import scipy.interpolate as interp

    g = data.g.values
    b_r = data.b_r.values
    Mg = data.Mg.values
    isochrones_return = []
    for isochrone in isochrones:
        g = isochrone['Gmag'].values
        y = g
        b = isochrone['G_BP'].values
        r = isochrone['G_RP'].values
        x = (b - r)
        x = x[np.where((x>=min(b_r)) & (x<=max(b_r)))]
        y = y[np.where((x>=min(b_r)) & (x<max(b_r)))]
        age = isochrone['logAge'].values
        age = age[np.where((x>=min(b_r)) & (x<max(b_r)))]
        x = x[np.where(y>=min(Mg))]
        y = y[np.where(y>=min(Mg))]
        age = age[np.where(y>=min(Mg))]
        f_y = interp.interp1d(x,y,fill_value='extrapolate')
        new_y = f_y(b_r)

        error = np.sum((Mg - new_y)**2)
        isochrone['square_residuals'] = error
        isochrone = pd.DataFrame({'Gmag':y,'B_R':x,'logAge':age, 'square_residuals':error})
        isochrones_return.append(isochrone)
    return isochrones_return

def best_fit(data, isochrones):
    isochrones = calc_square_error(data,isochrones)
    ages = []
    errors = []
    for iso in isochrones:
        ages.append(iso['logAge'].iloc[0])
        errors.append(iso['square_residuals'].iloc[0])
    return ages[np.argmin(errors)]

def plot_best_fit(data,isochrones):
    isochrones = calc_square_error(data,isochrones)
    ages = []
    errors = []
    for iso in isochrones:
        ages.append(iso['logAge'].iloc[0])
        errors.append(iso['square_residuals'].iloc[0])
    isochrones = pd.concat(isochrones)
    age = isochrones.loc[isochrones['logAge'] == ages[np.argmin(errors)]]['logAge'].iloc[0]
    isochrone_best_fit = isochrones.loc[isochrones['logAge'] == ages[np.argmin(errors)]]
    g = data.g.values
    b_r = data.b_r.values
    Mg = data.Mg.values

    g_model = isochrone_best_fit['Gmag'].values
    b_r_model = isochrone_best_fit['B_R'].values

    xy = np.vstack([b_r,Mg])
    numb_density = gaussian_kde(xy)(xy)

    plt.figure()
    plt.scatter(b_r,Mg,c=numb_density,s=1,alpha=.4)
    plt.gca().invert_yaxis()

    plt.plot(b_r_model,g_model,linewidth=.5,c='black',label = r'$\tau = $' + str(round((10**age)/1e6,3)) + 'Myr')
    plt.legend()
    plt.show()
#     plt.savefig('/Users/admin/Desktop/orion/orion_populations/ages_plots/' + str(data.label.iloc[0]),overwrite = True,dpi = 300)

def one_cluster(cluster,isochrone_data):
    isochrones = isochrone_grid(isochrone_data,'logAge')
    my_data_tabulated, df_bad = tabulate_my_data(cluster)
#     age_best_fit = best_fit(my_data_tabulated,isochrones)
#     print(age_best_fit)
    plot_best_fit(my_data_tabulated,isochrones)
    
#def all_clusters(df):
#    clusters = df.groupby(['label'])
#    for name,cluster in clusters:
#        isochrones = isochrone_grid(isochrone_data,'logAge')
#        my_data_tabulated, df_bad = tabulate_my_data(cluster)
#    #     age_best_fit = best_fit(my_data_tabulated,isochrones)
#    #     print(age_best_fit)
#        plot_best_fit(my_data_tabulated,isochrones)




