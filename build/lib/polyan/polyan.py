"""polyan provides functions for simulating polysome profiles from other
types of data.

It requires pandas and numpy to be installed within the Python 
environment.

fp2poly: models polysome peak info from ribosome footrinting data

plot_poly: generates plotable x and y coordinates of traces from fp2poly output

compare_profiles: compares two ribosome footprinting datasets based on
mRNA shifts between peaks in corresponding polysome profiles
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import Data  # relative-import the *package* containing the data

def fp2poly(df, df_columns=['ORF', 'Ribo_Prints', 'RNA_Prints'], parset = 'Scer', include_idle = True, poly_limit=30, remove_spurious_RNAs=True):
    """Calculates peak volumes of a polysome profile from a 
    footprinting dataset provided as a dataframe.
    
    The dataframe must be formatted to contain one column with gene 
    names, one column with Ribo-Seq data, and one column with matching
    RNA-Seq data (optional).
    
    Parameters
    ----------
    df: the pandas dataframe that contains the footprinting data
    
    df_columns: list including the column names for columns that contain 
    the     gene names, Ribo-Seq data and RNA-Seq data in that order. 
    Default is ['ORF','Ribo_Prints','RNA_Prints'].
    
    parset: either a keyword specifying a pre-defined parameterset 
    (currently 'Scer' or 'HEK', or a dictionary specifying 'Species':
    species name, 'RNA_content':number of RNAs per cell, 'Ribo_content':
    number of ribosomes per cell, 'frac_act': fraction of actively 
    translating ribosomes, 'frac_split':     fraction of inactive 
    ribosomes that are split into separate subunits, 'gene_reference':
    path to a file listing all genes of the organism with their lengths
    and (optionally, if no column with RNA abundances is contained in df)
    'RNA_reference': path to a file listing RNA abundances .
    
    include_idle: bool, if True (default) idle ribosomes are included.
    
    poly_limit:int. ignore any mRNAs containing more ribosomes than 
    specified here. Deault is 30.
    
    remove_spurious_RNAs: bool, if True (default) ignores RNAs 
    expressed at less than 0.05 RNAs per cell.
    
    Returns
    -------
    numpy array, containing the volume of individual peaks of the 
    modelled polysome.
    """

    #read in the parameter set to be used
    with pkg_resources.open_text(Data, 'parameters.json') as read_file:
        parameterset = json.load(read_file)
    if type(parset)== str:
        if parset in parameterset.keys():
            parameters = parameterset[parset]
            print('Using parameterset for ' + parameters['Species'])
    elif type(parset)==dict and all(elem in parset.keys() for elem in ['Species','RNA_content','Ribo_content','frac_act','frac_split','gene_reference']):
        parameters=parset
        print('Using parameterset for ' + parameters['Species'])
    else:
        print('\'parset\' needs to be either a keyword specifying a pre-defined parameterset, or a dictionary with the keys Species, RNA_content,Ribo_content,frac_act, frac_split, gene_reference and (optional) RNA_reference.')
        return
    
    #####assemble a processible dataset
    
    #if non-default columns are specified, check whether columns exist:
    if df_columns != ['ORF', 'Ribo_Prints', 'RNA_Prints']:
        if set(df_columns).issubset(df.columns):
            pass
        else:
            print('Specified columns do not exist in data frame. Available columns include: \n' +
                  df.columns)
        return
    elif 'RNA_Prints' in df.columns:
        dats = df[df_columns]
    else:
        dats = df[df_columns[:2]]
        

    #check whether reference RNA-Seq data need to be used
    if len(dats.columns) == 2:
        if 'RNA_reference' in parameters.keys():
            with pkg_resources.open_text(Data, parameters['RNA_reference']) as read_file:
                RNA_ref = pd.read_csv(read_file)
        else:
            print('No column for RNA data or RNA reference set have been specified.')
            return
        dats.columns = ['ORF', 'Ribo_Prints']
        dats = dats.merge(RNA_ref, how='inner', on='ORF')
        print('No mRNA data specified, using reference RNA-Seq data.')
    elif len(df_columns) == 3:
        dats.columns = ['ORF', 'Ribo_Prints', 'RNA_Prints']
    else:
        print('Incorrect column number. \nfp2poly requires two or three columns in this order: \n' +
              '1) gene names, 2) Ribo-Seq counts, 3) (optional) RNA-Seq counts')
        return

    #remove rows where either the RNA prints or Ribo prints are 0
    # (these do not contribute information to the profile)
    dats = dats.loc[(dats['RNA_Prints'] > 0) & (dats['Ribo_Prints'] > 0)]

    ######convert Ribo-Seq reads to RPKM


    #combine input dataset with gene length information
    with pkg_resources.open_text(Data, parameters['gene_reference']) as read_file:
                genes = pd.read_csv(read_file)
    dats = dats.merge(genes[['name', 'length']],
                      how='inner', left_on='ORF', right_on='name')[
                          ['ORF', 'RNA_Prints', 'Ribo_Prints', 'length']]
    #calculate RPKM
    dats['RNA_RPKM'] = dats['RNA_Prints']/(dats['length']/1000)

    #####sort genes into polysome peaks according to RPKM info

    #determine conversion factor from Ribo_Prints to no of Ribosomes
    RiboPrints2Ribos = (parameters['Ribo_content'] * parameters['frac_act']) / sum(dats['Ribo_Prints'])
    dats['Ribos_bound'] = dats['Ribo_Prints'] * RiboPrints2Ribos
    #determine conversion factor from RNA_RPKM to no of RNAs
    RNARPKM2RNAs = parameters['RNA_content'] / sum(dats['RNA_RPKM'])
    dats['RNAs_per_cell'] = dats['RNA_RPKM'] * RNARPKM2RNAs
    #calculate the ribosome load per RNA (RperR)
    dats['RperR'] = dats['Ribos_bound'] / dats['RNAs_per_cell']
    dats = dats.dropna()
    #remove rows where the number of ribosomes per RNA is > poly_limit
    dats = dats.loc[dats['RperR'] <= poly_limit]
    #remove spurious RNAs (< than 0.05 RNAs per cell)
    if remove_spurious_RNAs:
        dats = dats.loc[dats['RNAs_per_cell'] > 0.05]

    ######assign RNAs into polysome peaks

    #make an array to hold the relative weights for each polysome class
    poly_array = np.zeros(poly_limit+2)
    #if indicated, assign idle ribosomes to the first three peaks,
    # based on the fraction of active ribosomes and the fraction of split inactive ribosomes
    if include_idle:
        idle_ribos = (1 - parameters['frac_act']) * parameters['Ribo_content']
        poly_array[0] += idle_ribos * parameters['frac_split'] * 0.34
        poly_array[1] += idle_ribos * parameters['frac_split'] * 0.66
        poly_array[2] += idle_ribos * (1-parameters['frac_split'])
    #go through each row of dats and add ribosomes to the appropriate peak
    for _, row in dats.iterrows():
        this_RperR = row['RperR']
        these_Ribos_bound = row['Ribos_bound']
        #if the number of ribos per RNA is an exact integer,
        # assign the ribos to the corrresponding peak
        floor_val = int(this_RperR)
        if  float(floor_val) == this_RperR:
            poly_array[floor_val + 1] += these_Ribos_bound
        #if the number of ribos is between two integers,
        # split the ribos proportionally between the two adjacent peaks
        #for example, for 5.6 Ribos per RNA 60% of ribosomes go to the 6-some peak,
        #40% to the 5-some peak
        else:
            ceil_weight = (this_RperR-floor_val)
            floor_weight = 1 - ceil_weight
            if floor_val != 0:
                poly_array[floor_val + 1] += these_Ribos_bound * floor_weight
            poly_array[floor_val + 2] += these_Ribos_bound * ceil_weight
    #normalise values to total signal
    poly_array = poly_array/sum(poly_array)
    return poly_array

def plot_poly(peak_vols):
    """
    Returns x and y coordinates of a polysome profile from a list of peak volumes
    computed by fp2poly
    """

    #define the function for returning a normal distribution centred around mu with variance
    # (=peak width) sigma
    def normpdf(x, mu, sigma):
        return 1/(np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(-1 * ((x - mu) ** 2 / (2 * sigma ** 2)))

    #define the relative molecular weights corresponding to the different
    # peaks (40S, 60S, 1-,2-,3-,... some)
    peak_ids = [0.34, 0.66] + list(range(1, len(peak_vols)-1))
    #calculate a series of peak locations based on a typical polysome profile
    peak_locs = [0.15*np.log(19 * x -0.6) for x in peak_ids]
    #calculate a series of peak widths based on a typical polysome profile
    peak_widths = [-0.0004 * x + 0.015 for x in peak_ids]
    #definine the oD drift and the initial peak based on a typical polysome profile
    x = np.linspace(0, 1, num=400)
    drift = -0.07 + 0.19 * x + 0.05
    debris = np.exp((-x + 0.085) * 20) * 1.2
    #construct the plot of the polysome profile
    sum_trace = drift + debris
    for peak_no in range(len(peak_locs)):
        this_peak = normpdf(x, peak_locs[peak_no], peak_widths[peak_no]) * peak_vols[peak_no]
        sum_trace += this_peak
    return x, sum_trace




def compare_profiles(dats1, dats2, dats1_columns=['ORF', 'RNA_Prints', 'Ribo_Prints'],
                     dats2_columns=['ORF', 'RNA_Prints', 'Ribo_Prints'], parset = 'Scer',
                     colors=['steelblue', 'orange'], conditions=['Cond. 1', 'Cond. 2'],
                     return_df=False):

    """
    Computes and displays the predominant movements of transcripts between polysome peaks
    for two conditions.
    """


    def counts_to_RPKM(dats):
        #read in the gene length info dataset
        with pkg_resources.open_text(Data, parameters['gene_reference']) as read_file:
                genes = pd.read_csv(read_file)
        #combine input dataset with gene length information
        dats = dats.merge(genes[['name', 'length']], how='inner', left_on='ORF', right_on='name')[
            ['ORF', 'RNA_Prints', 'Ribo_Prints', 'length']]
        #calculate RPKM
        dats['RNA_RPKM'] = dats['RNA_Prints']/(dats['length']/1000)
        dats = dats.drop(['RNA_Prints', 'length'], axis=1)
        return dats

    def calc_RperR(Rdats,poly_limit=30):
        #determine conversion factor from Ribo_Prints to no of Ribosomes
        RiboPrints2Ribos = (parameters['Ribo_content'] * parameters['frac_act']) / sum(Rdats['Ribo_Prints'])
        Rdats['Ribos_bound'] = Rdats['Ribo_Prints']*RiboPrints2Ribos
        #determine conversion factor from RNA_RPKM to no of RNAs
        RNARPKM2RNAs = parameters['RNA_content'] / sum(Rdats['RNA_RPKM'])
        Rdats['RNAs_per_cell'] = Rdats['RNA_RPKM']*RNARPKM2RNAs
        #calculate the ribosome load per RNA (RperR)
        Rdats['RperR'] = np.round(Rdats['Ribos_bound'] / Rdats['RNAs_per_cell'])
        Rdats = Rdats.dropna()
        #remove rows where the number of ribosomes per RNA is > poly_limit
        Rdats = Rdats.loc[Rdats['RperR'] <= poly_limit]
        return Rdats
    
    #read in the parameter set to be used
    with pkg_resources.open_text(Data, 'parameters.json') as read_file:
        parameterset = json.load(read_file)
    if type(parset)== str:
        if parset in parameterset.keys():
            parameters = parameterset[parset]
            print('Using parameterset for ' + parameters['Species'])
    elif type(parset)==dict and all(elem in parset.keys() for elem in ['Species','RNA_content','Ribo_content','frac_act','frac_split','gene_reference']):
        parameters=parset
        print('Using parameterset for ' + parameters['Species'])
    else:
        print('\'parset\' needs to be either a keyword specifying a pre-defined parameterset, or a dictionary with the keys Species, RNA_content,Ribo_content,frac_act, frac_split, gene_reference and (optional) RNA_reference.')
        return
    
    #prepare datasets and compute RperR values via RNA RPKM values
    dats1 = dats1[dats1_columns]
    dats1.columns = ['ORF', 'RNA_Prints', 'Ribo_Prints']
    dats1 = counts_to_RPKM(dats1)
    dats1 = calc_RperR(dats1)
    dats2 = dats2[dats2_columns]
    dats2.columns = ['ORF', 'RNA_Prints', 'Ribo_Prints']
    dats2 = counts_to_RPKM(dats2)
    dats2 = calc_RperR(dats2)

    #compute column movements for individual transcripts
    fromto = dats1[['ORF', 'RperR']].merge(dats2[['ORF', 'RperR']], how='inner', on='ORF')
    fromto.columns = ['ORF', 'from', 'to']

    #calculate main destinations for each origin peak
    from_unique = fromto['from'].unique()
    from_unique.sort()
    from_vec, to_vec, dir_vec, alpha_vec = [], [], [], []
    for from_idx in range(len(from_unique)):
        this_from = fromto.loc[fromto['from'] == from_unique[from_idx]]
        this_to_unique = this_from['to'].unique()
        for to_idx in range(len(this_to_unique)):
            from_vec.append(from_unique[from_idx])
            to_vec.append(this_to_unique[to_idx])
            if from_unique[from_idx] < this_to_unique[to_idx]:
                dir_vec.append('up')
            elif from_unique[from_idx] > this_to_unique[to_idx]:
                dir_vec.append('down')
            else:
                dir_vec.append('nc')
            alpha_vec.append(np.log(
                this_from.loc[this_from['to'] == this_to_unique[to_idx]].shape[0]))
    alpha_vec_scaled = [(alpha)/(max(alpha_vec)*0.7) for alpha in alpha_vec]

    df = pd.DataFrame({'From':from_vec, 'To':to_vec, 'Direction':dir_vec, 'Alpha':alpha_vec_scaled})
    up_df = df.loc[df['Direction'] == 'up']
    down_df = df.loc[df['Direction'] == 'down']

    #prepare main figure
    fig,axs = plt.subplots(2,1,constrained_layout=True,figsize = (4,3))

    for idx in range(up_df.shape[0]):
        axs[0].plot([up_df.iloc[idx]['From'], up_df.iloc[idx]['To']], [3, 1], linewidth=4,
                 color=colors[0], alpha=up_df.iloc[idx]['Alpha'])
    axs[0].set_ylim(0.9, 3.1)
    axs[0].set_yticks([3, 1])
    axs[0].set_yticklabels([conditions[0], conditions[1]])
    axs[0].set_xlim((0, 30))


    for idx in range(down_df.shape[0]):
        axs[1].plot([down_df.iloc[idx]['From'], down_df.iloc[idx]['To']], [3, 1],
                 linewidth=4, color=colors[1], alpha=down_df.iloc[idx]['Alpha'])
    axs[1].set_ylim(0.9, 3.1)
    axs[1].set_yticks([3, 1])
    axs[1].set_yticklabels([conditions[0], conditions[1]])
    axs[1].set_xlabel('Polysome number')
    axs[1].set_xlim((0, 30))

    if return_df:
        return fig, fig.axes, df
    else:
        return fig, fig.axes


