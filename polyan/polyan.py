"""polyan provides functions for simulating polysome profiles from other
types of data.

It requires pandas and numpy to be installed within the Python
environment. For compare_profiles pyplot is also required.

fp2poly: models polysome peak info from ribosome footrinting data

plot_poly: generates plotable x and y coordinates of traces from fp2poly output

compare_profiles: compares two ribosome footprinting datasets based on
mRNA shifts between peaks in corresponding polysome profiles
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.stats import norm

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import Data  # relative-import the *package* containing the data


def fp2poly(df, has_RNA=True, has_length=False, parset='Scer',
            include_idle=True, poly_limit=30, remove_spurious_RNAs=True):
    """Calculates peak volumes of a polysome profile from a
    footprinting dataset provided as a dataframe.

    The dataframe must be formatted to contain one column with gene
    names, one column with Ribo-Seq data, and one column with matching
    RNA-Seq data (optional).

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        A pandas dataframe that contains the footprinting data

    has_RNA : bool
        True indicates that the third column contains RNA Seq counts

    has_length : bool
        True indicates that the third (if no RNA Seq data) or fourth
        (if with RNA Seq data) column contains gene length data for calculating
        RPK values

    parset : dict
        Either a keyword specifying a pre-defined parameterset
        (currently 'Scer' or 'HEK'), or a dictionary of the form {'Species':
        species name, 'RNA_content':number of RNAs per cell, 'Ribo_content':
        number of ribosomes per cell, 'frac_act': fraction of actively
        translating ribosomes, 'frac_split':     fraction of inactive
        ribosomes that are split into separate subunits, 'gene_reference':
        path to a file listing all genes of the organism with their lengths
        and (optionally, if no column with RNA abundances is contained in df)
        'RNA_reference': path to a file listing RNA abundances}

    Returns
    -------
    numpy.ndarray
        Volumes of individual peaks of the modelled polysome.
    """

    # ####prepare parset

    # read in the parameter set to be used
    with pkg_resources.open_text(Data, 'parameters.json') as read_file:
        parameterset = json.load(read_file)
    # determine what kind of parameter was assigned to parset
    # and process accordingly
    if type(parset) == str:
        if parset in parameterset.keys():
            parameters = parameterset[parset]
    elif type(parset) == dict and all(elem in parset.keys() for elem in [
                         'Species', 'RNA_content', 'Ribo_content', 'frac_act',
                         'frac_split', 'gene_reference']):
        parameters = parset
    else:
        print('\'parset\' needs to be either a keyword specifying a pre- \
              defined parameterset, or a dictionary with the keys Species, \
              RNA_content, Ribo_content,frac_act, frac_split, gene_reference \
              and (optional) RNA_reference.')
        return

    # ####assemble a processible dataset
    # test that the passed dataframe has a correect number of columns
    if df.shape[1] > 2 + has_RNA + has_length:
        print('Warning: too many colums. Columns ' +
              str(2 + has_RNA + has_length + 1) + ' to ' + str(df.shape[1]) +
              ' will be ignored.')
        df = df.iloc[:, :2 + has_RNA + has_length]
    elif df.shape[1] < 2 + has_RNA + has_length:
        raise ValueError('The dataframe does not have sufficient columns for \
                         the specified data')
        return
    colnames = ['ORF', 'Ribo_Prints']
    if has_RNA:
        colnames.append('RNA_Prints')
    if has_length:
        colnames.append('length')
    df.columns = colnames

    # combine data with reference for gene length and transcript if required)
    if (not has_RNA) or (not has_length):
        # load reference dataset specified in parameterset
        with pkg_resources.open_text(Data, parameters['gene_reference']) as \
             read_file:
            genes = pd.read_csv(read_file)
        # check which column of reference file contains relevant descriptors
        if parset == 'Scer':
            # do >95% of gene names start with Y -> systematic names ?
            # the 95% is to catch typos (if using non-systematic gene names)
            # the proportion not starting in Y should be much higher
            if sum(genes.systematic.str[0] == 'Y') > (0.95 * df.shape[0]):
                add_columns = ['systematic']
            else:
                add_columns = ['gene']
        elif parset == 'HEK':
            if df.ORF[0][:4] == 'ENSG':
                add_columns = ['name']
            elif df.ORF[0][:4] == 'ENST':
                add_columns = ['transcript']
            elif df.ORF[0][:3] == 'NM_':
                add_columns = ['refseq']
            else:
                add_columns = ['gene']
        # combine dataframe with reference data (RNA abundance and ORF lengths)
        if not has_RNA:
            add_columns.append('RNA_Prints')
        if not has_length:
            add_columns.append('length')
        dats = df.merge(genes[add_columns],
                        how='inner', left_on='ORF', right_on=add_columns[0])[
                        ['ORF', 'Ribo_Prints', 'RNA_Prints', 'length']]
    else:
        dats = df[['ORF', 'Ribo_Prints', 'RNA_Prints', 'length']]
    # remove rows where either the RNA prints or Ribo prints are 0
    # (these do not contribute information to the profile)
    dats = dats.loc[(dats['RNA_Prints'] > 0) & (dats['Ribo_Prints'] > 0)]

    # #####convert RNA-Seq reads to RPK
    dats['RNA_RPK'] = dats['RNA_Prints']/(dats['length']/1000)

    # ####sort genes into polysome peaks according to RPKM info

    # determine conversion factor from Ribo_Prints to no of Ribosomes
    RiboPrints2Ribos = ((parameters['Ribo_content'] *
                         parameters['frac_act']) / sum(dats['Ribo_Prints']))
    # calculate no of Ribos bound to each mRNA population
    dats['Ribos_bound'] = dats['Ribo_Prints'] * RiboPrints2Ribos
    # determine conversion factor from RNA_RPK to no of RNAs
    RNARPK2RNAs = parameters['RNA_content'] / sum(dats['RNA_RPK'])
    dats['RNAs_per_cell'] = dats['RNA_RPK'] * RNARPK2RNAs
    # calculate the ribosome load per RNA (RperR)
    dats['RperR'] = dats['Ribos_bound'] / dats['RNAs_per_cell']
    # remove any nans
    dats = dats.dropna()
    # remove rows where the number of ribosomes per RNA is > poly_limit
    dats = dats.loc[dats['RperR'] <= poly_limit]
    # remove spurious RNAs (< than 0.05 RNAs per cell)
    if remove_spurious_RNAs:
        dats = dats.loc[dats['RNAs_per_cell'] > 0.05]

    # #####assign RNAs into polysome peaks

    # make an array to hold the relative weights for each polysome class
    poly_array = np.zeros(poly_limit+2)
    # if indicated, assign idle ribosomes to the first three peaks,
    # based on the fractions of active split inactive ribosomes
    if include_idle:
        idle_ribos = (1 - parameters['frac_act']) * parameters['Ribo_content']
        poly_array[0] += idle_ribos * parameters['frac_split'] * 0.34
        poly_array[1] += idle_ribos * parameters['frac_split'] * 0.66
        poly_array[2] += idle_ribos * (1-parameters['frac_split'])
    # go through each row of dats and add ribosomes to the appropriate peak
    for _, row in dats.iterrows():
        this_RperR = row['RperR']
        these_Ribos_bound = row['Ribos_bound']
        # if the number of ribos per RNA is an exact integer,
        # assign the ribos to the corrresponding peak
        floor_val = int(this_RperR)
        if float(floor_val) == this_RperR:
            poly_array[floor_val + 1] += these_Ribos_bound
        # if the number of ribos is between two integers,
        # split the ribos proportionally between the two adjacent peaks
        # for example, for 5.6 Ribos per RNA 60% of ribosomes go to the 6-some
        # peak, 40% to the 5-some peak
        else:
            ceil_weight = (this_RperR-floor_val)
            floor_weight = 1 - ceil_weight
            if floor_val != 0:
                poly_array[floor_val + 1] += these_Ribos_bound * floor_weight
            poly_array[floor_val + 2] += these_Ribos_bound * ceil_weight
    # normalise values to total signal
    poly_array = poly_array/sum(poly_array)
    return poly_array


def plot_poly(peak_vols):

    """
    Returns x and y coordinates of a polysome profile from a list of peak
    volumes computed by fp2poly.

    Parameters
    ----------
    peak_vols : numpy.ndarray
        Output from fp2poly containing peak volumes for a dataset.

    Returns
    -------
    tuple
        the x and y parameters for the modelled polysome profile
    """

    # define the function for returning a normal distribution centred around mu
    # with variance (=peak width) sigma
    def normpdf(x, mu, sigma):
        return 1/(np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(-1 * ((x - mu) **
                                                            2 / (2 * sigma
                                                                 ** 2)))

    # define the relative molecular weights corresponding to the different
    # peaks (40S, 60S, 1-,2-,3-,... some)
    peak_ids = [0.34, 0.66] + list(range(1, len(peak_vols)-1))
    # calculate a series of peak locations based on a typical polysome profile
    peak_locs = [0.15 * np.log(19 * x - 0.6) for x in peak_ids]
    # calculate a series of peak widths based on a typical polysome profile
    peak_widths = [-0.0004 * x + 0.015 for x in peak_ids]
    # definine oD drift and initial peak based on a typical polysome profile
    x = np.linspace(0, 1, num=400)
    drift = -0.07 + 0.19 * x + 0.05
    debris = np.exp((-x + 0.085) * 20) * 1.2
    # construct the plot of the polysome profile
    sum_trace = drift + debris
    for peak_no in range(len(peak_locs)):
        this_peak = normpdf(x, peak_locs[peak_no], peak_widths[peak_no]) * \
                    peak_vols[peak_no]
        sum_trace += this_peak
    return x, sum_trace


def compare_profiles(dats1, dats2, dats1_columns=['ORF','Ribo_Prints',
                                                  'RNA_Prints'],
                      dats2_columns=['ORF', 'RNA_Prints', 'Ribo_Prints'],
                      parset='Scer', colors=['steelblue', 'orange'],
                      conditions=['Cond. 1', 'Cond. 2'],
                      return_df=False):

    """
    Computes and displays the predominant movements of transcripts between
    polysome peaks for two conditions.

    Parameters
    ----------
    dats1 : pandas.core.frame.DataFrame
    dats2 : pandas.core.frame.DataFrame
        pandas dataframes containing the footprinting data

    dats1_columns : list
    dats2_columns : list
        lists containing the column names with ORF, Ribo-Seq and RNA-Seq data

    parset: dict
        the context data for the organism (Scer or HEK)

    Returns
    -------
    (matplotlib.figure.Figure,matplotlib.axes._subplots.AxesSubplot)
        A tuple with the Figure and Axes objects of the plot.
    """

    def counts_to_RPKM(dats):
        # read in the gene length info dataset
        with pkg_resources.open_text(Data, parameters['gene_reference']) as \
             read_file:
            genes = pd.read_csv(read_file)
        # combine input dataset with gene length information
        dats = dats.merge(genes[['systematic', 'length']], how='inner',
                          left_on='ORF', right_on='systematic')[
                          ['ORF', 'RNA_Prints', 'Ribo_Prints', 'length']]
        # calculate RPK
        dats['RNA_RPKM'] = dats['RNA_Prints']/(dats['length']/1000)
        dats = dats.drop(['RNA_Prints', 'length'], axis=1)
        return dats

    def calc_RperR(Rdats, poly_limit=30):
        # determine conversion factor from Ribo_Prints to no of Ribosomes
        RiboPrints2Ribos = (parameters['Ribo_content'] *
                            parameters['frac_act']) / sum(Rdats['Ribo_Prints'])
        Rdats['Ribos_bound'] = Rdats['Ribo_Prints']*RiboPrints2Ribos
        # determine conversion factor from RNA_RPKM to no of RNAs
        RNARPKM2RNAs = parameters['RNA_content'] / sum(Rdats['RNA_RPKM'])
        Rdats['RNAs_per_cell'] = Rdats['RNA_RPKM']*RNARPKM2RNAs
        # calculate the ribosome load per RNA (RperR)
        Rdats['RperR'] = np.round(Rdats['Ribos_bound'] /
                                  Rdats['RNAs_per_cell'])
        Rdats = Rdats.dropna()
        # remove rows where the number of ribosomes per RNA is > poly_limit
        Rdats = Rdats.loc[Rdats['RperR'] <= poly_limit]
        return Rdats

    # read in the parameter set to be used
    with pkg_resources.open_text(Data, 'parameters.json') as read_file:
        parameterset = json.load(read_file)
    if type(parset) == str:
        if parset in parameterset.keys():
            parameters = parameterset[parset]
            print('Using parameterset for ' + parameters['Species'])
    elif type(parset) == dict and all(elem in parset.keys() for elem in [
                         'Species', 'RNA_content', 'Ribo_content', 'frac_act',
                         'frac_split', 'gene_reference']):
        parameters = parset
        print('Using parameterset for ' + parameters['Species'])
    else:
        print('\'parset\' needs to be either a keyword specifying a '
              'pre-defined parameterset, or a dictionary with the keys Species\
              , RNA_content,Ribo_content,frac_act, frac_split, gene_reference \
              and (optional) RNA_reference.')
        return

    # prepare datasets and compute RperR values via RNA RPKM values
    dats1 = dats1[dats1_columns]
    dats1.columns = ['ORF', 'RNA_Prints', 'Ribo_Prints']
    dats1 = counts_to_RPKM(dats1)
    dats1 = calc_RperR(dats1)
    dats2 = dats2[dats2_columns]
    dats2.columns = ['ORF', 'RNA_Prints', 'Ribo_Prints']
    dats2 = counts_to_RPKM(dats2)
    dats2 = calc_RperR(dats2)

    # compute column movements for individual transcripts
    fromto = dats1[['ORF', 'RperR']].merge(dats2[['ORF', 'RperR']],
                                           how='inner', on='ORF')
    fromto.columns = ['ORF', 'from', 'to']

    # calculate main destinations for each origin peak
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
                this_from.loc[this_from['to'] ==
                              this_to_unique[to_idx]].shape[0]))
    alpha_vec_scaled = [(alpha)/(max(alpha_vec)*0.7) for alpha in alpha_vec]

    df = pd.DataFrame({'From': from_vec, 'To': to_vec, 'Direction': dir_vec,
                       'Alpha': alpha_vec_scaled})
    up_df = df.loc[df['Direction'] == 'up']
    down_df = df.loc[df['Direction'] == 'down']

    # prepare main figure
    fig, axs = plt.subplots(2, 1, constrained_layout=True, figsize=(4, 3))

    for idx in range(up_df.shape[0]):
        axs[0].plot([up_df.iloc[idx]['From'], up_df.iloc[idx]['To']], [3, 1],
                    linewidth=4, color=colors[0],
                    alpha=up_df.iloc[idx]['Alpha'])
    axs[0].set_ylim(0.9, 3.1)
    axs[0].set_yticks([3, 1])
    axs[0].set_yticklabels([conditions[0], conditions[1]])
    axs[0].set_xlim((0, 30))

    for idx in range(down_df.shape[0]):
        axs[1].plot([down_df.iloc[idx]['From'], down_df.iloc[idx]['To']],
                    [3, 1], linewidth=4, color=colors[1],
                    alpha=down_df.iloc[idx]['Alpha'])
    axs[1].set_ylim(0.9, 3.1)
    axs[1].set_yticks([3, 1])
    axs[1].set_yticklabels([conditions[0], conditions[1]])
    axs[1].set_xlabel('Polysome number')
    axs[1].set_xlim((0, 30))

    if return_df:
        return fig, fig.axes, df
    else:
        return fig, fig.axes


# ===========================================================================

def rmsd_profile(peakvols1, peakvols2, parset='Scer'):

    '''Calculates the Root Mean Square Deviation between
    two polysome profiles simulated using polyan.fp2poly().

    Parameters
    ----------
    peakvols1 : numpy.ndarray
    peakvols2 : numpy.ndarray
        arrays containing peak volumes for two datasets calculated by fp2 poly

    parset : str
        The name of the reference peak volumes to be used if peakvols2 is "ref"
        This parameter is only evaluated if peakvols2 is ref.

    Returns
    -------
    float
        the calculated RMSD between the two peak volume arrays.
    '''

    if type(peakvols2) == str and peakvols2 == 'ref':
        if parset == 'Scer':
            ref_peaks = [0.003, 0.0051, 0.026, 0.018, 0.0205, 0.0155, 0.012,
                         0.0085, 0.008, 0.006, 0.005, 0.0035, 0.002, 0.0015,
                         0.001, 0.0005, 0.0002, 5e-05]
            peakvols2 = ref_peaks / np.sum(ref_peaks)
        elif parset == 'HEK':
            ref_peaks = [0.0048, 0.008, 0.0033, 0.0019, 0.0032, 0.004, 0.004,
                         0.0032, 0.003, 0.003, 0.00255, 0.00265, 0.002,
                         0.00195, 0.0018, 0.0017, 0.0013, 0.00097, 0.00079,
                         0.00043, 0.00035, 0.00036, 0.00036]
            peakvols2 = ref_peaks / np.sum(ref_peaks)

    # ensure that both peakvols are of equal length - pad with zeros
    if len(peakvols1) < len(peakvols2):
        temp = np.zeros_like(peakvols2)
        temp[:len(peakvols1)] = peakvols1
        peakvols1 = temp.copy()
    elif len(peakvols1) > len(peakvols2):
        temp = np.zeros_like(peakvols1)
        temp[:len(peakvols2)] = peakvols2
        peakvols2 = temp.copy()

    rmsd = np.sqrt(sum(np.square(np.subtract(peakvols1, peakvols2))) /
                   len(peakvols1))

    return rmsd

# =============================================================================


def nrmsd_profile(peakvols1, peakvols2, parset='Scer'):

    '''Calculates the normalised Root Mean Square Deviation
    between two polysome profiles simulated using
    polyan.fp2poly(). Normalisation is to the average
    RMSD for all pairwise comparison of datasets used in
    the original publication.

    Parameters
    ----------
    peakvols1 : numpy.ndarray
    peakvols2 : numpy.ndarray
        arrays containing peak volumes for two datasets calculated by fp2 poly

    parset : str
        The name of the reference peak set to be used for the normalisation
        Possibe values are 'Scer' or 'HEK'

    Returns
    -------
    float
        the calculated normalised RMSD between the two peak volume arrays.
    '''

    # define the reference for normalisation, ie the average
    # rmsd for all pairwise comparions for the datasets used
    # in the original paper
    if parset == 'Scer':
        normval = 0.018049851598865935
    elif parset == 'HEK':
        normval = 0.01264714889273241
    else: 
        return -1

    return rmsd_profile(peakvols1, peakvols2) / normval

# =============================================================================


def prmsd_profile(peakvols1, peakvols2, parset='Scer'):

    '''Calculates the p-value for observing a particular
    Root Mean Square Deviation between two polysome 
    profiles simulated using polyan.fp2poly(). The 
    p-value is calculated based on the distribution of 
    individual RMSD values for for all pairwise comparisons
    of "known good" datasets used in the original publication.

    Parameters
    ----------
    peakvols1 : numpy.ndarray
    peakvols2 : numpy.ndarray
        arrays containing peak volumes for two datasets 
        calculated by fp2poly

    parset : str
        The name of the reference peak set to be used if 
        prmsd_profile is called with the "peakvols2 = ref" 
        option. Possibe values are 'Scer' or 'HEK'.

    Returns
    -------
    float
        The probability of observing an RMSD with similar magnitude as that 
        for peakvols1 and peakvols2 in the reference RMSD values.
    '''

    # return the probability of observing the RMSD value for the
    #specified comparison, based on the parameters of a normal
    #distribution fitted to the reference RMSD values.

    return 1-norm(0.014, 0.005).cdf(rmsd_profile(peakvols1, peakvols2))