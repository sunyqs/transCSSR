# This code is based on David Darmon's transCSSR code https://github.com/ddarmon/transCSSR written for Python 3.7
# It derives epsilon-machines (eM) from one string of sequence data (univariate) via an implementation of Cosma Shalizi's Causal State Splitting Reconstruction (CSSR) http://bactra.org/CSSR/
# Ihis file should be placed in the main folder of the transCSSR download, 'transCSSR-master' and adds the following:
    # reading from csv: automatically reading temporal sequences coded as columns (with headers) from a csv file
    # alphabet detection: takes the symbols of the involved alphabet directly, without the need to input them manually in the script
    # output file: calculates and prints measures for (1) e-machine
        # The output is:
            # Cmu	H[X_{0}]	hmu  	E
            # the output results file will be stored in a new folder in transCSSR-master, together with a new sub-folder with the .dot and .dat results
##### this script was coded by Timothy Zhang, at UC Davis, C^2-Lab https://github.com/3tz/transCSSR ####
    
# to-do:    
    # define l_max: the maximum word-length to be tested for
    # to read data, CREATE a folder called csv directly in main folder where this code is stored
      # store data in csv format, first row header of sequence, columns with data. All columns must have the same length
      # the code will read in the alphabet automatically. The categorical alphabet must be coded as one-letter per event
      # define filename of the csv's in last code-block below

from transCSSR_bc import *
from os.path import join
import os
import csv
from utils import csv_to_list, get_uniques_from_2d_list


def cssr(string_y, ays, yt_name, l_max):
    """
    Perform
    Arguments:
        string_y: str
            Individual outcomes for each trial concatenated together.
            Ex: '01001100' for 8 trials of coin flipping.
        ays: list
            List of unique outcomes for output
            Ex: ['0', '1'] for coin flipping.
        yt_name: string
            Name that will be used for the output files.
            Ex: coin
        l_max: int, default 4
            The maximum history length to look back when inferring predictive
            distributions.
    """
    string_x = '0' * len(string_y)
    xt_name = ''
    axs = ['0']  # the input alphabets
    # All of the possible pairs of emission symbols for (x, y)
    e_symbols = list(itertools.product(axs, ays))
    # the significance level associated with CSSR's hypothesis tests.
    alpha = 0.001

    inf_alg = 'transCSSR'
    Tx = len(string_x)
    Ty = len(string_y)
    assert Tx == Ty, 'The two time series must have the same length.'
    word_lookup_marg, word_lookup_fut = estimate_predictive_distributions(
        string_x, string_y, l_max
    )
    epsilon, inv_epsilon, morph_by_state = run_transCSSR(
        word_lookup_marg, word_lookup_fut, l_max, axs, ays, e_symbols, xt_name,
        yt_name, alpha=alpha)

    print('The epsilon-transducer has {} states.'.format(len(inv_epsilon)))
    print_morph_by_states(morph_by_state, axs, ays, e_symbols)
    filtered_states, filtered_probs, string_y_pred = filter_and_predict(
        string_x, string_y, epsilon, inv_epsilon, morph_by_state, axs, ays,
        e_symbols, l_max)
    machine_fname = join('transCSSR_results', '+.dot')
    transducer_fname = join('transCSSR_results', '+{}.dot'.format(yt_name))
    pred_probs_by_time, cur_states_by_time = filter_and_pred_probs(
        string_x, string_y, machine_fname, transducer_fname, axs, ays, inf_alg,
        verbose_filtering_errors=True)
    pred_probs_by_time_break, cur_states_by_time_break = \
        filter_and_pred_probs_breakforbidden(
            string_x, string_y, machine_fname, transducer_fname, axs, ays,
            inf_alg
        )

    for t in range(30):
        print((
            t, string_y[t], filtered_probs[t], pred_probs_by_time_break[t, 1],
            pred_probs_by_time[t, 1]
        ))


def main(path_outputyt, prefix_outdir, max_l):
    """ Run CSSR on the given input files and generate output files and
        conditional measures in the given directory.
    Arguments:
        path_outputyt: str
            Path to the outputYt csv file. E.g.: 'csv/outputYt.csv'
        prefix_outdir: str
            Prefix naming of the output directory. Actual output directory will
            be appended with @max_l.
            E.g.: 'output_trans' with @max_l=1 will have outputs under
                './output_trans_L1'
        max_l: int
            The maximum L value for computing the transducer.
    Returns
    """
    # create the final output directory
    dir_out = prefix_outdir + '_L' + str(max_l)
    PATH_DOT_RESULTS = os.path.join(dir_out, 'dot_results')
    os.makedirs(PATH_DOT_RESULTS, exist_ok=True)
    # load in the .csv as column-major 2d-lists
    cols_y = csv_to_list(path_outputyt)
    # Find a set of unique outcomes for each of x & y
    ays = get_uniques_from_2d_list(cols_y)

    # List to save C_X & h_X for each pair so we can write to results.csv later
    results = []
    # Headers of CSV:
    header = [
        'machine_name', 'machine_H[X_{0}]', "machine_E", 'machine_C_mu',
        'machine_h_mu',
    ]
    # Loop through each pair of columns in the CSVs
    for col_y in cols_y:
        yt_name = col_y[0]
        name_machine = '+%s.dot' % yt_name

        stringY = ''.join(col_y[1:])
        cssr(stringY, ays, yt_name, max_l)
        HLs, _, h_mu_mach, _, E, C_mu_mach, _ = compute_ict_measures(
            join('transCSSR_results', name_machine),
            ays, 'transCSSR', L_max=max_l
        )
        # Also move away the cssr. e.g.: +Y1.dot & +Y1.dat_results
        name_machine_dat = name_machine[:-4] + '.dat_results'
        os.rename(
            join('transCSSR_results', name_machine),
            join(
                PATH_DOT_RESULTS, 'mach_' + name_machine.replace('+', '_')
            )
        )
        os.rename(
            join('transCSSR_results', name_machine_dat),
            join(
                PATH_DOT_RESULTS, 'mach_' + name_machine_dat.replace('+', '_')
            )
        )
        results.append([
            'mach_' + name_machine.replace('+', '_'),
            HLs[0], E, C_mu_mach, h_mu_mach,
        ])

    # Now that we have moved all of the .dot & .dat_results, let's create the
    #   results.csv that has the C_mu & h_mu
    with open(join(dir_out, 'results.csv'), 'w') as f:
        writer = csv.writer(f, quoting=2)
        writer.writerow(header)
        for row in results:
            writer.writerow(row)


if __name__ == '__main__':
    main('csv/Yt_test.csv', 'output_cssr_yt_test', 1)
