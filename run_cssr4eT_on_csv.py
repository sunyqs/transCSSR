import csv
import itertools
import transCSSR as tC
import transCSSR_bc as tC_bc
from os.path import join
import os
from run_cssr_on_csv import cssr
from utils import csv_to_list, get_uniques_from_2d_list
import time


def trans(string_x, string_y, axs, ays, xt_name, yt_name, l_max):
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # Set the parameters and associated quantities:
    # 	axs, ays -- the input / output alphabets
    # 	alpha    -- the significance level associated with
    # 	            CSSR's hypothesis tests.
    # 	L        -- The maximum history length to look
    #               back when inferring predictive
    #               distributions.
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e_symbols = list(
        itertools.product(axs, ays))  # All of the possible pairs of emission
    alpha = 0.001

    # The two time series must have the same length.
    assert (len(string_x) == len(string_y))

    word_lookup_marg, word_lookup_fut = tC.estimate_predictive_distributions(
        string_x, string_y, l_max)

    epsilon, invepsilon, morph_by_state = tC.run_transCSSR(
        word_lookup_marg, word_lookup_fut, l_max, axs, ays, e_symbols, xt_name,
        yt_name, alpha=alpha)

    print('The epsilon-transducer has {} states.'.format(len(invepsilon)))
    tC.print_morph_by_states(morph_by_state)


def trans_bc(string_x, string_y, axs, ays, xt_name, yt_name, l_max):
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # Set the parameters and associated quantities:
    # 	axs, ays -- the input / output alphabets
    # 	alpha    -- the significance level associated with
    # 	            CSSR's hypothesis tests.
    # 	L        -- The maximum history length to look
    #               back when inferring predictive
    #               distributions.
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # All of the possible pairs of emission
    e_symbols = list(
        itertools.product(axs, ays))
    alpha = 0.001

    counting_method = 1
    # The two time series must have the same length.
    assert (len(string_x) == len(string_y))

    startTime = time.time()
    word_lookup_marg, word_lookup_fut = tC_bc.estimate_predictive_distributions(
        string_x, string_y, l_max,
        counting_method=counting_method, axs=axs, ays=ays
    )
    print('The transCSSR counting took {0} seconds...'.format(
        time.time() - startTime))

    epsilon, invepsilon, morph_by_state = tC_bc.run_transCSSR(
        word_lookup_marg, word_lookup_fut, l_max, axs, ays, e_symbols, xt_name,
        yt_name, alpha=alpha, verbose=False)

    print('The epsilon-transducer has {} states.'.format(len(invepsilon)))
    tC_bc.print_morph_by_states(morph_by_state, axs, ays, e_symbols)


def main(path_inputxt, path_outputyt, prefix_outdir, max_l):
    """ Run CSSR, transCSSR, and transCSSR_bc on the given input files and
        generate output files and conditional measures in the given directory.
    Arguments:
        path_inputxt: str
            Path to the inputXt csv file. E.g.: 'csv/inputXt.csv'
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
    cols_x = csv_to_list(path_inputxt)
    cols_y = csv_to_list(path_outputyt)
    # Find a set of unique outcomes for each of x & y
    axs = get_uniques_from_2d_list(cols_x)
    ays = get_uniques_from_2d_list(cols_y)

    # List to save C_X & h_X for each pair so we can write to results.csv later
    results = []
    # Headers of CSV:
    header = [
        'machine_name', 'machine_H[X_{0}]', "machine_E", 'machine_C_mu',
        'machine_h_mu',
        'transducer_name', 'transducer_C_mu', 'transducer_h_mu',
        'transducer_bc_name', 'transducer_bc_C_mu', 'transducer_bc_h_mu'

    ]

    # Loop through each pair of columns in the CSVs
    for col_x, col_y in zip(cols_x, cols_y):
        xt_name, yt_name = col_x[0], col_y[0]
        name_machine = '+%s.dot' % xt_name
        name_transducer = '%s+%s.dot' % (xt_name, yt_name)

        stringX, stringY = ''.join(col_x[1:]), ''.join(col_y[1:])
        cssr(stringX, axs, xt_name, max_l)
        HLs, _, h_mu_mach, _, E, C_mu_mach, _ = tC_bc.compute_ict_measures(
            join('transCSSR_results', name_machine),
            axs, 'transCSSR', L_max=max_l
        )

        # generate output files for transducer without BC
        trans(stringX, stringY, axs, ays, xt_name, yt_name, max_l)
        C_mu, h_mu = tC.compute_conditional_measures(
            join('transCSSR_results', name_machine),
            join('transCSSR_results', name_transducer),
            axs, ays, inf_alg='transCSSR')

        # Move transducer output files out of the way so we can generate
        #   trans_bc.
        # Transducer output files: 'X1+Y1.dot' & 'X1+Y1.dat_results'
        # name_transducer has 'X1+Y1.dot' so just need the latter
        name_transducer_dat = name_transducer[:-4] + '.dat_results'
        os.rename(
            join('transCSSR_results', name_transducer),
            join(
                PATH_DOT_RESULTS,
                'trans_' + name_transducer.replace('+', '_')
            )
        )
        os.rename(
            join('transCSSR_results', name_transducer_dat),
            join(
                PATH_DOT_RESULTS,
                'trans_' + name_transducer_dat.replace('+', '_')
            )
        )

        # Now that trans outputs are moved,let's generate trans_bc
        trans_bc(stringX, stringY, axs, ays, xt_name, yt_name, max_l)
        C_mu_bc, h_mu_bc = tC.compute_conditional_measures(
            join('transCSSR_results', name_machine),
            join('transCSSR_results', name_transducer),
            axs, ays, inf_alg='transCSSR')

        # Move them away too but append a '_bc'
        os.rename(
            join('transCSSR_results', name_transducer),
            join(
                PATH_DOT_RESULTS,
                'trans_bc_' + name_transducer.replace('+', '_')
            )
        )
        os.rename(
            join('transCSSR_results', name_transducer_dat),
            join(
                PATH_DOT_RESULTS,
                'trans_bc_' + name_transducer_dat.replace('+', '_')
            )
        )
        # Also move away the cssr. e.g.: +X1.dot & +X1.dat_results
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
            'trans_' + name_transducer.replace('+', '_'), C_mu, h_mu,
            'trans_bc_' + name_transducer.replace('+', '_'), C_mu_bc, h_mu_bc
        ])

    # Now that we have moved all of the .dot & .dat_results, let's create the
    #   results.csv that has the C_mu & h_mu
    with open(join(dir_out, 'results.csv'), 'w') as f:
        writer = csv.writer(f, quoting=2)
        writer.writerow(header)
        for row in results:
            writer.writerow(row)


if __name__ == '__main__':
    main(
        path_inputxt='csv/full_Xt.csv',
        path_outputyt='csv/full_Yt.csv',
        prefix_outdir='output_trans',
        max_l=1,
    )
