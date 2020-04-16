import glob
from transCSSR_bc import *
from os.path import join
import os
import re
import csv
from utils import csv_to_list, sort_nicely, get_uniques_from_2d_list


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


if __name__ == '__main__':
    # parameters
    # specify data source here: CREATE csv directly in main folder!
    MAX_L = 1
    # Directory Prefix to hold the final output files.
    #   The final directory will be appended with MAX_L.
    #   Ex: 'output_cssr' => 'output_cssr_L1'
    PATH_OUTPUT = 'output_cssr'

    # create the final output directory
    PATH_OUTPUT += '_L' + str(MAX_L)
    PATH_DOT_RESULTS = os.path.join(PATH_OUTPUT, 'dot_results')
    os.makedirs(PATH_DOT_RESULTS, exist_ok=True)
    # cols = csv_to_list('csv/RawDataCategorical_halfs.csv')
    cols = csv_to_list('csv/Yt_test.csv')
    unique_outcomes = get_uniques_from_2d_list(cols)

    for col in cols:
        header = col[0]
        stringY = ''.join(col[1:])
        cssr(stringY, unique_outcomes, header, MAX_L)

    # move all of the generated files to the new designated directory
    pattern = r'\+W[0-9]*\.[1-2]\.(dat_results|dot)'
    paths = [
        f for f in glob.glob('transCSSR_results/*') if re.search(pattern, f)]
    for path in paths:
        fname = os.path.split(path)[-1]
        os.rename(
            path, os.path.join(PATH_DOT_RESULTS, fname.replace('+', '_'))
        )

    # Sort the files in a natural order and generate the 4var .csv
    dot_paths = glob.glob(os.path.join(PATH_DOT_RESULTS, '*.dot'))
    sort_nicely(dot_paths)
    with open(os.path.join(PATH_OUTPUT, 'results.csv'), 'w') as f:
        writer = csv.writer(f, quoting=2)
        header = ['dot_filename', 'Cmu', 'H[X_{0}]', 'hmu', 'E']
        writer.writerow(header)
        for p in dot_paths:
            print(p)
            HLs, hLs, hmu, ELs, E, Cmu, etas_matrix = compute_ict_measures(
                p, unique_outcomes, 'transCSSR', L_max=MAX_L)
            fn = os.path.split(p)[-1]
            writer.writerow([fn, Cmu, HLs[0], hmu, E])
