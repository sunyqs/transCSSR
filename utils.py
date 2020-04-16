import csv
import re


def csv_to_list(path):
    """ Convert csv @path into a column-major list """
    # read in csv
    with open(path) as f:
        reader = csv.reader(f, delimiter=',')
        data = [row for row in reader]
    # transpose to column major
    return list(zip(*data))


def get_uniques_from_2d_list(data):
    """ get a list of all unique values from a 2d list sorted alphabetically"""
    # set of all possible outcomes
    uniques = set()
    for column in data:
        # exclude header
        uniques.update(column[1:])
    uniques = list(uniques)
    uniques.sort()
    return uniques


def try_int(s):
    """ https://stackoverflow.com/a/4623518 """
    try:
        return int(s)
    except ValueError:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
        https://stackoverflow.com/a/4623518
    """
    return [try_int(c) for c in re.split('([0-9]+)', s)]


def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
        https://stackoverflow.com/a/4623518
    """
    l.sort(key=alphanum_key)


def rdc_to_xtyt_csv(
        path_rdc='csv/RawDataCategorical_halfs.csv',
        prefix_out='csv/full_'):
    data = csv_to_list(path_rdc)
    xt = []
    yt = []

    for i, col in enumerate(data):
        # even col; xt
        if i % 2 == 0:
            xt.append(('X' + str(len(xt)+1),) + col[1:])
        # odd col; yt
        else:
            yt.append(('Y' + str(len(yt)+1),) + col[1:])

    with open(prefix_out+'Xt.csv', 'w') as f:
        writer = csv.writer(f)
        # transposed
        for row in list(zip(*xt)):
            writer.writerow(row)

    with open(prefix_out+'Yt.csv', 'w') as f:
        writer = csv.writer(f)
        # transposed
        for row in list(zip(*yt)):
            writer.writerow(row)
