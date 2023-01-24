#!/usr/bin/env python

import sys
import csv
import re
import argparse
import matplotlib.pyplot as plt
from matplotlib import rcParams

from dark.dimension import dimensionalIterator


DIGITS = re.compile(r"\d+")


def getData(filename):
    """
    Read the CSV data

    @return: A C{dict} of C{dict}s keyed by C{str} name of the mutation (e.g.,
        "ORF1:Y2435K") and C{str} sample name, with value the C{float} percentage
        of the mutation.
    """
    result = {}

    with open(filename) as fp:
        reader = csv.reader(fp, delimiter=",")
        header = next(reader)
        assert header == ["Mutation", "sample", "%Mutation", "%Rest"]

        for row in reader:
            name, approach, percent, _ = row
            # name, approach, percent = row[:3]

            if name not in result:
                result[name] = {}

            if approach in result[name]:
                print(f"Approach name {approach!r} found twice for {name!r}.",
                      file=sys.stderr)
                sys.exit(1)

            if percent == "":
                result[name][approach] = None
            else:
                result[name][approach] = float(percent)

    return result


def checkApproachNames(data):
    """
    Check that the same set of approach names is present for each name
    in the data.

    @return: the C{int} number of approaches.
    """
    first = True

    for name in data:
        if first:
            approaches = set(data[name])
            first = False
        else:
            assert set(data[name]) == approaches

    return approaches


def plotPie(percent, ax, row, col, approach, name):
    """
    Make a pie chart with the given percentage, drawn in ax.
    """
    if row == 0:
        ax.set_title(name)

    if col == 0:
        ax.set_ylabel(approach)

    if percent is None:
        # Don't show this subplot, as coverage was too low to have a percent.
        ax.axis("off")
    else:
        # labels = "Mutation", "No mutation"
        sizes = percent, 100.0 - percent
        ax.pie(sizes)


def rowKey(name):
    """
    name is a value like "Serum native".
    """
    tissue, approach = name.split()
    approachNumber = ["native", "cap"].index(approach)

    # Old-fashioned verbose way to do it.
    # if tissue == "Serum":
    #     tissueNumber = 0
    # else:
    #     tissueNumber = 1

    # New way to do it.
    # tissueNumber = 0 if tissue == "Serum" else 1

    # How to do it if you like logic and can convert a Boolean to a number.
    # Look up George Boole on Wikipedia.
    tissueNumber = int(tissue != "Serum")

    return tissueNumber, tissue, approachNumber


def columnKey(name):
    """
    name is a value like "ORF1:V345K". We want to return the ORF and then
    the integer site (345 in this case).
    """
    orf, mutation = name.split(':')
    orfNumber = int(DIGITS.search(orf).group())
    position = int(DIGITS.search(mutation).group())

    return orfNumber, position


def makePlot(data, output):
    """
    Make the plot for Jules.
    """
    approaches = checkApproachNames(data)
    rows = len(approaches)
    cols = len(data)

    figure, ax = plt.subplots(
        rows, cols, squeeze=False, sharex=True, figsize=(3 * cols, 1.5 * rows)
    )

    coords = dimensionalIterator((rows, cols))

    sortedApproaches = sorted(approaches, key=rowKey)

    for approach in sortedApproaches:
        for name in sorted(data, key=columnKey):
            row, col = next(coords)
            plotPie(data[name][approach], ax[row][col], row, col, approach, name)

    plt.savefig(output, bbox_inches="tight")


def main():
    parser = argparse.ArgumentParser(description="Plot mutations for Jules")

    parser.add_argument(
        "file",
        help="The CSV file to read.",
    )

    parser.add_argument(
        "--output", required=True,
        help="The output file to write.",
    )

    args = parser.parse_args()

    rcParams.update(
        {
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "figure.titlesize": 16,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "legend.fontsize": 12,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )

    data = getData(args.file)
    makePlot(data, args.output)


if __name__ == "__main__":
    main()
