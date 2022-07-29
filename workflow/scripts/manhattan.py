#!/usr/bin/env python

import re
import sys
import click
import numpy as np
import pandas as pd
from pathlib import Path
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype


@click.command()
@click.option(
    "-l",
    "--linear",
    multiple=True,
    type=click.Path(exists=True, path_type=Path),
    default=[Path("-")],
    show_default="stdin",
    help="PLINK2 .linear files containing the assocation results",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A PNG file containing a Manhattan plot of the results",
)
def main(linear=[sys.stdin], output=sys.stdout):
    """
    Create a manhattan plot from the results of a PLINK2 GWAS
    """
    plink_cols = {
        "#CHROM": "chromosome",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "A1": "allele",
        "TEST": "test",
        "OBS_CT": "num_samps",
        "BETA": "beta",
        "SE": "se",
        "T_STAT": "tstat",
        "P": "pval",
    }
    keep_cols = ["chromosome", "pos", "beta", "se", "pval"]

    # create the plot
    fig, ax = plt.subplots(1, len(linear), figsize=(14, 8), sharex=True, sharey=True)

    # parse the .linear files
    dfs = {}
    max_pval = -1
    for idx, linear_fname in enumerate(linear):
        df = pd.read_csv(
            linear_fname,
            sep="\t",
            header=0,
            names=plink_cols.values(),
            usecols=keep_cols,
        )
        # throw out anything outside of a 1 Mbp window around the peak
        width = 1e6 / 3
        max_idx = min(df.index, key=lambda i: df.iloc[i]["pval"])
        max_pos = df.iloc[max_idx]["pos"]
        df = df[(df["pos"] > max_pos-width) & (df["pos"] < max_pos+width)]
        # -log_10(pvalue)
        df['-log10(p)'] = -np.log(df["pval"])
        df.chromosome = df.chromosome.astype('category')
        df.chromosome = df.chromosome.astype(
            CategoricalDtype(sorted(map(int, df.chromosome.dtype.categories)), ordered=True)
        )
        df = df.sort_values('chromosome')
        # create the plot using pandas and add it to the figure
        df.plot(kind='scatter', x='pos', y='-log10(p)', ax=ax[idx])
        ax_name = Path(Path(linear_fname.stem).stem).stem
        ax[idx].title.set_text(ax_name.replace('-', " + "))
        ax[idx].set(xlabel=None, ylabel=None)
        dfs[ax_name] = df
        max_val = ax[idx].get_ylim()[1]
        if max_pval < max_val:
            max_pval = max_val
    df = pd.concat(dfs, ignore_index=True)

    # set the y-axis limit so that both axes have the same limit
    print(max_pval)
    for idx, linear_fname in enumerate(linear):
        ax[idx].set_ylim(top=max_pval)

    # save the graph
    fig.supxlabel('Chromosomal Position')
    fig.supylabel('-log10(p value)')
    plt.savefig(output)


if __name__ == "__main__":
    main()
