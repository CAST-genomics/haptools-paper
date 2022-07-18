#!/usr/bin/env python

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
    type=click.Path(exists=True, path_type=Path),
    default=Path("-"),
    show_default="stdin",
    help="A PLINK2 .linear file containing the assocation results",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A PNG file containing a Manhattan plot of the results",
)
def main(linear=sys.stdin, output=sys.stdout):
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

    # parse the .linear file
    df = pd.read_csv(
        linear,
        sep="\t",
        header=0,
        names=plink_cols.values(),
        usecols=keep_cols,
    )
    # -log_10(pvalue)
    df['-10log10(p)'] = -10*np.log(df["pval"])
    df.chromosome = df.chromosome.astype('category')
    df.chromosome = df.chromosome.astype(
        CategoricalDtype(sorted(map(int, df.chromosome.dtype.categories)), ordered=True)
    )
    df = df.sort_values('chromosome')

    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    df_grouped = df.groupby(('chromosome'))
    fig, ax = plt.subplots(1, 2, figsize=(14, 8))

    # qqplot
    sm.qqplot(df['pval'], ax=ax[0], line='45')

    # manhattan plot
    x_labels = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='pos', y='-10log10(p)', ax=ax[1])
        x_labels.append(name)
    ax[1].set_xticklabels(x_labels)
    # set axis limits
    ax[1].set_xlim([0, max(df['pos'])])
    # ax[1].set_ylim([0, 3.5])
    # x axis label
    ax[1].set_xlabel('Chromosome')

    # show the graph
    plt.savefig(output)


if __name__ == "__main__":
    main()
