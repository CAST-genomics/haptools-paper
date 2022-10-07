#!/usr/bin/env python

from pathlib import Path

import click
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt


def read_pheno(path: Path) -> pd.Series:
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
        "ERRCODE": "error",
    }
    keep_cols = ["chromosome", "pos", "id", "beta", "se", "pval"]
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        names=plink_cols.values(),
        usecols=keep_cols,
    ).sort_values('pos')
    df["pval"] = df["pval"].fillna(np.inf)
    df['-log10(p)'] = -np.log10(df["pval"])
    # replace -infinity values with 0
    df['-log10(p)'].replace([-np.inf], 0, inplace=True)
    return df

@click.command()
@click.argument("haptools", type=click.Path(exists=True, path_type=Path))
@click.argument("gcta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A PNG file containing the results",
)
def main(haptools: Path, gcta: Path, output: Path):
    """
    Compare two distributions of phenotypes (from haptools vs from GCTA)

    \f
    Parameters
    ----------
    haptools: Path
        The path to a .pheno file output by haptools simphenotype
    gcta: Path
        The path to a .phen file output by gcta --simu-qt
    output: Path
        The path to a PNG file summarizing the results
    """
    phens = pd.concat(
        {'haptools': read_pheno(haptools), 'gcta': read_pheno(gcta)},
        ignore_index=False, names=("method", "sample"),
    ).reset_index()
    haptools = phens[phens["method"] == "haptools"]
    gcta = phens[phens["method"] == "gcta"]

    fig, ax = plt.subplots(1, 2, figsize=(14, 8))
    haptools_color, gcta_color = tuple([plt.cm.tab10(i) for i in range(2)])

    sns.histplot(
        data=haptools,
        ax=ax[0],
        x="-log10(p)",
        color=haptools_color,
        label="haptools",
        element="step",
        alpha=0.25,
    )
    sns.histplot(
        data=gcta,
        ax=ax[0],
        x="-log10(p)",
        color=gcta_color,
        label="gcta",
        element="step",
        alpha=0.25,
    )
    ax[0].legend()

    x0, x1 = ax[0].get_xlim()  # extract the endpoints for the x-axis
    x_pdf = np.linspace(x0, x1, 1000)
    y_pdf = stats.norm.pdf(x_pdf)
    ax[0].plot(x_pdf, y_pdf, 'r', lw=2, label='pdf')

    sm.qqplot(
        data=haptools["-log10(p)"],
        line="45",
        ax=ax[1],
        marker=".",
        label="haptools",
        markerfacecolor=haptools_color,
        markeredgecolor=haptools_color,
        markersize=10,
    )
    sm.qqplot(
        data=gcta["-log10(p)"],
        line="45",
        ax=ax[1],
        marker=".",
        label="gcta",
        markerfacecolor=gcta_color,
        markeredgecolor=gcta_color,
        markersize=10,
    )
    ax[1].legend()

    plt.savefig(output)


if __name__ == "__main__":
    main()
