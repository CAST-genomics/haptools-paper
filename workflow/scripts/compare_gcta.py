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
    df = pd.read_csv(
        path, header=None, delim_whitespace=True, index_col=0, names=("sample", "Phenotypes"), comment="#"
    )
    df["Phenotypes"] = (df["Phenotypes"]-df["Phenotypes"].mean())/df["Phenotypes"].std()
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

    fig, ax = plt.subplots(1, 3, figsize=(14, 8))

    sns.histplot(data=phens, x="Phenotypes", hue="method", ax=ax[0], element="step")

    x0, x1 = ax[0].get_xlim()  # extract the endpoints for the x-axis
    x_pdf = np.linspace(x0, x1, 1000)
    y_pdf = stats.norm.pdf(x_pdf)
    ax[0].plot(x_pdf, y_pdf, 'r', lw=2, label='pdf')

    sm.qqplot(data=haptools["Phenotypes"], line="s", ax=ax[1])
    sm.qqplot(data=gcta["Phenotypes"], line="s", ax=ax[2])

    plt.savefig(output)


if __name__ == "__main__":
    main()
