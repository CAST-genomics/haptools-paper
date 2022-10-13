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
        path,
        header=None,
        delim_whitespace=True,
        index_col=0,
        names=("sample", "Phenotypes"),
        comment="#",
    )
    # df["Phenotypes"] = (df["Phenotypes"]-df["Phenotypes"].mean())/df["Phenotypes"].std()
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
    # create a DataFrame where each row is indexed by its sample ID and
    # the columns are ('haptools', 'gcta')
    phens = pd.concat(
        {'haptools': read_pheno(haptools), 'gcta': read_pheno(gcta)},
        ignore_index=False, names=("method", "sample"),
    ).unstack(level=0)["Phenotypes"]

    sns.scatterplot(data=phens, x="gcta", y="haptools")
    plt.axline([0, 0], [1, 1])

    plt.savefig(output)


if __name__ == "__main__":
    main()
