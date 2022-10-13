#!/usr/bin/env python

import re
import sys
import click
import numpy as np
import pandas as pd
from glob import glob
from pathlib import Path
import matplotlib.pyplot as plt


def get_paths(path: Path, binary: bool = False) -> dict:
    """
    Retrieve the paths to each of the PLINK2 GWAS files

    Parameters
    ----------
    path: Path
        The path to a directory containing all of the PLINK2 GWAS files

        The directory path itself must contain a "beta" wildcard
    binary: bool, optional
        True if the phenotype was simulated for a case/control trait

    Returns
    -------
    dict
        A dict mapping matched beta values to a list of PLINK2 GWAS files
    """
    paths = glob(str(path).format(beta="*"))
    rgx = re.compile(str(path).format(beta="(.*)"))
    ending = ".glm.logistic.hybrid" if binary else ".glm.linear"
    return {
        float(rgx.search(p)[1]): list(Path(p).glob("*"+ending)) for p in paths
    }


def check_significant(snp: str, gwas: Path, alpha: float = 5e-8, binary: bool=False) -> bool:
    """
    Check whether a SNP is significant in a GWAS result

    Parameters
    ----------
    snp: str
        The ID of the SNP
    gwas: Path
        The path to a PLINK2 ".glm.linear" file containing the results of a GWAS
    alpha: float, optional
        The significance threshold
    binary: bool, optional
        True if the phenotype was simulated for a case/control trait

    Returns
    -------
    bool
        True if the SNP is significant under the chosen threshold and False otherwise
    """
    if binary:
        plink_cols = {
            "#CHROM": "chromosome",
            "POS": "pos",
            "ID": "id",
            "REF": "ref",
            "ALT": "alt",
            "A1": "allele",
            "FIRTH?": "firth",
            "TEST": "test",
            "OBS_CT": "num_samps",
            "OR": "or",
            "LOG(OR)_SE": "se",
            "Z_STAT": "stat",
            "P": "pval",
            "ERRCODE": "error",
        }
    else:
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
            "T_STAT": "stat",
            "P": "pval",
            "ERRCODE": "error",
        }
    keep_cols = ["id", "pval"]
    df = pd.read_csv(
        gwas,
        sep="\t",
        header=0,
        names=plink_cols.values(),
        usecols=keep_cols,
    )
    return df[df["id"] == snp]["pval"].values[0] < alpha


@click.command()
@click.argument("snp", type=str)
@click.argument("ancestry", type=click.Path(path_type=Path))
@click.argument("normal", type=click.Path(path_type=Path))
@click.option(
    "-b",
    "--binary",
    is_flag=True,
    default=False,
    show_default=True,
    help="True if the simulated phenotype was binary",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A PNG file containing two power plots comparing the results",
)
def main(snp: str, ancestry: Path, normal: Path, binary: bool=False, output: str=sys.stdout):
    """
    Create power plots from the results of a PLINK2 GWAS with an ancestral effect
    vs a normal effect

    \f
    Parameters
    ----------
    snp: str
        The ID of the simulated causal variant
    ancestry: Path
        The path to a directory containing the GWAS on a simulated ancestral effect

        The path should contain a wildcard "beta" for each of the beta values
    normal: Path
        The path to a directory containing the GWAS on a simulated normal effect

        The path should contain a wildcard "beta" for each of the beta values
    output: Path
        The path to a PNG file to which to write the resulting power plots
    """
    # open each path and compute a power summary statistic:
    #   the number of times that the chosen causal SNP is significant
    ancestry = {
        beta: sum(
            check_significant(snp, path, binary=binary)
            for path in paths
        )/len(paths) for beta, paths in get_paths(ancestry, binary=binary).items()
    }
    normal = {
        beta: sum(
            check_significant(snp, path, binary=binary)
            for path in paths
        )/len(paths) for beta, paths in get_paths(normal, binary=binary).items()
    }
    # sort the values so they are in ascending order
    ancestry = dict(sorted(ancestry.items(), key=lambda i: i[0]))
    normal = dict(sorted(normal.items(), key=lambda i: i[0]))
    # now, plot these values
    plt.figure(figsize=(14,8))
    #plt.gca().set_xticks(ancestry.keys())
    plt.plot(ancestry.keys(), ancestry.values(), linestyle="-", marker="o", label="ancestry")
    plt.plot(normal.keys(), normal.values(), linestyle="--", marker="s", label="normal")
    plt.ylabel("Power")
    plt.xlabel("Effect size")
    plt.legend()
    plt.savefig(output)


if __name__ == "__main__":
    main()
