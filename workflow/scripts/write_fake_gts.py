#!/usr/bin/env python

import sys
import logging
from pathlib import Path
from time import process_time

import click
import numpy as np

from haptools.data import GenotypesPLINK, GenotypesRefAlt


def getLogger(name: str = None, level: str = "ERROR"):
    """
    Retrieve a Logger object

    Parameters
    ----------
    name : str, optional
        The name of the logging object
    level : str, optional
        The level of verbosity for the logger
    """
    if name is None:
        name = ""
    else:
        name = "." + name

    # create logger
    logger = logging.getLogger("haptools" + name)
    logger.setLevel(level)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # create formatter
    db_time = "|%(asctime)s" if level == "DEBUG" else ""
    formatter = logging.Formatter(
        fmt="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        datefmt="%H:%M:%S",
    )

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    return logger


@click.command()
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "-v",
    "--variants",
    type=int,
    default=30000,
    show_default=True,
    help="The number of variants to write",
)
@click.option(
    "-s",
    "--samples",
    type=int,
    default=500000,
    show_default=True,
    help="The number of samples to write",
)
@click.option(
    "--plink",
    is_flag=True,
    default=False,
    show_default="vcf",
    help="Whether to write PLINK files or a VCF",
)
def main(
    output: Path,
    variants: int = 30000,
    samples: int = 500000,
    plink: bool = False,
):
    """
    Create a random set of genotypes and write them to a genotype file
    """
    log = getLogger("test", "DEBUG")
    if plink:
        gts = GenotypesPLINK(output, log=log)
    else:
        gts = GenotypesRefAlt(output, log=log)
    log.info("Creating variants and samples")
    if 'aaf' in gts.variants.dtype.names:
        gts.variants = np.array([
            (f"SNP{i}", "1", i, 0, "A", "T") for i in range(1, variants+1)
        ], dtype=gts.variants.dtype)
    else:
        gts.variants = np.array([
            (f"SNP{i}", "1", i, "A", "T") for i in range(1, variants+1)
        ], dtype=gts.variants.dtype)
    gts.samples = tuple(
        f"Sample_{i}" for i in range(1, samples+1)
    )
    log.info("Generating fake genotype matrix")
    gts.data = np.random.randint(2, size=variants*samples*2, dtype=np.uint8).reshape(
        (samples, variants, 2)
    )
    log.info("Calling the write method")
    start = process_time()
    gts.write()
    end = process_time()
    print(end - start)


if __name__ == "__main__":
    main()
