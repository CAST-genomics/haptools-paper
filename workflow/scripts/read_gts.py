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
@click.argument("gt_file", type=click.Path(path_type=Path))
@click.option(
    "-v",
    "--max-variants",
    type=int,
    default=None,
    show_default=True,
    help="The max number of variants to read",
)
def main(
    gt_file: Path,
    max_variants: int = 30000,
):
    """
    Read a genotype file and see how long it takes
    """
    log = getLogger("test", "DEBUG")
    if gt_file.suffix == ".pgen":
        gts = GenotypesPLINK(gt_file, log=log)
    else:
        gts = GenotypesRefAlt(gt_file, log=log)
    start = process_time()
    gts.read(max_variants=max_variants)
    end = process_time()
    print(end - start)


if __name__ == "__main__":
    main()
