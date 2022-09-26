#!/usr/bin/env python

import click
import numpy as np
from pathlib import Path
from haptools.data import GenotypesPLINK


@click.command()
@click.argument("file1", type=click.Path(exists=True, path_type=Path))
@click.argument("file2", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
def main(file1: Path, file2: Path, output: Path):
    """
    Merge variants from two PGEN files that have the same set of samples
    """
    gts1 = GenotypesPLINK(file1)
    gts2 = GenotypesPLINK(file2)
    out = GenotypesPLINK(output)

    gts1.read()
    gts2.read()

    assert gts1.samples == gts2.samples

    out.samples = gts1.samples
    out.variants = np.concatenate((gts1.variants, gts2.variants))
    out.data = np.concatenate((gts1.data, gts2.data), axis=1)

    out.write()
    

if __name__ == "__main__":
    main()