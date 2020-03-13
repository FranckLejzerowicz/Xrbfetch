# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import biom
import pandas as pd
from biom.util import biom_open
from os.path import dirname, isdir


def write_outputs(
        o_biom_file: str,
        o_metadata_file: str,
        biom_updated: biom.table,
        metadata_edit_best: pd.DataFrame,
        redbiom_output: str,
        redbiom_samples: str) -> None:
    """
    Write the metadata and the biom table outputs.

    Parameters
    ----------
    o_metadata_file : str
        Path to the output metadata table file.
    o_biom_file : str
        Path to the output biom table file.
    biom_updated : biom.table
        The biom table without ambiguous samples,
        with a min number of reads per sample, and
        without duplicated sample per host, and
        with samples re-named as per AGP system.
    metadata_edit_best : pd.DataFrame
        Corresponding metadata table.
    redbiom_output : str
        Path to the biom table returned by redbiom.
    redbiom_samples : str
        Path to the file containing the samples used for fetching.
    """
    print(' - Write files:')
    if not isdir(dirname(o_metadata_file)):
        os.makedirs(dirname(o_metadata_file))
    metadata_edit_best.to_csv(o_metadata_file, index=False, sep='\t')
    print('   *', o_metadata_file)

    if not isdir(dirname(o_biom_file)):
        os.makedirs(dirname(o_biom_file))
    with biom_open(o_biom_file, 'w') as f:
        biom_updated.to_hdf5(f, 'Xrbfetch')
    print('   *', o_biom_file)
    os.remove(redbiom_output)
    os.remove(redbiom_samples)


def read_meta_pd(metadata_file: str) -> pd.DataFrame:
    """
    Read metadata with first column as index.

    Parameters
    ----------
    metadata_file : str
        Path to metadata file containing the samples ot fetch and filter.

    Returns
    -------
    metadata : pd.DataFrame
        Metadata table.
    """
    with open(metadata_file) as f:
        for line in f:
            break
    header = line.strip()
    for sep in ['\t', ';', ',']:
        if len(header.split(sep))>1 and len(header.split(sep)) == (header.count(sep)+1):
            first_col = line.split(sep)[0]
            break
    else:
        print('no separator found among: "<tab>", ",", ";"\nExiting')
        sys.exit(1)
    metadata = pd.read_csv(metadata_file, header=0, sep=sep,
                          dtype={first_col: str}, low_memory=False)
    metadata.rename(columns={first_col: 'sample_name'}, inplace=True)
    metadata.columns = [x.lower() for x in metadata.columns]
    # remove NaN only columns
    metadata = metadata.loc[:,~metadata.isna().all()]
    return metadata