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
from os.path import dirname, isdir, splitext


def write_outputs(
        o_biom_file: str,
        o_metadata_file: str,
        biom_updated: biom.table,
        metadata_edit_best: pd.DataFrame,
        redbiom_output: str,
        redbiom_samples: str,
        dim: bool) -> None:
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
    dim : bool
        Whether to add the number of samples in the final biom file name before extension or not.
    """
    if dim:
        n_sams = metadata_edit_best.shape[0]
        o_metadata_file = '%s_%ss%s' % (
            splitext(o_metadata_file)[0],
            n_sams,
            splitext(o_metadata_file)[1]
        )
        o_biom_file = '%s_%ss%s' % (
            splitext(o_biom_file)[0],
            n_sams,
            splitext(o_biom_file)[1]
        )

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


def read_biom(redbiom_output: str) -> tuple:
    """
    Read biom file

    Parameters
    ----------
    redbiom_output : str
        The biom table returned by redbiom.

    Returns
    -------
    biom_tab : biom.table
        Feature table retrieved from redbiom.
    biom_tab_sams : list
        Samples of the feature table.
    """
    print(' - Load biom table... ', end='')
    biom_tab = biom.load_table(redbiom_output)
    biom_tab_sams = biom_tab.ids(axis='sample').tolist()
    print('Done -> %s samples' % len(biom_tab_sams))
    return biom_tab, biom_tab_sams


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
            first_line = line.split(sep)
            break
    else:
        print('no separator found among: "<tab>", ",", ";"\nExiting')
        sys.exit(1)

    strs = {first_line[0]: 'str'}
    strs_up = {}
    for col in first_line:
        if 'qiita_prep_id' in col:
            strs_up[col] = 'str'
        elif 'sample_name' in col:
            strs_up[col] = 'str'
    strs.update(strs_up)

    metadata = pd.read_csv(metadata_file, header=0, sep=sep,
                          dtype=strs, low_memory=False)
    metadata.rename(columns={first_line[0]: 'sample_name'}, inplace=True)
    metadata.columns = [x.lower() for x in metadata.columns]
    # remove NaN only columns
    metadata = metadata.loc[:,~metadata.isna().all()]
    return metadata