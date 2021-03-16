# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import json
import biom
import subprocess
import pandas as pd

from biom.util import biom_open
from os.path import dirname, isdir, isfile, splitext


def get_outputs(metadata_edit_best: pd.DataFrame,
                o_metadata_file: str, o_biom_file: str) -> tuple:
    """
    Write the metadata and the biom table outputs.

    Parameters
    ----------
    metadata_edit_best : pd.DataFrame
        Corresponding metadata table.
    o_metadata_file : str
        Path to the output metadata table file.
    o_biom_file : str
        Path to the output biom table file.

    Returns
    -------
    o_metadata_file : str
        Path to the output metadata table file.
    o_biom_file : str
        Path to the output biom table file.
    """
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
    return o_metadata_file, o_biom_file


def write_outputs(
        o_biom_file: str,
        o_metadata_file: str,
        biom_updated: biom.Table,
        metadata_edit_best: pd.DataFrame,
        dim: bool = False) -> None:
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
    dim : bool
        Whether to add the number of samples in the
        final biom file name before extension or not.
    """
    if dim:
        o_metadata_file, o_biom_file = get_outputs(
            metadata_edit_best, o_metadata_file, o_biom_file)

    if biom_updated.shape[0]:
        print('Outputs:')
        if o_metadata_file[0] == '/':
            if not isdir(dirname(o_metadata_file)):
                os.makedirs(dirname(o_metadata_file))
        metadata_edit_best.to_csv(o_metadata_file, index=False, sep='\t')
        print(o_metadata_file)

        if not isdir(dirname(o_biom_file)):
            os.makedirs(dirname(o_biom_file))
        with biom_open(o_biom_file, 'w') as f:
            biom_updated.to_hdf5(f, 'Xrbfetch')
        print(o_biom_file)


def read_biom(redbiom_output: str) -> biom.Table:
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
    print('- Load biom table... ', end='')
    biom_tab = biom.load_table(redbiom_output)
    print('Done -> %s samples (all preps)' % biom_tab.shape[0])
    return biom_tab


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
    metadata = pd.read_table(metadata_file, dtype=str)
    metadata.rename(columns={metadata.columns.tolist()[0]: 'sample_name'}, inplace=True)
    metadata.columns = [x.lower() for x in metadata.columns]
    # remove NaN only columns
    metadata = metadata.loc[:, ~metadata.isna().all()]
    return metadata


def make_samples_list_tmp(m_metadata_file: str,
                          metadata: pd.DataFrame,
                          context: str,
                          timetoken: str) -> str:
    """
    Create the file passed to redbiom and
    containing the samples in one list.

    Parameters
    ----------
    m_metadata_file : str
        Path to metadata file containing
        the samples to fetch and filter.
    metadata : pd.DataFrame
        Metadata table.
    context : str
        Redbiom context for fetching 16S data from Qiita.
    timetoken : str
        Time stamp for unique file naming.

    Returns
    -------
    redbiom_samples : str
        Path to the file passed to redbiom and
        containing the samples in one list.
    """
    redbiom_samples = '%s_redbiom_sams_%s_%s.tmp' % (
        splitext(m_metadata_file)[0], context, timetoken)
    with open(redbiom_samples, 'w') as o:
        for feature_name in metadata.sample_name:
            o.write('%s\n' % feature_name)
    return redbiom_samples


def run_fetch(
        redbiom_samples: str,
        context: str,
        redbiom_output: str) -> None:
    """
    Delete intermediate metadata files.

    Parameters
    ----------
    redbiom_samples : str
        Path to metadata file containing
        the samples ot fetch and filter.
    context : str
        redbiom context name.
    redbiom_output : str
        Path to the fetched biom file.
    """
    cmd = [
        'redbiom', 'fetch', 'samples',
        '--from', redbiom_samples,
        '--context', context,
        '--output', redbiom_output
    ]
    subprocess.call(cmd)
    print(' '.join(cmd))


def delete_files(
        redbiom_output: str,
        redbiom_samples: str) -> None:
    """
    Delete intermediate metadata files.

    Parameters
    ----------
    redbiom_output : str
        Path to fetched biom file.
    redbiom_samples : str
        Path to samples to fetch.
    """
    if isfile(redbiom_output):
        os.remove(redbiom_output)
    redbiom_output_amb = '%s.ambiguities' % splitext(redbiom_output)[0]
    if isfile(redbiom_output_amb):
        os.remove(redbiom_output_amb)
    if isfile(redbiom_samples):
        os.remove(redbiom_samples)


def write_summary(
        o_summary_file: str,
        summary: list) -> None:
    """
    Write a two-columns file.

    Parameters
    ----------
    o_summary_file : str
        Path to summary file.
    summary : list
        List of two-items lists
        [[str, int], [str, int], ...]
    """
    with open(o_summary_file, 'w') as o:
        o.write('steps\tsamples\n')
        for step, num in summary:
            o.write('%s\t%s\n' % (step, num))
    print(o_summary_file)


def read_json_ambiguities_file(
        redbiom_output: str) -> dict:
    """
    Read a json file.

    Parameters
    ----------
    redbiom_output : str
        Path to json file.
    """
    json_ambi = {}
    json_ambi_fp = '%s.ambiguities' % redbiom_output
    if isfile(json_ambi_fp):
        with open(json_ambi_fp) as f:
            json_ambi = json.load(f)
    return json_ambi
