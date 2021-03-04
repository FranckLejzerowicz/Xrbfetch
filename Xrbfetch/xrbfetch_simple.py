# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from redbiom import fetch
import pandas as pd
import pkg_resources

from os.path import abspath, isfile

RESOURCES = pkg_resources.resource_filename('Xrbfetch', 'resources')


def run_simple(
        m_metadata_file: str,
        o_metadata_file: str,
        o_biom_file: str,
        p_redbiom_context: str,
        p_bloom_sequences: str,
        p_reads_filter: int,
        unique: bool,
        force: bool):
    """
    Main script for fetching a metadata's samples on redbiom,
    solving the ambiguities by keep the prep's sample with the most reads.
    Also possibly activates the keeping of the "best" host sample.
    This is the exact code from Daniel, which is much more concise.

    Parameters
    ----------
    m_metadata_file : str
        Path to metadata file containing the
        samples ot fetch and filter.
    o_metadata_file : str
        Path to the output metadata table file.
    o_biom_file : str
        Path to the output biom table file.
    p_redbiom_context : str
        Redbiom context for fetching 16S data from Qiita.
    p_bloom_sequences : str
        Fasta file containing the sequences known
        to bloom in fecal samples.
    p_reads_filter : int
        Minimum number of reads per sample.
    unique : bool
        Whether to keep a unique sample per host or not.
    force : bool
        Whether to overwrite an already generated output.
    """
    # Read input metadata
    metadata = pd.read_table(m_metadata_file, dtype = str)
    metadata.set_index(metadata.columns.tolist()[0], inplace=True)
    # get the sample IDs
    ids = set(metadata.index)

    # get the default or passed bloom sequences file path
    bloom_sequences_fp = '%s/newblooms.all.fasta' % RESOURCES
    if p_bloom_sequences:
        bloom_sequences_fp = abspath(p_bloom_sequences)
    # read the actual bloom sequences
    blooms = {s.strip() for s in open(bloom_sequences_fp) if not s.startswith('>')}
    # parse length of sequences to fetch from context and trim blooms accordingly
    if p_redbiom_context.split('-')[-2].endswith('nt'):
        length = int(p_redbiom_context.split('-')[-2][:-2])
        blooms = {x[:length] for x in blooms}

    # perform the fetch using the redbiom API
    tab, am = fetch.data_from_samples(p_redbiom_context, ids)
    # ambiguity resolution drops the redbiom tag
    tab = fetch._ambiguity_keep_most_reads(tab, am)

    # filter bloom sequences if present on the sample
    if set(tab.ids(axis='observation')) & blooms:
        tab.filter(set(tab.ids(axis='observation')) - blooms, axis='observation')

    # filter sample having too few reads in total
    tab.filter(lambda v, i, m: v.sum() > p_reads_filter)
    # get the read counts per sample
    reads = dict(zip(tab.ids(axis='sample'), tab.sum(axis='sample')))

    # subset the metadata to the remaining samples
    metadata = metadata.loc[set(tab.ids())]

    if unique:
        # only keep one sample per host for both metadata and data
        metadata['temporary_read_counts'] = [reads[x] for x in metadata.index]
        metadata.sort_values('temporary_read_counts', ascending=False, inplace=True)
        metadata = metadata[~metadata.host_subject_id.duplicated()]
        metadata.drop(columns='temporary_read_counts', inplace=True)
        tab.filter(set(metadata.index)).remove_empty()

    if not isfile(o_biom_file) or force:
        with biom.util.biom_open(o_biom_file, 'w') as o_biom_file_handle:
            tab.to_hdf5(o_biom_file_handle, 'custom')

    if not isfile(o_metadata_file) or force:
        metadata.to_csv(o_metadata_file, index=True, sep='\t')
