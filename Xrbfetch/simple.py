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
from Xrbfetch.checks import check_fetched_samples, check_replicates_amount
from Xrbfetch.io import write_summary

RESOURCES = pkg_resources.resource_filename('Xrbfetch', 'resources')


def run_simple(
        m_metadata_file: str,
        o_metadata_file: str,
        o_summary_file: str,
        o_biom_file: str,
        p_redbiom_context: str,
        p_bloom_sequences: str,
        p_reads_filter: int,
        unique: bool,
        force: bool,
        verbose: bool
):
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
    o_summary_file : str
        Path to the output summary file.
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
    verbose : bool
        Whether to show missing, non-fetched samples
        and duplicates or not.
    """
    # -----------
    # Read inputs
    # -----------
    # metadata
    metadata = pd.read_table(m_metadata_file, dtype=str)
    metadata.set_index(metadata.columns.tolist()[0], inplace=True)
    # get the sample IDs
    ids = set(metadata.index)
    # bloom sequences
    bloom_sequences_fp = '%s/newblooms.all.fasta' % RESOURCES
    if p_bloom_sequences:
        bloom_sequences_fp = abspath(p_bloom_sequences)
    # read the actual bloom sequences
    blooms = {s.strip() for s in open(bloom_sequences_fp) if not s.startswith('>')}
    # parse length of sequences to fetch from context and trim blooms accordingly
    if p_redbiom_context.split('-')[-2].endswith('nt'):
        length = int(p_redbiom_context.split('-')[-2][:-2])
        blooms = {x[:length] for x in blooms}

    # ----------------
    # start processing
    # ----------------
    # perform the fetch using the redbiom API
    if verbose:
        print('Fetching %s samples from redbiom... ' % len(ids), end='')
    tab, am = fetch.data_from_samples(p_redbiom_context, ids)
    if verbose:
        print('Done')
    # init summary object with number of samples to fetch
    summary = [['Fetching samples from redbiom', len(ids)]]

    if verbose:
        print('- Samples (all preps) in the fetched biom table:', tab.shape[1])
        check_fetched_samples(list(metadata.index), tab)
        check_replicates_amount(tab)
    summary.append(['Samples in the fetched biom table (all preps)', tab.shape[1]])

    # ambiguity resolution drops the redbiom tag
    if verbose:
        print('- Get most-reads sample from ambiguous redbiom-fetched samples... ', end='')
    tab = fetch._ambiguity_keep_most_reads(tab, am)
    if verbose:
        print('Done -> %s samples' % tab.shape[1])
    summary.append(['Get most-reads sample from ambiguous redbiom-fetched samples', tab.shape[1]])

    # filter bloom sequences if present on the sample
    blooms_in = set(tab.ids(axis='observation')) & blooms
    if blooms_in:
        if verbose:
            print('- Filter blooms... ', end='')
        tab.filter(set(tab.ids(axis='observation')) - blooms, axis='observation')
        if verbose:
            print('Done -> %s samples' % tab.shape[1])
        summary.append(['Filtered %s blooms sequences' % len(blooms_in), tab.shape[1]])
    elif verbose:
        print('- No bloom sequence to filter')
        summary.append(['No bloom sequence to filter', tab.shape[1]])

    # filter sample having too few reads in total
    if verbose:
        print('- Filter biom for min %s reads per sample... ' % p_reads_filter, end='')
    tab.filter(lambda v, i, m: v.sum() > p_reads_filter)
    if verbose:
        print('Done -> %s samples' % tab.shape[1])
    summary.append(['Filter biom for min %s reads per sample' % p_reads_filter, tab.shape[1]])

    # get the read counts per sample
    reads = dict(zip(tab.ids(axis='sample'), tab.sum(axis='sample')))

    # subset the metadata to the remaining samples
    metadata = metadata.loc[set(tab.ids())]

    if unique:
        # only keep one sample per host for both metadata and data
        if verbose:
            print('- Keep the best sample per "host_subject_id"... ', end='')
        metadata['temporary_read_counts'] = [reads[x] for x in metadata.index]
        metadata.sort_values('temporary_read_counts', ascending=False, inplace=True)
        metadata = metadata[~metadata.host_subject_id.duplicated()]
        metadata.drop(columns='temporary_read_counts', inplace=True)
        tab.filter(set(metadata.index)).remove_empty()
        if verbose:
            print('Done -> %s samples' % tab.shape[1])
        summary.append(['Keep the best sample per "host_subject_id"', tab.shape[1]])

    if not isfile(o_biom_file) or force:
        with biom.util.biom_open(o_biom_file, 'w') as o_biom_file_handle:
            tab.to_hdf5(o_biom_file_handle, 'custom')

    if not isfile(o_metadata_file) or force:
        metadata.to_csv(o_metadata_file, index=True, sep='\t')

    write_summary(o_summary_file, summary)

    print('Outputs:')
    print(o_biom_file)
    print(o_metadata_file)
    print(o_summary_file)
