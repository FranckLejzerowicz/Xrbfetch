# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import abspath, splitext
import pkg_resources

from Xrbfetch.simple import run_simple
from Xrbfetch.io import read_meta_pd, read_biom, write_outputs
from Xrbfetch.data import (
    run_redbiom_fetch, remove_blooms, get_reads_features_counts, solve_ambiguous_preps,
    merge_read_features_counts, filter_reads, update_sample_name
)
from Xrbfetch.checks import check_fetched_samples, check_replicates_amount, potential_stop
from Xrbfetch.duplicates import remove_duplicates

RESOURCES = pkg_resources.resource_filename('Xrbfetch', 'resources')


def xrbfetch(
        m_metadata_file: str,
        o_metadata_file: str,
        o_biom_file: str,
        p_redbiom_context: str,
        p_bloom_sequences: str,
        p_reads_filter: int,
        unique: bool,
        update: bool,
        dim: bool,
        force: bool,
        simple: bool,
        verbose: bool) -> None:
    """
    Main script for fetching a metadata's samples on redbiom and then,
    filtering the retrieved samples to keep only the "best" in terms of
    numbers of reads and in for ties, in terms of number of features.
    Also includes a filter for keepin such "best" samples for each host.

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
    update : bool
        Update the sample names to remove Qiita-prep info.
    dim : bool
        Whether to add the number of samples in the final
        biom file name before extension or not.
    force : bool
        Re-fetch and not use an 'already-fetched' table.
    simple : bool
        Whether to perform the steps using one-liners or not.
    verbose : bool
        Whether to show missing, non-fetched samples
        and duplicates or not.
    """

    # Read metadata with first column as index.
    m_metadata_file = abspath(m_metadata_file)
    o_metadata_file = abspath(o_metadata_file)
    o_summary_file = '%s_summary.tsv' % splitext(o_metadata_file)[0]
    o_biom_file = abspath(o_biom_file)

    if simple:
        # added as per suggestions from Daniel
        run_simple(
            m_metadata_file, o_metadata_file, o_summary_file,
            o_biom_file, p_redbiom_context, p_bloom_sequences,
            p_reads_filter, unique, force, verbose
        )
    else:
        metadata = read_meta_pd(m_metadata_file)
        summary = [['Fetching samples from redbiom', metadata['sample_name'].nunique()]]
        # Fetch the samples using redbiom
        redbiom_output, redbiom_samples = run_redbiom_fetch(
            metadata, m_metadata_file, o_metadata_file, p_redbiom_context, force)

        # Read biom file and show non fetched samples and replication amount
        biom_tab = read_biom(redbiom_output)
        potential_stop(biom_tab, o_summary_file, summary, redbiom_output, redbiom_samples)

        summary.append(['Samples in the fetched biom table (all preps)', biom_tab.shape[0]])
        if verbose:
            check_fetched_samples(list(metadata['sample_name']), biom_tab)
            check_replicates_amount(biom_tab)

        # Remove the bloom sequences from the fetched samples.
        if p_bloom_sequences != 'no':
            length = int(p_redbiom_context.split('-')[-2][:-2])
            biom_tab = remove_blooms(biom_tab, length, p_bloom_sequences)
            summary.append(['Filtered blooms sequences', biom_tab.shape[1]])
        potential_stop(biom_tab, o_summary_file, summary, redbiom_output, redbiom_samples)

        read_counts, feat_counts = get_reads_features_counts(biom_tab)

        # Pick the best prep per sample to solve the ambiguous fetched results.
        biom_tab_no_ambi = solve_ambiguous_preps(
            redbiom_output, biom_tab, read_counts, feat_counts)
        summary.append(['Get most-reads (or most-features) sample from ambiguous samples',
                        biom_tab_no_ambi.shape[1]])
        potential_stop(biom_tab_no_ambi, o_summary_file, summary, redbiom_output, redbiom_samples)

        # Merge reads and features count to metadata
        metadata_counts = merge_read_features_counts(
            metadata, biom_tab_no_ambi, read_counts, feat_counts)

        # Filter to a minimum number of reads.
        biom_tab_filt, metadata_filt = filter_reads(
            metadata_counts, biom_tab_no_ambi, p_reads_filter)
        summary.append(['Filter biom for min %s reads per sample' % p_reads_filter, biom_tab_filt.shape[1]])
        potential_stop(biom_tab_filt, o_summary_file, summary, redbiom_output, redbiom_samples)

        # Remove duplicates (host and sample preps).
        biom_nodup, metadata_edit_best = remove_duplicates(
            biom_tab_filt, metadata_filt, unique)
        summary.append(['Keep the best sample per "host_subject_id"', biom_nodup.shape[1]])
        potential_stop(biom_nodup, o_summary_file, summary, redbiom_output, redbiom_samples)

        # Update the biom sample name.
        biom_updated = update_sample_name(update, biom_nodup)

        # write the metadata and the biom table outputs
        write_outputs(
            o_biom_file, o_metadata_file, biom_updated,
            metadata_edit_best, dim)

        potential_stop(biom_tab, o_summary_file, summary, redbiom_output, redbiom_samples, True)
