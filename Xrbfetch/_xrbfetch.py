# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from Xrbfetch._xrbfetch_io import read_meta_pd, write_outputs
from Xrbfetch._xrbfetch_data import (
    run_redbiom_fetch,
    remove_blooms,
    solve_ambiguous_preps,
    filter_reads,
    update_sample_name
)
from Xrbfetch._xrbfetch_duplicates import remove_duplicates

RESOURCES = pkg_resources.resource_filename('Xrbfetch', 'resources')


def xrbfetch(
        m_metadata_file: str,
        o_metadata_file: str,
        o_biom_file: str,
        p_redbiom_context: str,
        p_bloom_sequences: str,
        p_reads_filter : int,
        unique: bool,
        update: bool) -> None:
    """
    Main script for fetching a metadata's samples on redbiom and then,
    filtering the retrieved samples to keep only the "best" in terms of
    numbers of reads and in for ties, in terms of number of features.
    Also includes a filter for keepin such "best" samples for each host.

    Parameters
    ----------
    m_metadata_file : str
        Path to metadata file containing the samples ot fetch and filter.
    o_metadata_file : str
        Path to the output metadata table file.
    o_biom_file : str
        Path to the output biom table file.
    p_redbiom_context : str
        Redbiom context for fetching 16S data from Qiita.
    p_bloom_sequences : str
        Fasta file containing the sequences known to bloom in fecal samples.
    p_reads_filter : int
        Minimum number of reads per sample.
    unique : bool
        Whether to keep a unique sample per host or not.
    update : bool
        Update the sample names to remove Qiita-prep info.
    """
    # Read metadata with first column as index.
    metadata = read_meta_pd(m_metadata_file)
    # Fetch the samples using RedBiom
    redbiom_output, redbiom_samples = run_redbiom_fetch(
        metadata, m_metadata_file, p_redbiom_context)
    # Remove the bloom sequences from the fetched samples.
    biom_tab, biom_tab_removed_ids = remove_blooms(redbiom_output, p_bloom_sequences)
    # Pick the best prep sample per sample to solve the ambiguous fetched results.
    biom_tab_no_ambi, ids_read_counts, ids_feat_counts = solve_ambiguous_preps(
        redbiom_output, biom_tab, list(biom_tab_removed_ids))
    # Filter to a minimum number of reads.
    biom_tab_filt, metadata_filt = filter_reads(
        metadata, biom_tab_no_ambi, p_reads_filter,
        ids_read_counts, ids_feat_counts
    )
    # Remove duplicates (host and sample preps).
    biom_nodup, metadata_edit_best = remove_duplicates(biom_tab_filt, metadata_filt)
    # Update the biom sample name.
    biom_updated = update_sample_name(update, biom_nodup)
    # write the metadata and the biom table outputs
    write_outputs(o_biom_file, o_metadata_file, biom_updated, metadata_edit_best)
