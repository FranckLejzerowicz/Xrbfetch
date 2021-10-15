# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import datetime
import pandas as pd
import pkg_resources

from os.path import abspath, isfile, splitext
from Xrbfetch.io import make_samples_list_tmp, run_fetch, read_json_ambiguities_file


RESOURCES = pkg_resources.resource_filename('Xrbfetch', 'resources')


def update_sample_name(
        update: bool,
        biom_nodup: biom.Table) -> biom.Table:
    """
    Update the biom sample name.
    For edge cases I spotted (confirmed by Daniel).

    Parameters
    ----------
    update : bool
        Update the sample names to remove Qiita-prep info.
    biom_nodup : biom.table
        The biom table without ambiguous samples,
        with a min number of reads per sample, and
        without duplicated sample per host.

    Returns
    -------
    biom_updated : biom.table
        The biom table without ambiguous samples,
        with a min number of reads per sample, and
        without duplicated sample per host, and
        with samples re-named as per AGP system.
    """
    if update:
        ids_map = {}
        ids_res = []
        for sam in biom_nodup.ids(axis='sample'):
            updated_sam = sam.rsplit('.', 1)[0]
            ids_res.append([sam, updated_sam])
            ids_map[sam] = updated_sam
        print('- Update sample name to remove prep file info... ', end='')
        biom_nodup.update_ids(id_map=ids_map, axis='sample', inplace=True)
        print('Done')
    return biom_nodup


def merge_read_features_counts(
        metadata: pd.DataFrame,
        biom_tab_no_ambi: biom.Table,
        read_counts: dict,
        feat_counts: dict) -> pd.DataFrame:
    """
    Filter to a minimum number of reads.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for the included samples only.
    biom_tab_no_ambi : biom.table
        The biom table without ambiguous samples.
    read_counts : dict
        Number of reads per sample.
    feat_counts : dict
        Number of features per sample.

    Returns
    -------
    metadata_counts : pd.DataFrame
        Metadata table with reads and features counts
    """
    # get the sample IDs (with prep)
    ids = biom_tab_no_ambi.ids(axis='sample')
    # get sample IDs and reads/features counts
    # removed from the metadata if present
    columns = ['sample_name']
    for col in ['qiita_prep_id', 'read_count',
                'feature_count', 'orig_sample_name']:
        if col in metadata.columns:
            metadata.drop(columns=col, inplace=True)
        columns.append(col)

    # re-create these IDs and reads/features counts afresh
    ids_read_feat_counts_pd = pd.DataFrame([
        [
            ID.rsplit('.', 1)[0],  # sample ID without prep number
            ID.split('.')[-1],     # prep number
            read_counts[ID],       # read number
            feat_counts[ID],       # feature number
            ID                     # full sample ID
        ] for ID in ids],
        columns=columns
    )
    # merge these fresh values back into the metadata
    metadata_counts = metadata.merge(
        ids_read_feat_counts_pd, on='sample_name', how='right')
    return metadata_counts


def filter_reads(
        metadata_counts: pd.DataFrame,
        biom_tab_no_ambi: biom.Table,
        reads_filter: int) -> tuple:
    """
    Filter to a minimum number of reads.

    Parameters
    ----------
    metadata_counts : pd.DataFrame
        Metadata for the included samples only.
    biom_tab_no_ambi : biom.table
        The biom table without ambiguous samples.
    reads_filter : int
        Minimum number of reads per sample.

    Returns
    -------
    biom_tab_filt : biom.table
        The biom table without ambiguous samples
        and with a min number of reads per sample.
    metadata_filt : pd.DataFrame
        Corresponding metadata table.
    """

    # Filter to keep only the samples with min number reads
    print('- Filter biom for min %s reads per sample... ' % reads_filter, end='')
    metadata_filt = metadata_counts.loc[
        metadata_counts['read_count'] >= reads_filter].copy()
    biom_tab_filt = biom_tab_no_ambi.filter(
        ids_to_keep=metadata_filt['orig_sample_name'].tolist(),
        axis='sample'
    ).copy()
    biom_tab_filt.remove_empty(
        axis='observation', inplace=True
    )
    print('Done -> %s samples' % biom_tab_filt.shape[1])
    return biom_tab_filt, metadata_filt


def run_redbiom_fetch(
        metadata: pd.DataFrame,
        m_metadata_file: str,
        o_metadata_file: str,
        p_redbiom_context: str,
        force: bool) -> tuple:
    """
    Fetch the samples using RedBiom.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for the included samples only.
    m_metadata_file : str
        Path to input metadata for the included samples only.
    o_metadata_file : str
        Path to output metadata for the included samples only.
    p_redbiom_context : str
        Redbiom context for fetching 16S data from Qiita.

    Returns
    -------
    redbiom_output : str
        Path to the biom table returned by redbiom.
    """
    # define temporary files: for sample_ids file and redbiom output
    timetoken = str(datetime.datetime.now()).split('.')[0].replace(' ', '-').replace(':', '-')
    redbiom_output = '%s_redbiom_%s_%s.biom' % (splitext(o_metadata_file)[0], p_redbiom_context, timetoken)
    redbiom_samples = ''
    if not force and isfile(redbiom_output):
        print('Using an already generated biom file for this metadata file:\n '
              '-> %s (%s samples) <-' % (redbiom_output, metadata.shape[0]))
    else:
        context = 'Deblur-Illumina-16S-V4-150nt-780653'
        if p_redbiom_context:
            context = p_redbiom_context

        redbiom_samples = make_samples_list_tmp(m_metadata_file, metadata, context, timetoken)
        print('Fetching %s samples from redbiom... ' % metadata.shape[0], end='')
        run_fetch(redbiom_samples, context, redbiom_output)
        print('Done')

    return redbiom_output, redbiom_samples


def remove_blooms(
        biom_tab: biom.Table,
        p_bloom_sequences: str) -> biom.Table:
    """
    Remove the bloom sequences from the fetched samples.

    Parameters
    ----------
    biom_tab : biom.table
        Feature table retrieved from redbiom.
    p_bloom_sequences : str
        Fasta file containing the sequences known to bloom in fecal samples.

    Returns
    -------
    biom_tab : biom.table
        Feature table retrieved from redbiom, without blooms.
    biom_tab_removed_ids : list
        Samples ids removed as bloom dropouts.
    """

    bloom_sequences_fp = '%s/newblooms.all.fasta' % RESOURCES
    if p_bloom_sequences:
        bloom_sequences_fp = abspath(p_bloom_sequences)
    bloom_seqs = set([x.strip() for x in open(bloom_sequences_fp).readlines()
                      if x[0] != '>'])
    print(bloom_seqs)

    in_out = {'in': [], 'out': []}
    for feature in biom_tab.ids(axis='observation'):
        if feature in bloom_seqs:
            print(feature)
            in_out['in'].append(feature)
        else:
            in_out['out'].append(feature)

    print('- Filter blooms... ', end='')
    biom_tab.filter(
        ids_to_keep=in_out['out'],
        axis='observation',
        inplace=True
    )
    biom_tab.remove_empty(axis='sample', inplace=True)
    print('Done -> %s samples' % biom_tab.shape[1])
    return biom_tab


def get_reads_features_counts(
        biom_tab: biom.Table) -> tuple:
    """
    Count the number of reads and features
    for each prep's sample.

    Parameters
    ----------
    biom_tab : biom.Table
        Current data retrieved from redbiom,
        without blooms.

    Returns
    -------
    ids_read_counts : dict
        Number of reads per sample.
    ids_feat_counts : dict
        Number of features per sample.
    """
    ids = biom_tab.ids(axis='sample')
    read_counts = biom_tab.sum(axis='sample')
    ids_read_counts = dict((ids[x], read_counts[x]) for x in range(len(ids)))
    feat_counts = biom_tab.nonzero_counts(axis='sample', binary=True)
    ids_feat_counts = dict((ids[x], feat_counts[x]) for x in range(len(ids)))
    return ids_read_counts, ids_feat_counts


def reverse_json_ambi(
        json_ambi: dict,
        ids: set) -> dict:
    """
    Reverse a dict (values as keys, keys as list values).

    Parameters
    ----------
    json_ambi : dict
        A dict having one string value per key.
    ids : set
        sample IDs from the biom table.

    Returns
    -------
    json_ambi_rev : dict
        A dict having the input key(s) in lists
        as values to the input values as keys.
    """
    json_ambi_rev = {}
    for k, v in json_ambi.items():
        if k in ids:
            json_ambi_rev.setdefault(v, []).append(k)
    return json_ambi_rev


def get_most_read_or_features(
        amb_sams: list,
        counts: dict) -> list:
    """
    Get the samples that have the most
    counts (reads of features).

    Parameters
    ----------
    amb_sams : list
        Sample IDs.
    counts : dict
        Count per Sample ID.

    Returns
    -------
    cur_best_samples : list
        Selected samples (with max reads/features).
    """
    cur_counts = dict((amb_sam, counts[amb_sam]) for amb_sam in amb_sams)
    count_max = max(cur_counts.values())
    cur_best_samples = [amb_sam for amb_sam in amb_sams
                        if cur_counts[amb_sam] == count_max]
    return cur_best_samples


def solve_ambiguous_preps(
        redbiom_output: str,
        biom_tab: biom.Table,
        read_counts: dict,
        feat_counts: dict) -> biom.Table:
    """
    Pick the best prep sample per sample to
    solve the ambiguous fetched results.

    Parameters
    ----------
    redbiom_output : str
        The biom table returned by redbiom.
    biom_tab : biom.Table
        Current data retrieved from redbiom, without blooms.
    read_counts : dict
        Number of reads per sample.
    feat_counts : dict
        Number of features per sample.

    Returns
    -------
    biom_tab_no_ambi : biom.Table
        The biom table without ambiguous samples.
    """
    # read the ambiguities
    json_ambi = read_json_ambiguities_file(redbiom_output)
    ids = set(biom_tab.ids(axis='sample'))
    json_ambi_rev = reverse_json_ambi(json_ambi, ids)

    best_samples = []
    if json_ambi_rev:
        print('- Get best samples from ambiguous redbiom results... ', end='')
        # for each sample and its list of preps version
        for sam, amb_sams in json_ambi_rev.items():
            cur_best_samples = get_most_read_or_features(amb_sams, read_counts)
            if len(cur_best_samples) > 1:
                # if >1 samples have same amount of
                # reads, take that with the most features
                cur_best_samples = get_most_read_or_features(cur_best_samples, feat_counts)
            best_samples.append(cur_best_samples[0])
        print('Done -> %s ambiguous samples' % len(best_samples))

    best_samples_noprep = set([best_sample.rsplit('.', 1)[0] for best_sample in best_samples])
    non_ambiguous = [i for i in ids if i.rsplit('.', 1)[0] not in best_samples_noprep]

    ids_to_keep = best_samples + non_ambiguous
    print('- Non-ambiguous + best of ambiguous ', end='')
    biom_tab.filter(ids_to_keep=ids_to_keep, axis='sample', inplace=True)
    biom_tab.remove_empty(axis='observation', inplace=True)
    print('Done -> %s samples' % len(ids_to_keep))
    return biom_tab
