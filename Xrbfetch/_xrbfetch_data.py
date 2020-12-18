# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import json
import subprocess
import pandas as pd
import pkg_resources

from os.path import abspath, isfile, splitext

RESOURCES = pkg_resources.resource_filename('Xrbfetch', 'resources')


def update_sample_name(update: bool, biom_nodup: biom.table) -> biom.table:
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
        for sam in biom_nodup.ids(axis='sample'):
            if str(sam.split('.')[1]).isdigit():
                updated_sam = '.'.join(sam.split('.')[:2])
            else:
                updated_sam = '.'.join(sam.split('.')[:2])[:-1]
            ids_map[sam] = updated_sam
        print(' - Update sample name to remove prep file info... ', end='')
        biom_nodup.update_ids(id_map=ids_map, axis='sample', inplace=True)
        print('Done')
    return biom_nodup


def filter_reads(
        metadata: pd.DataFrame,
        biom_tab_no_ambi: biom.table,
        reads_filter: int,
        ids_read_counts: dict,
        ids_feat_counts: dict) -> tuple:
    """
    Filter to a minimum number of reads.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for the included samples only.
    biom_tab_no_ambi : biom.table
        The biom table without ambiguous samples.
    reads_filter : int
        Minimum number of reads per sample.
    ids_read_counts : dict
        Number of reads per sample.
    ids_feat_counts : dict
        Number of features per sample.

    Returns
    -------
    biom_tab_filt : biom.table
        The biom table without ambiguous samples
        and with a min number of reads per sample.
    metadata_filt : pd.DataFrame
        Corresponding metadata table.
    """
    print(' - Merge reads and features counts to metadata... ', end='')
    ids = biom_tab_no_ambi.ids(axis = 'sample')
    columns = ['sample_name']
    for col in ['qiita_prep_id',
                'read_count',
                'feature_count',
                'orig_sample_name']:
        if col in metadata.columns:
            metadata.drop(columns=col, inplace=True)
        columns.append(col)

    ids_read_feat_counts_pd = pd.DataFrame([
        [
            '.'.join(ID.split('.')[:-1]),
            ID.split('.')[-1],
            ids_read_counts[ID],
            ids_feat_counts[ID],
            ID
        ] for ID in ids],
        columns = columns
    )
    metadata_no_ambi = metadata.merge(ids_read_feat_counts_pd, on='sample_name', how='right')
    print('Done -> %s samples in merged metadata' % metadata_no_ambi.shape[0])

    # Filter to keep only the samples with min number reads
    print(' - Filter biom for min %s reads per sample... ' % reads_filter, end='')
    metadata_filt = metadata_no_ambi.loc[metadata_no_ambi['read_count'] >= reads_filter].copy()
    biom_tab_filt = biom_tab_no_ambi.filter(
        ids_to_keep = metadata_filt['orig_sample_name'].tolist(),
        axis = 'sample'
    ).copy()
    biom_tab_filt.remove_empty(
        axis = 'observation', inplace=True
    )
    print('Done -> %s samples' % biom_tab_filt.shape[1])
    return biom_tab_filt, metadata_filt


def run_redbiom_fetch(
        metadata: pd.DataFrame,
        m_metadata_file: str,
        p_redbiom_context: str,
        force: bool) -> tuple:
    """
    Fetch the samples using RedBiom.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata for the included samples only.
    m_metadata_file : str
        Path to output metadata for the included samples only.
    p_redbiom_context : str
        Redbiom context for fetching 16S data from Qiita.

    Returns
    -------
    redbiom_output : str
        Path to the biom table returned by redbiom.
    redbiom_samples : str
        Path to the file containing the samples used for fetching.
    """
    redbiom_samples = '%s_redbiom_sams.tmp' % splitext(m_metadata_file)[0]
    redbiom_output = '%s_redbiom.biom' % splitext(m_metadata_file)[0]
    if not force and isfile(redbiom_output):
        print('Using an already generated biom file for this metadata file:\n '
              '-> %s (%s samples) <-' % (redbiom_output, metadata.shape[0]))
    else:
        context = 'Deblur-Illumina-16S-V4-150nt-780653'
        if p_redbiom_context:
            context = p_redbiom_context

        with open(redbiom_samples, 'w') as o:
            for feature_name in metadata.sample_name:
                o.write('%s\n' % feature_name)

        cmd = [
            'redbiom', 'fetch', 'samples',
            '--from', redbiom_samples,
            '--context', context,
            '--output', redbiom_output
        ]
        print('Fetching %s samples from RedBiom... ' % metadata.shape[0], end='')
        print(' '.join(cmd))
        subprocess.call(cmd)
        print('Done')
        subprocess.call(['rm', redbiom_samples])
    return redbiom_output, redbiom_samples


def remove_blooms(biom_tab: biom.table, biom_tab_sams: list,
                  p_bloom_sequences: str) -> tuple:
    """
    Remove the bloom sequences from the fetched samples.

    Parameters
    ----------
    biom_tab : biom.table
        Feature table retrieved from redbiom.
    biom_tab_sams : list
        Samples of the feature table.
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

    bloom_seqs = set([x.strip() for x in open(bloom_sequences_fp).readlines() if x[0] != '>'])

    in_out = {'in': [], 'out': []}
    for feature in biom_tab.ids(axis='observation'):
        if feature in bloom_seqs:
            in_out['in'].append(feature)
        else:
            in_out['out'].append(feature)

    print(' - filter blooms... ', end='')
    biom_tab.filter(
        ids_to_keep=in_out['out'],
        axis='observation',
        inplace=True
    )
    biom_tab.remove_empty(axis='sample', inplace=True)
    biom_tab_removed_ids = list(set(biom_tab_sams) ^ set(biom_tab.ids(axis='sample')))
    print('Done -> %s samples (%s dropouts)' % (biom_tab.shape[1], len(biom_tab_removed_ids)))
    return biom_tab, biom_tab_removed_ids


def solve_ambiguous_preps(
        redbiom_output: str,
        biom_tab: biom.table,
        biom_tab_removed_ids: list) -> tuple:
    """
    Pick the best prep sample per sample to solve the ambiguous fetched results.

    Parameters
    ----------
    redbiom_output : str
        The biom table returned by redbiom.
    biom_tab : biom.table
        Current data retrieved from redbiom, without blooms.
    biom_tab_removed_ids : list
        Samples ids removed as bloom dropouts.

    Returns
    -------
    biom_tab_no_ambi : biom.table
        The biom table without ambiguous samples.
    ids_read_counts : dict
        Number of reads per sample.
    ids_feat_counts : dict
        Number of features per sample.
    """
    ids = biom_tab.ids(axis='sample')
    print(' - Count reads per sample... ', end='')
    read_counts = biom_tab.sum(axis='sample')
    print('Done')

    print(' - Count features per sample... ', end='')
    feat_counts = biom_tab.nonzero_counts(axis='sample', binary=True)
    print('Done')

    ids_read_counts = dict([ids[x], read_counts[x]] for x in range(len(ids)))
    ids_feat_counts = dict([ids[x], feat_counts[x]] for x in range(len(ids)))

    json_ambi_fp = '%s.ambiguities' % redbiom_output
    if isfile(json_ambi_fp):
        with open(json_ambi_fp) as f:
            json_ambi = json.load(f)

    ambi_selection = []
    if json_ambi:
        print(' - Get best samples from ambiguous redbiom results... ', end='')
        for sam, amb_sams in json_ambi.items():
            if isinstance(amb_sams, list):
                dict_amb_sams_read_counts = {}
                dict_amb_sams_feat_counts = {}
                amb_sams_fixes = ['%s.%s' % (x.split('_')[1], x.split('_')[0]) for x in amb_sams]
                amb_sams_fixes_kept = [x for x in amb_sams_fixes if x not in biom_tab_removed_ids]
                if len(amb_sams_fixes_kept):
                    for amb_sam_id in amb_sams_fixes_kept:
                        dict_amb_sams_read_counts[amb_sam_id] = ids_read_counts[amb_sam_id]
                        dict_amb_sams_feat_counts[amb_sam_id] = ids_feat_counts[amb_sam_id]

                    read_count_max = max([x for x in dict_amb_sams_read_counts.values()])
                    read_count_max_samples = [x for x, n in dict_amb_sams_read_counts.items() if n == read_count_max]
                    if len(read_count_max_samples) > 1:
                        feat_count_max = max([dict_amb_sams_feat_counts[feat] for feat in read_count_max_samples])
                        feat_count_max_samples = [x for x, n in dict_amb_sams_feat_counts.items() if
                                                  n == feat_count_max and x in read_count_max_samples]
                        ambi_selection.append(feat_count_max_samples[0])
                    else:
                        ambi_selection.append(read_count_max_samples[0])
        print('Done -> %s ambiguous samples' % len(ambi_selection))

    ambi_selection_noPrep = set(['.'.join(x.split('.')[:-1]) for x in ambi_selection])
    non_ambi = [x for x in ids if '.'.join(x.split('.')[:-1]) not in ambi_selection_noPrep]
    ids_to_keep = ambi_selection + non_ambi
    print(' - Non-ambiguous + best of ambiguous ', end='')
    biom_tab.filter(ids_to_keep=ids_to_keep, axis='sample', inplace=True)
    biom_tab.remove_empty(axis='observation', inplace=True)
    print('= %s samples' % len(ids_to_keep))
    return biom_tab, ids_read_counts, ids_feat_counts

