# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd


def add_edit_and_working_sample_name(metadata_filt: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters
    ----------
    metadata_filt : pd.DataFrame
        Metadata table for the non-ambiguous samples
        and with a min number of reads per sample.

    Returns
    -------
    metadata_edit : pd.DataFrame
        Metadata table for the non-ambiguous samples
        and with a min number of reads per sample and
        new, renamed samples names (for checking
        duplicates).
    """
    metadata_edit = metadata_filt.copy()
    metadata_edit['edit_sample_name'] = [
        x if str(x.split('.')[1]).isdigit() else
        '%s.%s.%s' % (
            x.split('.')[0],
            x.split('.')[1][:-1],
            x.split('.')[2]
        ) for x in metadata_edit['orig_sample_name']
    ]
    metadata_edit['working_sample_name'] = [
        str(int(x.split('.')[1])) if
        str(x.split('.')[1]).isdigit() else
        str(int(x.split('.')[1][:-1])
            ) for x in metadata_edit['orig_sample_name']
    ]
    return metadata_edit


def keep_the_best_root_sample_name_sample(metadata_edit_host: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters
    ----------
    metadata_edit_host : pd.DataFrame
        Metadata table for the non-ambiguous samples
        and with a min number of reads per sample and
        with only one sample per host.

    Returns
    -------
    metadata_edit_best : pd.DataFrame
        Metadata table for the non-ambiguous samples
        and with a min number of reads per sample and
        with only one sample per host and per sample
    """
    metadata_edit_best = metadata_edit_host.copy()
    best_samples = []
    print('   * Keep the best working_sample_name per sample... ', end='')
    for working_sample_name, subtab in metadata_edit_best.groupby('working_sample_name'):
        curr_subtab = subtab[['read_count', 'feature_count', 'edit_sample_name']]
        max_read_count_sample = curr_subtab.loc[curr_subtab.read_count == max(curr_subtab.read_count)]
        if max_read_count_sample.shape[0] > 1:
            max_feat_count_sample = max_read_count_sample.loc[
                max_read_count_sample.feature_count == max(max_read_count_sample.feature_count)]
            max_feat_count_sample.sort_values('feature_count', ascending=False, inplace=True)
            best_sample = str(next(iter(max_feat_count_sample.edit_sample_name), 'no match'))
        else:
            best_sample = str(next(iter(max_read_count_sample.edit_sample_name), 'no match'))
        best_samples.append(best_sample)

    metadata_edit_best_out = metadata_edit_best.loc[metadata_edit_best['edit_sample_name'].isin(best_samples)]
    print('Done -> %s samples' % metadata_edit_best_out.shape[0])
    return metadata_edit_best_out


def keep_the_best_host_subject_id_sample(metadata_edit: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters
    ----------
    metadata_edit : pd.DataFrame
        Metadata table for the non-ambiguous samples
        and with a min number of reads per sample and
        new, renamed samples names (for checking
        duplicates).

    Returns
    -------
    metadata_edit_host : pd.DataFrame
        Metadata table for the non-ambiguous samples
        and with a min number of reads per sample and
        with only one sample per host.
    """
    best_samples = []
    print('   * Keep the best host_subject_id per sample... ', end='')
    metadata_edit_host = metadata_edit.copy()
    for host_subject_id, subtab in metadata_edit_host.groupby('host_subject_id'):
        curr_subtab = subtab[['sample_name', 'read_count', 'feature_count']]
        max_read_count_sample = curr_subtab.loc[curr_subtab.read_count == max(curr_subtab.read_count)]
        if max_read_count_sample.shape[0] > 1:
            max_feat_count_sample = max_read_count_sample.loc[
                max_read_count_sample.feature_count == max(max_read_count_sample.feature_count)]
            max_feat_count_sample.sort_values('feature_count', ascending=False, inplace=True)
            best_sample = str(next(iter(max_feat_count_sample.sample_name), 'no match'))
        else:
            best_sample = str(next(iter(max_read_count_sample.sample_name), 'no match'))
        best_samples.append(best_sample)

    metadata_edit_host_out = metadata_edit_host.loc[metadata_edit_host['sample_name'].isin(best_samples)]
    print('Done -> %s samples' % metadata_edit_host_out.shape[0])
    return metadata_edit_host_out


def remove_duplicates(biom_tab_filt: biom.table,
                      metadata_filt: pd.DataFrame,
                      unique: bool) -> tuple:
    """
    Remove duplicates (host and sample preps).

    Parameters
    ----------
    biom_tab_filt : biom.table
        The biom table without ambiguous samples
        and with a min number of reads per sample.
    metadata_filt : pd.DataFrame
        Corresponding metadata table.
    unique : bool
        Whether to keep a unique sample per host or not.

    Returns
    -------
    biom_nodup : biom.table
        The biom table without ambiguous samples,
        with a min number of reads per sample, and
        without duplicated sample per host.
    metadata_edit_best : pd.DataFrame
        Corresponding metadata table.
    """
    if unique:
        print(' - Remove duplicates:')
        metadata_edit = add_edit_and_working_sample_name(metadata_filt)
        metadata_edit_host = keep_the_best_host_subject_id_sample(metadata_edit)
        metadata_edit_best = keep_the_best_root_sample_name_sample(metadata_edit_host)
        biom_nodup = biom_tab_filt.filter(
            ids_to_keep = metadata_edit_best['orig_sample_name'].tolist(),
            axis = 'sample').copy()
        biom_nodup.remove_empty(axis='observation', inplace=True)
        return biom_nodup, metadata_edit_best
    else:
        return biom_tab_filt, metadata_filt


def check_replicates_amount(biom_tab_sams: list) -> None:
    """
    Parameters
    ----------
    biom_tab_sams : list
        Samples of the feature table.
    """
    dups = {}
    for sam in biom_tab_sams:
        dups.setdefault('.'.join(sam.split('.')[:-1]), []).append(sam)

    dups_to_print = {}
    for sam, sams in dups.items():
        if len(sams) > 1:
            dups_to_print.setdefault(len(sams), []).append(sam)

    if dups_to_print:
        print(' * Sample ID duplication (without prep number):')
        for n, sams in sorted(dups_to_print.items(), key=lambda x: len(x[1])):
            print('     --> %s samples have %s replicates' % (len(sams), n))


def check_fetched_samples(metadata: pd.DataFrame, biom_tab_sams: list) -> None:
    """
    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata table.
    biom_tab_sams : list
        Samples of the feature table.
    """

    biom_sams = ['.'.join(sam.split('.')[:-1]) for sam in biom_tab_sams]
    common_sams = set(metadata.sample_name) & set(biom_sams)
    only_meta = set(metadata.sample_name) ^ common_sams
    if len(only_meta):
        print(' [Warning] %s samples not fecthed in any prep:' % len(only_meta))
        print(', '.join(only_meta))