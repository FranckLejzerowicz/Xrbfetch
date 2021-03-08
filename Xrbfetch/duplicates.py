# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd


def add_edit_and_working_sample_name(
        metadata_filt: pd.DataFrame) -> pd.DataFrame:
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

    edit_sample_names = []
    working_sample_names = []
    for sam in metadata_edit['orig_sample_name']:
        sam_agp = sam.split('.')[1]
        if sam_agp.isdigit():
            edit_sample_names.append(sam)
            working_sample_names.append(str(int(sam_agp)))
        else:
            sam_edit = '%s.%s.%s' % (sam.split('.')[0],
                                     sam_agp[:-1],
                                     sam.split('.')[-1])
            edit_sample_names.append(sam_edit)
            working_sample_names.append(str(int(sam_agp[:-1])))

    metadata_edit['edit_sample_name'] = edit_sample_names
    metadata_edit['working_sample_name'] = working_sample_names
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
    print(' - Keep the best working_sample_name per sample... ', end='')
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


def keep_the_best_host_subject_id_sample(
        metadata_edit: pd.DataFrame) -> pd.DataFrame:
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
    metadata_edit_host = metadata_edit.copy()
    if max(metadata_edit_host.host_subject_id.value_counts()) == 1:
        return metadata_edit_host
    print('- Keep the best sample per host_subject_id... ', end='')
    metadata_edit_host.sort_values(
        ['read_count', 'feature_count'],
        ascending=False,
        inplace=True)
    metadata_edit_host = metadata_edit_host[
        ~metadata_edit_host.host_subject_id.duplicated()]
    print('Done -> %s samples' % metadata_edit_host.shape[0])
    return metadata_edit_host


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
    metadata_filt = add_edit_and_working_sample_name(metadata_filt)
    if unique:
        metadata_filt = keep_the_best_host_subject_id_sample(metadata_filt)
    metadata_edit = keep_the_best_root_sample_name_sample(metadata_filt)
    biom_nodup = biom_tab_filt.filter(
        ids_to_keep=metadata_edit['orig_sample_name'].tolist(),
        axis='sample').copy()
    biom_nodup.remove_empty(axis='observation', inplace=True)
    return biom_nodup, metadata_edit
