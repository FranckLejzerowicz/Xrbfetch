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
        x if
        str(x.split('.')[1]).isdigit() else
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
    for working_sample_name, subtab in metadata_edit_best.groupby('working_sample_name'):
        curr_subtab = subtab[['read_count', 'feature_count', 'edit_sample_name']]
        max_read_count_sample = curr_subtab.loc[curr_subtab.read_count == max(curr_subtab.read_count)]
        if max_read_count_sample.shape[0] > 1:
            curr_sample_tab = max_read_count_sample.loc[
                max_read_count_sample.feature_count == max(max_read_count_sample.feature_count), 'edit_sample_name']
            if curr_sample_tab.shape[0]:
                best_samples.append(curr_sample_tab.tolist()[0])
            else:
                best_samples.append(curr_sample_tab.item()[0])
        else:
            best_samples.append(max_read_count_sample.edit_sample_name.item())
    print('number of best samples:', len(best_samples))
    metadata_edit_best_out = metadata_edit_best.loc[metadata_edit_best['edit_sample_name'].isin(best_samples)]
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
    metadata_edit_host = metadata_edit.copy()
    for host_subject_id, subtab in metadata_edit_host.groupby('host_subject_id'):
        curr_subtab = subtab[['sample_name', 'read_count', 'feature_count']]
        max_read_count_sample = curr_subtab.loc[curr_subtab.read_count == max(curr_subtab.read_count)]
        if max_read_count_sample.shape[0] > 1:
            curr_sample_tab = max_read_count_sample.loc[
                max_read_count_sample.feature_count == max(max_read_count_sample.feature_count), 'sample_name']
            if curr_sample_tab.shape[0]:
                best_samples.append(curr_sample_tab.tolist()[0])
            else:
                best_samples.append(curr_sample_tab.item()[0])
        else:
            best_samples.append(max_read_count_sample['sample_name'].item())

    print('number of best samples:', len(best_samples))
    metadata_edit_host_out = metadata_edit_host.loc[metadata_edit_host['sample_name'].isin(best_samples)]
    return metadata_edit_host_out


def remove_duplicates(biom_tab_filt: biom.table, metadata_filt: pd.DataFrame) -> tuple:
    """
    Remove duplicates (host and sample preps).

    Parameters
    ----------
    biom_tab_filt : biom.table
        The biom table without ambiguous samples
        and with a min number of reads per sample.
    metadata_filt : pd.DataFrame
        Corresponding metadata table.

    Returns
    -------
    biom_nodup : biom.table
        The biom table without ambiguous samples,
        with a min number of reads per sample, and
        without duplicated sample per host.
    metadata_edit_best : pd.DataFrame
        Corresponding metadata table.
    """
    print(' - Remove duplicates...', end='')
    metadata_edit = add_edit_and_working_sample_name(metadata_filt)
    metadata_edit_host = keep_the_best_host_subject_id_sample(metadata_edit)
    metadata_edit_best = keep_the_best_root_sample_name_sample(metadata_edit_host)
    biom_nodup = biom_tab_filt.filter(
        ids_to_keep=metadata_edit_best['orig_sample_name'].tolist(),
        axis='sample').copy()
    biom_nodup.remove_empty(axis='observation', inplace=True)
    print('Done')
    return biom_nodup, metadata_edit_best

