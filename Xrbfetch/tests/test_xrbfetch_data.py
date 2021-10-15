# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest
import pkg_resources
import pandas as pd
import numpy as np
from biom.table import Table

from pandas.testing import assert_frame_equal

from Xrbfetch.data import (
    merge_read_features_counts,
    run_redbiom_fetch,
    filter_reads,
    update_sample_name,
    remove_blooms,
    get_reads_features_counts,
    reverse_json_ambi,
    get_most_read_or_features,
    solve_ambiguous_preps
)

ROOT = pkg_resources.resource_filename('Xrbfetch', 'tests')


class TestData(unittest.TestCase):

    def setUp(self):

        self.metadata = pd.DataFrame({'sample_name': ['10317.000001778']})
        self.m_metadata_file = '%s/samples/test.tsv' % ROOT
        self.metadata.to_csv(self.m_metadata_file, index=False, sep='\t')
        self.redbiom_output = '%s/samples/test_redbiom.biom' % ROOT
        self.redbiom_output_amb = '%s/samples/test_redbiom.biom.ambiguities' % ROOT
        self.bloom_sequences = '%s/samples/seqs.fasta' % ROOT
        self.context = 'Deblur-Illumina-16S-V4-150nt-780653'

        self.metadata_counts = pd.DataFrame({
            'sample_name': ['a.1', 'b.2', 'c.3'],
            'orig_sample_name': ['a.1.x', 'b.2.y', 'c.3.z'],
            'read_count': [999, 1001, 1001]})
        self.metadata_counts_merge = pd.DataFrame({
            'sample_name': ['a.1', 'b.2', 'c.3'],
            'qiita_prep_id': ['x', 'y', 'z'],
            'read_count': [999, 1001, 1001],
            'feature_count': [3, 2, 2],
            'orig_sample_name': ['a.1.x', 'b.2.y', 'c.3.z']
        })
        self.read_counts = {'a.1.x': 999, 'b.2.y': 1001, 'c.3.z': 1001}
        self.feat_counts = {'a.1.x': 3, 'b.2.y': 2, 'c.3.z': 2}
        self.metadata_counts_out = self.metadata_counts.iloc[1:, :].copy()
        self.biom_tab_no_ambi = Table(
            np.array([[1, 0, 0], [1, 1, 1], [1, 1, 1]]),
            ['AA', 'CC', 'GG'], ['a.1.x', 'b.2.y', 'c.3.z'])
        self.biom_tab_no_ambi_out = Table(
            np.array([[1, 1], [1, 1]]),
            ['CC', 'GG'], ['b.2.y', 'c.3.z'])
        self.biom_tab_no_ambi_out_updated = Table(
            np.array([[1, 1], [1, 1]]),
            ['CC', 'GG'], ['b.2', 'c.3'])

        self.biom_with_bloom_1nt = Table(
            np.array([[1, 1, 1], [1, 1, 0]]),
            ['CC', 'TT'], ['a', 'b', 'c'])
        self.biom_with_bloom_2nt = Table(
            np.array([[1, 1, 1], [1, 1, 0]]),
            ['C', 'T'], ['a', 'b', 'c'])
        self.biom_with_nobloom = Table(
            np.array([[1, 1]]),
            ['TT'], ['a', 'b'])
        self.bloom_sammple = ['c']

        self.reads_filter = 1000

    def test_run_redbiom_fetch(self):
        redbiom_output = run_redbiom_fetch(self.metadata, self.m_metadata_file, self.context, False)
        self.assertEqual(self.redbiom_output, redbiom_output)
        self.assertTrue(os.path.isfile(redbiom_output))
        redbiom_output = run_redbiom_fetch(self.metadata, self.m_metadata_file, self.context, True)
        self.assertEqual(self.redbiom_output, redbiom_output)
        self.assertTrue(os.path.isfile(redbiom_output))
        redbiom_output = run_redbiom_fetch(self.metadata, self.m_metadata_file, self.context, False)
        self.assertEqual(self.redbiom_output, redbiom_output)
        self.assertTrue(os.path.isfile(redbiom_output))

    def test_filter_reads(self):
        biom_tab_filt, metadata_filt = filter_reads(
            self.metadata_counts, self.biom_tab_no_ambi, self.reads_filter)
        self.assertEqual(biom_tab_filt, self.biom_tab_no_ambi_out)
        assert_frame_equal(metadata_filt, self.metadata_counts_out)

    def test_merge_read_features_counts(self):
        metadata_counts = merge_read_features_counts(
            self.metadata_counts, self.biom_tab_no_ambi,
            self.read_counts, self.feat_counts
        )
        assert_frame_equal(self.metadata_counts_merge, metadata_counts)

    def test_update_sample_name(self):
        biom_tab_no_ambi_out_updated = update_sample_name(False, self.biom_tab_no_ambi_out)
        self.assertEqual(biom_tab_no_ambi_out_updated, self.biom_tab_no_ambi_out)
        biom_tab_no_ambi_out_updated = update_sample_name(True, self.biom_tab_no_ambi_out)
        self.assertEqual(biom_tab_no_ambi_out_updated, self.biom_tab_no_ambi_out_updated)

    def test_remove_blooms(self):
        biom_tab = remove_blooms(
            2, self.biom_with_nobloom, self.bloom_sequences)
        self.assertEqual(self.biom_with_nobloom, biom_tab)
        biom_tab = remove_blooms(
            2, self.biom_with_bloom_2nt, self.bloom_sequences)
        self.assertEqual(self.biom_with_nobloom, biom_tab)
        biom_tab = remove_blooms(
            1, self.biom_with_bloom_1nt, self.bloom_sequences)
        self.assertEqual(self.biom_with_nobloom, biom_tab)

    def test_get_reads_features_counts(self):
        # because the counts are dummy counts set to 1; reads counts = features counts...
        read_counts, feat_counts = get_reads_features_counts(self.biom_tab_no_ambi)
        self.assertEqual(read_counts, {'a.1.x': 3, 'b.2.y': 2, 'c.3.z': 2})
        self.assertEqual(feat_counts, {'a.1.x': 3, 'b.2.y': 2, 'c.3.z': 2})
        read_counts, feat_counts = get_reads_features_counts(self.biom_tab_no_ambi_out)
        self.assertEqual(read_counts, {'b.2.y': 2, 'c.3.z': 2})
        self.assertEqual(feat_counts, {'b.2.y': 2, 'c.3.z': 2})
        read_counts, feat_counts = get_reads_features_counts(self.biom_tab_no_ambi_out_updated)
        self.assertEqual(read_counts, {'b.2': 2, 'c.3': 2})
        self.assertEqual(feat_counts, {'b.2': 2, 'c.3': 2})
        read_counts, feat_counts = get_reads_features_counts(self.biom_with_bloom)
        self.assertEqual(read_counts, {'a': 2, 'b': 2, 'c': 1})
        self.assertEqual(feat_counts, {'a': 2, 'b': 2, 'c': 1})
        read_counts, feat_counts = get_reads_features_counts(self.biom_with_nobloom)
        self.assertEqual(read_counts, {'a': 1, 'b': 1})
        self.assertEqual(feat_counts, {'a': 1, 'b': 1})

    def test_reverse_json_ambi(self):

        rev_json = reverse_json_ambi({'a': '1', 'b': '1', 'c': '2'}, set())
        self.assertEqual(rev_json, {'1': ['a', 'b'], '2': ['c']})
        rev_json = reverse_json_ambi({'a': '1', 'b': '2', 'c': '3'}, set())
        self.assertEqual(rev_json, {'1': ['a'], '2': ['b'], '3': ['c']})

        rev_json = reverse_json_ambi({'a': '1', 'b': '1', 'c': '2'}, {'a'})
        self.assertEqual(rev_json, {'1': ['b'], '2': ['c']})
        rev_json = reverse_json_ambi({'a': '1', 'b': '2', 'c': '3'}, {'b'})
        self.assertEqual(rev_json, {'1': ['a'], '3': ['c']})
        rev_json = reverse_json_ambi({'a': '1', 'b': '2', 'c': '3'}, {'a', 'b'})
        self.assertEqual(rev_json, {'3': ['c']})

    def test_get_most_read_or_features(self):
        sams = get_most_read_or_features(['a.1', 'a.2'],  {'a.1': 1, 'a.2': 2})
        self.assertEqual(sams, ['a.2'])
        sams = get_most_read_or_features(['a.1', 'a.2', 'a.3'], {'a.1': 1, 'a.2': 2, 'a.3': 2})
        self.assertEqual(sams, ['a.2', 'a.3'])

    def tearDown(self):
        if os.path.isfile(self.redbiom_output):
            os.remove(self.redbiom_output)
        if os.path.isfile(self.redbiom_output_amb):
            os.remove(self.redbiom_output_amb)
        if os.path.isfile(self.m_metadata_file):
            os.remove(self.m_metadata_file)


class TestData(unittest.TestCase):

    def setUp(self):
        self.ambiguities_json = '%s/ambiguities/json' % ROOT
        self.read_counts_diff = {"10317.000001778.57016": 1,
                                 "10317.000002860.57016": 2,
                                 "10317.000002860.58862": 1}
        self.read_counts_equal = {"10317.000001778.57016": 1,
                                  "10317.000002860.57016": 1,
                                  "10317.000002860.58862": 1}
        self.feat_counts = {"10317.000001778.57016": 1,
                            "10317.000002860.57016": 1,
                            "10317.000002860.58862": 2}
        self.biom = Table(np.array([[1, 1, 1], [1, 1, 1]]), ['sp1', 'sp2'],
                          ["10317.000001778.57016",
                           "10317.000002860.57016",
                           "10317.000002860.58862"])
        self.biom_diff = Table(np.array([[1, 1], [1, 1]]), ['sp1', 'sp2'],
                               ["10317.000001778.57016", "10317.000002860.57016"])
        self.biom_equal = Table(np.array([[1, 1], [1, 1]]), ['sp1', 'sp2'],
                                ["10317.000001778.57016", "10317.000002860.58862"])

    def test_solve_ambiguous_preps(self):
        biom_noprep = solve_ambiguous_preps(
            self.ambiguities_json, self.biom.copy(),
            self.read_counts_diff, self.feat_counts)
        self.assertEqual(biom_noprep, self.biom_diff)

        biom_noprep = solve_ambiguous_preps(
            self.ambiguities_json, self.biom.copy(),
            self.read_counts_equal, self.feat_counts)
        self.assertEqual(biom_noprep, self.biom_equal)


if __name__ == '__main__':
    unittest.main()
