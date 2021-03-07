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
from biom.table import Table

from pandas.testing import assert_frame_equal

from Xrbfetch.duplicates import (
    add_edit_and_working_sample_name,
    keep_the_best_host_subject_id_sample
)

ROOT = pkg_resources.resource_filename('Xrbfetch', 'tests')


class TestData(unittest.TestCase):
    def setUp(self):
        self.metadata_digit = pd.DataFrame({'orig_sample_name': ['1.01.1']})
        self.metadata_nodigit = pd.DataFrame({'orig_sample_name': ['1.01A.1']})

        self.metadata_hosts_reads = pd.DataFrame({'sample_name': ['1', '2', '3'],
                                                  'host_subject_id': ['a', 'a', 'b'],
                                                  'read_count': ['1', '2', '1'],
                                                  'feature_count': ['1', '1', '1']})
        self.metadata_hosts_reads_out = self.metadata_hosts_reads.drop(index=[0])
        self.metadata_hosts_feats = pd.DataFrame({'sample_name': ['1', '2', '3'],
                                                  'host_subject_id': ['a', 'a', 'b'],
                                                  'read_count': ['1', '1', '1'],
                                                  'feature_count': ['2', '1', '1']})
        self.metadata_hosts_feats_out = self.metadata_hosts_feats.drop(index=[1])

    def test_add_edit_and_working_sample_name(self):
        metadata_digit = add_edit_and_working_sample_name(self.metadata_digit)
        assert_frame_equal(metadata_digit, pd.DataFrame({'orig_sample_name': ['1.01.1'],
                                                         'edit_sample_name': ['1.01.1'],
                                                         'working_sample_name': ['1']}))
        metadata_nodigit = add_edit_and_working_sample_name(self.metadata_nodigit)
        assert_frame_equal(metadata_nodigit, pd.DataFrame({'orig_sample_name': ['1.01A.1'],
                                                           'edit_sample_name': ['1.01.1'],
                                                           'working_sample_name': ['1']}))

    def test_keep_the_best_host_subject_id_sample(self):
        metadata_hosts_reads = keep_the_best_host_subject_id_sample(self.metadata_hosts_reads)
        metadata_hosts_feats = keep_the_best_host_subject_id_sample(self.metadata_hosts_feats)
        assert_frame_equal(metadata_hosts_reads, self.metadata_hosts_reads)
        assert_frame_equal(metadata_hosts_feats, self.metadata_hosts_feats)


if __name__ == '__main__':
    unittest.main()
