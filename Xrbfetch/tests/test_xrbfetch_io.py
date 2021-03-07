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
from Xrbfetch.io import (
    read_meta_pd,
    run_fetch,
    delete_files,
    write_summary,
    make_samples_list_tmp,
    get_outputs,
    write_outputs
)

ROOT = pkg_resources.resource_filename('Xrbfetch', 'tests')


class TestMd(unittest.TestCase):

    def setUp(self):
        self.md_res = pd.DataFrame({
            'sample_name': ['A', 'B', 'C'],
            'col1': ['A1', 'B1', 'C1'],
            'col2': ['A2', 'B2', 'C2'],
            'col3': ['A3', 'B3', 'C3']
        })

    def test_read_meta_pd(self):
        md_tabs_fp = '%s/metadata/test_md/md_tabs.tsv' % ROOT
        md_tabs = read_meta_pd(md_tabs_fp)
        assert_frame_equal(md_tabs, self.md_res)
        md_upper_fp = '%s/metadata/test_md/md_upper.tsv' % ROOT
        md_upper = read_meta_pd(md_upper_fp)
        assert_frame_equal(md_upper, self.md_res)
        md_missing_fp = '%s/metadata/test_md/md_missing.tsv' % ROOT
        md_missing = read_meta_pd(md_missing_fp)
        assert_frame_equal(md_missing, self.md_res)


class TestRun(unittest.TestCase):

    def setUp(self):
        self.redbiom_samples = '%s/samples/redbiom_samples.txt' % ROOT
        self.redbiom_output = '%s/samples/redbiom_samples.biom' % ROOT
        self.redbiom_output_amb = '%s/samples/redbiom_samples.biom.ambiguities' % ROOT
        self.context = 'Deblur-Illumina-16S-V4-150nt-780653'

    def test_run_fetch(self):
        with open(self.redbiom_samples, 'w') as o:
            o.write('sample_name\n10317.000001778\n')
        run_fetch(self.redbiom_samples, self.context, self.redbiom_output)
        self.assertTrue(os.path.isfile(self.redbiom_output))
        self.assertFalse(os.path.isfile(self.redbiom_samples))

    def tearDown(self):
        os.remove(self.redbiom_output)
        if os.path.isfile(self.redbiom_output_amb):
            os.remove(self.redbiom_output_amb)


class TestIO(unittest.TestCase):

    def setUp(self):

        self.redbiom_samples = '%s/samples/redbiom_samples.tsv' % ROOT
        self.redbiom_output = '%s/samples/redbiom_samples_redbiom.biom' % ROOT
        self.redbiom_output_amb = '%s/samples/redbiom_samples_redbiom.biom.ambiguities' % ROOT
        open(self.redbiom_output, 'w').close()
        open(self.redbiom_output_amb, 'w').close()

        self.o_summary_file = '%s/summary/summary.txt' % ROOT
        self.summary = ['steps\tsamples\n', 'a\t1\n', 'b\t2\n']

        self.redbiom_samples_temp = '%s/samples/redbiom_samples_redbiom_sams.tmp' % ROOT

        self.metadata = pd.DataFrame({'sample_name': ['a', 'b']})

        self.o_metadata_file = '%s/samples/redbiom_samples.tsv' % ROOT
        self.o_metadata_file_dim = '%s/samples/redbiom_samples_2s.tsv' % ROOT
        self.o_biom_file = '%s/samples/redbiom_samples_redbiom.biom' % ROOT
        self.o_biom_file_dim = '%s/samples/redbiom_samples_redbiom_2s.biom' % ROOT

        self.biom_table = Table(np.array([[1, 2], [1, 2]]), ['sp1', 'sp2'], ['s1', 's2'])

    def test_write_summary(self):
        write_summary(self.o_summary_file, [['a', 1], ['b', 2]])
        summary = open(self.o_summary_file).readlines()
        self.assertEqual(summary, self.summary)

    def test_delete_files(self):
        delete_files(self.redbiom_samples)
        self.assertFalse(os.path.exists(self.redbiom_output))
        self.assertFalse(os.path.exists(self.redbiom_output_amb))

    def test_make_samples_list_tmp(self):
        redbiom_samples = make_samples_list_tmp(self.redbiom_samples, self.metadata)
        redbiom_samples_temp = open(self.redbiom_samples_temp).readlines()
        self.assertEqual(redbiom_samples_temp, ['a\n', 'b\n'])
        self.assertEqual(redbiom_samples, self.redbiom_samples_temp)

    def test_get_outputs(self):
        o_metadata_file, o_biom_file = get_outputs(
            self.metadata, self.redbiom_samples, self.redbiom_output)
        self.assertEqual(o_metadata_file, self.o_metadata_file_dim)
        self.assertEqual(o_biom_file, self.o_biom_file_dim)

    def test_write_outputs(self):
        write_outputs(self.o_biom_file, self.o_metadata_file, self.biom_table, self.metadata, False)
        self.assertTrue(os.path.isfile(self.o_biom_file))
        self.assertTrue(os.path.isfile(self.o_metadata_file))
        write_outputs(self.o_biom_file, self.o_metadata_file, self.biom_table, self.metadata, True)
        self.assertTrue(os.path.isfile(self.o_biom_file_dim))
        self.assertTrue(os.path.isfile(self.o_metadata_file_dim))

    def tearDown(self):
        if os.path.isfile(self.o_biom_file):
            os.remove(self.o_biom_file)
        if os.path.isfile(self.o_biom_file_dim):
            os.remove(self.o_biom_file_dim)
        if os.path.isfile(self.o_metadata_file):
            os.remove(self.o_metadata_file)
        if os.path.isfile(self.o_metadata_file_dim):
            os.remove(self.o_metadata_file_dim)
        if os.path.isfile(self.redbiom_samples_temp):
            os.remove(self.redbiom_samples_temp)
        if os.path.isfile(self.redbiom_samples):
            os.remove(self.redbiom_samples)
        if os.path.isfile(self.redbiom_output):
            os.remove(self.redbiom_output)
        if os.path.isfile(self.redbiom_output_amb):
            os.remove(self.redbiom_output_amb)


if __name__ == '__main__':
    unittest.main()
