# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import biom
from Xrbfetch.io import write_summary, delete_files


def check_fetched_samples(metadata_sams: list, biom_tab: biom.Table) -> None:
    """
    Parameters
    ----------
    metadata_sams : list
        Samples of the metadata table.
    biom_tab : biom.Table
        Feature table.
    """
    # get the sample names without the prep info number
    sams = set(['.'.join(biom_sam.split('.')[:-1]) for biom_sam in biom_tab.ids(axis='sample')])
    # get sample only present in metadata, i.e., not fetched
    only_meta = set(metadata_sams).difference(sams)
    if len(only_meta):
        print(' [Warning] %s samples not fetched in any prep:' % len(only_meta))
        print(', '.join(only_meta))


def check_replicates_amount(biom_tab: biom.Table) -> None:
    """
    Parameters
    ----------
    biom_tab : biom.Table
        Biom table.
    """
    dups = {}
    for sam in biom_tab.ids(axis='sample'):
        dups.setdefault('.'.join(sam.split('.')[:-1]), []).append(sam)

    dups_to_print = {}
    for sam, sams in dups.items():
        if len(sams) > 1:
            dups_to_print.setdefault(len(sams), []).append(sam)

    if dups_to_print:
        print(' * Sample ID duplication (without prep number):')
        for n, sams in sorted(dups_to_print.items(), key=lambda x: len(x[1])):
            print('     --> %s samples have %s replicates' % (len(sams), n))


def potential_stop(
        biom_tab: biom.Table,
        o_summary_file: str,
        summary: list,
        redbiom_output: str,
        redbiom_samples: str,
        do_stop: bool = False):

    if not biom_tab.shape[0] or do_stop:
        # write the metadata and the biom table outputs
        write_summary(o_summary_file, summary)
        # delete intermediate files
        delete_files(redbiom_output, redbiom_samples)
        sys.exit(0)


