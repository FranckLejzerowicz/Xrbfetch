# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from Xrbfetch._xrbfetch import xrbfetch
from Xrbfetch import __version__


@click.command()
@click.option(
    "-m", "--m-metadata-file", required=True,
    help="Path to metadata file containing the samples to fetch "
         "and filter (must be in the first column)."
)
@click.option(
    "-o", "--o-metadata-file", required=True,
    help="Path to the output metadata table file."
)
@click.option(
    "-b", "--o-biom-file", required=True,
    help="Path to the output biom table file."
)
@click.option(
    "-r", "--p-redbiom-context", required=False,
    default="Deblur-Illumina-16S-V4-150nt-780653", show_default=True,
    help="Redbiom context for fetching 16S data from Qiita."
)
@click.option(
    "-s", "--p-bloom-sequences", required=False, default=None, show_default=False,
    help="Fasta file containing the sequences known to bloom in fecal samples "
         "(defaults to 'newblooms.all.fasta' file from package's folder 'resources')."
)
@click.option(
    "-f", "--p-reads-filter", default=1500, show_default=True, type=int,
    help="Minimum number of reads per sample."
)
@click.option(
    "--unique/--no-unique", default=True, show_default=True,
    help="Keep a unique sample per host (most read, or most features)."
)
@click.option(
    "--update/--no-update", default=True, show_default=True,
    help="Update the sample names to remove Qiita-prep info."
)
@click.option(
    "--dim/--no-dim", default=True, show_default=True,
    help="Add the number of samples in the final biom file name before "
         "extension (e.g. for '-b out.biom' it becomes 'out_1000s.biom')."
)
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Re-fetch and not use an 'already-fetched' table."
)
@click.option(
    "--simple/--no-simple", default=True, show_default=True,
    help="Perform the steps using one-liners (less checks)."
)
@click.option(
    "--verbose/--no-verbose", default=True, show_default=True,
    help="Show missing, non-fetched samples and duplicates."
)
@click.version_option(__version__, prog_name="Xrbfetch")


def standalone_xrbfetch(
        m_metadata_file,
        o_metadata_file,
        o_biom_file,
        p_redbiom_context,
        p_bloom_sequences,
        p_reads_filter,
        unique, update,
        dim, force,
        simple,
        verbose
):

    xrbfetch(
        m_metadata_file,
        o_metadata_file,
        o_biom_file,
        p_redbiom_context,
        p_bloom_sequences,
        p_reads_filter,
        unique, update,
        dim, force,
        simple,
        verbose
    )


if __name__ == "__main__":
    standalone_xrbfetch()