# Xrbfetch

Fetch microbiome tables from Qiita using RedBiom based on a metadata file and filter:
- blooms,
- read-sum based
- duplicate hosts and preps

## Installation

For the first install:
```
pip install git+https://github.com/FranckLejzerowicz/Xrbfetch.git
```

and then if there are updates...
```
pip install --upgrade git+https://github.com/FranckLejzerowicz/Xrbfetch.git
```

*_Note that python and pip should be python3_

## Input

`-m` option **[required]**: Must be an existing path to a metadata file containing the samples 
to fetch and filter (the samples must be listed in the first column). 

## Outputs

`-o` option **[required]**: Path to the output metadata table file. 
Will contain the same info as the input but only for the samples that pass filters. 
Only for samples named according to the AGP standard (i.e. `<study-id>.<9-digits-sample-id>`)

`-b` option **[required]**: Path to the output biom table file. Will contain only the samples that pass filters. 

## Example

- Running:

    ```
    Xrbfetch \
        -m ./Xrbfetch/tests/metadata/metadata.tsv \
        -o ./Xrbfetch/tests/output/metadata_out.tsv \
        -b ./Xrbfetch/tests/output/featuredata.biom
    ```
    Would return `metadata_out_<number>s.tsv` and `featuredata_<number>s.tsv`

- Running:

    ```
    Xrbfetch \
        -m ./Xrbfetch/tests/metadata/metadata.tsv \
        -o ./Xrbfetch/tests/output/metadata_out.tsv \
        -b ./Xrbfetch/tests/output/featuredata.biom
        --no-dim
    ```
    Would return `metadata_out.tsv` and `featuredata.tsv`


### Optional arguments

```
  -m, --m-metadata-file TEXT    Path to metadata file containing the samples
                                to fetch and filter (must be in the first
                                column).  [required]
  -o, --o-metadata-file TEXT    Path to the output metadata table file.
                                [required]
  -b, --o-biom-file TEXT        Path to the output biom table file.
                                [required]
  -r, --p-redbiom-context TEXT  Redbiom context for fetching 16S data from
                                Qiita.  [default: Deblur-
                                Illumina-16S-V4-150nt-780653]
  -s, --p-bloom-sequences TEXT  Fasta file containing the sequences known to
                                bloom in fecal samples (defaults to
                                'newblooms.all.fasta' file from package's
                                folder 'resources').
  -f, --p-reads-filter INTEGER  Minimum number of reads per sample.  [default:
                                1500]
  --unique / --no-unique        Keep a unique sample per host (most read, or
                                most features).  [default: True]
  --update / --no-update        Update the sample names to remove Qiita-prep
                                info.  [default: True]
  --dim / --no-dim              Add the number of samples in the final biom
                                file name before extension (e.g. for '-b
                                out.biom' it becomes 'out_1000s.biom').
                                [default: True]
  --verbose / --no-verbose      Show missing, non-fetched samples and
                                duplicates.  [default: True]
  --version                     Show the version and exit.
  --help                        Show this message and exit.
```



### Bug Reports

contact `flejzerowicz@health.ucsd.edu`