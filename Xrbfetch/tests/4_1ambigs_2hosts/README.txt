This test data is simple (only four sample from AGP):
10317.000001778
10317.000002860
10317.000020106
10317.000027859

But it is useful:

One sample has few reads: test the option `--p-reads-filter` (and only one prep):
  - 10317.000001778.57016: 33 reads

One sample other has two ambiguous preps: test filtering to keep the most-reads sample:
  - 10317.000002860.57016: 15877 reads
  - 10317.000002860.58862: 31707 reads

The two others samples are assigned to the same host: test filtering to keep best sample per host:
  - 10317.000001778 -> host `A`
  - 10317.000002860 -> host `B`
  - 10317.000027859 -> host `C`
  - 10317.000020106 -> host `C`


The command to test on this dataset is:

# run the simple workflow
Xrbfetch \
  -m Xrbfetch/tests/4_1ambigs_2hosts/4_1ambigs_2hosts.tsv \
  -o Xrbfetch/tests/4_1ambigs_2hosts/4_1ambigs_2hosts_out.tsv \
  -b Xrbfetch/tests/4_1ambigs_2hosts/4_1ambigs_2hosts_out.biom \
  --simple

# run the more complex workflow (which includes bloom filtering prior to most-reads ambiguity filtering)
Xrbfetch \
  -m Xrbfetch/tests/4_1ambigs_2hosts/4_1ambigs_2hosts.tsv \
  -o Xrbfetch/tests/4_1ambigs_2hosts/4_1ambigs_2hosts_out.tsv \
  -b Xrbfetch/tests/4_1ambigs_2hosts/4_1ambigs_2hosts_out.biom


