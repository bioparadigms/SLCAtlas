# Running the SLC discovery pipeline

Data is stored partly in MySQL databases and partly in files. MySQL tables can be "ephemeral", which means each time they are regenerated, the old table is cleared and its data replaced by new data. Tables can also be "persistent", in which case they have a history mechanism. Changes in fields are recorded each time a modification (change, deletion) is performed. More about this later.

The commands to run the pipeline are described in the sections below. Python was used under Conda, and environment descriptions can be found in `environment.yml` (full) and `environment-from-history.yml` (only explicitly installed packages). The environment might be referred to in scripts as `pyweb`. Note that some scripts have been run in a HPC environment under SLURM and thus will require modification for your specific environment.

## Generate HMMs of TCDB families

* `cd align-tcdb-subfamilies`

The FASTA files of selected SLC-like (sub)families from the TCDB are available in `align-tcdb-subfamilies/input`. The alignments were performed using T-Coffee on a computing cluster. The submit scripts are available for reference.

Output alignments are generated in `align-tcdb-subfamilies/output`.

* `cd make-tcdb-hmms`
* `python runme.py`

Generates HMMs using `hmmbuild` from HMMER based on the alignments generated in the previous step.

Output HMMs are generated in `make-tcdb-hmms/output/*.hmm`.

## Run the initial HMM search

* `python all-search.py >all-search.out 2>all-search.err`

Retrieves all TCDB families with `status > 0` (cross-check with TCDB sequences having `status > 0`), also Pfam domains with `status > 0`, which will be searched for.

HMM files are read from `hmms` subfolder, which contains the merged contents of folders `make-tcdb-hmms/output/*.hmm` and selected SLC-like Pfam models stored in `pfam-hmms/*.hmm`.

Then some sanity checks are performed:
* Check that a HMM file exists for each family/Pfam domain.
* Check that the HMM file for each TCDB family contains the same sequences as in the database.
* Check if all HMM files have a corresponding selected family.

Afterwards, runs `hmmsearch` on each family on the sequences in `uniprot-all.fasta`, output will be placed in the `output` subfolder.

## Upload HMM hits into the MySQL database

* `python upload-hmm-hits.py >upload-hmm-hits.sql 2>upload-hmm-hits.err`
* upload `upload-hmm-hits.sql` to the MySQL DB (ephemeral table).

Reads the `output/*.domtblout.txt` files and selects best (accession, hmm) pairs.

Outputs SQL for table `hmm_hits`. Selects hits (`status = 1`) with bit score > 50.

You have to manually upload the resulting `upload-hmm-hits.sql` into the database.

## Sync UniProt info for HMM hits

* `bash sync-uniprot-info-all.sh >sync-uniprot-info-all.sql 2>sync-uniprot-info-all.err`
* upload `sync-uniprot-info-all.sql` to the MySQL DB (changes to permanent table).

Internally runs `sync-uniprot-info.py`. Processes given `.fasta` and `.txt.gz` files (in UniProt native format). Searches both FASTA and TXT files for proteins existing in the `hmm_hits` and `slc_like` tables with `status > 0`. Proteins where data are new or have been changed compared to those in the `uniprot_proteins` and `uniprot_dbrefs` tables will be updated.

Outputs SQL to update `uniprot_proteins` and `uniprot_dbrefs` tables.

You have to manually upload the resulting `sync-uniprot-info-all.sql` into the database.

## Generate FASTA files of HMM hits for BLAST search

* `python make-hits-fasta.py >make-hits-fasta.out 2>make-hits-fasta.err`

Reads selected HMM hits from `hmm_hits` (`status > 0`) and the corresponding sequences from `uniprot_proteins`.

Outputs files `hits.*.fasta`, where the asterisk stands for the NCBI taxonomy ID of the organism.

## Run all-against-all BLAST search for HMM hits

* `python run-blast.py >run-blast.out 2>run-blast.err`

For each organism (`hits.*.fasta` files), runs `makeblastdb` to create a BLAST DB of all HMM hits, followed by a `blastp` of all hits against this DB. The `makeblastdb` run generates `phr`, `pin` and `psq` files, and each `blastp` run generates a `hits.*.blast.xml` file. These files are then processed to retrieve all hits in a tabular (TSV) format that is written into `run-blast.out` (`stdout`).

## Cluster HMM hits

* `python cluster-hits.py >cluster-hits.sql 2>cluster-hits.err`
* upload `cluster-hits.sql` to the MySQL DB (ephemeral table).

Hits are taken from `hmm_hits` (`status > 0`) and their annotations from `uniprot_proteins`. A graph of sequences is built where sequences are nodes and they are connected if they share any of the annotations ("dbrefs") or their assigned gene symbols.

Contradicting annotations (e.g., different values with the same dbref) are reported. UniGene annotations are dropped from the protein if contradicting UniGene annotations are found, since these have been found very prone to problems.

Sequences are also connected if they share a BLAST hit with sequence identity > 95% from the all-against-all search.

Split constraints are implemented that prevent certain sequence pairs from being joined into the same cluster. This is currently only used to counteract incorrect sequences annotations in UniProt.

One representative sequence is selected from each connected component. If there is a Swiss-Prot (verified) sequence in a component, it is selected as representative. Otherwise the longest sequence is selected. Components with multiple Swiss-Prot sequence are an error and are reported.

Outputs SQL for the `hmm_hits_clusters` table (ephemeral table).

You have to manually upload the resulting `cluster-hits.sql` into the database.

## Sync SLC-like proteins based on HMM hit clusters

* `python sync-slc-like.py >sync-slc-like.sql 2>sync-slc-like.err`
* upload `sync-slc-like.sql` to the MySQL DB (changes to permanent table)

The cluster representatives generated in the last step are actually the found SLC-like proteins. However, the `slc_like` table is permanent and may include manual modifications. Therefore a sync mechanism is used to update it.

Synchronizes the set of representative sequences in `hmm_hits_clusters` with `slc_like`. Inserts new proteins and deletes ones not present, unless they were manually added.

# Generation of SLC-like families

* `python cluster-slc-like.py >cluster-slc-like.sql 2>cluster-slc-like.err`
* upload `cluster-slc-like.sql` to the MySQL DB (ephemeral table).

Performs HMM fingerprint-based clustering of selected proteins from `slc_like` (`status >= 0`). HMM fingerprints are built based on data in `hmm_hits`. Clustering constraints are taken into account from `slc_like_cluster_constraints`. This script also generates flat clusters at threshold value 0.7, while the other script `cluster-slc-like-variable-threshold.py` can use a user-specified threshold.

Outputs SQL for ephemeral tables `slc_like_clusters` and `slc_like_clusters_linkage`.

You have to manually upload the resulting `cluster-slc-like.sql` into the database.

# Generating circular dendrogram

* `python circular-dendrogram.py >circular-dendrogram.out 2>circular-dendrogram.err`

Draws a circular dendrogram of the clustering of SLC-like proteins based on `slc_like`, `slc_like_clusters` and `slc_like_clusters_linkage`. Removes non-human proteins by pruning the tree.

Generates `circular-dendrogam.svg`.

# Phylogenetic trees

* `python make-slc-family-fastas.py >make-slc-family-fastas.out 2>make-slc-family-fastas.err`

Generates FASTA files for each family cluster based on `slc_like_clusters`, for proteins from `slc_like` having `status >= 0`.

FASTA files are generated in the `branches-fastas` folder and the file name contains a hash key of member accession codes. Groups of proteins with the same members should give the same phylogenetic tree.

* `cd branches-phylip ; bash all-trees.sh ; cd ..`

Runs ClustalO and PhyML to generate per-family multiple alignments and (unrooted) phylogenetic trees.

Output is generated in the folder `branches-phylip`.

* `cd branches-notung ; bash all-trees-notung.sh ; cd ..`

Runs NOTUNG to root the phylogenetic trees generated in the previous step. Rooting, rearrangement and reconciliation are performed.

## Generate phylogenetic tree SVG and PDF files

* `cd branches-notung`
* `bash all-trees-svg.sh`

Generates phylogenetic trees for reconciled trees in the `branch_*` subfolders in SVG format. Uses the `make-tree-svg.py` script. Generates output files `make-tree-svg.out` (SVG text file) and `make-tree-svg.err` in each `branch_*` subfolders.

* `python make-slc-family-trees.py >make-slc-family-trees.out 2>make-slc-family-trees.err`

Goes through each family according to `slc_like_clusters` and converts phylogenetic tree images generated in the previous step to PDF files using Inkscape. The SVG files generated in the previous step are modified to include family name annotations based on the `slc_like_clusters` table and stored in `branch_*/branch_*.svg`. The PDF files are generated as `branch_*/branch_*.pdf`. Then, PDF file images for all families are merged using GhostScript into the single file `slc-all-trees.pdf`.

## Generating figure with orthologs and paralogs

* `cd branches-notung`
* `python get-homologues.py >get-homologues.svg 2>get-homologues.err`

Goes through each family according to `slc_like_clusters` and finds orthologs of human proteins in other organisms. Generates a shaded SVG table to show results.

Outputs SVG text.

# Search for Pfam models in TCDB families

* `cd tcdb-pfam-search`
* `python all-fastas.py >all-fastas.fasta 2>all-fastas.err`
* `bash all-fastas-search.sh >all-fastas-search.out 2>all-fastas-search.err`


