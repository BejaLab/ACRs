ACRs
====

This repository hosts the bioinformatic analysis worklow used in [Oppermann et al, **Robust Optogenetic Inhibition with Red-light-sensitive Anion-conducting Channelrhodopsins**, *eLife*](https://doi.org/10.7554/eLife.90100).

The [snakemake](https://snakemake.readthedocs.io/)-based code for the workflow is under `workflow/`. The dependencies are under conda control (`snakemake --use-conda`), see `workflow/envs`. The analysis files are under `analysis/`, in particular `ACRs/analysis/pdb` contains the output files for the structural alignment (`sequences.aln` is the main output) and `analysis/all` hosts the various files generated for the phylogenetic analysis:

* `sequences.fasta`: unaligned sequences 
* `cdhit.fasta`: non-redundant set of sequences
* `mafft.fasta`: mafft alignment
* `trimal.fasta`: trimal trimmed alignment
* `iqtree.treefile`: the phylogenetic tree from iqtree with ultrafasta bootstrap support values

The metadata for the sequences can be found in `metadata/Channelrhodopsins_Updated_List.xlsx` which is a snapshot of the [Catalogue of Natural Channelrhodopsins](https://doi.org/10.5281/zenodo.5749640).
