# G2/Jockey-3 Phylogenetic Analysis

## Transposon Annotation and Alignment

### Annotation

The genomic sequences of each *Drosophila* species were annotated with a custom transposable element library in a fasta format with Repeatmasker v4.1.0. Example below with an expected GFF output file:

```
RepeatMasker -no_is -gff -a -inv -pa 18 -lib specieslib_mod2_centromere_v4.fasta -div 30 dsim_scaffold2_2019.fasta
```

### Extracting sequences

The GFF file with transposon annotations were used to extract the sequences. G2/Jockey-3 and other related transposon sequences were extracted with Bedtools v2.29.2 with an example bash script below:

```
## simulans
fastaFile="dsim_scaffold2_2019.fasta"
gffFile="dsim_scaffold2_2019.fasta.out.gff"
species="simulans"
speciesAb="Dsim" ## abreviated species name
tes=("Jockey-3_DM" "Jockey-3_DSim" "Jockey-1_DSe" "Jockey-7_DYa")

for te in ${tes[@]}; do
	grep -w ${te} $gffFile > ${species}.${te}.gff
	bedtools getfasta -fi $fastaFile -bed ${species}.${te}.gff -s -fo ${species}.${te}.fasta
	# cleanup of sequence names
	sed -i '' 's/(/|/g' ${species}.${te}.fasta 
	sed -i '' 's/)/|/g' ${species}.${te}.fasta
	sed -i '' 's/:/|/g' ${species}.${te}.fasta
	sed -i '' "s/>/>${speciesAb}|${te}|/g" ${species}.${te}.fasta
done
```

### Sequence alignment

Fasta file sequences were imported into Geneious v8.1.6. Sequences were aligned with MAFFT and errors were manually corrected. The open reading frames of the original transposons with [NCBI OrfFinder](https://www.ncbi.nlm.nih.gov/orffinder/). Sequences aligning to ORF2 were extracted and exported for phylogenetic analysis. Consensus sequences were generated for each clade of G2/Jockey-3 by isolating the specific species or clade sequences from the master alignment file and using the Geneious "Generate Consensus Sequence" tool.


## Phylogenetic analysis

The resulting alignment files were used for phylogenetic tree construction using RAxML v8.2.11. An example bash script to generate the AUTOMRE files, in NEWICK format, is below:

```
alignment="alignment_Jockey-3_melsimyak_400_ORF2_mafft_Jockey-1_Dse"
raxmlHPC-PTHREADS -s "${alignment}.fasta" -n "${alignment}.automre" -m GTRGAMMA -T 24 -d -p 12345 -# autoMRE -k -x 12345 -f a
```

Two newick files were generated, RAxML\_bipartitions.alignment\_Jockey-3\_melsimyak\_400\_ORF2\_mafft\_Jockey-1\_Dse.automre for all G2/Jockey-3 fragments and RAxML\_bipartitions.alignment\_jockey\_family\_Final.automre for the consensus sequences.


## Tree illustration

The phylogenetic trees were illustrated with the ggtree package, a phylogenetic illustration package extension for ggplot. Additional packages including ape were used to import the NEWICK phylogenetic tree files and other tidyverse packages. The specific R scripts are listed below and can be run from the command line, changing the path of the AUTOMRE files as needed within the R files.

```
## All fragments
Rscript Jockey-3_melsimyak_ORF2_PhylogeneticTree.R

## Consensus sequences
Rscript Jockey_family_PhylogeneticTree.R
``` 


