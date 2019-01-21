# Make a bedfile of the yeast proteome
R script that provides the codon position for every amino acid in the yeast proteome. This comes handy when say you have a list of transcription start sites in the yeast genome and want to find the closest start codon to each transcription start site (TSS), calculate the distance between the TSS and the start codon, get the length of the resulting coding region etc.

## How the script works
The script uses a FASTA file - orf_trans.fasta, of all protein sequences in *S. cerevisiae* as the input. This file was obtained from [SGD](https://www.yeastgenome.org/). The first line in the FASTA file for each protein sequence contains the systematic name for the gene along with the chromosomal coordinates for the exons that code for the protein. The script makes use of this information to assign each amino acid in the protein sequence to the chromosomal coordinates of its respective codon.

If you want to try running the script yourself, just make sure you run the script from the folder where it is saved and that the FASTA file is also present in the same folder. It is recommended to run the script on a HPC cluster.

## The output file
The output file of the script has also been provided in this repository - cerevisiae_proteome.bedgraph. The sorted version of the same file is also provided - The format of the output file is as follows:
- Column 1: Chromosome number
- Column 2: 5' position of codon
- Column 3: 3' position of codon
- Column 4: strand (+ is Watson strand, - is Crick strand)
- Column 5: Amino acid
- Column 6: Gene name
```
chrI    1807    1809    -       *       YAL068C
chrI    1810    1812    -       N       YAL068C
chrI    1813    1815    -       A       YAL068C
chrI    1816    1818    -       I       YAL068C
chrI    1819    1821    -       T       YAL068C
chrI    1822    1824    -       Y       YAL068C
```
Since this example is of a gene on the negative strand, the 5' most gene coordinate will be a stop codon (\*) and the 3' most gene coordinate will be the start codon.

For genes that have a codon spanning the exon-exon junction, the codon is represented as follows:
```
chrII   462205  462207  +       E       YBR111W-A
chrII   462208  462209  +       L       YBR111W-A
chrII   462290  462290  +       L       YBR111W-A
chrII   462291  462293  +       I       YBR111W-A
```
The codon coding for the leucine spans the exon-exon junction in this gene and is represented twice in the file. The region from 462210 - 462289 contains the intron sequence.
