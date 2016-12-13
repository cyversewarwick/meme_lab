# MEME-LaB

## The Purpose of the Algorithm

MEME-LaB is a *de novo* motif overrepresentation algorithm, finding the top overrepresented motifs in provided gene groups. In contrast to the hypergeometric motif test tool (HMT), also offered on CyVerse, MEME-LaB does not use any predefined motif sequences, instead identifying novel potential binding sites. The search itself is performed using [MEME][meme]. The output is formatted in a user-friendly manner, allowing at-a-glance information on the quality of the detected motifs and more detailed reports on positioning within the promoter and exact sequence matches.

## Test Run

If you want to try MEME-LaB out and get a feel for its output formatting, you can find demonstration data at `iplantcollaborative/example_data/cyverseuk/meme-lab_testdata` under Community Data. Parameters can be left unaltered for the demonstration run.

## Input

### Genome

The app requires a genome sequence on input so that promoter sequences can be easily automatically extracted. Please provide the genome as a single FASTA file.

### GFF3 Annotation

The GFF3 annotation needs to be compatible with the provided genome, and is used to locate the genes for promoter extraction. In the second app, the GFF3 annotation is used to create a complete universe to use as the background in hypergeometric testing.

### GFF3 Gene ID Attribute

A GFF3 file can carry a lot of information about an organism's genes, whilst the program is after the very basics - a distinct and discernible gene ID style that is also used in the input of gene groups for the second app, along with corresponding positioning. The GFF3 file is filtered to the lines that contain gene information, but the script subsequently needs information on which of the information fields to use as the identifier. For example, the Arabidopsis test data provided (`iplantcollaborative/example_data/cyverseuk/meme-lab_testdata/annot.gff3` under Community Data) has the `gene_id=` field correspond to AGI identifiers, which are the widely accepted locus code nomenclature for Arabidopsis.

### Promoter Length

How many base pairs upstream of the transcription start site to take as the promoter sequence for analysis. Some common values include 200, 500 and 1000 base pairs.

### Gene Group Input File

A tab-delimited file with two columns, with the first column being the number of the gene group and the second column being a single gene ID in that gene group. For example, if you wanted to analyse two gene groups, with the first one containing 10 genes and the second one containing 15 genes, the input file would be 25 lines long, with 10 lines of gene group 1 IDs and 15 lines of gene group 2 IDs. Consult an example Arabidopsis file at `iplantcollaborative/example_data/cyverseuk/meme-lab_testdata/input.txt` under Community Data.

### Number Of Motifs Per Gene Group

The number of top motifs to report per gene group. The motifs will be sorted on e-value, i.e. an estimate of the number of motifs at least as interesting as the found motif that would be found by chance if the input sequence bases were randomly shuffled, but no cutoff will be placed on the e-value itself.

### Minimum Motif Width

The minimum number of base pairs of length a motif has to have to be reported.

### Maximum Motif Width

The maximum number of base pairs of length a motif has to have to be reported.

## Output

### `output/results.html`

The basic MEME-LaB report, featuring a condensed version of the information featured in the MEME motif discovery output. The webapp shows the top motifs for each of the provided clusters, with the provided measures (information content, e-value) being helpful in identifying the most relevant motifs. Information content rewards the motifs for having well defined conservation (i.e. little variation in possible occurring bases) at each of the positions, the higher the number the better. The e-value tracks how likely it would be to spot motifs at least as interesting as the found one if the bases of the input sequences were to be randomly shuffled. The lower the e-value the better, with an e-value higher than 0.01 possibly implying the motif as an artifact. For details on information content and e-values, consult [this protocol][protocol]. The positioning of the found motif within the promoters of the genes in the group can be viewed in a pop-up window opened by pressing the blue `more...` hyperlink.

### `output/meme/`

The MEME analysis output for each of the gene groups, featuring more detailed information on the top motifs. A particular gene group's analysis output is housed in a subfolder named after the gene group name provided in the input file. Within said subfolder, the detailed results can be accessed in `meme.html`. Pressing the blue down arrow in the More column of the motif of interest will reveal the exact sequence matches of the motif in the promoters of the genes in the group, along with several additional MEME measures of motif quality. This output is a more advanced form of the basic report and only needs to be investigated in specific in-depth cases. The folder also features the MEME output in a number of other formats along with `.eps` files of the plots of the motifs' position weight matrices.

[meme]: http://nar.oxfordjournals.org/content/early/2009/05/20/nar.gkp335.short
[protocol]: http://www.sdsc.edu/~tbailey/MEME-protocol-draft2/protocols.html