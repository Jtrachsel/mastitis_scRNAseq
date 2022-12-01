## John Lippolis mastitis milk and blood single cell RNAseq project

## Updates 12-1-2022  

- Jayne's investigations revealed that blood data is whole blood, not PBMC  
  - large proportion of myleloid cells is expected because of this
- Many cells may in fact be neutrophils (large majority in the milk, many in blood)  
- Because of the above discoveries, the weird cell proportions are not due to poor quality data and low UMI/gene counts are expected for these cell types.  
  - We can relax the QC filtering cutoffs to include more cells.  
- Jayne added new marker genes and rough cell identities.  
- Jayne started work collecting data sets to map our data onto.  This will help with cell type assignments.  
  - One milk dataset of interest uses non ensembl geneIDs.  Need to figure out which database these IDs are from and map to ensembl.  
  
## TODO 12-1-2022  
- Re-run QC with Jayne's relaxed cutoffs  
- Use Jayne's markers to assign cell types and make a nice looking report.  
- Figure out weird gene IDs in human milk data and map to ensembl.  
  
  
  



## reference genome info
from ensembl.org/Bos_taurus/Info/Annotation

Assembly
The ARS-UCD1.2 assembly was submitted by Usda Ars on April 2018. The assembly is on chromosome level, consisting of 2,597 contigs assembled into 2,211 scaffolds. From these sequences, 30 chromosomes have been built. The N50 size is the length such that 50% of the assembled genome lies in blocks of the N50 size or longer. The N50 length for the contigs is 25,896,116 while the scaffold N50 is 103,308,737.


From: https://www.ebi.ac.uk/ena/browser/view/GCA_002263795.2 
Assembly: GCA_002263795.2
This is a de novo assembly of Bos taurus using long reads for assembly and short reads for scaffolding and polishing. The assembly was scaffolded using a combination of in vitro chromosome conformation capture sequencing, optical maps, and a recombination map. Subsequent scaffolding was confirmed using genetic linkage maps and a Radiation Hybrid map. This material is based on work supported by the National Institute of Food and Agriculture, U.S. Department of Agriculture, National Research Support Project NRSP8 (Cattle Genome Coordination)
Comment
In May 2022 254 sequences were suppressed because they were found to be contaminants.

Organism: Bos taurus
Accession:GCA_002263795
Assembly Name:ARS-UCD1.2
Assembly Level:chromosome
Genome Representation:full
Total Length:2,711,209,833
Ungapped Length:2,711,181,671
N50:103,308,737
Spanned Gaps:386
Unspanned Gaps:0
Scaffold Count:1,957
Count Contig:2,343
Contig N50:25,896,116
Contig L50:32
Contig N75:14,624,078
Contig N90:4,439,585
Scaf L50:12
Scaf N75:69,862,954
Scaf N90:51,992,305
Replicon Count:31
Count Non Chromosome Replicon:1
Count Alt Loci Units:0
Count Regions:0
Count Patches:0
ENA-LAST-UPDATED:2022-05-20





Gene annotation
The gene annotation process was carried out using a combination of protein-to-genome alignments, annotation mapping from a suitable reference species and RNA-seq alignments (where RNA-seq data with appropriate meta data were publicly available). For each candidate gene region, a selection process was applied to choose the most appropriate set of transcripts based on evolutionary distance, experimental evidence for the source data and quality of the alignments. Small ncRNAs were obtained using a combination of BLAST and Infernal/RNAfold. Pseudogenes were calculated by looking at genes with a large percentage of non-biological introns (introns of <10bp), where the gene was covered in repeats, or where the gene was single exon and evidence of a functional multi-exon paralog was found elsewhere in the genome. lincRNAs were generated via RNA-seq data where no evidence of protein homology or protein domains could be found in the transcript.


Statistics
Summary
Assembly	ARS-UCD1.2, INSDC Assembly GCA_002263795.2, Apr 2018
Base Pairs	2,715,853,792
Golden Path Length	2,715,853,792
Annotation provider	Ensembl
Annotation method	Full genebuild
Genebuild started	Sep 2018
Genebuild released	Dec 2018
Genebuild last updated/patched	Nov 2018
Database version	107.12
Gene counts
Coding genes	21,880
Non coding genes	5,235
Small non coding genes	3,375
Long non coding genes	1,488
Misc non coding genes	372
Pseudogenes	492
Gene transcripts	43,984
Other
Genscan gene predictions	46,441
Short Variants	98,844,739
Structural variants	18,942

# mastitis_scRNAseq
