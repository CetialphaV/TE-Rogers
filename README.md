# TE-ROGERS: Transposable Element - Regulatory Overview and Gene Environment Relationship Surveyor


##Overview

TE-ROGERS is a toolkit designed to analyze transposable elements (TEs) and their surrounding genomic environments. It helps identify 5' UTR regions of TEs, find nearby genes, and integrate differential expression data. The toolkit also calculates distances between TEs and nearby genes to assess the regulatory landscape. BEDTools is required to process genomic intervals.

##Prerequisites

Ensure the following dependencies are installed:

Python 3.x
Required Python libraries: pandas, os, subprocess, Bio (for SeqIO)
BEDTools
Genome reference files in FASTA, BED, and GTF format
TE sequence files
Differential expression data (CSV)
Necessary Base Files

BED file format for TEs (Required for Step 1): The BED file used in Step 1 must contain the following columns:
Chromosome (chr)
Start position
End position
TE name
TE unique ID (mandatory, not optional)
Strand (+ or -)
TE class (e.g., LINE/L1, Simple_repeat)
TE length (last column)
Example BED file content:
chr1    3000001    3000097    L1MdFanc_I    0    -    LINE/L1    96
chr1    3000098    3000123    (T)n    1    +    Simple_repeat    25

Genome Reference File (FASTA): Path to the genome reference file (e.g., GRCm38).
GTF Annotation File: GTF file for gene annotations (e.g., gencode.v46.basic.annotation.gtf).
Pipeline Steps

##Step 1: Create TE Reference Data
In this step, a reference FASTA file is created from the BED coordinates of transposable elements (TEs) using the reference genome file. The BED file must follow the format mentioned above, with the last column representing the TE length and the fifth column containing the TE unique ID.

The function convertBEDToFasta will use BEDTools to extract the sequences of TEs from the genome reference based on their coordinates in the BED file.

Convert BED File to FASTA using BEDTools

referenceCreator = TE_ReferenceCreator()
referenceCreator.convertBEDToFasta(bed_path="testBEDFileForRogers.bed",
                                   genomeRefPath="/path/to/GRCm38.primary_assembly.genome.fa",
                                   outputName="TESequenceReference.fasta")
This process includes:

Ensuring the BED file is correctly formatted.
Running the bedtools getfasta command to generate the FASTA file of TE sequences.
Enhancing the headers in the resulting FASTA file to include the TE name, class, family, chromosome location, start and stop coordinates, strand orientation, unique ID, and length.
Example FASTA Header:
>TE_Class@TE_ID|TE_Family|strand|chromosome:start-end|TE_Length
AGCTAGCTAGCTAGCTAG...

The final FASTA file is saved and used in subsequent steps for further analysis.

##Step 2: Find 5' UTR Regions in TEs
Use the UTR_Finder class to identify the 5' UTR regions of TEs by performing a BLAST search. The BLAST database is created from the TE reference file, and the first translatable protein is used for the search.

 
 
UTRFinder = UTR_Finder(TE_Family="LINE/L1",
                       pathToFirstProteinSequence=pathToFirstProteinSequence,
                       pathToTERef=pathToTERef)
                       
UTRFinder.createBLASTDB()

UTRFinder.performBLASTSearch(pval="1e-7")

UTRFinder.findUTR_TEs(utrLengthThreshold=500, startCodonDistance=None,
                      outputFasta="LINE_L1_5PrimeUTRs.fa", startCodons=[])
                      
##Step 3: Create Gene Neighborhood Data

Using the NeighborhoodCreator, you can generate a neighborhood file around TEs and calculate the distances between TEs and nearby genes using BEDTools.

1. Create Neighborhood BED File

This step creates a BED file around each TE with the specified neighborhood size.

 
 
neighborHoodCreator = NeighborhoodCreator(pathToUTRFasta="LINE_L1_5PrimeUTRs.fa", 
                                          pathToGTPFile="gencode.v46.basic.annotation.gtf",
                                          baseName="TE_analysis")
neighborHoodCreator.createNeighborhoodBED(neighborhoodSize=5000000)  # Adjust size as needed
2. Create Gene Reference BED File

Convert the GTF annotation file into a BED file for use in BEDTools intersection.

 
 
neighborHoodCreator.createGeneRefBED()
3. Determine Overlapping Genes in TE Neighborhood

Find overlapping genes within the TE neighborhood by intersecting the TE BED file and the gene BED file using BEDTools.

 
 
neighborHoodCreator.determineOverlappingFeatures()
Step 4: Associate Genes with the Closest TE
In this step, an output file is created where each gene is associated with the closest TE based on calculated distances. The final result includes TE-gene associations and the distances between them.

 
 
neighborHoodCreator.createFinalDistanceFile(outputName="final_neighborhoodSummaryFile.tsv")
The final file contains the following fields:

TE_chr: Chromosome of the TE

TE_start, TE_stop: Coordinates of the TE

TE_strand: Strand of the TE

TE_Class: Classification of the TE

Gene_ID, Gene_name: Gene identifier and name

Gene_strand: Strand of the gene

Gene_class: Classification of the gene

Distance_TE_Gene: Calculated distance between the TE and the gene

Unique_Gene_Count: Number of unique genes associated with the TE

TE_Number: TE identifier number

Step 5: Process Differential Expression Data

Merge TE and gene differential expression data into a final summary file for visualization. The differential expression data for TEs and genes is incorporated using fold change and p-value cutoffs.

Extract differential expression data:
 
 
-RNA differential data
proba_not_de_df.to_csv(os.path.join(pathToOutput, "RNA_pval.csv"), sep="\t")
lfc_mean_df.to_csv(os.path.join(pathToOutput, "RNA_LFC.csv"), sep="\t")

-TE differential data
proba_not_de_df.to_csv(os.path.join(pathToOutput, "TE_pval.csv"), sep="\t")
lfc_mean_df.to_csv(os.path.join(pathToOutput, "TE_LFC.csv"), sep="\t")
Create the final summary file:
 
 
teGeneProcessor = TEGeneProcessor(TEFoldPath="scRNA_cleanedDiffData/TE_LFC.csv",
                                  TEPvalPath="scRNA_cleanedDiffData/TE_pval.csv",
                                  GeneFoldPath="scRNA_cleanedDiffData/RNA_LFC.csv",
                                  GenePvalPath="scRNA_cleanedDiffData/RNA_pval.csv",
                                  TE_gene_neighborhoodPath="TE_analysis_gene_neighborhood.tsv")
teGeneProcessor.createSummaryFile("TE_Line_neighborhoodSummaryFile.tsv")

##Step 6: Visualization
Use the TE_Gene_ResultVisualizer to create visualizations of TE and gene relationships.


visualizer = TE_Gene_ResultVisualizer(summaryFileName="TE_Line_neighborhoodSummaryFile.tsv")
visualizer.plotResults()
##Output Files

TESequenceReference.fasta: The final TE reference FASTA file.
LINE_L1_5PrimeUTRs.fa: FASTA file of 5' UTR sequences.
TE_analysis_gene_neighborhood.tsv: Gene neighborhood file with distances to TEs.
closestDistance_neighborhoodSummaryFile.tsv: File associating each gene with the closest TE and the calculated distance.
TE_Line_neighborhoodSummaryFile.tsv: Summary file combining TE and gene expression data.
Required Tools


