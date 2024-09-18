from Bio import SeqIO
import pandas as pd
import subprocess

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


class NeighborhoodCreator():
    def __init__(self, pathToUTRFasta, pathToGTPFile, baseName="TE_analysis"):
        self.pathToUTRFasta = pathToUTRFasta
        self.pathToGTFFile = pathToGTPFile


        self.pathToGTF_BED = f"{baseName}_geneReferenceFile.bed"
        self.neighborhood_BED = f"{baseName}_neighborCoords.bed"
        self.pathToOverlapFile = f"{baseName}_neighborhoodGenes.bed"


    def createNeighborhoodBED(self, neighborhoodSize=10000):
        print("Creating Neighborhood TE BED File")
        print(self.pathToUTRFasta)
        with open(self.pathToUTRFasta, 'r') as TE_fasta, open(self.neighborhood_BED, "w") as neighborBED:
            for record in SeqIO.parse(TE_fasta, "fasta"):
                header = record.id
                headerSplit = header.split("|")
                headerPos = headerSplit[-1]
                chromosome, coordinates = headerPos.split(":")
                coordinates = coordinates.split("_")[0]
                startPos, endPos = coordinates.split("-")
                startPos = int(startPos) - 1 # BED format  uses 0 based start positions
                endPos = int(endPos)

                startPosNew = max([(startPos - neighborhoodSize), 0])
                endPosNew = endPos + neighborhoodSize

                headerName = "_".join(headerSplit[:2])
                headerStrand = headerSplit[2]

                if headerStrand == ".":
                    headerStrand = "+"
                elif headerStrand == "C":
                    headerStrand = "-"

                lineForBed = f"{chromosome}\t{startPosNew}\t{endPosNew}\t{headerName}\t0\t{headerStrand}\t{startPos}\t{endPos}\n"
                neighborBED.write(lineForBed)

        print("Finished creating neighborhood TE file")


    def _extract_feature(self, attributes, attribute_name):
        for attribute in attributes.split(';'):
            if attribute_name in attribute:
                return attribute.split('"')[1]  # Extract the value within quotes
        return None




    def createGeneRefBED(self):
        print("Loading GTF annotation file")
        GTFFrame = pd.read_csv(self.pathToGTFFile, sep="\t", skiprows=[0, 1, 2, 3, 4], header=None)
        print("Finished loading GTF annotation file. Beginning creation of BED file")
        gtf_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        GTFFrame.columns = gtf_columns
        GTFFrame = GTFFrame[GTFFrame["feature"] == "gene"]
        GTFFrame = GTFFrame.drop(columns=['source', 'score', 'frame'])

        GTFFrame['gene_id'] = GTFFrame['attribute'].apply(lambda x: self._extract_feature(x, 'gene_id'))
        GTFFrame['gene_type'] = GTFFrame['attribute'].apply(lambda x: self._extract_feature(x, 'gene_type'))
        GTFFrame['gene_name'] = GTFFrame['attribute'].apply(lambda x: self._extract_feature(x, 'gene_name'))

        GTFFrame = GTFFrame.drop(columns=['attribute', 'feature'])
        GTFFrame = GTFFrame[['seqname', 'start', 'end', 'gene_id', 'strand', 'gene_type', 'gene_name']]

        GTFFrame['start'] = GTFFrame['start'] - 1

        GTFFrame.to_csv(self.pathToGTF_BED, sep='\t', header=False, index=False)
        print(f"Bed annotation file created successfully and saved to {self.pathToGTF_BED}")


    def determineOverlappingFeatures(self):

        print("Finding Genes in Neighborhood")
        command = [
            "bedtools", "intersect",
            "-a", self.neighborhood_BED,
            "-b", self.pathToGTF_BED,
            "-wa",  # Write the original entry in A for each overlap
            "-wb"  # Write the original entry in B for each overlap
        ]
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            print("Bedtools intersect executed successfully.")
            print(f"Writing results to {self.pathToOverlapFile}")
            with open(self.pathToOverlapFile, "w") as outputFile:
                for line in result.stdout.split("\n"):
                    outputFile.write(f"{line}\n")
            print("Results saved successfully")
        except subprocess.CalledProcessError as e:
            print("Error executing bedtools intersect:", e)
            print("Error output:", e.stderr)
        print("Genes in Neighborhood Have Been Found!")


    def _findDistanceBetweenFeatures(self, featureOneCoords, featureTwoCoords):
        if featureOneCoords[0] > featureTwoCoords[1]: #Check if the TE is after the end of the gene
            distance = featureTwoCoords[1] - featureOneCoords[0]
            #If the gene is upsteam so distance is negative (end cord of gene - start cord of TE)
        elif featureOneCoords[1] < featureTwoCoords[0]: #Check if the end of TE is before the start of the gene
            distance = featureTwoCoords[0] - featureOneCoords[1]
            # If the gene is downstream so distance is postive (Beginning of gene - end of TE)
        else:
            distance = 0 #If they overlapp then distance is zero
        return distance


    def createFinalDistanceFile(self, outputName):
        print("Calculating Gene Distances and creating final result file")
        overlapBEDFrame = pd.read_csv(self.pathToOverlapFile, sep="\t", header=None)
        columnNames = ["TE_chr", "neighbor_start", "neighbor_stop", "TE_Class", "TE_score", "TE_strand", "TE_start", "TE_stop",
                       "Gene_chr", "Gene_start", "Gene_end", "Gene_ID", "Gene_strand", "Gene_class", "Gene_name"]
        overlapBEDFrame.columns = columnNames

        overlapBEDFrame['Distance_TE_Gene'] = overlapBEDFrame.apply(lambda row: self._findDistanceBetweenFeatures((row['TE_start'], row['TE_stop']),
                                                                                  (row['Gene_start'], row['Gene_end'])), axis=1)

        overlapBEDFrame = overlapBEDFrame.drop(columns=['neighbor_start', 'neighbor_stop', 'TE_score',
                                                        'Gene_chr', 'Gene_start', 'Gene_end'])

        overlapBEDFrame = overlapBEDFrame[['TE_chr', 'TE_start', 'TE_stop', 'TE_strand', 'TE_Class',
                 'Gene_ID', 'Gene_name', 'Gene_strand', 'Gene_class', 'Distance_TE_Gene']]

        overlapBEDFrame.sort_values(by=['TE_chr', 'TE_start', 'TE_stop'], inplace=True)
        overlapBEDFrame['Unique_Gene_Count'] = overlapBEDFrame.groupby(['TE_start', 'TE_stop', 'TE_chr'])['Gene_ID'].transform('nunique')
        print(overlapBEDFrame.head())
        overlapBEDFrame['TE_Number'] = overlapBEDFrame['TE_Class'].str.extract(r'@(\d+)_')
        overlapBEDFrame['TE_Class'] = overlapBEDFrame['TE_Class'].str.replace(r'@\d+', '', regex=True)
        print(overlapBEDFrame['TE_Number'].unique())
        overlapBEDFrame.to_csv(outputName, sep="\t")
        print("Finished! Process Complete")




if __name__ == "__main__":
    pathToUTRFasta = "LINE_L1_5PrimeUTRs.fa"
    pathToGTPFile = "gencode.v46.basic.annotation.gtf"

    neighborHoodCreator = NeighborhoodCreator(pathToUTRFasta, pathToGTPFile)


    finalGeneNeighborhoodName = "TE_gene_neighborhood.tsv"
    #
    neighborHoodCreator.createGeneRefBED()
    neighborHoodCreator.createNeighborhoodBED(neighborhoodSize=5000000)
    neighborHoodCreator.determineOverlappingFeatures()
    neighborHoodCreator.createFinalDistanceFile(finalGeneNeighborhoodName)


