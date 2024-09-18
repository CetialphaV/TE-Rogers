import pandas as pd
import subprocess
import os
import csv

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)




class TE_ReferenceCreator:
    import os
    import subprocess

    def _dataframe_to_fasta(self, frameToConvert, fasta_filename):
        print("Beginning conversion of dataframe to fasta file")
        with open(fasta_filename, 'w') as fasta_file:
            for index, row in frameToConvert.iterrows():
                header = f">{row[0]}"  # Add '>' before the header (column 0)
                sequence = row[1]  # Sequence (column 1)
                fasta_file.write(f"{header}\n")
                fasta_file.write(f"{sequence}\n")

    def convertBEDToFasta(self, bed_path, genomeRefPath, outputName="TESequenceReference"):
        # Check if the BED file is in the current working directory, and if so, make the path absolute
        if os.path.isfile(bed_path) and not os.path.isabs(bed_path):
            bed_path = os.path.abspath(bed_path)

        # Check if the genome reference path is in the current working directory, and if so, make the path absolute
        if os.path.isfile(genomeRefPath) and not os.path.isabs(genomeRefPath):
            genomeRefPath = os.path.abspath(genomeRefPath)

        # Construct the output fasta file path


        # Construct the bedtools command
        command = [
            "bedtools", "getfasta",
            "-fi", genomeRefPath,
            "-bed", bed_path,
            "-fo", outputName,
            "-tab"
        ]

        # Run the bedtools command
        try:
            subprocess.run(command, check=True)
            print(f"FASTA file created: {outputName}")
        except subprocess.CalledProcessError as e:
            print(f"Error running bedtools: {e}")

        referenceBedFrame = pd.read_csv(bed_path, sep="\t", header=None)
        referenceBedFrame[5] = referenceBedFrame[5].replace({'+': '.', '-': 'C'})
        part1 = referenceBedFrame[3].astype(str) + "@" + referenceBedFrame[4].astype(str) + '|' + referenceBedFrame[6].astype(str) + '|' + referenceBedFrame[5].astype(str) + '|' + referenceBedFrame[0].astype(str)
        part2 = referenceBedFrame[1].astype(str) + '-' + referenceBedFrame[2].astype(str)
        combined = part1 + ':' + part2 + '\\'
        referenceBedFrame['Combined'] = combined


        outputFrameWithSequence = pd.read_csv(outputName, sep="\t", header=None)
        outputFrameWithSequence[0] = referenceBedFrame['Combined']

        self._dataframe_to_fasta(outputFrameWithSequence, outputName)
        print("Conversion complete")




    def parseAlignFile(self, pathToAlignFile, outName="TE_Sequence_Reference.tsv"):
        print("Starting to parse repeatmasker alignment file")
        TEInfo = []
        with open(pathToAlignFile) as alignFile:
            header = ""
            lastLine = ""
            entrySequence = ""
            for lineNum, line in enumerate(alignFile):

                # if lineNum % 1000000 == 0:
                #     print(f"Starting {int(lineNum / 1000000)}th millionth line")

                #Entries end with Gap_init in the line
                if "Gap_init" in lastLine:

                    #Add last entry as a line in list with sequence
                    headerLine = header.split()

                    if headerLine[8] != "C":
                        headerLine.insert(8, ".")


                    headerLine.append(entrySequence)

                    TEInfo.append(headerLine)


                    #reset header, lastline, and sequence
                    header = ""
                    lastLine = ""
                    entrySequence = ""

                    continue

                lastLine = line

                strippedLine = line.strip()

                # Skip blank lines
                if strippedLine == "\n" or (not strippedLine):
                    continue

                #If header not set, assign first line of new entry as header
                if not header:
                    header = strippedLine
                    continue


                #Sequence will be in lines beginning in chr
                if "chr" in strippedLine:
                    #Sequence is in 3rd position (index 2)
                    sequence = strippedLine.split()[2]
                    #Remove any non-nucleotides
                    sequence = sequence.replace('-', '')
                    sequence = sequence.replace('X', 'N')
                    entrySequence += sequence

        print("Finished parsing from Repeatmasker alignment file. Converting to DataFrame")
        TE_sequenceFrame = pd.DataFrame(TEInfo)
        print("Saving DataFrame to CSV")
        TE_sequenceFrame.to_csv(outName, sep="\t")


        return TE_sequenceFrame


if __name__ == "__main__":
    pathToAlignFile = "hg38.fa.align"
    referenceCreator = TE_ReferenceCreator()
    # TE_sequenceFrame = referenceCreator.parseAlignFile(pathToAlignFile)
    # print(TE_sequenceFrame.head())
    pathToGenomeReference = "/Users/loganrivera/PycharmProjects/TE_neighborAnalysis/GRCm38.primary_assembly.genome.fa"
    pathToTEBedFile = "/Users/loganrivera/PycharmProjects/TE_neighborAnalysis/testBEDFileForRogers.bed"
    outPrefix = "IAmATestOut.bed"
    referenceCreator.convertBEDToFasta(bed_path=pathToTEBedFile,
                                       genomeRefPath=pathToGenomeReference,
                                       outputName=outPrefix)
