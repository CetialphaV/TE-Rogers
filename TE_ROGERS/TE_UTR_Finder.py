import pandas as pd
import os
import csv
import subprocess
from Bio import SeqIO


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

class UTR_Finder:

    def __init__(self, TE_Family, pathToFirstProteinSequence, pathToTERef):
        self.TE_Family = TE_Family
        self.pathToFirstProteinSequence = pathToFirstProteinSequence

        TERefEnding = pathToTERef.split(".")[-1]
        if TERefEnding == "tsv":
            self.TERefIsTSV = True
        elif TERefEnding == "fa" or TERefEnding == "fasta":
            self.TERefIsTSV = False
        else:
            raise ValueError(f"Invalid file extension '{TERefEnding}'. Allowed extensions are: .tsv, .fa, .fasta")

        self.pathToTEReference = pathToTERef

        self.correctedTEFamilyName = TE_Family.replace('/', '_')
        self.outDataBaseName = f"{self.correctedTEFamilyName}_TE_Ref"
        self.fastaRefFile = self.outDataBaseName + ".fasta"
        self.blastOutFilePath = f"BLAST_{self.correctedTEFamilyName}/resultsFor_{self.pathToFirstProteinSequence.split('.')[0]}"
        self.fastaRefFile = f"{self.correctedTEFamilyName}_TE_Ref.fasta"

    def _loadTERefDataFrame(self):
        TERefFrameFull = pd.read_csv(self.pathToTEReference, sep="\t", index_col=0)
        TERefFrameFull[['TE_Motif', 'TE_Family']] = TERefFrameFull['9'].str.split('#', expand=True)
        print(f"Possible Families: {TERefFrameFull['TE_Family'].unique()}")

        # Filter for SelectedFamily
        TERefFrameFull = TERefFrameFull[TERefFrameFull["TE_Family"] == self.TE_Family]


        return TERefFrameFull


    def _convertFilteredDataFrameToFasta(self, filteredTERef, output_fasta_name):
        # Construct header
        filteredTERef['header'] = '>' + filteredTERef['TE_Motif'] + \
                                  '|' + filteredTERef['TE_Family'] + \
                                  '|' + filteredTERef["8"] + \
                                  '|' + filteredTERef["4"] + ':' + \
                                  filteredTERef["5"].astype(str) + \
                                  '-' + filteredTERef["6"].astype(str)

        fasta_lines = filteredTERef['header'] + '\n' + filteredTERef['15']
        print(output_fasta_name)
        fasta_lines.to_csv(
            output_fasta_name,
            index=False,
            header=False,
            sep='\n',
            quoting=csv.QUOTE_NONE,
            escapechar='\\'
        )

    def _filterTEFasta(self, outputName):
        filtered_sequences = []
        
        # Read the FASTA file and filter sequences by the header
        for record in SeqIO.parse(self.pathToTEReference, "fasta"):
            if self.TE_Family in record.description:  # Check if TE_Family is in the header
                filtered_sequences.append(record)
        
        # Check if any sequences matched the filter_value
        if not filtered_sequences:
            raise ValueError(f"No sequences found with '{self.TE_Family}' in the header.")
        
        # Write the filtered sequences to the output FASTA file with each sequence on one line
        with open(outputName, "w") as output_handle:
            for record in filtered_sequences:
                output_handle.write(f">{record.id}\n{str(record.seq)}\n")
        
        print(f"Filtered sequences saved to {outputName}.")
       
    def createBLASTDB(self):
        if self.TERefIsTSV:
            filteredTERef = self._loadTERefDataFrame()
            self._convertFilteredDataFrameToFasta(filteredTERef, f"{self.outDataBaseName}.fasta")
        else:
            self._filterTEFasta(outputName=f"{self.outDataBaseName}.fasta")
        os.makedirs(f"{self.outDataBaseName}_database", exist_ok=True)
        os.system(
            f"makeblastdb -in {self.outDataBaseName}.fasta -dbtype nucl -out {self.outDataBaseName}_database/{self.outDataBaseName}")


    def performBLASTSearch(self, pval="1e-7"):
        # Ensure the output directory exists
        output_dir = os.path.dirname(self.blastOutFilePath)
        os.makedirs(output_dir, exist_ok=True)

        # Construct the tBLASTn command
        tblastn_command = [
            "tblastn",
            "-query", self.pathToFirstProteinSequence,
            "-db", f"{self.outDataBaseName}_database/{self.outDataBaseName}",
            "-out", self.blastOutFilePath,
            "-evalue", f'{pval}',
            "-outfmt", "6", # Tabular format
            "-max_target_seqs", "100000000"
        ]

        print(f"Running tBLASTn with command: {' '.join(tblastn_command)}")

        try:
            result = subprocess.run(tblastn_command, check=True, capture_output=True, text=True)
            print("tBLASTn search completed successfully.")
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e.stderr}")

    def _split_into_codons(self, nucleotidesToCheckForStart):
        # Calculate remainder to find extra nucleotides
        remainder = len(nucleotidesToCheckForStart) % 3

        # Create the first substring with the extra nucleotides (if any)
        if remainder != 0:
            codonsToCheck = [nucleotidesToCheckForStart[:remainder]]
            start_index = remainder
        else:
            codonsToCheck = []
            start_index = 0

        # Append the remaining substrings of length 3
        codonsToCheck += [nucleotidesToCheckForStart[i:i + 3] for i in range(start_index, len(nucleotidesToCheckForStart), 3)]

        return codonsToCheck

    def _checkForValidStart(self, codonsToCheck, startCodons):
        validCodons = {}

        for codonNum, codon in enumerate(codonsToCheck):
            if codon in startCodons:
                if codon in validCodons:
                    validCodons[codon].append(codonNum)
                else:
                    validCodons[codon] = [codonNum]
        print(codonsToCheck)
        print(startCodons)
        print(validCodons)
        return validCodons


    def _reverse_complement(self, dna_sequence):
        # Nucleotide Dict
        complement = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
            'N': 'N'  # If N in dict
        }

        # Generate the complement for each nucleotide
        complimentSeq = [complement[base] for base in dna_sequence]

        # Reverse compliment
        reverseCompSeq = complimentSeq[::-1]
        return ''.join(reverseCompSeq)

    def _writeUTRs_to_fasta(self, sequence_dict, fasta_filename):
        with open(fasta_filename, 'w') as fasta_file:
            for header, sequenceAndPos in sequence_dict.items():
                # Write the header and fix Coordinates to be UTR sequence position
                headerComps = header.split("|")
                locationComps = headerComps[-1].split(":")
                headerCoords = locationComps[-1]
                cleanedCoords = ""
                for char in str(headerCoords):
                    if char.isdigit() or char == "-":
                        cleanedCoords += char

                startCord, endCord = [int(coord) for coord in cleanedCoords.split("-")]
                if not sequenceAndPos[2]:
                    endCord = startCord + (int(sequenceAndPos[1]) - 1)
                else:
                    startCord = endCord - (int(sequenceAndPos[1]) - 1)
                newPos = f"{locationComps[0]}:{startCord}-{endCord}"
                newHeader = "|".join(headerComps[:-1])
                newHeader = newHeader + "|" + newPos



                fasta_file.write(f'>{newHeader}_&{sequenceAndPos[1]}\n')
                # Write the sequence, wrapping lines to a max width of 80 characters
                for i in range(0, len(sequenceAndPos[0]), 80):
                    fasta_file.write(sequenceAndPos[0][i:i+80] + '\n')

        print(f"FASTA file '{fasta_filename}' has been created.")

    def _get_full_upstream_sequence(self, sequence, s_start, s_end):

        if s_start < s_end:
            # Forward alignment: Take everything upstream of s_start
            upstream_sequence = sequence[:s_start - 1]  # s_start is 1-based, so subtract 1 for 0-based indexing
            is_reverse_complement = False
        else:
            # Reverse alignment: Take everything downstream of s_start, then reverse complement
            upstream_sequence = sequence[s_start:]  # s_start is the beginning of the alignment
            upstream_sequence = self._reverse_complement(upstream_sequence)
            is_reverse_complement = True

        return upstream_sequence, is_reverse_complement

    def findUTR_TEs(self, utrLengthThreshold=900, startCodonDistance=50,
                    startCodons=[], outputFasta="5Prime_UTRs.fa"):
        #Load BLAST Results as Dataframe
        blastColNames = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                         "gap_opens", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"]
        print("Reading in Blast Results")
        blastResultFrame = pd.read_csv(self.blastOutFilePath, sep="\t")
        blastResultFrame.columns = blastColNames




        sequenceNamesToKeep = list(set(blastResultFrame["subject_id"]))


        #Get upstream nucleotides of match for each sequence from fasta reference
        fastaDict = {}
        print("Getting sequences for BLAST Results")
        with open(self.fastaRefFile) as refFasta:
            currentSequenceName = ""
            sequenceNum = 0

            for line in refFasta:


                if ">" in line:
                    sequenceName = line.strip(">").strip()

                    if sequenceName in sequenceNamesToKeep:

                        fastaDict[sequenceName] = ""
                        currentSequenceName = sequenceNum
                    else:
                        currentSequenceName = ""

                    sequenceNum += 1

                    if sequenceNum % 100000 == 0:
                        #blastResultFrame = blastResultFrame[blastResultFrame['subject_id'].isin(fastaDict.keys())]
                        #break
                        print(f"Looked through {sequenceNum} sequences so far in the reference File")
                elif currentSequenceName:
                    fastaDict[sequenceName] = line.strip()


        if startCodons:
            print("Checking for valid start codons")
            #Keep track of codon statistics
            codonCountDict = {codon: 0 for codon in startCodons}
            codonCountDict["No_Codon"] = 0
            #Check for valid Start Codons in the correct frame

        validUTRRows = {}
        numRows = len(list(blastResultFrame.index))
        for rowIndex, resultRow in blastResultFrame.iterrows():
            sequenceName = resultRow["subject_id"]
            sequence = fastaDict[sequenceName]
            sequenceStart = resultRow["s_start"]
            sequenceEnd = resultRow["s_end"]
            upstreamRegion, isReverse = self._get_full_upstream_sequence(sequence, sequenceStart, sequenceEnd)





            if startCodons:
                nucleotidesToCheckForStart = upstreamRegion[-startCodonDistance:]
                codonsToCheck = self._split_into_codons(nucleotidesToCheckForStart)
                validCodons = self._checkForValidStart(codonsToCheck, startCodons)
                bestStartCodonPosition = None

                if not validCodons:
                    codonCountDict["No_Codon"] += 1
                    continue

                #Choose the best start codon position with preference to ATG
                for codon, codonPositions in validCodons.items():
                    codonCountDict[codon] += 1
                    closestStartCodon = max(codonPositions)

                    if codon == "ATG":
                        bestStartCodonPosition = closestStartCodon
                        break
                    else:
                        if (not bestStartCodonPosition) or (closestStartCodon > bestStartCodonPosition):
                            bestStartCodonPosition = closestStartCodon



                #Convert start codon position back to read position
                bestStartCodonPositionInRead = (len(upstreamRegion) - startCodonDistance) + sum([len(codon) for codon in codonsToCheck[:bestStartCodonPosition]])

                #Ensure that the potential UTR still meets size threshold
                potentialUTR = upstreamRegion[:bestStartCodonPositionInRead]
            else:

                potentialUTR = upstreamRegion
                bestStartCodonPositionInRead = len(potentialUTR)

            if len(potentialUTR) > utrLengthThreshold:
                #If multiple hits with the same sequence name have a valid UTR, choose the one with the longest UTR
                if sequenceName in validUTRRows:
                    if len(potentialUTR) > len(validUTRRows[sequenceName][0]):
                        if isReverse:
                            potentialUTR = self._reverse_complement(potentialUTR)
                        validUTRRows[sequenceName] = [potentialUTR, bestStartCodonPositionInRead, isReverse]
                else:
                    if isReverse:
                        potentialUTR = self._reverse_complement(potentialUTR)
                    validUTRRows[sequenceName] = [potentialUTR, bestStartCodonPositionInRead, isReverse]

        print("Best UTRs have been choosen")
        print()
        if startCodons:
            for codonType, codonNum in codonCountDict.items():
                print(f"Among the TE elements: {codonType} appreared roughly {codonNum/numRows}% of the time (out of all BLAST hits not sequence IDs)")
            print("Note: This is just a very rough estimate since some UTRs may not have been long enough after accounting for the stop codon and were discarded")
            print()

        print(f"{len(validUTRRows)} valid UTRs were found from {len(sequenceNamesToKeep)} sequences with valid BLAST hits {len(validUTRRows)/len(sequenceNamesToKeep)}% of the sequences had a valid UTR")
        if startCodons:
            print(f"All valid UTRs had a start codon {startCodons} within {startCodonDistance} nucleotides of the beginning of the BLAST hit")
        print(f"And all UTRs have a length of at least {utrLengthThreshold} nucleotides")
        print()
        print()
        print("Writing 5' UTRs to FASTA file")
        self._writeUTRs_to_fasta(validUTRRows, outputFasta)
        print("Process Completed Successfully")



if __name__ == "__main__":
    import sys
    TE_Family = "LINE/L1"
    pathToTERef = "TE_Sequence_Reference.tsv"
    pathToFirstProteinSequence = "ORF1.fasta"

    pathToTERef = "/Users/loganrivera/PycharmProjects/TE_neighborAnalysis/TE_ROGERS/IAmATestOut.bed.fasta"
    utrFinder = UTR_Finder(TE_Family=TE_Family,
                           pathToFirstProteinSequence=pathToFirstProteinSequence,
                           pathToTERef=pathToTERef)

    #Create the BLAST DataBase
    utrFinder.createBLASTDB()

    sys.exit()
    # Perform the BLAST search
    utrFinder.performBLASTSearch()

    #alternateStartCodons = ["ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"]
    # alternateStartCodons = ["ATG"]

    outputUTRFile = f"LINE_L1_5PrimeUTRs.fa"
    #Get 5' UTR Regions

    utrFinder.findUTR_TEs(outputFasta=outputUTRFile,
                          startCodons=[])

