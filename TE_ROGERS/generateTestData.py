import random
import pandas as pd
import os


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)





class TestDataGenerator():

    def __init__(self, numGeneIDs, numTEIDs, pathToTEneighborhoodFile):
        self.numGeneIDs = numGeneIDs
        self.numTEIDs = numTEIDs
        self.pathToTENeighborhoodFile = pathToTEneighborhoodFile

        self.geneNeighborFrame = pd.read_csv(self.pathToTENeighborhoodFile, sep="\t")
        self.geneIDList = self.geneNeighborFrame["Gene_ID"].unique()
        self.TEIDList = self.geneNeighborFrame["TE_Number"].unique()

        self.geneIDToUse = self.geneIDList[:self.numGeneIDs]
        self.TEIDToUse = self.TEIDList[:self.numTEIDs]

    def _createTestDataFile(self, featureNames, outputPath, randomBarcodes, featureRange=(0, 1)):

        barcodeFeatureDict = {}
        for barcode in randomBarcodes:
            barcodeFeatureDict[barcode] = []
            for _ in range(len(featureNames)):
                featureVal = random.uniform(featureRange[0], featureRange[1])
                barcodeFeatureDict[barcode].append(featureVal)

        testDataFrame = pd.DataFrame(barcodeFeatureDict)
        testDataFrame.index = featureNames

        testDataFrame.to_csv(outputPath, sep="\t")

    def generateAllTestData(self, barcodeLength, numBarcodes,  outputDir, featureRanges=((-2, 2), (0, 1), (-2, 2), (0, 1))):
        if not os.path.exists(outputDir):
            # Create the directory
            os.makedirs(outputDir)
            print(f"Directory '{outputDir}' was created.")
        else:
            print(f"Directory '{outputDir}' already exists.")

        randomBarcodes = testDataGenerator._generateRandomBarcodes(barcodeLength=barcodeLength, numBarcodes=numBarcodes)
        print("Creating Gene fold data")
        self._createTestDataFile(self.geneIDToUse, f"{outputDir}/Gene_Fold", randomBarcodes,
                                             featureRange=featureRanges[0])
        print("Creating Gene P-val data")
        self._createTestDataFile(self.geneIDToUse, f"{outputDir}/Gene_Pval", randomBarcodes,
                                             featureRange=featureRanges[1])
        print("Creating TE fold data")
        self._createTestDataFile(self.TEIDToUse, f"{outputDir}/TE_Fold", randomBarcodes,
                                             featureRange=featureRanges[2])
        print("Creating TE P-val data")
        self._createTestDataFile(self.TEIDToUse, f"{outputDir}/TE_Pval", randomBarcodes,
                                             featureRange=featureRanges[3])
        print("All Test data created successfully")


    def _generateRandomBarcodes(self, numBarcodes, barcodeLength, possibleNucleotides = ("A", "T", "C", "G")):
        randomBarcodes = [''.join(random.choices(possibleNucleotides, k=barcodeLength)) for _ in range(numBarcodes)]
        return randomBarcodes



if __name__ == "__main__":
    testDataGenerator = TestDataGenerator(numGeneIDs=3000, numTEIDs=2000, pathToTEneighborhoodFile="TE_gene_neighborhood.tsv")

    testDataGenerator.generateAllTestData(barcodeLength=16, numBarcodes=3000, outputDir='testData42Dir')
