import pandas as pd
import numpy as np
import time
import os



pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

class TEGeneProcessor:

    def __init__(self, TEFoldPath, TEPvalPath, GeneFoldPath, GenePvalPath, TE_gene_neighborhoodPath):
        self.TEFoldPath = TEFoldPath
        self.TEPvalPath = TEPvalPath
        self.GeneFoldPath = GeneFoldPath
        self.GenePvalPath = GenePvalPath

        self.TE_gene_neighborhoodPath = TE_gene_neighborhoodPath

        self.TEFoldFrame = None
        self.TEPvalFrame = None
        self.GeneFoldFrame = None
        self.GenePvalFrame = None

        self.neighborhoodFrame = None

        self.TECombinedFrame = None
        self.GeneCombinedFrame = None





    def _loadDataFrames(self):
        print("Loading and initializing reference dataframes")
        self.TEFoldFrame = pd.read_csv(self.TEFoldPath, sep="\t", index_col=0).T
        self.TEPvalFrame = pd.read_csv(self.TEPvalPath, sep="\t", index_col=0).T
        self.GeneFoldFrame = pd.read_csv(self.GeneFoldPath, sep="\t", index_col=0).T
        self.GenePvalFrame = pd.read_csv(self.GenePvalPath, sep="\t", index_col=0).T

        self.neighborhoodFrame = pd.read_csv(self.TE_gene_neighborhoodPath, sep="\t")


    def filterTEsAndGenes(self, TEFoldCutoff, TEPvalCutoff, GeneFoldCutoff, GenePvalCutoff):
        self._loadDataFrames()
        print("beginning filtering of TEs and genes")
        dataFrameToFilter = [self.TEFoldFrame, self.TEPvalFrame, self.GeneFoldFrame, self.GenePvalFrame]
        cutoffs = [TEFoldCutoff, TEPvalCutoff, GeneFoldCutoff, GenePvalCutoff]
        greaterOrLess = ["<", ">", "<", ">"] # direction to mark values to filter
        filteredFrames = []

        for currentFrame, cutoff, direction in zip(dataFrameToFilter, cutoffs, greaterOrLess):

            if direction == "<":
                currentFrame = currentFrame.mask(abs(currentFrame) < cutoff, other=np.nan)
            elif direction == ">":
                currentFrame = currentFrame.mask(abs(currentFrame) > cutoff, other=np.nan)

            filteredFrames.append(currentFrame)

        TEFoldFrame, TEPvalFrame, GeneFoldFrame, GenePvalFrame = filteredFrames
        TECombinedMask = TEFoldFrame.isna() | TEPvalFrame.isna()
        TECombinedFrame = pd.DataFrame(1, index=TEFoldFrame.index, columns=TEFoldFrame.columns)
        TECombinedFrame[TECombinedMask] = np.nan

        GeneCombinedMask = GeneFoldFrame.isna() | GenePvalFrame.isna()
        GeneCombinedFrame = pd.DataFrame(1, index=GeneFoldFrame.index, columns=GeneFoldFrame.columns)
        GeneCombinedFrame[GeneCombinedMask] = np.nan

        print("finished filtering of TEs and genes")

        self.TECombinedFrame = TECombinedFrame
        self.GeneCombinedFrame = GeneCombinedFrame

    def createSummaryFile(self, outFileName):
        print("Beginning to create summary file!")
        start_time = time.time()
        barcodes = list(self.TEFoldFrame.index)
        

        te_combined = self.TEFoldFrame.astype(str) + '@' + self.TEPvalFrame.astype(str)
        te_combined.columns = self.TEFoldFrame.columns  # Preserve TE identifiers as columns
        te_combined['barcode'] = te_combined.index  # Add barcode as a column


        # Combine Gene fold changes and p-values into one DataFrame with '@' separator
        gene_combined = self.GeneFoldFrame.astype(str) + '@' + self.GenePvalFrame.astype(str)
        gene_combined.columns = self.GeneFoldFrame.columns # Preserve gene identifiers as columns
        gene_combined['barcode'] = gene_combined.index  # Add barcode as a column

        results = []
        for barcode in barcodes:
            # Determine valid TEs and genes for the current barcode
            valid_TE = self.TECombinedFrame.loc[barcode].dropna().index
            valid_genes = self.GeneCombinedFrame.loc[barcode].dropna().index
            

            # Filter neighborhoodFrame based on valid TEs and genes for the current barcode
            filtered_neighborhood = self.neighborhoodFrame[
                (self.neighborhoodFrame['TE_Number'].isin(valid_TE)) &
                (self.neighborhoodFrame['Gene_name'].isin(valid_genes))
                ]

            if filtered_neighborhood.empty:
                print("we empty")
                continue

            te_data = te_combined.loc[barcode, valid_TE].reset_index().rename(
                columns={'index': 'TE_Number', barcode: 'TE_combined'})
            gene_data = gene_combined.loc[barcode, valid_genes].reset_index().rename(
                columns={'index': 'Gene_name', barcode: 'Gene_combined'})


            # Merge with filtered_neighborhood on TE and Gene IDs
            
            merged_data = pd.merge(filtered_neighborhood, te_data, on='TE_Number', how='inner')
            merged_data = pd.merge(merged_data, gene_data, on='Gene_name', how='inner')

            merged_data['barcode'] = barcode

            results.extend(merged_data.to_dict('records'))

        summary_df = pd.DataFrame(results)

        summary_df[['TE_fold_change', 'TE_p_value']] = summary_df['TE_combined'].str.split('@', expand=True)
        summary_df[['Gene_fold_change', 'Gene_p_value']] = summary_df['Gene_combined'].str.split('@', expand=True)
        summary_df.drop(columns=['TE_combined', 'Gene_combined'], inplace=True)
        try:
            column_order = ["barcode", "TE_chr", "TE_start", "TE_stop",
                            "TE_strand", "TE_Class", "TE_Number", "TE_fold_change",
                            "TE_p_value", "Gene_ID", "Gene_name", "Gene_strand",
                            "Gene_class", "Gene_fold_change", "Gene_p_value",
                            "Distance_TE_Gene", "Unique_Gene_Count", "Unnamed: 0"]
            summary_df = summary_df[column_order]
        except:
            pass

        summary_df.to_csv(outFileName, sep="\t")
        end_time = time.time()  # End the timer
        print(f"Function took {end_time - start_time:.2f} seconds to complete.")

    def getClosestTEDistanceToGene(self, pathToTE_SummaryFile, outputPath):
        print("Began Looking for TEs closest to each gene")

        # Load the data
        summaryFrame = pd.read_csv(pathToTE_SummaryFile, sep="\t", index_col=0)

        # Create a column with the absolute distance
        summaryFrame['Absolute_Distance'] = summaryFrame['Distance_TE_Gene'].abs()

        # Find the minimum absolute distance for each barcode and gene combination
        closest_TE = summaryFrame.groupby(['barcode', 'Gene_ID'])['Absolute_Distance'].min().reset_index()

        # Merge to get rows corresponding to these minimum distances
        closestDistanceFrame = pd.merge(
            summaryFrame,
            closest_TE,
            left_on=['barcode', 'Gene_ID', 'Absolute_Distance'],
            right_on=['barcode', 'Gene_ID', 'Absolute_Distance']
        )

        # Drop the temporary 'Absolute_Distance' column as it's no longer needed
        closestDistanceFrame = closestDistanceFrame.drop(columns=['Absolute_Distance'])

        # Remove duplicates (if any)
        closestDistanceFrame = closestDistanceFrame.drop_duplicates(subset=['barcode', 'Gene_ID', 'Distance_TE_Gene'])

        # Sort by barcode for better readability
        closestDistanceFrame = closestDistanceFrame.sort_values(by='barcode')

        # Save the final output to a file
        closestDistanceFrame.to_csv(outputPath, sep="\t", index=False)

        print("Found TEs closest to each gene")




if __name__ == "__main__":
    pathToTestData = "testSingleCellData"
    TEFoldPath = os.path.join(pathToTestData, "TE_Fold")
    TEPvalPath = os.path.join(pathToTestData, "TE_Pval")
    GeneFoldPath = os.path.join(pathToTestData, "Gene_Fold")
    GenePvalPath = os.path.join(pathToTestData, "Gene_Pval")



    TE_gene_neighborhoodPath = "TE_gene_neighborhood.tsv"

    teGeneProcessor = TEGeneProcessor(TEFoldPath,
                                      TEPvalPath,
                                      GeneFoldPath,
                                      GenePvalPath,
                                      TE_gene_neighborhoodPath)



    teGeneProcessor.filterTEsAndGenes(TEFoldCutoff=1.2, TEPvalCutoff=0.05, GeneFoldCutoff=0, GenePvalCutoff=0.05)

    outFileName = "neighborhoodSummaryFile_ForPlotting.tsv"
    teGeneProcessor.createSummaryFile(outFileName)

    outFileNameDistanceFile = "closestDistance_neighborhoodSummaryFileForPlotting.tsv"
    teGeneProcessor.getClosestTEDistanceToGene(outFileName, outFileNameDistanceFile)
