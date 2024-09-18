import seaborn
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd


pd.set_option('display.max_columns', None)
class TE_Gene_ResultVisualizer():
    def __init__(self, pathToFinalSummaryFile, pathToClosestDistanceFile):
        self.pathToFinalSummaryFile = pathToFinalSummaryFile
        self.pathToClosestDistanceFile = pathToClosestDistanceFile




    def createClosestDistancePlot(self, pvalFilter):
        closestDistFrame = pd.read_csv(self.pathToClosestDistanceFile, sep="\t")
        columnsToKeep = ["Gene_fold_change", "Gene_p_value", "Distance_TE_Gene", "Unique_Gene_Count"]
        closestDistFrame = closestDistFrame[columnsToKeep]
        closestDistFrame = closestDistFrame[closestDistFrame["Gene_p_value"] < pvalFilter]
        seaborn.scatterplot(closestDistFrame, x="Distance_TE_Gene", y="Gene_fold_change")
        plt.show()

    def createFoldCorrelationPlot(self, genePValFilter=0.01, TEPValFilter=0.01):
        mainSummaryFile = pd.read_csv(self.pathToFinalSummaryFile, sep="\t")
        columnsToKeep = ["TE_fold_change", "TE_p_value", "Gene_fold_change", "Gene_p_value", "Distance_TE_Gene"]
        mainSummaryFile = mainSummaryFile[columnsToKeep]
        mainSummaryFile["Distance_TE_Gene"] = (mainSummaryFile["Distance_TE_Gene"].abs() / 1000)
        mainSummaryFile = mainSummaryFile[mainSummaryFile["Gene_p_value"] < genePValFilter]
        mainSummaryFile = mainSummaryFile[mainSummaryFile["TE_p_value"] < TEPValFilter]
        seaborn.scatterplot(mainSummaryFile, x="Gene_fold_change", y="TE_fold_change", size="Distance_TE_Gene",
                            sizes=(20, 200))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Distance (kb)")
    
    # Show the plot
        plt.show()









if __name__ == "__main__":
    pathToClostestDistancefile = "closestDistance_neighborhoodSummaryFileForPlotting.tsv"
    pathToFinalSummaryFile = "neighborhoodSummaryFile_ForPlotting.tsv"

    visualizer = TE_Gene_ResultVisualizer(pathToFinalSummaryFile=pathToFinalSummaryFile,
                                          pathToClosestDistanceFile=pathToClostestDistancefile)

    visualizer.createClosestDistancePlot(pvalFilter=0.001)
    visualizer.createFoldCorrelationPlot()
