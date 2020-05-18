#=====================================================================================
#                                                                [CMSC 636] Drake Wong
#  WGCNA Script
#
#=====================================================================================
library(WGCNA)
options(stringsAsFactors = FALSE)


# Step 1: Read in data:
D = read.csv("concated_small_data.csv", check.names=FALSE)
dim(D)
names(D)
datExpr0 = as.data.frame((D[,-1]))


# Check for excessive missing values and outlier samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#Cluster the samples
sampleTree = hclust(dist(datExpr0), method = "average")


#=====================================================================================
#
#  Pick Power
#
#=====================================================================================
powers = c(1:10)
# bicor is for robust correlation, Pearson is for regular correl. For practical purposes "signed hybrid" vs. "signed" are almost identical, some prefer "signed hybrid" as this provides correlation that's positive and otherwise zero (it is called hybrid because it is a hybrid of weighted and unweighted networks). The similarity is raised to a power which is usually twice as big in signed networks than it is in hybrid networks. For correlations near 1, ( (1+cor)/2 )^(2*beta) approximately equals cor^beta, and for low correlations both adjacencies are near zero.
sft = pickSoftThreshold(datExpr0,corFnc = "bicor", corOptions=list(maxPOutliers=0.1),networkType = "signed hybrid", powerVector = powers, verbose = 5)
# Plot the results as scale-free topology fit index as a function of the soft-thresholding power
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.6
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="l",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="l",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="darkgreen")


#Save plot
dev.print(pdf, 'Report Main: Bicor Power Plot.pdf')

softPower = 8


#=====================================================================================
#
#  Model
#
#=====================================================================================
# Adjustable Params
param_adjacency_type = c("unsigned", "signed")    #P1
param_minModuleSize = c(10,30, 50)                #P2
param_cutreeDynamic_method = c("tree", "hybrid")  #P3
param_cutreeDynamic_deepSplit = c(0, 4)           #P4
param_MEDissThres = c(0.2, 0.6)                   #P5

cat("Adjacency","MinModuleSize", "CutMethod", "DeepSplit", "MergeThres", "NumModules", "GreySize","\n", file="report.csv",sep=",",append=F)

for (p_1 in param_adjacency_type){
  for (p_2 in param_minModuleSize){
    for (p_3 in param_cutreeDynamic_method){
      for (p_4 in param_cutreeDynamic_deepSplit){
        for (p_5 in param_MEDissThres){
          adjacency = adjacency(datExpr0, power = softPower, type = p_1)
          
          # Turn adjacency into topological overlap matrix
          TOM = TOMsimilarity(adjacency, verbose = 5)
          dissTOM = 1-TOM
          
          # Call the hierarchical clustering function
          geneTree = hclust(as.dist(dissTOM), method = "average")
          
          # Set the minimum module size
          minModuleSize = p_2
          
          # Module identification using dynamic tree cut:
          dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                      deepSplit = p_4, method = p_3, pamStage = TRUE, pamRespectsDendro = TRUE, useMedoids = FALSE,
                                      respectSmallClusters = TRUE, minClusterSize = minModuleSize, verbose = 5)
          
          # Convert numeric labels into colors
          dynamicColors = labels2colors(dynamicMods)
          table(dynamicColors)
          
          # Calculate eigengenes
          MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
          MEs = MEList$eigengenes
          
          # Calculate dissimilarity of module eigengenes
          MEDiss = 1-cor(MEs)
          
          # Cluster module eigengenes
          METree = hclust(as.dist(MEDiss), method = "average")
          
          #Highcut is determined by the clustering of the modules eigengenes modules
          MEDissThres = p_5
          
          # Call an automatic merging function
          merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 5)
          
          # The merged module colors and eigengenes of the new merged modules
          mergedColors = merge$colors
          mergedMEs = merge$newMEs
          
          
          # Plot Dendrogram
          sizeGrWindow(8,6)
          plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), "\n\n\nUnmerged\n\n Merged",
                              dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05,
                            main = "Gene dendrogram and module colors")
          #Save plot
          dev.print(pdf, paste0('Report 2.pdf'))
          
          num_of_modules = nrow(table(mergedColors))
          num_within_grey = sum(mergedColors == "grey")
          
          cat(p_1,p_2,p_3,p_4,p_5,num_of_modules,num_within_grey,"\n",file="report.csv",sep=",",append=T)
        }
      }
    }
  }
}











