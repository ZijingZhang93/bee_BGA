#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/apple/Desktop/rnaseq_brain_zzj/8WGCNA/data/";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
datExpr <- data.frame(datExpr)
# Load network data saved in the second part.
lnames = load(file = "./02-networkConstruction-stepByStep.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 5);
# Read in the annotation file
annot = read.csv(file = "MetaboliteAnnotation.csv");
# Select module probes
probes = names(data.frame(datExpr))

# Select modules
modNames = substring(names(MEs), 3)
for (modules in modNames) {
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modGenes = annot$Name_des[match(modProbes, annot$name)];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.01,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
  
  # Export the attribute files of aim gene
  modClass = annot$Class.English.[match(modProbes, annot$name)]
  modAbund = annot$abundance[match(modProbes, annot$name)]
  modAttr = data.frame(modProbes,modGenes,modClass,modAbund)
  write.table(modAttr,sep = '\t',file = paste("CytoscapeInput-attr-", paste(modules, collapse="-"), ".txt", sep=""),row.names = F,quote = F)
}
