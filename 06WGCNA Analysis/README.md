#=====================================================================================
#  01-dataInput.R
#=====================================================================================
This is the ﬁrst step of any network analysis.
1)We load typical expression data
2)Pre-process them into a format suitable for network analysis
3)Clean the data by removing obvious outlier samples as well as metabolites and samples with excessive numbers of missing entries

#=====================================================================================
#  02-networkConstr-man.R
#=====================================================================================
The construction of the metabolites network and identiﬁcation of modules was used Step-by-step strategy
1)We first Choose the soft-thresholding power
2)The co-expression similarity and adjacency was calculated
3)To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity
4)Clustering using TOM
5)Merging of modules whose expression proﬁles are very similar

#=====================================================================================
#  03-relateModsToExt.R
#=====================================================================================
In this script, we do the following things:
1)Quantifying module–trait associations
2)Gene relationship to trait and important modules: Gene Signiﬁcance and Module Membership
3)Intramodular analysis: identifying genes with high GS and MM
4)Summary output of network analysis results

#=====================================================================================
#  04-ExportNetwork.R
#=====================================================================================
This script is used to export a Cytoscape format

