#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

rm(list=ls())
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/apple/Desktop/rnaseq_brain_zzj/8WGCNA/data/";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
library(ggplot2)
library(ggrepel)
library(grid)
library(RColorBrewer)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "02-networkConstruction-stepByStep.RData");
lnames

datExpr <- data.frame(datExpr)
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
color_use = rev(colorRampPalette(brewer.pal(10, "RdBu"))(20))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf('../Plots/04Module-trait_relationships_blue.pdf',width = 6,height = 5)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = color_use,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# names of the traits
traitNames = names(datTraits)
# names (colors) of the modules
modNames = substring(names(MEs), 3)
# names of the gene
annot = read.csv(file = "MetaboliteAnnotation.csv");
probes = names(datExpr)
probes2annot = match(probes, annot$name)
geneName = annot$Name_des[probes2annot]
geneNum = strsplit(probes,'_')
geneNum = sapply(geneNum,function(x){return(x[2])})

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", traitNames, sep="");
names(GSPvalue) = paste("p.GS.", traitNames, sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
#Using the GS and MM measures, we can identify genes that have a high signiﬁcance for weight 
#as well as high module membership in interesting modules. 
vplayout <- function(x,y){
  viewport(layout.pos.row=x,layout.pos.col=y)
}
pdf('../Plots/04_3MM-GS_annot.pdf',width = 5*8,height = 5*9)
grid.newpage()
pushViewport(viewport(layout=grid.layout(9,8)))
for (module in modNames) {
  for (trait in traitNames) {
    column = match(module, modNames);
    columnTrait = match(trait, traitNames);
    moduleGenes = (moduleColors==module);
    gene = geneName[moduleGenes]
    num = geneNum[moduleGenes]
    
    x = abs(geneModuleMembership[moduleGenes, column])
    y = abs(geneTraitSignificance[moduleGenes, columnTrait])
    data = data.frame(x,y,gene,num)
    textAnnote = data[data$x>0.8,]
    
    fit = lm(y~x)
    corr = cor(y,x,use = "p")
    p = corPvalueStudent(as.matrix(corr), sum(moduleGenes)+1)
    if(p < 0.0001){
      p = format(p,digits = 2,scientific = T)
    }else{
      p = signif(p,2)
    }
    
    pic<- ggplot(data = data,aes(x = x,y = y))+
      geom_smooth(method = 'lm',formula = y~x,color = 'gray',se = F,size = 0.5)+
      geom_vline(xintercept = 0.8,color = 'gray')+
      geom_point(color = module,size = 2,alpha = 0.8)+
      xlab(paste("Module Membership in", module, "module"))+
      ylab(paste("Gene significance for",trait))+
      ggtitle(paste("Module membership vs. gene significance\ncor =",round(corr,2),', p =',p))+
      geom_text_repel(data = textAnnote,aes(x,y,label = gene),size = 1.2,segment.size = 0.1,force = T)+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5))
    
    print(pic,vp = vplayout(column,columnTrait))
  }
}
dev.off()

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

datExpr = as.data.frame(datExpr)
names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


names(datExpr)[moduleColors=="brown"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = "MetaboliteAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(name = probes,
                       Name_des = annot$Name_des[probes2annot],
                       Chinese.name = annot$Chinese.name[probes2annot],
                       Formula = annot$Formula[probes2annot],
                       Molecular.Weight = annot$Molecular.Weight[probes2annot],
                       RT.min = annot$RT..min.[probes2annot],
                       Class.English = annot$Class.English.[probes2annot],
                       Class.Chinese = annot$Class.Chinese.[probes2annot],
                       abundance = annot$abundance[probes2annot],
                       moduleColor = moduleColors)

# names (colors) of the modules
modNames = substring(names(MEs), 3)
for (module in modNames) {
  column = match(module, modNames);
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, column], 
                         MMPvalue[, column]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[column], sep=""),
                       paste("p.MM.", modNames[column], sep=""))
  }
# Order the genes in the geneInfo variable first by module color
geneOrder = order(geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "../Plots/04MetaboliteInfo_MM.csv",row.names = FALSE)
#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================
# Create the starting data frame
geneInfo0 = data.frame(name = probes,
                       Name_des = annot$Name_des[probes2annot],
                       Chinese.name = annot$Chinese.name[probes2annot],
                       Formula = annot$Formula[probes2annot],
                       Molecular.Weight = annot$Molecular.Weight[probes2annot],
                       RT.min = annot$RT..min.[probes2annot],
                       Class.English = annot$Class.English.[probes2annot],
                       Class.Chinese = annot$Class.Chinese.[probes2annot],
                       abundance = annot$abundance[probes2annot],
                       moduleColor = moduleColors)

# names of the traits
traitNames = names(datTraits)
for (trait in traitNames) {
  column = match(trait, traitNames);
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[, column], 
                         MMPvalue[, column]);
  names(geneInfo0) = c(oldNames, paste("GS.", traitNames[column], sep=""),
                       paste("p.GS.", traitNames[column], sep=""))
}
# Order the genes in the geneInfo variable first by module color
geneOrder = order(geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "../Plots/04MetaboliteInfo_GS.csv",row.names = FALSE)
