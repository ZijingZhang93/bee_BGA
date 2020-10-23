#=====================================================================================
#  Load package
#=====================================================================================
library("clusterProfiler")
library('grid')
library('readr')
library('ggplot2')

setwd("C:\\Users\\apple\\Desktop\\rnaseq_brain_zzj\\5KO analysis差异基因\\1单菌-富集分析-level3")

#Import annotation file
name_term_gene <- read_csv("./name_term_gene_3.csv",col_types = cols(term = col_character()))
term2gene <- name_term_gene[,c(2,3)]
term2name <- name_term_gene[,c(2,1)]

#file name, the file including gene list
sp_up <- c(paste(c('Ba','Bi','F5','Gi','Sn','CV','F4'),rep('up',7),sep = '_'))
sp_down <- c(paste(c('Ba','Bi','F5','Gi','Sn','CV','F4'),rep('down',7),sep = '_'))
len <- length(sp_up)

png('./kegg_function_1.png',width = 10*len*100,height = 8*100,res = 120)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,7)))
for (i0 in 1:len){
  #import gene list
  gene <- read.delim(paste0('./',sp_up[i0],'.txt'),header = FALSE)
  gene <- gene$V1
  #functional enrichment analysis
  x <- enricher(gene,TERM2GENE = term2gene,TERM2NAME = term2name,
                pvalueCutoff = 1,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
  write.csv(x,paste0('./result_csv/',sp_up[i0],'.csv'))
  #d <- heatplot(x)+ggtitle(sp_up[i0])
  d <- dotplot(x,color = 'pvalue',title = sp_up[i0])
  print(d,vp = vplayout(1,i0))
}
for (i0 in 1:(len)){
  #import gene list
  gene <- read.delim(paste0('./',sp_down[i0],'.txt'),header = FALSE)
  gene <- gene$V1
  #functional enrichment analysis
  x <- enricher(gene,TERM2GENE = term2gene,TERM2NAME = term2name,
                pvalueCutoff = 1,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
  write.csv(x,paste0('./result_csv/',sp_down[i0],'.csv'))
  #d <- heatplot(x)+ggtitle(sp_down[i0])
  d <- dotplot(x,color = 'pvalue',title = sp_down[i0])
  print(d,vp = vplayout(2,i0))
}
dev.off()

#=====================================================================================
#  function
#=====================================================================================
#vplayout
vplayout <- function(x,y){
  viewport(layout.pos.row=x,layout.pos.col=y)
}
