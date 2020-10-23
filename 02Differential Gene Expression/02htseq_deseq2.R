#=====================================================================================
#  Load package
#=====================================================================================
library(DESeq2)
library(ggplot2)
library(plotly)
library(grid)

#=====================================================================================
#  Import data
#=====================================================================================
path <- 'C:\\Users\\apple\\Desktop\\rnaseq_brain_zzj\\diff\\old_htseq_gff-padj'
count_data <- read.table(paste0(path,"\\merge_htcount.txt"),header = T, sep = "\t", row.names = 1)
colnames(count_data) <- c("Ba1", "Ba2", "Ba3", 
                          "Bi1", "Bi2", "Bi3", 
                          "CV1", "CV2", "CV3",
                          "F41", "F42", "F43", 
                          "F51", "F52", "F53", 
                          "GF1", "GF2", "GF3",
                          "Gi1", "Gi2", "Gi3",
                          "Sn1", "Sn2", "Sn3")

count_data <- as.matrix(count_data)
count_data[is.na(count_data)] <- 0

#=====================================================================================
#  DESeq2 analysis
#=====================================================================================
condition <- factor(c("Ba", "Ba", "Ba", 
                      "Bi", "Bi", "Bi",
                      "CV", "CV", "CV",
                      "F4", "F4", "F4", 
                      "F5", "F5", "F5",
                      "GF", "GF", "GF",
                      "Gi", "Gi", "Gi", 
                      "Sn", "Sn", "Sn"))
dds <- DESeqDataSetFromMatrix(count_data, DataFrame(condition), design= ~condition )
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

#normlization matrix
rld <- rlogTransformation(dds)
expr_data <- assay(rld) 
write.csv(expr_data,paste0(path,'\\','expr.csv'))

all <- unique(as.character(condition))
control <- 'GF'
test <- all[-which(all == control)]

#=====================================================================================
#  Drawing volcano map
#=====================================================================================
pdf(paste0(path,'\\','adj_pic.pdf'),width = 7,height = 5*length(test))
grid.newpage()
pushViewport(viewport(layout=grid.layout(length(test),1)))
for(i0 in 1:length(test)){
  print(i0)
  res <- results(dds,contrast=c("condition",test[i0],control))
  
  diff_stat <- as.data.frame(res)
  diff_stat[which(diff_stat$padj < 0.05 & diff_stat$log2FoldChange >= 1),'diff'] <- 'up'
  diff_stat[which(diff_stat$padj < 0.05 & diff_stat$log2FoldChange <= -1),'diff'] <- 'down'
  diff_stat[!(diff_stat$diff %in% c('up', 'down')),'diff'] <- 'no'
  diff_stat$id <- row.names(diff_stat)
  
  write.csv(diff_stat,paste0(path,'\\',test[i0],'-',control,'.csv'))
  
  p <- deseq2_plot(diff_stat)
  p <- p+ggtitle(paste0(test[i0],'-',control))
  print(p,vp = vplayout(i0,1))
  p <- ggplotly(p)
  htmlwidgets::saveWidget(as.widget(p), paste0(path,'\\',test[i0],'-',control,'.html'))
}
dev.off()

#=====================================================================================
#  function
#=====================================================================================
deseq2_plot <- function(diff_stat){
  mytheme <- 
    theme_bw()+
    theme(legend.text = element_text(color = 'gray40',size = rel(0.9)),
          legend.title = element_text(color = 'gray20',size = rel(1)))+
    theme(panel.border = element_rect(color='gray20'),
          panel.grid.major.y = element_blank(), #去掉横的线
          panel.grid.minor.y = element_blank(), #去掉横的线
          panel.grid.major.x = element_blank(), #去掉竖线
          panel.grid.minor.x = element_blank(), #去掉竖线
          panel.background = element_rect(color = 'gray',fill = 'transparent'))
  
  ggplot(diff_stat, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(text = id,color = diff), size = 1, alpha = 0.6) +
    scale_colour_manual(limits = c('up', 'down', 'no'), values = c('red', 'blue', 'gray40')) +
    labs(x = 'log2 Fold Change', y = '-log10 padj-value') +
    geom_vline(xintercept = c(-1, 1), color = 'gray70', size = 0.5) + 
    geom_hline(yintercept = -log10(0.05), color = 'gray70', size = 0.5)+
    mytheme
}

vplayout <- function(x,y){
  viewport(layout.pos.row=x,layout.pos.col=y)
}



