library(optparse,quietly = T)
library(yaml,quietly = T)        
library(tidyverse, quietly = T)
library(dplyr, quietly = T)
library(ggplot2,quietly = T)
library(ggridges,quietly = T)
library(ggrepel, quietly = T)

pval_gene <- function(gene_in,time_point_in)
{
  in_gene = gene_in
  tp = time_point_in 
  sel_gene <- survival %>% filter(gene == in_gene) %>% 
    filter(STUDY_ID %in% tx)
  sel_gene = sel_gene$STUDY_ID[which(sel_gene$STUDY_ID %in% tx)]
  
  sel_gene_ir = ir[sel_gene,sel_gene,]
  
  sel_gene_tbl = as.data.frame(timex_1)
  sel_gene_tbl$timex_1 = as.character(sel_gene_tbl$timex_1)
  
#  lapply(1:NROW(sel_gene_tbl), FUN = function(x) sum(ir[sel_gene,sel_gene,x], na.rm = T) ) %>% unlist -> sel_gene_tbl$sum 
#  lapply(1:NROW(sel_gene_tbl), FUN = function(x) mean(ir[sel_gene,sel_gene,x], na.rm = T) ) %>% unlist -> sel_gene_tbl$mean
#  lapply(1:NROW(sel_gene_tbl), FUN = function(x) median(ir[sel_gene,sel_gene,x], na.rm = T) ) %>% unlist -> sel_gene_tbl$median  
  
  sel_gene_tbl_1 = as.data.frame(timex_1)
  sel_gene_tbl_1$timex_1 = as.character(sel_gene_tbl_1$timex_1)
  
  arr_2_25_stx <- ir[sel_gene,sel_gene,tp]
  not_gene <- c(tx[!tx %in% sel_gene])
  arr_2_25_other <- ir[not_gene, not_gene,tp]
  
  gene_present <- reshape2::melt(arr_2_25_stx) %>% filter(!is.na(value))
  gene_absent <- reshape2::melt(arr_2_25_other) %>% filter(!is.na(value))


  m_present <- mean(gene_present$value)
  m_absent <- mean(gene_absent$value)  
  ret <- NA
  if (NROW(gene_present) > 0 & NROW(gene_absent) > 0)
  {
    resx <- wilcox.test(gene_present$value,gene_absent$value,alternative = "greater")
    resx$p.value -> ret[1]
  }
  ret[2] <- m_present
  ret[3] <-  m_absent
  ret %>% return

}



survival$gene %>% table %>% as.data.frame() %>% arrange(Freq %>% desc) -> gene_list
names(gene_list) <- c('gene','freq')
gene_list <- gene_list %>% filter(gene != 'NA' & freq > 1)

plotx_sum <- matrix(ncol=2,nrow=1) %>% as.data.frame()
names(plotx_sum) <- c('time_point','pval')
plotx_sum$gene <- NA
plotx_sum$log <- NA
plotx_sum <- plotx_sum[0,]

pvals <- seq(1,303,by=3)
pres <- seq(2,303,by=3)
abs <- seq(3,303,by=3) 

for(numx in 1:NROW(gene_list)) {
  genexx = gene_list$gene[numx] %>% as.character
  print(genexx)
  print('--------------------')
  
  plotx <- matrix(ncol=2,nrow=length(timex_1)) %>% as.data.frame()
  names(plotx) <- c('time_point','pval')
  plotx$time_point <- timex
  plotx$gene <- genexx


  setx = lapply(1:NROW(timex_1), FUN = function(x) pval_gene(genexx, timex_1[x]) ) %>% unlist
  plotx$pval = setx[c(pvals)]
   plotx$mean_present = setx[c(pres)]
   plotx$mean_absent = setx[c(abs)]

  plotx$log <- -log10(plotx$pval)
  plotx_sum <- rbind(plotx_sum,plotx)
  
}

plotx_sum$log_mod = plotx_sum$log - (-log10(0.05))
plotx_sum$log_mod[plotx_sum$log_mod < 0] <- 0
plotx_sum$log_mod[is.na(plotx_sum$log_mod)] <- 0

plotx_sum$log_mod[ which(is.infinite(plotx_sum$log_mod))] = 210
plotx_sum$mean_ratio = plotx_sum$mean_present/plotx_sum$mean_absent

write_csv(plotx_sum,paste0(input.yaml$outputDir,'/results/phensim_scores.csv'))

if(input.yaml$plot == T ){
pdf(paste0(input.yaml$outputDir,"/plots/ridge_plot.pdf"),width = 16, height = 9  )
print(ggplot(plotx_sum, aes(x=time_point, y=gene,
                       height = log_mod, fill=gene)) +
  geom_ridgeline(show.legend = F,size=0,alpha=0.3) +
  xlab("Age") +
  ylab('-log(p-value) if p-value is < 0.05') +
  ggtitle("Wilcoxon Test for significant similarities"))
  
dev.off()  
}



