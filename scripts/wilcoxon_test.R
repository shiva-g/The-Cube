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
  
  lapply(1:NROW(sel_gene_tbl), FUN = function(x) sum(ir[sel_gene,sel_gene,x], na.rm = T) ) %>% unlist -> sel_gene_tbl$sum 
  lapply(1:NROW(sel_gene_tbl), FUN = function(x) mean(ir[sel_gene,sel_gene,x], na.rm = T) ) %>% unlist -> sel_gene_tbl$mean
  lapply(1:NROW(sel_gene_tbl), FUN = function(x) median(ir[sel_gene,sel_gene,x], na.rm = T) ) %>% unlist -> sel_gene_tbl$median  
  
  sel_gene_tbl_1 = as.data.frame(timex_1)
  sel_gene_tbl_1$timex_1 = as.character(sel_gene_tbl_1$timex_1)
  
  arr_2_25_stx <- ir[sel_gene,sel_gene,tp]
  not_gene <- c(tx[!tx %in% sel_gene])
  arr_2_25_other <- ir[not_gene, not_gene,tp]
  
  gene_present <- reshape2::melt(arr_2_25_stx) %>% filter(!is.na(value))
  gene_absent <- reshape2::melt(arr_2_25_other) %>% filter(!is.na(value))
  ret <- NA
  if (NROW(gene_present) > 0 & NROW(gene_absent) > 0)
  {
    resx <- wilcox.test(gene_present$value,gene_absent$value,alternative = "two.sided")
    resx$p.value -> ret
  }
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

for(numx in 1:NROW(gene_list)) {
  genexx = gene_list$gene[numx] %>% as.character
  print(genexx)
  print('--------------------')
  
  plotx <- matrix(ncol=2,nrow=length(timex_1)) %>% as.data.frame()
  names(plotx) <- c('time_point','pval')
  plotx$time_point <- timex
  plotx$gene <- genexx

  plotx$pval = lapply(1:NROW(timex_1), FUN = function(x) pval_gene(genexx, timex_1[x]) ) %>% unlist
  
  plotx$log <- -log10(plotx$pval)
  plotx_sum <- rbind(plotx_sum,plotx)
  
}

plotx_sum$log_mod = plotx_sum$log - (-log10(0.05))
plotx_sum$log_mod[plotx_sum$log_mod < 0] <- 0
plotx_sum$log_mod[is.na(plotx_sum$log_mod)] <- 0

write_csv(t_prime,paste0(input.yaml$outputDir,'/results/phensim_scores.csv'))

if(input.yaml$plot == T ){
pdf(paste0(input.yaml$outputDir,"/plots/ridge_plot.pdf"),width = 16, height = 9  )
print(ggplot(plotx_sum, aes(x=time_point, y=gene,
                       height = log_mod, fill=gene)) +
  geom_ridgeline(show.legend = F,size=0,alpha=0.3) +
  xlab("Age") +
  ylab('-log(p-value) if p-value is < 0.05') +
  ggtitle("Wilcoxon Test for significant similarities between all genetic epilepies in EGRP"))
  
dev.off()  
}



