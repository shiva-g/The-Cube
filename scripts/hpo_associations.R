library(optparse,quietly = T)
library(yaml,quietly = T)        
library(tidyverse, quietly = T)
library(dplyr, quietly = T)
library(ggplot2,quietly = T)
library(ggridges,quietly = T)
library(ggrepel, quietly = T)

dir.create(paste0(input.yaml$outputDir,"/plots"))
survival %>% count(gene) %>% filter(n >1 & gene != "NA") %>% arrange(n %>% desc) -> gene_x

names(gene_x) = c("gene", "Freq")
t_prime = NULL  

fish_test <- function(tp,hpo_id){
  a_1 <- ar_present[,tp,hpo_id]
  a_2 <- ar_absent[,tp,hpo_id]
  a_1 <- a_1[!is.na(a_1)]
  a_2 <- a_2[!is.na(a_2)]
  pval = NA
  if( (a_1[a_1 == 1] %>% length()) > 0)
  {
    a_1x <- table(a_1) %>% as.data.frame() #gene_present
    a_2x <- table(a_2) %>% as.data.frame() #gene_absent
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('gene_present','gene_absent')
    
    #fish['gene_present','hpo_present'] <- a_1x$Freq[a_1x$a_1 == 1] 
    fish['gene_present','hpo_present'] <- 0 
    if( (a_1x$Freq[a_1x$a_1 == 1] %>% length) > 0)
    {
      fish['gene_present','hpo_present'] <- a_1x$Freq[a_1x$a_1 == 1] 
    }
    
    #fish['gene_present','hpo_absent'] <- a_1x$Freq[a_1x$a_1 == 0]
    fish['gene_present','hpo_absent'] <- 0 
    if( (a_1x$Freq[a_1x$a_1 == 0] %>% length) > 0)
    {
      fish['gene_present','hpo_absent'] <- a_1x$Freq[a_1x$a_1 == 0] 
    }
    
    #fish['gene_absent','hpo_present'] <- a_2x$Freq[a_2x$a_2 == 1]
    fish['gene_absent','hpo_present'] <- 0 
    if( (a_2x$Freq[a_2x$a_2 == 1] %>% length) > 0)
    {
      fish['gene_absent','hpo_present'] <- a_2x$Freq[a_2x$a_2 == 1] 
    }
    
    #fish['gene_absent','hpo_absent'] <- a_2x$Freq[a_2x$a_2 == 0]
    fish['gene_absent','hpo_absent'] <- 0 
    if( (a_2x$Freq[a_2x$a_2 == 0] %>% length) > 0)
    {
      fish['gene_absent','hpo_absent'] <- a_2x$Freq[a_2x$a_2 == 0]  
    }
    
    # using fisher's test to calculate pval
    pval <- fisher.test(fish)$p.value  
  }
  logx = -log10(pval)
  return(logx)
}

for (gene_count in 1:NROW(gene_x)) {
  in_gene <- gene_x$gene[gene_count]
  sel_gene <- survival %>% filter(gene == in_gene) %>%   filter(STUDY_ID %in% tx)
  sel_gene = sel_gene$STUDY_ID[which(sel_gene$STUDY_ID %in% tx)]
  not_gene <- tx[!tx %in% sel_gene]
  
  ar_present <- ar[sel_gene,,]
  ar_absent <- ar[not_gene,,]
  # hpo_freq is the matrix of all time points x hpo names
  
  hpo_terms <- dimnames(ar)[[3]]
  hpo_freq <- matrix(nrow=length(timex), ncol=length(hpo_terms)+1) %>% as.data.frame
  names(hpo_freq) <- c('time_point',hpo_terms)
  hpo_freq$time_point <- timex
  
   print(in_gene)
  
  tmp <- unlist( lapply(2:NCOL(hpo_freq), function(x){ 
    lapply(1:NROW(timex_1), function(y){
      fish_test(timex_1[y],names(hpo_freq)[x])
    }) }) )
  
  hpo_freq[,2:NCOL(hpo_freq)] <- tmp
  
  txx <- hpo_freq %>% select(-time_point)
  
  tyy <- matrix(ncol =3, nrow = ncol(txx)) %>% as.data.frame
  tyy[,1] <- names(txx)
  for(i in 1:NROW(tyy))
  {
    temp <- txx[,i]
    max_temp <- NA
    time_temp <- NA
    if (length(temp[!is.na(temp)]) > 0)
    {
      max_temp <- max(temp,na.rm = T)
      time_temp <- timex[which(temp == max(temp,na.rm = T))]
      if(length(time_temp) > 1) {time_temp <- time_temp[1]}
    }
    tyy[i,2] <- max_temp 
    tyy[i,3] <- time_temp
  }
  
  names(tyy) <- c('hpo_id',"max_log","time_point")
  hpo_dict <- hpo_ancs %>% select(term, def)
  tzz <- merge(tyy,hpo_dict, 
               by.x = 'hpo_id',
               by.y = 'term', 
               all.x = T, all.y = F)
  tzz$gene <- in_gene
  tzz$max_log[tzz$max_log == '-Inf'] <- NA
  tzz <- tzz %>% filter(!is.na(max_log))
  
  t_prime <- rbind(t_prime,tzz)
} 
write_csv(t_prime,paste0(input.yaml$outputDir,'/results/hpo_associations.csv'))

if(input.yaml$plot == T ){
  for (gn in 1:NROW(gene_x)) {
    spec_gene = gene_x$gene[gn]
    print(spec_gene)
    t_prime_gene <- t_prime %>%
      arrange(max_log %>% desc) %>%
      filter(max_log > -log10(0.05)) %>%
      filter(gene == spec_gene)
    
    titlex <- paste('Association of HPO terms and',spec_gene,'over time')
    
    pdf(paste0(input.yaml$outputDir,"/plots/Asociation_of_HPO_terms_and_",spec_gene,'_over_time.pdf'),width = 16, height = 9  )
    
    print(ggplot(t_prime_gene,aes(x=time_point,
                            y=max_log,
                            size= max_log,
                            color= max_log,
                            label = def
    )) +
      geom_point() +
      geom_label_repel(size = 4) +
      theme(legend.position = 'none') +
      ggtitle(titlex) +
      xlab('Age') +
      ylab('-log10(pvalue) at time point with most significant association'))
    dev.off()
    # ggsave(paste0(input.yaml$outputDir,"/plots/Asociation_of_HPO_terms_and_",spec_gene,'_over_time.jpeg'),  dpi = 750, width = 16, height = 9)
  }
  
}
