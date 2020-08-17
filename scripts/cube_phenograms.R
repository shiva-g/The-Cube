####################################
##### Initialize
####################################

library(data.table)
library(ggrepel)
library(ggplot2)
library(heatmaply)
library(av)

####################################
##### Read files
####################################

## phensim score results
plotx_sum <- fread.csv("../results/phensim_scores.csv",sep)

## Phenotypic associations results
t_prime <- read_csv("../results/hpo_associations_all.csv")

## 
survival = read_csv("../files/survival_v2_clean.csv")
merge_count = read_csv("../files/merge_count_v2_clean.csv")
dx_shifted_hpo = read_csv("../files/dx_shifted_hpo_v2_clean.csv")
sr = readRDS('../results/sr_mat.rds')
ar = readRDS('../results/ar_mat.rds')

## hpo file with id and def
pheno_def = read_csv("/Volumes/helbig_lab/projects/The_Cube/1 - old/dev/pheno_abnrm.csv") ## file not iin git
hpo_tree = read_csv("../files/hpo_is.a_tree.csv")
hpo_ancs = read_csv("../files/hpo_ancestors.csv")

####################################
##### Phenogram
####################################

##  filter Modifier terms from phenotypic associations results
t_prime_pheno = t_prime[t_prime$hpo %in% pheno_def$def,]

### List the time_point for each gene with max sim

plotx_sum %>% group_by(gene) %>% filter(log == max(log,na.rm = T)) %>% filter(gene != "no_gene") -> max_phensim

## Filter phenotypic associatoin file for gene and all associations at that time point 
hpo_assoc <- do.call( rbind, 
                  lapply(1:nrow(max_phensim), function(xx){
                    t_prime_pheno %>% filter( gene == max_phensim$gene[xx], time_point == max_phensim$time_point[xx]) -> tmp
                    tmp  %>% left_join(pheno_def,by=c("hpo"="def")) -> tmp
                    return(tmp) }
                    ))

## write results
fwrite(hpo_assoc,"../results/phensim_phenogram_Aug2020.csv",row.names = F,sep = ",")

## save each file for gene 
# lapply(1:length(unique(hpo_assoc$gene)), function(xx) write.csv(hpo_assoc %>% filter(gene %in% unique(hpo_assoc$gene)[xx]),paste0("../results/",unique(hpo_assoc$gene)[xx],"_phensim_phenogram_Aug2020.csv"),row.names = F) )

## 
tx = dimnames(ar)[[1]]
timex_1 = dimnames(ar)[[2]]
as.numeric(unlist(strsplit(timex_1,split = "_"))) -> tmp
timex = tmp[c(F,T)]
rm(tmp)

hpo_assoc %>% count(gene)  %>% arrange(n %>% desc) -> gene_x
names(gene_x) = c("gene", "Freq")

## Frequencies
frequencies <- lapply(1:nrow(hpo_assoc), function(gene_count) {
  in_gene <- hpo_assoc$gene[gene_count]
  sel_gene <- survival %>% filter(gene == in_gene) %>% 
    filter(STUDY_ID %in% tx)
  sel_gene = sel_gene$STUDY_ID[which(sel_gene$STUDY_ID %in% tx)]
  not_gene <- tx[!tx %in% sel_gene]
  
  ar_present <- ar[sel_gene,,]
  ar_absent <- ar[not_gene,,]
  
  print(gene_count)
  a_1 <- ar_present[,paste0("t_",hpo_assoc$time_point[gene_count]),hpo_assoc$hpo[gene_count]]
  a_2 <- ar_absent[,paste0("t_",hpo_assoc$time_point[gene_count]),hpo_assoc$hpo[gene_count]]
  a_1 <- a_1[!is.na(a_1)]
  a_2 <- a_2[!is.na(a_2)]
  pval = NA
  if( (a_1[a_1 == 1] %>% length()) > 0)
  {
    a_1x <- table(a_1) %>% as.data.frame() #gene_present
    ## freq_gene 
    #hpo_assoc$freq_gene[gene_count] <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
    result_1 <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
    
    if( (a_2[a_2 == 1] %>% length()) > 0){
      a_2x <- table(a_2) %>% as.data.frame() #gene_absent
      ## freq_no_gene
      result_2 <- a_2x$Freq[a_2x$a_2 ==1]/sum(a_2x$Freq)
    } else {
      result_2 <- 1/658
    }
  }
 return(list(result_1,result_2)) 
}
)

## 
hpo_assoc$freq_gene <- unlist(frequencies)[seq(1,1588,2)]
hpo_assoc$freq_nogene <- unlist(frequencies)[seq(2,1588,2)]
  

## data for plotting

hpo_assoc_fil <- hpo_assoc %>% filter(freq_gene > (.2)) %>% 
  mutate(color.code = case_when(logp >= -log10(0.05) ~ "significant",
                                TRUE ~ "not significant")) 

## Separate plots 
plots <- list()
for(i in 1:length(unique(hpo_assoc_fil$gene))){
  
plots[[i]] <- ggplot(hpo_assoc_fil %>% filter(gene == unique(hpo_assoc_fil$gene)[i]), aes(x = freq_nogene,y = freq_gene, 
                  label = definition, color = logp)) + # color = log_pval,
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  # scale_color_gradient(low = "steelblue", high = "steelblue4") +
  geom_point(aes(size = logp, color = factor(color.code))) +
  scale_color_manual(values = c("steelblue", "red")) +
  geom_label_repel(size = 4.5, point.padding = 0.15,box.padding =0.2, show.legend = F,segment.size = 0.2, max.iter = 3500,colour = "black", force=28, family = "Helvetica") +
  theme_classic() +theme(text = element_text(size=12),legend.position = "none") + 
  ggtitle(paste0("Phenogram of ",unique(hpo_assoc_fil$gene)[i]," at ",unique(hpo_assoc_fil$time_point[hpo_assoc_fil$gene == unique(hpo_assoc_fil$gene)[i]]),"years"))
}
pdf("../results/cube_phenograms.pdf",width = 12,height = 9)
sapply(plots,print)
dev.off()



### Movie function

phenogram_freq <- function(gene_name){
## subset data
  t_prime_pheno %>% filter( gene == gene_name) -> tmp_gene
  tmp_gene %>% left_join(pheno_def,by=c("hpo"="def")) -> tmp_gene_1
  tmp_gene_1 %>% mutate(color.code = case_when(logp >= -log10(0.05) ~ "significant",
                                                TRUE ~ "not significant")) -> tmp_gene_1
  ## Frequencies  
  gene_freq <- lapply( 1:nrow(tmp_gene_1),
                       function(gene_count){
                         in_gene <- tmp_gene_1$gene[gene_count]
                         sel_gene <- survival %>% filter(gene == in_gene) %>% 
                           filter(STUDY_ID %in% tx)
                         sel_gene = sel_gene$STUDY_ID[which(sel_gene$STUDY_ID %in% tx)]
                         not_gene <- tx[!tx %in% sel_gene]
                         
                         ar_present <- ar[sel_gene,,]
                         ar_absent <- ar[not_gene,,]
                         print(gene_count)
                         
                         a_1 <- ar_present[,paste0("t_",tmp_gene_1$time_point[gene_count]),tmp_gene_1$hpo[gene_count]]
                         a_2 <- ar_absent[,paste0("t_",tmp_gene_1$time_point[gene_count]),tmp_gene_1$hpo[gene_count]]
                         a_1 <- a_1[!is.na(a_1)]
                         a_2 <- a_2[!is.na(a_2)]
                         pval = NA
                         if( (a_1[a_1 == 1] %>% length()) > 0)
                         {
                           a_1x <- table(a_1) %>% as.data.frame() #gene_present
                           ## freq_gene 
                          # tmp_gene_1$freq_gene[gene_count] <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
                           result_1 <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
                           if( (a_2[a_2 == 1] %>% length()) > 0){
                             a_2x <- table(a_2) %>% as.data.frame() #gene_absent
                             ## freq_no_gene
                           #  tmp_gene_1$freq_nogene[gene_count] <- 
                             result_2 <-  a_2x$Freq[a_2x$a_2 ==1]/sum(a_2x$Freq)
                           } else {
                            # tmp_gene_1$freq_nogene[gene_count] <- 1/658
                             result_2 <-  1/658
                           }
                         }
                         return(list(result_1,result_2))
                       }
    )
  tmp_gene_1$freq_gene <- unlist(gene_freq)[seq(1,length(unlist(gene_freq)),2)]
    tmp_gene_1$freq_nogene <- unlist(gene_freq)[seq(2,length(unlist(gene_freq)),2)]
    rm(gene_freq)
    print(tmp_gene_1)
    return(tmp_gene_1)
}


## scn1a example - will be generalizable for all genes

scn1a <- phenogram_freq("SCN1A")
scn1a %>% arrange(desc(logp)) %>% pull(definition) %>% unique -> sig_terms_scn1a
scn1a %>% mutate(color.code.2 = case_when(definition %in% sig_terms_scn1a[1:10] ~ definition,
                                              TRUE ~ "not significant")) -> scn1a

write.csv(scn1a, "~/Desktop/Cube_Files/Cytoscape_phenogram/scn1a_movie_data.csv",row.names = F)

## function

lapply(1:nrow(gene), function(gn) { 
tmp <- phenogram_freq("SCN1A")
tmp %>% arrange(desc(logp)) %>% pull(definition) %>% unique -> sig_terms
tmp %>% mutate(color.code.2 = case_when(definition %in% sig_terms[1:10] ~ definition,
                                              TRUE ~ "not significant")) -> tmp

write.csv(tmp, paste0("../results/,gn,"_movie_data.csv"),row.names = F)

})


#### phenograms movie functtion
makeplot <- function(gn){ 
  tmp_gene_1 <- read.csv(paste0("../results/",gn,"_movie_data.csv"),stringsAsFactors = F)
  for(i in 1:length(unique(tmp_gene_1$time_point))){
    
    plot <- ggplot(tmp_gene_1 %>% filter(time_point == unique(tmp_gene_1$time_point)[i],freq_gene > (.2),color.code == "significant",color.code.2 != "not significant"), aes(x = freq_nogene,y = freq_gene, 
                                                                              label = definition, color = logp)) + 
      geom_point(aes(size = logp, color = factor(color.code))) +
      # color = log_pval,
      scale_x_continuous(limits = c(0,1)) +
      scale_y_continuous(limits = c(0,1)) +
      geom_abline(intercept = 0, slope = 1) +
      # scale_color_gradient(low = "steelblue", high = "steelblue4") +
      geom_point(aes(size = logp, color = factor(color.code.2))) +
     geom_label_repel(size = 3, point.padding = 0.15,box.padding =0.2, show.legend = F,segment.size = 0.2, max.iter = 3500,colour = "black", force=28, family = "Helvetica") +
      geom_point(data = tmp_gene_1 %>% filter(time_point == unique(tmp_gene_1$time_point)[i],freq_gene > (.2),color.code.2=="not significant"),mapping= aes(size = logp, color = factor(color.code))) +
      # scale_color_manual(values = c(
      #   "not significant" = "gray",
      #   "significant" = "red")) +
      scale_color_manual(values = c( "not significant" = "gray",
                  "Seizures" = "#E56C50",
                  "Generalized seizures" = "#BE8900",
                  "Status epilepticus" = "#909A00",
                  "Neurodevelopmental abnormality" = "#00AA00",
                  "Intellectual disability" = "#00B272",
                  "Neurodevelopmental delay" = "#00B2AE",
                  "Generalized tonic-clonic seizures" = "#00A4E8",
                  "Migraine" = "#927EFC",
                  "Migraine without aura" = "#D65EE9",
                  "Intellectual disability, severe" = "#F25799")) +
      theme_classic() +theme(text = element_text(size=9),legend.position = "none") + 
        ggtitle(paste0("Phenogram of ",unique(tmp_gene_1$gene[tmp_gene_1$time_point == unique(tmp_gene_1$time_point)[i]])," at ",unique(tmp_gene_1$time_point)[i],"years"))
    print(plot)
    }
  
}

## stop-motion movie
av_capture_graphics(makeplot("scn1a"),'~/Desktop/Cube_Files/Cytoscape_phenogram/SCN1A_phenogram_overtime_v1.mp4',res = 180, width = 1080, height = 720, vfilter = 'fps=2',framerate = .25)
# av_capture_graphics(makeplot(),'~/Desktop/Cube_Files/Cytoscape_phenogram/SCN1A_phenogram_overtime_v6.mp4',res = 180, width = 1080, height = 720, vfilter = 'fps=2',framerate = .25)


### HEATMAP
 
 scn1a_heatmap <- function(gn){
   tmp_gene <- read.csv(paste0("~/Desktop/Cube_Files/Cytoscape_phenogram/",gn,"_movie_data.csv"),stringsAsFactors = F)

   
 }
 ## usiing ggplot
 scn1a$freq_ratio <- scn1a$freq_gene / scn1a$freq_nogene
 ggplot(data = scn1a %>% filter(color.code.2 != "not significant"), aes(x = time_point,y=definition,fill = log(freq_ratio))) + geom_tile() + theme_classic() + theme(axis.text.y = element_text(angle = 45),axis.title.y = element_blank())+xlab('Age in years')+scale_fill_gradientn(colours = c("turquoise2", "#CCFF00FF","#FF758F","#FF7578"),na.value = "azure2")
 
 ## usiing heatmap 
 scn1a_mat <- scn1a %>% select(1,9,10)  %>% spread(scn1a,key = color.code.2,value = freq_ratio)
 row.names(scn1a_mat) <- scn1a_mat
 heatmaply()
 
 
 
 
 
 
 
 
 ####################################
##### Initialize
####################################

library(data.table)
library(ggrepel)
library(ggplot2)
library(heatmaply)
library(av)

####################################
##### Read files
####################################

## phensim score results
plotx_sum <- fread.csv("/Volumes/helbig_lab/Users/ganesans/The-Cube/results_3month/phensim_scores.csv",sep)

## hpo file with id and def
pheno_def = read_csv("/Volumes/helbig_lab/projects/The_Cube/1 - old/dev/pheno_abnrm.csv")

## Phenotypic associations results
t_prime <- read_csv("/Volumes/helbig_lab/Users/ganesans/The-Cube/results_3month/hpo_associations_all.csv")

## 
survival = read_csv("/Volumes/helbig_lab/Users/ganesans/The-Cube/files/survival_v2_clean.csv")
merge_count = read_csv("/Volumes/helbig_lab/Users/ganesans/The-Cube/files/merge_count_v2_clean.csv")
dx_shifted_hpo = read_csv("/Volumes/helbig_lab/Users/ganesans/The-Cube/files/dx_shifted_hpo_v2_clean.csv")
sr = readRDS('/Volumes/helbig_lab/Users/ganesans/The-Cube/results_3month/sr_mat.rds')
ar = readRDS('/Volumes/helbig_lab/Users/ganesans/The-Cube/results_3month/ar_mat.rds')

pheno_def = read_csv("/Volumes/helbig_lab/projects/The_Cube/1 - old/dev/pheno_abnrm.csv")
hpo_tree = read_csv("files/hpo_is.a_tree.csv")
hpo_ancs = read_csv("files/hpo_ancestors.csv")

## clean the file

dx_shifted_hpo %>% filter(! NEURO_DX %>% is.na)  -> tmp

merge_count %>% filter(lower < 25) %>%
  filter(STUDY_ID %in% tmp$STUDY_ID) -> merge_count

dx_shifted_hpo %>% filter(STUDY_ID %in% merge_count$STUDY_ID & AGE < 25 & !(is.na(NEURO_DX))) -> dx_shifted_hpo

survival %>% filter(STUDY_ID %in%  merge_count$STUDY_ID ) -> survival

merge_count$upper[which(merge_count$upper > 25)] = 25

# survival %>% count(gene) %>% filter(n >1 & gene != "NA") %>% arrange(n %>% desc) -> gene_x
# names(gene_x) = c("gene", "Freq")

####################################
##### Phenogram
####################################

##  filter Modifier terms from phenotypic associations results
t_prime_pheno = t_prime[t_prime$hpo %in% pheno_def$def,]

### List the time_point for each gene with max sim

plotx_sum %>% group_by(gene) %>% filter(log == max(log,na.rm = T)) %>% filter(gene != "no_gene") -> max_phensim

## Filter phenotypic associatoin file for gene and all associations at that time point 
tmp_1 <- do.call( rbind, 
                  lapply(1:nrow(max_phensim), function(xx){
                    t_prime_pheno %>% filter( gene == max_phensim$gene[xx], time_point == max_phensim$time_point[xx]) -> tmp
                    tmp  %>% left_join(pheno_def,by=c("hpo"="def")) -> tmp
                    return(tmp) }
                    ))

## write results
fwrite(tmp_1,"phensim_phenogram_Aug2020.csv",row.names = F,sep = ",")

## save each file for gene 
# lapply(1:length(unique(tmp_1$gene)), function(xx) write.csv(tmp_1 %>% filter(gene %in% unique(tmp_1$gene)[xx]),paste0("~/Desktop/Cube_Files/Cytoscape_phenogram/",unique(tmp_1$gene)[xx],"_phensim_phenogram_Aug2020.csv"),row.names = F) )

## 
tx = dimnames(ar)[[1]]
timex_1 = dimnames(ar)[[2]]
as.numeric(unlist(strsplit(timex_1,split = "_"))) -> tmp
timex = tmp[c(F,T)]

tmp_1 %>% count(gene)  %>% arrange(n %>% desc) -> gene_x
names(gene_x) = c("gene", "Freq")

## Frequencies
frequencies <- lapply(1:nrow(tmp_1), function(gene_count) {
  in_gene <- tmp_1$gene[gene_count]
  sel_gene <- survival %>% filter(gene == in_gene) %>% 
    filter(STUDY_ID %in% tx)
  sel_gene = sel_gene$STUDY_ID[which(sel_gene$STUDY_ID %in% tx)]
  not_gene <- tx[!tx %in% sel_gene]
  
  ar_present <- ar[sel_gene,,]
  ar_absent <- ar[not_gene,,]
  
  print(gene_count)
  a_1 <- ar_present[,paste0("t_",tmp_1$time_point[gene_count]),tmp_1$hpo[gene_count]]
  a_2 <- ar_absent[,paste0("t_",tmp_1$time_point[gene_count]),tmp_1$hpo[gene_count]]
  a_1 <- a_1[!is.na(a_1)]
  a_2 <- a_2[!is.na(a_2)]
  pval = NA
  if( (a_1[a_1 == 1] %>% length()) > 0)
  {
    a_1x <- table(a_1) %>% as.data.frame() #gene_present
    ## freq_gene 
    #tmp_1$freq_gene[gene_count] <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
    result_1 <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
    
    if( (a_2[a_2 == 1] %>% length()) > 0){
      a_2x <- table(a_2) %>% as.data.frame() #gene_absent
      ## freq_no_gene
      result_2 <- a_2x$Freq[a_2x$a_2 ==1]/sum(a_2x$Freq)
    } else {
      result_2 <- 1/658
    }
  }
 return(list(result_1,result_2)) 
}
)

## 
tmp_1$freq_gene <- unlist(frequencies)[seq(1,1588,2)]
tmp_1$freq_nogene <- unlist(frequencies)[seq(2,1588,2)]
  

## data for plotting

tmp_2 <- tmp_1 %>% filter(freq_gene > (.2)) %>% 
  mutate(color.code = case_when(logp >= -log10(0.05) ~ "significant",
                                TRUE ~ "not significant")) 

## Separate plots 
plots <- list()
for(i in 1:length(unique(tmp_2$gene))){
  
plots[[i]] <- ggplot(tmp_2 %>% filter(gene == unique(tmp_2$gene)[i]), aes(x = freq_nogene,y = freq_gene, 
                  label = definition, color = logp)) + # color = log_pval,
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  # scale_color_gradient(low = "steelblue", high = "steelblue4") +
  geom_point(aes(size = logp, color = factor(color.code))) +
  scale_color_manual(values = c("steelblue", "red")) +
  geom_label_repel(size = 4.5, point.padding = 0.15,box.padding =0.2, show.legend = F,segment.size = 0.2, max.iter = 3500,colour = "black", force=28, family = "Helvetica") +
  theme_classic() +theme(text = element_text(size=12),legend.position = "none") + 
  ggtitle(paste0("Phenogram of ",unique(tmp_2$gene)[i]," at ",unique(tmp_2$time_point[tmp_2$gene == unique(tmp_2$gene)[i]]),"years"))
}
pdf("~/Desktop/Cube_Files/Cytoscape_phenogram/cube_phenograms.pdf",width = 12,height = 9)
sapply(plots,print)
dev.off()



### Movie... for SCN1A

phenogram_freq <- function(gene_name){
## subset data
  t_prime_pheno %>% filter( gene == gene_name) -> tmp_gene
  tmp_gene %>% left_join(pheno_def,by=c("hpo"="def")) -> tmp_gene_1
  tmp_gene_1 %>% mutate(color.code = case_when(logp >= -log10(0.05) ~ "significant",
                                                TRUE ~ "not significant")) -> tmp_gene_1
  ## Frequencies  
  gene_freq <- lapply( 1:nrow(tmp_gene_1),
                       function(gene_count){
                         in_gene <- tmp_gene_1$gene[gene_count]
                         sel_gene <- survival %>% filter(gene == in_gene) %>% 
                           filter(STUDY_ID %in% tx)
                         sel_gene = sel_gene$STUDY_ID[which(sel_gene$STUDY_ID %in% tx)]
                         not_gene <- tx[!tx %in% sel_gene]
                         
                         ar_present <- ar[sel_gene,,]
                         ar_absent <- ar[not_gene,,]
                         print(gene_count)
                         
                         a_1 <- ar_present[,paste0("t_",tmp_gene_1$time_point[gene_count]),tmp_gene_1$hpo[gene_count]]
                         a_2 <- ar_absent[,paste0("t_",tmp_gene_1$time_point[gene_count]),tmp_gene_1$hpo[gene_count]]
                         a_1 <- a_1[!is.na(a_1)]
                         a_2 <- a_2[!is.na(a_2)]
                         pval = NA
                         if( (a_1[a_1 == 1] %>% length()) > 0)
                         {
                           a_1x <- table(a_1) %>% as.data.frame() #gene_present
                           ## freq_gene 
                          # tmp_gene_1$freq_gene[gene_count] <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
                           result_1 <- a_1x$Freq[a_1x$a_1 ==1]/sum(a_1x$Freq)
                           if( (a_2[a_2 == 1] %>% length()) > 0){
                             a_2x <- table(a_2) %>% as.data.frame() #gene_absent
                             ## freq_no_gene
                           #  tmp_gene_1$freq_nogene[gene_count] <- 
                             result_2 <-  a_2x$Freq[a_2x$a_2 ==1]/sum(a_2x$Freq)
                           } else {
                            # tmp_gene_1$freq_nogene[gene_count] <- 1/658
                             result_2 <-  1/658
                           }
                         }
                         return(list(result_1,result_2))
                       }
    )
  tmp_gene_1$freq_gene <- unlist(gene_freq)[seq(1,length(unlist(gene_freq)),2)]
    tmp_gene_1$freq_nogene <- unlist(gene_freq)[seq(2,length(unlist(gene_freq)),2)]
    rm(gene_freq)
    print(tmp_gene_1)
    return(tmp_gene_1)
}


## scn1a example
scn1a <- phenogram_freq("SCN1A")
scn1a %>% arrange(desc(logp)) %>% pull(definition) %>% unique -> sig_terms_scn1a
scn1a %>% mutate(color.code.2 = case_when(definition %in% sig_terms_scn1a[1:10] ~ definition,
                                              TRUE ~ "not significant")) -> scn1a

write.csv(scn1a, "~/Desktop/Cube_Files/Cytoscape_phenogram/scn1a_movie_data.csv",row.names = F)


#### phenograms movie functtion
makeplot <- function(gn){ 
  tmp_gene_1 <- read.csv(paste0("~/Desktop/Cube_Files/Cytoscape_phenogram/",gn,"_movie_data.csv"),stringsAsFactors = F)
  for(i in 1:length(unique(tmp_gene_1$time_point))){
    
    plot <- ggplot(tmp_gene_1 %>% filter(time_point == unique(tmp_gene_1$time_point)[i],freq_gene > (.2),color.code == "significant",color.code.2 != "not significant"), aes(x = freq_nogene,y = freq_gene, 
                                                                              label = definition, color = logp)) + 
      geom_point(aes(size = logp, color = factor(color.code))) +
      # color = log_pval,
      scale_x_continuous(limits = c(0,1)) +
      scale_y_continuous(limits = c(0,1)) +
      geom_abline(intercept = 0, slope = 1) +
      # scale_color_gradient(low = "steelblue", high = "steelblue4") +
      geom_point(aes(size = logp, color = factor(color.code.2))) +
     geom_label_repel(size = 3, point.padding = 0.15,box.padding =0.2, show.legend = F,segment.size = 0.2, max.iter = 3500,colour = "black", force=28, family = "Helvetica") +
      geom_point(data = tmp_gene_1 %>% filter(time_point == unique(tmp_gene_1$time_point)[i],freq_gene > (.2),color.code.2=="not significant"),mapping= aes(size = logp, color = factor(color.code))) +
      # scale_color_manual(values = c(
      #   "not significant" = "gray",
      #   "significant" = "red")) +
      scale_color_manual(values = c( "not significant" = "gray",
                  "Seizures" = "#E56C50",
                  "Generalized seizures" = "#BE8900",
                  "Status epilepticus" = "#909A00",
                  "Neurodevelopmental abnormality" = "#00AA00",
                  "Intellectual disability" = "#00B272",
                  "Neurodevelopmental delay" = "#00B2AE",
                  "Generalized tonic-clonic seizures" = "#00A4E8",
                  "Migraine" = "#927EFC",
                  "Migraine without aura" = "#D65EE9",
                  "Intellectual disability, severe" = "#F25799")) +
      theme_classic() +theme(text = element_text(size=9),legend.position = "none") + 
        ggtitle(paste0("Phenogram of ",unique(tmp_gene_1$gene[tmp_gene_1$time_point == unique(tmp_gene_1$time_point)[i]])," at ",unique(tmp_gene_1$time_point)[i],"years"))
    print(plot)
    }
  
}

## stop-motion movie
av_capture_graphics(makeplot("scn1a"),'~/Desktop/Cube_Files/Cytoscape_phenogram/SCN1A_phenogram_overtime_v10.mp4',res = 180, width = 1080, height = 720, vfilter = 'fps=2',framerate = .25)
av_capture_graphics(makeplot(),'~/Desktop/Cube_Files/Cytoscape_phenogram/SCN1A_phenogram_overtime_v6.mp4',res = 180, width = 1080, height = 720, vfilter = 'fps=2',framerate = .25)

 av_capture_graphics(makeplot(),'~/Desktop/Cube_Files/Cytoscape_phenogram/SCN1A_phenogram_overtime.mp4',res = 180, width = 1080, height = 720, framerate = 2, vfilter = 'fps=35')


### HEATMAP

 ## usiing ggplot
 
gg_heatmap <- function(gn){
tmp_gene_1 <- read.csv(paste0("../results/",gn,"_movie_data.csv"),stringsAsFactors = F)
tmp_gene_1$freq_ratio <- tmp_gene_1$freq_gene / tmp_gene_1$freq_nogene
return(print( ggplot(data = scn1a %>% filter(color.code.2 != "not significant"), aes(x = time_point,y=definition,fill = log(freq_ratio))) + 
                   ggtitle(paste0("Heatmap of ",gn))+ geom_tile() + theme_classic() + theme(axis.text.y = element_text(angle = 45),axis.title.y = element_blank())+xlab('Age in years')+scale_fill_gradientn(colours = c("turquoise2", "#CCFF00FF","#FF758F","#FF7578"),na.value = "azure2")))
 }


 
