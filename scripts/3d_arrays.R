library(optparse,quietly = T)
library(yaml,quietly = T)        
library(tidyverse, quietly = T)
library(dplyr, quietly = T)
library(ggplot2,quietly = T)
library(ggridges,quietly = T)
library(ggrepel, quietly = T)

dir.create(paste0(input.yaml$outputDir,"/results"))

tx <- base_temp$STUDY_ID %>% unique()

max_age = 25
timex <- seq(0,max_age,by=0.25)
timex_1 = paste('t_',round(timex,2),sep='')
working_temp <- prop_temp %>% filter(AGE <= max_age) 
n_pat <- prop_temp %>% select(STUDY_ID) %>% unique() 
unique_hpo_prop <- prop_temp %>% select(HPO_ID_prop) %>% unique() 

tx_timex <- matrix(ncol = length(timex),nrow = length(tx)) %>% as.data.frame()
row.names(tx_timex) <- tx
names(tx_timex) <- paste('t_',round(timex,2),sep='')

# ar matrix  - patient, time, hpo - values with 1,0 and NA
ar <- array(dim=c(length(tx),
                  length(timex),
                  NROW(unique_hpo_prop)
))

dimnames(ar)[[1]] <- c(tx)
dimnames(ar)[[2]] <- paste('t_',round(timex,2),sep='')
dimnames(ar)[[3]] <- unique_hpo_prop %>% t %>% as.vector()

# sr matrix - summary from ar matrix 
sr <- array(dim=c(8,
                  length(timex),
                  NROW(unique_hpo_prop)
))

dimnames(sr)[[1]] <- c('all','emr_usage','na_count','present','absent','freq','IC','timex')
dimnames(sr)[[2]] <- paste('t_',round(timex,2),sep='')
dimnames(sr)[[3]] <- unique_hpo_prop %>% t %>% as.vector()

for(numx in 1:NROW(unique_hpo_prop))
{
  hpo_term <- unique_hpo_prop[numx,1] %>% as.character()
  
  tx_timex <- matrix(ncol = length(timex),nrow = length(tx)) %>% as.data.frame()
  row.names(tx_timex) <- tx
  names(tx_timex) <- paste('t_',round(timex,2),sep='')
  
  for(i in 1:NROW(tx_timex))
  {
    patID <- row.names(tx_timex)[i]
    
    emr <- merge_count %>% filter(STUDY_ID == patID) %>% 
      select(lower, upper) %>% 
    rename(min = lower, max = upper)
    emr_lower = which(timex < emr$min); if(length(emr_lower) == 0){emr_lower = 0}
    emr_upper = which(timex > emr$max); if(length(emr_upper) == 0){emr_upper = length(timex)}
    emr_ix <- c(max(emr_lower):min(emr_upper))
    tx_timex[i,emr_ix] <- 0 #write zeroes in all epochs covered by emr, then add layer with 1's
    
    hpo <- working_temp %>% filter(STUDY_ID == patID, HPO_ID_prop == hpo_term, AGE <= max_age)
    hpo <- hpo %>% group_by(STUDY_ID,HPO_ID_prop) %>% summarize(min = min(AGE),max = max(AGE))
    
    if(NROW(hpo) > 0)
    {
      hpo_lower = which(timex < hpo$min); if(length(hpo_lower) == 0){hpo_lower = 0}
      hpo_upper = which(timex > hpo$max); if(length(hpo_upper) == 0){hpo_upper = length(timex)}
      hpo_ix <- c(max(hpo_lower):min(hpo_upper))
      tx_timex[i,hpo_ix] <- 1
    }
    
  }
  
  ## create lr dataframe
  lr <- matrix(ncol = length(timex),nrow = 7) %>% as.data.frame()
  names(lr) <- paste('t_',round(timex,2),sep='')
  row.names(lr) <- c('all','emr_usage','na_count','present','absent','freq','IC')
  for(i in 1:NCOL(lr))
  {
    lr['present',i] <- (tx_timex[,i] == 1) %>% sum(na.rm = T) 
    lr['absent',i] <- (tx_timex[,i] == 0) %>% sum(na.rm = T) 
    lr['na_count',i] <- (tx_timex[,i]  %>% is.na) %>% sum(na.rm = T) 
    lr['emr_usage',i] = lr['present',i] + lr['absent',i]
    lr['all',i] = lr['emr_usage',i] + lr['na_count',i]
    lr['freq',i] = lr['present',i]/lr['emr_usage',i]
    lr['IC',i] = -log10(lr['freq',i])
    if(is.infinite(lr['IC',i])){ lr['IC',i] <- NA }
  }
  
  lr['timex',] <- timex
  lr %>% t %>% as.data.frame -> lrx 
  
  lrx_present <- lrx %>% select(timex,present) %>% mutate(hpo_term = '1 - present')
  lrx_absent <- lrx %>% select(timex,absent) %>% rename(present = absent) %>% mutate(hpo_term = '2 - absent')
  lrx <- rbind(lrx_present,lrx_absent)
  
  sumx <- prop_temp %>% filter(HPO_ID_prop == hpo_term) %>% 
    select(HPO_ID_prop,AGE,HPO_def_prop) %>% 
    group_by(HPO_ID_prop, HPO_def_prop) %>% 
    summarize(total = n(),
              min = min(AGE) %>% round(2),
              max=max(AGE)%>% round(2),
              median=median(AGE)%>% round(2))
  
  ar[,,numx] <- tx_timex %>% as.matrix
  sr[,,numx] <- lr %>% as.matrix
}


# ir matrix - information content of patient pairs . 
ir <- array(dim=c(length(tx),
                  length(tx),
                  length(timex)
))

dimnames(ir)[[1]] <- c(tx)
dimnames(ir)[[2]] <- c(tx)
dimnames(ir)[[3]] <- paste('t_',round(timex,2),sep='')


sim_score = function(pat1,pat2,time_point) {
  
  ar[c(pat1,pat2),time_point,] %>% t %>% as.data.frame -> merge_x
  merge_x$hpo_overlap <- NA 
  merge_x$hpo_overlap <- merge_x[,1] * merge_x[,2] 
  merge_x$hpo_id <- row.names(merge_x)
  
  merge_t <- sr[7,time_point,] %>% as.data.frame
  names(merge_t) <- c("IC")
  merge_t$hpo_id <- row.names(merge_t)
  
  merge_y <- merge_x %>% left_join(merge_t, by = "hpo_id") %>% filter(hpo_overlap == 1)
  
  return(sum(merge_y$IC))
}


comb_table = as.data.frame(t(combn(tx,2)) )
comb_table$V1 = as.character(comb_table$V1)
comb_table$V2 = as.character(comb_table$V2)
timex_1 = paste('t_',round(timex,2),sep='')

for(i in 1:NROW(comb_table)){
  
  for(j in 1:NROW(timex_1)){
  #  print(c(i,j))
    ir[comb_table$V1[i],comb_table$V2[i],timex_1[j]] = sim_score(comb_table$V1[i],comb_table$V2[i],timex_1[j])
  }
}

saveRDS(ar, paste0(input.yaml$outputDir,'/results/ar_mat.rds'))
saveRDS(sr,paste0(input.yaml$outputDir,'/results/sr_mat.rds'))
saveRDS(ir,paste0(input.yaml$outputDir,'/results/ir_mat.rds'))

