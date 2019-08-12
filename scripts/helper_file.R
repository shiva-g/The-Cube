library(optparse,quietly = T)
library(yaml,quietly = T)        
library(tidyverse, quietly = T)
library(dplyr, quietly = T)
library(ggplot2,quietly = T)
library(ggridges,quietly = T)
library(ggrepel, quietly = T)

if(exists("input.yaml")){
  input.yaml <- input.yaml
}else{
  message('Input YAML not found.\n')
  break;
}

dir.create(paste0(input.yaml$outputDir,"/temp"))

merge_count <- read_csv(input.yaml$aed.encounter)
dx_shifted_hpo_raw <- read_csv(input.yaml$diagnosis)
survival = read_csv(input.yaml$survival)

if(is.null(input.yaml$hpo_tree)){
hpo_tree = read_csv("files/hpo_is.a_tree.csv")
} else {
  hpo_tree = read_csv(input.yaml$hpo_tree)
}

if(is.null(input.yaml$hpo_ancs)){
  hpo_ancs = read_csv("files/hpo_ancestors.csv")
} else {
  hpo_ancs = read_csv(input.yaml$hpo_ancs)
}

# Filtering the datasets for patients with encounter before the age of 25. 
dx_shifted_hpo_raw %>% filter(! NEURO_DX %>% is.na)  -> tmp

merge_count %>% filter(lower < 25) %>% 
  filter(STUDY_ID %in% tmp$STUDY_ID) -> merge_count

dx_shifted_hpo_raw %>% filter(STUDY_ID %in% merge_count$STUDY_ID & AGE < 25 & !(is.na(NEURO_DX))) -> dx_shifted_hpo

survival %>% filter(STUDY_ID %in%  merge_count$STUDY_ID ) -> survival

merge_count$upper[which(merge_count$upper > 25)] = 25

rm(tmp)

length(unique(survival$STUDY_ID)) -> tmp 
length(unique(dx_shifted_hpo$STUDY_ID)) -> tmp_1
length(unique(merge_count$STUDY_ID)) -> tmp_2

if(tmp == tmp_1) {
  if(tmp == tmp_2){
    if(tmp_1 == tmp_2){
      indv <- tmp
    } else {
      message('Number of individuals in diagnosis and aed enconters are different')
      break;
    }
  } else {
    message('Number of individuals in survival and aed enconters are different')
    break;
  }
  } else {
    message('Number of individuals in survival and diagnosis are different')
    break;
  }
rm(tmp,tmp_1,tmp_2)



message(paste0('Numbers from the data provided :\n ',
               length(unique(survival$STUDY_ID)),'\tIndividuals \n',
               round(sum(merge_count$upper) - sum(merge_count$lower)) ,'\tPatient years.\n'))


######### base to prop v1 
# base_temp <- adding def of HPO terms to diagnosis table
base_temp <- dx_shifted_hpo %>%
  select(HPO_IMO_ID,AGE,STUDY_ID) %>%
  unique() %>% 
  separate_rows(HPO_IMO_ID,sep=';') %>% 
  left_join((hpo_tree %>% select(term, def)), by = c("HPO_IMO_ID" = "term")) %>% 
  rename(HPO_IMO_def = def) 

#all base terms propagate up to the parent term  
prop2 <- base_temp %>% 
  select(STUDY_ID,AGE,HPO_IMO_ID) %>%
  left_join(hpo_ancs,by = c("HPO_IMO_ID"="term")) %>%
  mutate(complete = grepl('HP:0000001',ancs)) %>%
  group_by(complete) %>% summarize(n = n())

prop2 <- base_temp %>% 
  select(STUDY_ID,AGE,HPO_IMO_ID) %>%
  left_join(hpo_ancs,by = c("HPO_IMO_ID"="term"))

prop3 <- prop2 %>%
  select(STUDY_ID,AGE,ancs) %>%
  separate_rows(ancs,sep=';') %>%
  rename(HPO_ID_prop = ancs)

prop4 <- prop3 %>% 
  unique()

prop4b <- prop4 %>%
  left_join(hpo_ancs,by = c("HPO_ID_prop"="term")) %>% 
  select(STUDY_ID,AGE,HPO_ID_prop,def)  %>% 
  rename(HPO_def_prop = def)

#Does every time point have 'HP:0000001'

prop5 <- prop4 %>% filter(HPO_ID_prop == 'HP:0000001') %>% group_by(AGE) %>% summarise(n = n())

if( prop5 %>% filter(n == 0) %>% NROW() != 0)
{
  message("Every time point did not propogate to HP:0000001, check the datasets.")
  break;
} else {
  rm(prop5,prop3,prop4,prop2)
  
  prop_temp <- prop4b; rm(prop4b)
  
}

write_csv(base_temp, paste0(input.yaml$outputDir,"/temp/base_hpo.csv"))
write_csv(base_temp, paste0(input.yaml$outputDir,"/temp/prop_hpo.csv"))

