# EMR based time-dependent genetic associations

This repository contains a set of scripts which analyzes electronic health records (EHR) records of a population accross a set time. Time-stamped problem list information from neurology encounters is extracted from the EHR and then translated into HPO (Human Phenotype Onotology) terms. For each patient, terms are split into 3 month time bins according to their time-stamp. Worth or Information Content (IC) is assigned to each term at each time point according to its prevalence within that 3 month time bin. Patients are then grouped according to their diagnosed genetic diagnosis and the similarity within each genetic-based group is compared to other patients at each particular bin. Consequently, we can find which genes at each time bins are signficantly more similar to other patients within their gene group than the rest of the cohort at the respective time bin. Furthermore, we can find which HPO terms are significantly more present in a gene group at each time bin.


## Scripts: ##


## Files: ##

hpo_is.a_tree.csv - This file contains the ontological information and definition for every single HPO term. The 'is.a' term is the parent term for each respective HPO term.
hpo_ancestors.csv -  This file contains the higher level terms of each HPO term. This is essential for calculating the MICA (most informative common ancestor) between 2 HPO terms which is one of the first steps to finding the similarity score between patients.

