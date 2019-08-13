# EMR based time-dependent genetic associations

This repository contains a set of scripts which analyzes electronic health records (EHR) records of a population accross a set time. Time-stamped problem list information from neurology encounters is extracted from the EHR and then translated into HPO (Human Phenotype Onotology) terms. For each patient, terms are split into 3 month time bins according to their time-stamp. Value or Information Content (IC) is assigned to each term at each time point according to its prevalence within that 3 month time bin. Patients are then grouped according to their diagnosed genetic diagnosis and the similarity within each genetic-based group is compared to other patients at each particular bin. Consequently, we can find which genes at each time bins are signficantly more similar to other patients within their gene group than the rest of the cohort at the respective time bin. Furthermore, we can find which HPO terms are significantly more present in a gene group at each time bin.


## Scripts: ##

This [wrapper](https://github.com/shiva-g/The-Cube/blob/master/wrapper.R) script needs to be submitted to run the entire pipeline.

  [Helper file](https://github.com/shiva-g/The-Cube/blob/master/scripts/helper_file.R)  - loads data and cleans it. Creates base and prop hpo files. 
  
  [3d Array creation](https://github.com/shiva-g/The-Cube/blob/master/scripts/3d_arrays.R) - creates the 3d matrices.
  
  [Fishers test](https://github.com/shiva-g/The-Cube/blob/master/scripts/hpo_associations.R) - hpo_associations.
  
  [Wilcoxon test](https://github.com/shiva-g/The-Cube/blob/master/scripts/wilcoxon_test.R) - ridge plot.
  

## Files: ##

[hpo_is.a_tree.csv](https://github.com/shiva-g/The-Cube/blob/master/files/hpo_is.a_tree.csv) - This file contains the ontological information and definition for every single HPO term. The 'is.a' term is the parent term for each respective HPO term.

[hpo_ancestors.csv](https://github.com/shiva-g/The-Cube/blob/master/files/hpo_ancestors.csv) -  This file contains the higher level terms of each HPO term. This is essential for calculating the MICA (most informative common ancestor) between 2 HPO terms which is one of the first steps to finding the similarity score between patients.


### Requirements:
  [R](https://www.r-project.org/) with the following packages:
    1. [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
    2. [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
    3. [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
    4. [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
    5. [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
    6. [ggridges](https://cran.r-project.org/web/packages/ggridges/index.html)
    7. [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)

### Submitting :

```
Rscript wrapper.R --input input.yml
```
