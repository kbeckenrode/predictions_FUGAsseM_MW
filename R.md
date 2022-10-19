title: "Predictions"
output: html_document
date: "2022-10-19"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##Look for other enzymes like this (i.e. mucus utilization) in Bacteroides/Parabacteroides
##Snip off top MetaWIBELE prioritizations in Bacteroides (etc.)
##Look for overlap of these with FUGAsseM predictions in glycan utilization terms
##Probably starting with GO:0006486 or GO:0009100, or wherever that overlaps with the FUGAsseM prediction “slice” set


## Step 1: Import FUGAsseM prediction file, METAWIBELE file, and list of FUGAsseM taxa into R project.
```{r}
library(readr)
HMP2_fugassem_GO_func_assignment_tsv <- read_delim("HMP2_fugassem_GO.func_assignment.tsv.gz", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
##View(HMP2_fugassem_GO_func_assignment_tsv)

library(readr)
HMP2_fugassem_string_taxa <- read_table("HMP2_fugassem_string_taxa.txt")
##View(HMP2_fugassem_string_taxa)

##load MetaWIBELE files
library(readr)
only_bacteroidetes <- read_delim("only_bacteroidetes.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
View(only_bacteroidetes)

library(readr)
only_proteobacteria <- read_delim("only_proteobacteria.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
View(only_proteobacteria)

library(readr)
only_firmicutes <- read_delim("only_firmicutes.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
View(only_firmicutes)

##load GO terms
library(readxl)
GO_terms_20220805 <- read_excel("GO_terms_20220805.xlsx")
View(GO_terms_20220805)

##load FUGAsseM
library(readr)
HMP2_fugassem_GO_func_assignment <- read_delim("HMP2_fugassem_GO.func_assignment.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
View(HMP2_fugassem_GO_func_assignment)

```

```{r}
## Step 2: Load libraries

library(dplyr)
library(stringr)

## adding  %>% filter(`Priority score` >0.5)  will pipe more filters as a TRUE/FALSE statement

```

```{r}
##Step 3: (only need to do once to make file) Get GO terms at and beneath parent terms 'GO:0006486' and 'GO:0009100'

library(GOfuncR)
x<-get_child_nodes(c('GO:0042710'))
##GO_272<-get_child_nodes(c('GO:0000272'))


```
```{r}

##only_bacteroidetes_filtered <- only_bacteroidetes %>% filter(`Priority score` >0.8)
##only_proteobacteria_filtered <- only_proteobacteria %>% filter(`Priority score` >0.9)
##only_firmicutes_filtered <- only_firmicutes %>% filter(`Priority score` >0.9)
```

## Step 4: Find Bacteroides and filter in FUGAsseM file. Also, re-arrange syntax of func column by adding a new fixed column. There was extra information attached to the GO term, so it could not intersect as is.

```{r}
##b <- HMP2_fugassem_GO_func_assignment_tsv %>% filter(str_detect(taxa_lineage, "p__Bacteroidetes"))
##new_b<-data.frame(str_split_fixed(b$func,"\\[",2 ))
##new_b3<-gsub(":","", new_b$X1)
##new_b4<-gsub("GO","GO:", new_b3)
##b$fixed<-new_b4

bb <- HMP2_fugassem_GO_func_assignment_tsv %>% filter(str_detect(taxa_lineage, "p__Bacteroidetes"))
new_b<-data.frame(str_split_fixed(bb$func,"\\[",2 ))
new_b3<-gsub(":","", new_b$X1)
new_b4<-gsub("GO","GO:", new_b3)
bb$fixed<-new_b4

```

```{r}
bp <- HMP2_fugassem_GO_func_assignment_tsv %>% filter(str_detect(taxa_lineage, "p__Proteobacteria"))
new_b<-data.frame(str_split_fixed(bp$func,"\\[",2 ))
new_b3<-gsub(":","", new_b$X1)
new_b4<-gsub("GO","GO:", new_b3)
bp$fixed<-new_b4

```

```{r}
bf <- HMP2_fugassem_GO_func_assignment_tsv %>% filter(str_detect(taxa_lineage, "p__Firmicutes"))
new_b<-data.frame(str_split_fixed(bf$func,"\\[",2 ))
new_b3<-gsub(":","", new_b$X1)
new_b4<-gsub("GO","GO:", new_b3)
bf$fixed<-new_b4

```
## Step 5: Intersect Bacteroides filtered MetaWIBELE priorities and FUGAsseM predictions with GO terms

```{r}
##Intersect MetaWIBELE with GO list
both_MW_go_b <-(data.frame(intersect(only_bacteroidetes_filtered$`GO(BP)`, x$child_go_id)))
both_MW_go_p <-(data.frame(intersect(only_proteobacteria_filtered$`GO(BP)`, x$child_go_id)))
both_MW_go_f <-(data.frame(intersect(only_firmicutes_filtered$`GO(BP)`, x$child_go_id)))


##Intersect FUGAsseM for GO terms
##both_FUGAsseM_go<-(data.frame(intersect(b$fixed, GO_expanded$child_go_id)))

FUGAsseM_go_b_biofilm<-(data.frame(intersect(bb$fixed, x$child_go_id)))
both_FUGAsseM_go_f<-(data.frame(intersect(bf$fixed, x$child_go_id)))
both_FUGAsseM_go_p<-(data.frame(intersect(bp$fixed, x$child_go_id)))

FUGAsseM_go_b_intersect<-(data.frame(intersect(bb2$feature, only_bacteroidetes_filtered$familyID)))

write_csv(both_FUGAsseM_go, "FUGAsseM_GO")

both_FUGAsseM_go_b_filtered <- bb2 %>% filter(`score` >0.95)


```
