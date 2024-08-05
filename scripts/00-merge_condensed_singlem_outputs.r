# This script merges large batch of SingleM condensed profiles
# It produces both the coverage table and the relative abundance/coverage table
# Script created and last updated by Yu Yang on Jul 31 2024

# This script requires an input file named as "output_otu_archive_condensed_defaultfilter.filepaths_nonemptyprofile"
# This file is a new line separated absolute paths of where the condensed profiles are located
# Please note to make sure only non-empty outputs should be included to run this script

library(dplyr)
library(data.table)
library(tibble)

############################################
### Merge singlem condensed result files ###
############################################

### obtain the file paths of condensed SingleM profiles
archive_condense_defaultfilter <- readLines("output_otu_archive_condensed_defaultfilter.filepaths_nonemptyprofile")

### remove preexisting file of name "OTU_table"
rm(OTU_table)

### code below recurssively merges condensed singlem output
for (file in archive_condense_defaultfilter){  
  
  # if the merged OTU_table doesn't exist, create it
  if (!exists("OTU_table")){
    OTU_table <- fread(file, header=TRUE, sep="\t") %>% 
      select(3,2,1)
    new_names <- c("taxonomy", unique(OTU_table$sample), "extra")
    setnames(OTU_table,names(OTU_table),new_names)
    OTU_table <- OTU_table %>% select(1,2)
  }
  
  # if the merged OTU_table does exist, append to it
  if (exists("OTU_table")){
    temp_dataset <-fread(file, header=TRUE, sep="\t")%>%
      select(3,2,1)
    temp_new_names <- c("taxonomy", unique(temp_dataset$sample), "extra")
    setnames(temp_dataset,names(temp_dataset),temp_new_names)
    temp_dataset<-temp_dataset %>% select(1,2)
    # all=T to display all taxon identified, not only taxon identified in ALL samples
    OTU_table<-merge(OTU_table, temp_dataset,all=T,by="taxonomy") 
    rm(temp_dataset)
  }
  
}

### convert NA values to 0 in the dataframe
OTU_table[is.na(OTU_table)] <- 0

### convert 1st to rowname, and remove the next first column (duplicate to the 2nd column)
OTU_table2 <- OTU_table %>%
  tibble::column_to_rownames(.,var="taxonomy")%>%
  select(c(-1))

# fix colnames by removing everything after "_R", run only if needed
colnames(OTU_table2) <- sub("_R.*$", "", colnames(OTU_table2))
OTU_table <- OTU_table2
rm(OTU_table2)

### convert to relative coverage abundances
OTU_table_rel <- OTU_table %>%
  mutate(across(where(is.numeric), ~ round(./sum(.)*100, digits = 8)))

### export merged relative abundance table
write.table(OTU_table_rel%>%
              rownames_to_column(.,var="Lineage"),
            row.names = F,col.names = T,sep="\t", quote = F,
            "./00-inputs/merged.condensed_table_rel_abund.csv")

### export merged table
write.table(OTU_table%>%
              rownames_to_column(.,var="Lineage"),
            row.names = F,col.names = T,sep="\t", quote = F,
              "./00-inputs/merged.condensed_table_abund.csv")

