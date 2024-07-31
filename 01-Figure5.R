# This script recreates Figure 5 in the big MfD paper 
# Script created and last updated by Yu Yang on Jul 31 2024

# This script includes
### novelty of metagenomes based on SingleM measured as the known species fraction
### novelty of short-read (SR) MFD MAGs compared to GTDB R220 

################################
##### Load needed packages ##### 
################################

library(readxl) 
library(data.table) 
library(dplyr)
library(stringr)
library(wesanderson)
library(tidyverse)
library(ggplot2)


################################
########## load data ###########
################################

### Load MFD metadata
mfd_metadata <- fread("00-inputs/merged_metadata_20240217.csv")

### Load MFD SingleM profile results in relative abundances using supplemented GTDB R214 metapackage without short-read (SR) MFD MAGs
OTU_table_rel<-fread("00-inputs/merged.condensed_table_rel_abund.csv")

### Load MFD SingleM profile results in abundances (coverages) using supplemented GTDB R214 metapackage without short-read (SR) MFD MAGs
OTU_table_cov<-fread(file="00-inputs/merged.condensed_table_abund.csv")

### Load MFD SingleM profile results in relative abundances using supplemented GTDB R214 metapackage with short-read (SR) MFD MAGs
OTU_table_rel_new<-fread(file="00-inputs/merged.condensed_table_rel_abund_withSRMAGs.csv")

### Load MFD SingleM profile results in abundances (coverages) using supplemented GTDB R214 metapackage with short-read (SR) MFD MAGs
OTU_table_cov_new<-fread(file="00-inputs/merged.condensed_table_abund_withSRMAGs.csv")

### Load NCBI datesets SingleM species fraction results using supplemented GTDB R214 metapackage 
# Results were provided by Ben J. Woodcroft, Queensland University of Technology, from his SingleM paper (doi: https://doi.org/10.1101/2024.01.30.578060)
# Data available at http://dx.doi.org/10.5281/zenodo.10547501
NCBI_species_anno<-fread("00-inputs/NCBI_withoutSRMAGs_species_frac_data_filtered.txt")


#################################
##### Create color palettes #####
#################################
# This was created by Thomas Bygh Nymann Jensen, Aalborg University
### Ontology palettes
sediment.palette <- colorRampPalette(c(wes_palette("IsleofDogs2")[3], wes_palette("FantasticFox1")[1]))
plot(rep(1, 5), col = sediment.palette(5), pch = 19, cex = 3)
soil.palette <- colorRampPalette(c(wes_palette("AsteroidCity1")[4], wes_palette("AsteroidCity1")[1]))
plot(rep(1, 9), col = soil.palette(9), pch = 19, cex = 3)
water.palette <- colorRampPalette(c(wes_palette("Darjeeling2")[2], wes_palette("Zissou1")[2]))
plot(rep(1, 4), col = water.palette(4), pch = 19, cex = 3)

### Sampletype palette
sampletype.palette <- c(soil.palette(1), sediment.palette(1), water.palette(1))
names(sampletype.palette) <- c("Soil", "Sediment", "Water")

sampletype.palette
plot(rep(1, 3), col = sampletype.palette, pch = 19, cex = 3)

### sampletype palette NCBI
sampletype.palette.NCBI <- c(soil.palette(1), sediment.palette(1), water.palette(1),"#e64b35")
names(sampletype.palette.NCBI) <- c("NCBI_Soil", "NCBI_Sediment", "NCBI_Water", "NCBI_Human")

### mfdo1 palette
mfdo1 <- read_excel('00-inputs/2024-02-13_mfd_db.xlsx') %>%
  select(mfd_sampletype:mfd_hab1) %>%
  filter(!is.na(mfd_hab1)) %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub"))) %>%
  distinct() %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  filter(!complex %in% c("Water, Subterranean, Freshwater",
                         "Water, Urban, Sandfilter",
                         "Soil, Urban, Other",
                         "Sediment, Urban, Other",
                         "Soil, Urban, Roadside",
                         "Soil, Subterranean, Urban")) %>%
  mutate(across(mfd_sampletype, ~factor(., levels = c("Soil", "Sediment", "Water"))),
         across(mfd_areatype, ~factor(., levels = c("Natural", "Subterranean", "Agriculture", "Urban")))) %>%
  arrange(mfd_sampletype, mfd_areatype, mfd_hab1)

mfdo1.palette <- c(soil.palette(9), sediment.palette(5), water.palette(4))
names(mfdo1.palette) <- mfdo1 %>% pull(complex)

mfdo1.palette
plot(rep(1, 18), col = mfdo1.palette, pch = 19, cex = 3)

### mfdo1 + NCBI palette
NCBI_MfD.palette<-c(sampletype.palette.NCBI,mfdo1.palette)

####################################################################################
##### For MfD samples, metapackage without SR MFD MAGs: species level coverage #####
####################################################################################

### apply filters: coverage (root coverage) 50, 1Gbp filter (the same as SingleM sandpiper paper: doi: https://doi.org/10.1101/2024.01.30.578060)
# 9820 samples remain
# count the total root coverage per sample
OTU_table_cov_sum<-OTU_table_cov%>%
  dplyr::select(-1)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  t()%>%
  as.data.frame()%>%
  rename(root_coverage = V1)
# count the total species coverage per sample
OTU_table_cov_species<-OTU_table_cov%>%
  filter(grepl("s__", Lineage, fixed = TRUE))%>%
  dplyr::select(-1)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  t()%>%
  as.data.frame()%>%
  rename(species_coverage = V1)
# merge root and species coverage and compute species_annotation: 
MfD_species_anno<-merge(OTU_table_cov_sum,OTU_table_cov_species,by=0,all=T)%>%
  column_to_rownames(.,var="Row.names")%>%
  rownames_to_column(.,var="flat_name")%>%
  mutate(species_annotation=100*(species_coverage/root_coverage))%>%
  merge(.,mfd_metadata%>%dplyr::select(c(2,14,36:40)),by="flat_name",all=F)%>%
  mutate(mbases=after_total_bases/1000000)%>%
  # filter: min. total root coverage of 50
  filter(root_coverage>=50) %>%
  # filter: min 1Gbp data
  filter(mbases>=1000) 

# obtain the ratio of overall root coverage and the overall species coverage
# Calculate the sum of "species_coverage"
MFD_sum_species_coverage <- sum(MfD_species_anno$species_coverage, na.rm = TRUE)
# Calculate the sum of "root_coverage"
MFD_sum_root_coverage <- sum(MfD_species_anno$root_coverage, na.rm = TRUE)
# Calculate the overall ratio: 7.8%
MFD_species_ratio <- MFD_sum_species_coverage / MFD_sum_root_coverage
# Calculate the average of species annotation ratio: 7.26%
MFD_species_ratio_average <-  mean(MfD_species_anno$species_annotation, na.rm = TRUE)

### extract data classified to the species level 
# Extract rows with "s__" in the "Lineage" column
OTU_table_rel_species <- OTU_table_rel %>%
  filter(str_detect(Lineage, "s__"))
# sum the total of relative abundances at the species level
OTU_table_rel_species_sum<-OTU_table_rel_species%>%
  dplyr::select(-1)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  t()%>%
  as.data.frame()%>%
  rename(speciesfrac = V1)%>%
  rownames_to_column(.,var="flat_name")

# keep the species level classification data for coverage and datasize filtered 9820 samples
OTU_table_rel_species_sum_filtered <- inner_join(OTU_table_rel_species_sum, 
                                                 MfD_species_anno, 
                                                 by = "flat_name")%>%
  dplyr::select(c(1:2))

# merge species frac with metdata
OTU_table_rel_species_sum_metadata<-merge(OTU_table_rel_species_sum_filtered,
                                          mfd_metadata,by="flat_name",all=F)%>%
  dplyr::select(c(1:3,5,15,37:41,45))%>%
  mutate(group="Metapackage without MfD MAGs")%>%
  # for sampletype_areatype_hab1
  mutate(Habitat = paste(mfd_sampletype,mfd_areatype,mfd_hab1, sep = "_"))%>%
  # for areatype_hab1
  mutate(Habitat2 = paste(mfd_areatype,mfd_hab1, sep = "_"))%>%
  # for sampletype_hab1
  mutate(Habitat3 = paste(mfd_sampletype,mfd_hab1, sep = "_"))

write.table(OTU_table_rel_species_sum_metadata,"00-inputs/MFD_withoutSRMAGs_species_frac_data_filtered.txt",
            sep="\t",quote=FALSE,row.names = FALSE)


#################################################################################
##### For MfD samples, metapackage with SR MFD MAGs: species level coverage #####
#################################################################################
# count the total root coverage per sample
OTU_table_cov_new_sum<-OTU_table_cov_new%>%
  dplyr::select(-1)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  t()%>%
  as.data.frame()%>%
  rename(root_coverage = V1)
# count the total species coverage per sample
OTU_table_cov_new_species<-OTU_table_cov_new%>%
  filter(grepl("s__", Lineage, fixed = TRUE))%>%
  dplyr::select(-1)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  t()%>%
  as.data.frame()%>%
  rename(species_coverage = V1)

# merge root and species coverage and compute species_annotation
# 9820 samples remained after applying the same cutoffs of 50 coverage and 1Gbp data
MfD_species_anno_new<-merge(OTU_table_cov_new_sum,OTU_table_cov_new_species,by=0,all=T)%>%
  column_to_rownames(.,var="Row.names")%>%
  rownames_to_column(.,var="flat_name")%>%
  mutate(species_annotation=100*(species_coverage/root_coverage))%>%
  merge(.,mfd_metadata%>%dplyr::select(c(2,14,36:40)),by="flat_name",all=F)%>%
  mutate(mbases=after_total_bases/1000000)%>%
  filter(root_coverage>=50) %>%
  filter(mbases>=1000) 

# obtain the ratio of overall root coverage and the overall species coverage
# Calculate the sum of "species_coverage"
MFD_sum_species_coverage_new <- sum(MfD_species_anno_new$species_coverage, na.rm = TRUE)
# Calculate the sum of "root_coverage"
MFD_sum_root_coverage_new <- sum(MfD_species_anno_new$root_coverage, na.rm = TRUE)
# Calculate the overall ratio: 12.2%
MFD_species_ratio_new <- MFD_sum_species_coverage_new / MFD_sum_root_coverage_new


### extract data classified to the species level 
# Extract rows with "s__" in the "Lineage" column
OTU_table_rel_new_species <- OTU_table_rel_new %>%
  filter(str_detect(Lineage, "s__"))
# sum the total of relative abundances at the species level
OTU_table_rel_species_new_sum<-OTU_table_rel_new_species%>%
  dplyr::select(-1)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  t()%>%
  as.data.frame()%>%
  rename(speciesfrac = V1)%>%
  rownames_to_column(.,var="flat_name")


OTU_table_rel_species_sum_filtered_new <- inner_join(OTU_table_rel_species_new_sum, 
                                                     MfD_species_anno_new, 
                                                     by = "flat_name")%>%
  dplyr::select(c(1:2))

# merge species frac with metdata
OTU_table_rel_species_sum_metadata_new<-merge(OTU_table_rel_species_sum_filtered_new,
                                              mfd_metadata,by="flat_name",all=F)%>%
  dplyr::select(c(1:3,5,15,37:41,45))%>%
  mutate(group="Metapackage with MfD MAGs")%>%
  # for sampletype_areatype_hab1
  mutate(Habitat = paste(mfd_sampletype,mfd_areatype,mfd_hab1, sep = "_"))%>%
  # for areatype_hab1
  mutate(Habitat2 = paste(mfd_areatype,mfd_hab1, sep = "_"))%>%
  # for sampletype_hab1
  mutate(Habitat3 = paste(mfd_sampletype,mfd_hab1, sep = "_"))

write.table(OTU_table_rel_species_sum_metadata_new,"00-inputs/MFD_withSRMAGs_species_frac_data_filtered.txt",
            sep="\t",quote=FALSE,row.names = FALSE)


########################################################################################################
##### For MfD samples, merged species level coverage for metapackages with and without SR MFD MAGs #####
########################################################################################################

# merge the know species fraction of results from with SR MFD MAGs and without SR MFD MAGs metapackages
OTU_table_rel_species_sum_metadata_merged<-rbind(OTU_table_rel_species_sum_metadata%>%
                                                   select(c(1:3,6:12)),
                                                 OTU_table_rel_species_sum_metadata_new%>%
                                                   select(c(1:3,6:12)))
write.table(OTU_table_rel_species_sum_metadata_merged,"00-inputs/MFD_withSRMAGs_withoutSRMAGs_long_species_frac_data_filtered.txt",
            sep="\t",quote=FALSE,row.names = FALSE)


# Calculate stats of known species fraction improvement in different samples for MfD
OTU_table_rel_species_sum_metadata_samplecounts_stats<-OTU_table_rel_species_sum_metadata_merged %>%
  mutate(Dataset="MfD")%>%
  group_by(Dataset,group,complex)%>%
  dplyr::summarise(across(
    .cols = speciesfrac,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))
MfD_species_anno_merged_counts_stats_cal<-OTU_table_rel_species_sum_metadata_samplecounts_stats %>%
  select(c(2,3,6))%>%
  pivot_wider(names_from = "group",values_from = "speciesfrac_MD")%>%
  mutate(improvement=`Metapackage with MfD MAGs` -`Metapackage without MfD MAGs`)%>%
  select(c(2,5))

####################################################################################
##### For NCBI samples, metapackage without SR MFD MAGs: species level coverage ####
####################################################################################
# select soil, sediment, marine, freshwater, and aquatic samples and human samples
keep_vector <- c("sediment","marine sediment","freshwater sediment",
                 "human",
                 "soil", 
                 "marine","freshwater","aquatic","drinking water","freshwater","groundwater","lake water","wastewater")

# changing aquatic, freshwater, and marine, and additional all to Water, same for sediments, use human as reference
NCBI_species_anno_select<-NCBI_species_anno%>%
  filter(organism %in% keep_vector)%>%
  mutate(group="NCBI")%>%
  mutate(organism = case_when(
    organism == "soil" ~ "NCBI_Soil",
    
    organism == "sediment" ~ "NCBI_Sediment", 
    organism == "marine sediment" ~ "NCBI_Sediment", 
    organism == "freshwater sediment" ~ "NCBI_Sediment",
    
    organism == "marine" ~ "NCBI_Water",
    organism == "freshwater" ~ "NCBI_Water",
    organism == "aquatic" ~ "NCBI_Water",
    organism == "drinking water" ~ "NCBI_Water",
    organism == "groundwater" ~ "NCBI_Water",
    organism == "lake water" ~ "NCBI_Water",
    organism == "wastewater" ~ "NCBI_Water",
    
    organism == "human" ~ "NCBI_Human",
    TRUE ~ organism
  ))

# obtain the ratio of overall root coverage and the overall species coverage
# Calculate the sum of "species_coverage"
NCBI_species_anno_select_nohuman<-NCBI_species_anno_select%>%
  filter(!grepl("Human", organism))
NCBI_sum_species_coverage <- sum(NCBI_species_anno_select_nohuman$species_coverage, na.rm = TRUE)
# Calculate the sum of "root_coverage"
NCBI_sum_root_coverage <- sum(NCBI_species_anno_select_nohuman$root_coverage, na.rm = TRUE)
# Calculate the overall species fraction: 40.5%
NCBI_species_ratio <- NCBI_sum_species_coverage / NCBI_sum_root_coverage


###########################################################################################################
###### combine species fraction data for MfD and NCBI based on supplemented GTDB R214 metapackage #########
###########################################################################################################
MfD_species_anno2<-MfD_species_anno%>%
  mutate(host_or_not="ecological")%>%
  dplyr::select(c(1:4,6,12,11))%>%
  mutate(group="MFD")%>%
  mutate(mfd_sampletype = case_when(
    mfd_sampletype == "Soil" ~ "MfD_Soil",
    mfd_sampletype == "Sediment" ~ "MfD_Sediment",
    mfd_sampletype == "Water" ~ "MfD_Water",
    TRUE ~ mfd_sampletype
  ))
colnames(MfD_species_anno2)<-c("acc","root_coverage","species_coverage",
                               "species_annotation","organism","host_or_not",
                               "mbases","group")
MfD_species_anno2_sampleID<-MfD_species_anno2%>%
  select(c("acc"))

MfD_NCBI_species_frac<-rbind(NCBI_species_anno_select,MfD_species_anno2)

write.table(MfD_NCBI_species_frac,"00-inputs/NCBI_MFD_withoutSRMAGs_long_species_frac_data_filtered.txt",
            sep="\t",quote=FALSE,row.names = FALSE)

# some statistics
MfD_NCBI_species_frac_stats<-MfD_NCBI_species_frac%>% 
  group_by(group,organism)%>%
  dplyr::summarise(across(
    .cols = species_annotation,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2), 
                 RSD=~round(100*(sd(., na.rm = T))/(mean(., na.rm = T)),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))%>%
  t()
MfD_NCBI_species_frac_stats_overall<-MfD_NCBI_species_frac%>% 
  # focus on only non-host ecological samples
  filter(organism!="NCBI_Human")%>%
  group_by(group)%>%
  dplyr::summarise(across(
    .cols = species_annotation,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 RSD=~round(100*(sd(., na.rm = T))/(mean(., na.rm = T)),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))%>%
  t()
MfD_NCBI_species_frac_stats_sample<-MfD_NCBI_species_frac%>% 
  group_by(group,organism)%>%
  dplyr::summarise(across(
    .cols = species_annotation,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na.rm = T),2),
                 RSD=~round(100*(sd(., na.rm = T))/(mean(., na.rm = T)),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))%>%
  t()


########################################################################################
#### Plot Figure 5a: jittered-boxplot per habitat/complex: sampletype_areatype_hab1 ####
########################################################################################

### combine MfD metagenomic datasets and NCBI datasets
OTU_table_rel_species_sum_metadata_merge_simple<-OTU_table_rel_species_sum_metadata_merged%>%
  select(c("flat_name","speciesfrac","complex","group"))%>%
  mutate(Dataset="MfD")
NCBI_species_anno_select_simple<-NCBI_species_anno_select%>%
  select(c("acc","species_annotation","organism"))%>%
  mutate(group=NA)%>%
  mutate(Dataset="NCBI")
colnames(NCBI_species_anno_select_simple)<-c("flat_name","speciesfrac","complex","group","Dataset")
OTU_NCBI_species_anno_select_simple<-rbind(OTU_table_rel_species_sum_metadata_merge_simple,
                                           NCBI_species_anno_select_simple)

### Calculate stats of known species fraction in different samples for MfD and NCBI
OTU_NCBI_table_rel_species_sum_metadata_samplecounts_stats_cal<-OTU_NCBI_species_anno_select_simple %>%
  group_by(Dataset,group,complex)%>% 
  dplyr::summarise(across(
    .cols = speciesfrac,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))

### Plot jittered boxplot for species fraction
plot_jitter_speciesfrac<-ggplot(OTU_NCBI_species_anno_select_simple ,
                                         aes(x = reorder(complex, speciesfrac,median), 
                                             y=speciesfrac,color=complex,shape=group)) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2,position = position_dodge(0.8))+
  geom_jitter(data=OTU_NCBI_species_anno_select_simple %>%
                filter(!grepl("NCBI", complex)),
              mapping=aes(x = reorder(complex, speciesfrac,median),
                          y=speciesfrac,color=complex,
                          shape=group,
              ),
              size=0.15, position = position_jitterdodge(0.15)) + 
  theme_bw()+  
  scale_fill_manual(values=NCBI_MfD.palette,guide="none")+
  scale_color_manual(values=NCBI_MfD.palette,guide="none")+
  scale_shape_manual(values=c("Metapackage without MfD MAGs"=16,
                              "Metapackage with MfD MAGs"=16),guide = "none")+
  labs(y="Known species fraction (%)", x="Habitats")+
  scale_y_continuous(limits = c(0,120))+
  geom_text(data = OTU_NCBI_table_rel_species_sum_metadata_samplecounts_stats_cal, 
            aes(complex, y = max(speciesfrac_MD)+20 ,
                label = paste("MD =", round(speciesfrac_MD, 1)), 
                fontface = "bold",family="Arial"), 
            position = position_dodge(width = 0.8), size = 2, vjust = 1,hjust=0.1)+ 
  geom_text(data = OTU_NCBI_table_rel_species_sum_metadata_samplecounts_stats_cal, 
            aes(complex, y = max(speciesfrac_MD)+20 ,
                label = paste("n =",speciesfrac_N), fontface = "bold",family="Arial"), 
            position = position_dodge(width = 0.8), size = 2, vjust = -0.5,hjust=0.1)+ 
  theme(axis.text.x = element_text(family="Arial",size=9),
  legend.position = c(0.9,0.8),legend.background=element_blank(),
  legend.text = element_text(family="Arial",size=9),
  legend.title = element_text(family="Arial",size=12),
  axis.title  = element_text(family="Arial",size=12),
  plot.title = element_text(family="Arial"),
  axis.text.y= element_text(family="Arial",size =9))+ 
  coord_flip()


ggsave(plot_jitter_speciesfrac,
       filename = "02-output/Figure5a.png",
       width = 5,height = 5)

################################################################
##### For Figure 5b: Novelty of MAGs in different habitats #####
################################################################

### Load MAG taxonomy nased on GTDB-tk, R220 
MAG_taxonomy<-fread("00-inputs/MAG_SR_MFD_taxonomy.tsv",
                    header=F)%>%
  separate(.,col="V1",sep="\t",into=c("GenomeID","Domain"))
colnames(MAG_taxonomy)<-c("GenomeID","Domain","Phylum","Class","Order","Family","Genus","Species")
MAG_taxonomy_sepcies<-MAG_taxonomy%>%
  mutate(species_annotation =  if_else(Species == "s__", "Not_Species", "Species"))%>%
  select(c(1,9))

### Load MAG quality
MAG_quality<-fread("00-inputs/MAG_SR_MFD_quality.tsv")%>%
  select(c(1,28))%>%
  mutate(GenomeID=bin)%>%
  select(c("GenomeID","MAG_status"))

### Compute MAG per habitat by assemblies 
MAG_sample<-MAG_quality%>%
  separate(.,col="GenomeID",into=c("flat_name","binID"),sep="\\.", remove = FALSE)%>%
  merge(.,mfd_metadata%>%
          select(c("flat_name","complex")),
        by="flat_name",all=F)%>%
  merge(.,MAG_taxonomy_sepcies,by="GenomeID",all=F)%>%
  group_by(complex,MAG_status,species_annotation)%>%
  count()%>%
  mutate(total_genome_habitat_assembly=n)%>%
  select(-c("n"))

# Find total number of MAGs recovered for each habitat
MAG_sample_summary_heatmap_totaln<-MAG_quality%>%
  separate(.,col="GenomeID",into=c("flat_name","binID"),sep="\\.", remove = FALSE)%>%
  merge(.,mfd_metadata%>%
          select(c("flat_name","complex")),
        by="flat_name",all=F)%>%
  merge(.,MAG_taxonomy_sepcies,by="GenomeID",all=F)%>%
  group_by(complex)%>%
  count()%>%
  mutate(total_MAG=n)%>%
  select(-c("n"))

# Find total number of HQ & MQ MAGs recovered for each habitat
MAG_sample_summary_heatmap_HQn<-MAG_quality%>%
  separate(.,col="GenomeID",into=c("flat_name","binID"),sep="\\.", remove = FALSE)%>%
  merge(.,mfd_metadata%>%
          select(c("flat_name","complex")),
        by="flat_name",all=F)%>%
  merge(.,MAG_taxonomy_sepcies,by="GenomeID",all=F)%>%
  group_by(complex,MAG_status)%>%
  count()%>%
  pivot_wider(.,names_from = "MAG_status",values_from = "n")

# count where recovered species representative SR MAGs are from in different habitats
rep_list<-fread("00-inputs/MAG_rep_sample.list",
                header=F)
colnames(rep_list)<-"flat_name"
rep_list_metadata<-merge(rep_list,mfd_metadata,by="flat_name",all=F)%>%
  group_by(complex)%>% 
  count()%>%
  # order complex
  mutate(complex = factor(complex,levels = level_complex, ordered = TRUE))
colnames(rep_list_metadata)<-c("complex","n_rep")

### barplot showing MAG quality for those recovered from each habitat
# show in heatmap
# MAG novelty at different taxonomic levels by gtdb-tk
MAG_sample_summary_taxonomiclevels<-MAG_quality%>%
  separate(.,col="GenomeID",into=c("flat_name","binID"),sep="\\.", remove = FALSE)%>%
  merge(.,mfd_metadata%>%
          select(c("flat_name","complex")),
        by="flat_name",all=F)%>%
  merge(.,MAG_taxonomy,by="GenomeID",all=F)%>%
  mutate(Species_end = ifelse(Species == "s__"& Genus != "g__", 0, 1))%>%
  mutate(Genus_end = ifelse(Genus == "g__"& Family != "f__", 0, 1))%>%
  # everything above the family level
  mutate(Family_end = ifelse(Family == "f__", 0, 1))%>% 
  group_by(complex)%>%
  summarise(Species = sum(Species_end),
            Genus = sum(Genus_end),
            Family= sum(Family_end))%>%
  merge(.,MAG_sample_summary_heatmap_totaln,by=1,all=T)%>%
  mutate(Not_species=total_MAG-Species,
         Not_genus=total_MAG-Genus,
         Not_family_above=total_MAG-Family)%>%
  select(c(1,6:8))


# define levels to align with Figure 5a
level_complex<-c("Soil, Natural, Bogs, mires and fens",
                 "Sediment, Urban, Freshwater",
                 "Soil, Urban, Greenspaces",
                 "Sediment, Natural, Freshwater",
                 "Soil, Natural, Coastal",
                 "Soil, Agriculture, Fields",
                 "Soil, Natural, Grassland formations",
                 "Soil, Natural, Dunes",
                 "Soil, Natural, Forests",
                 "Water, Urban, Drinking water",
                 "Soil, Natural, Rocky habitats and caves",
                 "Sediment, Natural, Saltwater",
                 "Soil, Natural, Temperate heath and scrub",
                 "Sediment, Urban, Saltwater",
                 "Sediment, Subterranean, Saltwater",
                 "Water, Natural, Saltwater",
                 "Water, Urban, Wastewater",
                 "Water, Urban, Biogas")
levels_columns<-c("improvement","n_rep","total_MAG","HQ","MQ","Not_family_above",
                  "Not_genus","Not_species")
MAG_sample_summary_heatmap_input<-merge(as.data.frame(MfD_species_anno_merged_counts_stats_cal)%>%
                                          select(-c(1)),
                                        MAG_sample_summary_heatmap_totaln,by="complex",all=T)%>%
  merge(.,rep_list_metadata,by=1,all=T)%>%
  merge(.,MAG_sample_summary_heatmap_HQn,by=1,all=T)%>%
  merge(.,MAG_sample_summary_taxonomiclevels,by=1,all=T)%>%
  pivot_longer(.,cols=5:9,names_to ="items",values_to = "n_MAGs" )%>%
  # order complex
  mutate(complex = factor(complex,levels = level_complex, ordered = TRUE))%>%
  # order columns
  mutate(items = factor(items,levels = levels_columns, ordered = TRUE))

MAG_sample_summary_heatmap_p<-ggplot(MAG_sample_summary_heatmap_input,
                                     aes(x=items, 
                                         y = factor(complex, level = level_complex))) +
  geom_tile(aes(fill = n_MAGs), color="black") +
  geom_text(aes(label = round(n_MAGs,0))) +
  scale_fill_gradient(low = "white", high = "white")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x=element_text(angle=90, hjust=1, vjust = -0.5),
        axis.text.y = element_text(),
        axis.text = element_text(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  # put x axis text to top
  scale_x_discrete(position = "top") 

ggsave(MAG_sample_summary_heatmap_p,
       filename = "02-output/Figure5b.png",
       width = 8,height = 5.5)
ggsave(MAG_sample_summary_heatmap_p,
       filename = "02-output/Figure5b.pdf",
       width = 8,height = 5.5)




save.image("01-Figure5.Rdata")

