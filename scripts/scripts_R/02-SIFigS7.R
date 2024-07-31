# This script recreates Supplementary Figure 7 in the big MfD paper 
# Script created and last updated by Yu Yang on Jul 31 2024

################################
##### Load needed packages ##### 
################################
library(data.table) 
library(ggplot2)
library(dplyr)
library(wesanderson)
library(tidyverse)
library("cowplot")

########################################
############## load data ###############
########################################
### Load MFD metadata
mfd_metadata <- fread("00-inputs/merged_metadata_20240217.csv")

### Load results of MAG coverage of micorbial community
community_cov<-fread("00-inputs/SIFig7_merged_singlem_appraise.txt")%>%
  mutate(group="Prokaryotes")
names(community_cov)[1] <- "flat_name"

### Merge bacteria and archaea results and metadata
community_cov_merged<-community_cov %>%
  group_by(flat_name,group)%>%
  summarise(total_num_not_found = sum(num_not_found, na.rm = TRUE),
            total_num_binned = sum(num_binned, na.rm = TRUE))%>%
  mutate(total=total_num_not_found+total_num_binned)%>%
  mutate(prokaryotic_cov=100*total_num_binned/total)
community_cov_merged_meta<-merge(mfd_metadata,community_cov_merged,
                                 by="flat_name",all=F)
#################################
##### Create color palettes #####
#################################
# This was created by Thomas Bygh Nymann Jensen, Aalborg University
## Ontology palettes
sediment.palette <- colorRampPalette(c(wes_palette("IsleofDogs2")[3], wes_palette("FantasticFox1")[1]))
plot(rep(1, 5), col = sediment.palette(5), pch = 19, cex = 3)
soil.palette <- colorRampPalette(c(wes_palette("AsteroidCity1")[4], wes_palette("AsteroidCity1")[1]))
plot(rep(1, 9), col = soil.palette(9), pch = 19, cex = 3)
water.palette <- colorRampPalette(c(wes_palette("Darjeeling2")[2], wes_palette("Zissou1")[2]))
plot(rep(1, 4), col = water.palette(4), pch = 19, cex = 3)

## Sampletype palette

sampletype.palette <- c(soil.palette(1), sediment.palette(1), water.palette(1))
names(sampletype.palette) <- c("Soil", "Sediment", "Water")

sampletype.palette
plot(rep(1, 3), col = sampletype.palette, pch = 19, cex = 3)


mfdo1 <- readxl::read_excel('/projects/microflora_danica/project_mfd_ARGs/01-Data/00-metadata/2024-02-13_mfd_db.xlsx') %>%
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

sample_area_palette<-c("Sediment_Natural"="#B6854D","Sediment_Subterranean"="#C9893B","Sediment_Urban"="#DD8D29",
                       "Soil_Agriculture"="#6C8645", "Soil_Natural"="#4B8E62", "Soil_Subterranean"="#2A967F",
                       "Soil_Urban"="#0A9F9D","Water_Natural"="#046C9A" ,"Water_Subterranean"="#3E91AF",
                       "Water_Urban"="#78B7C5")


############################################################
####################### Calculations #######################
############################################################
# overall
community_cov_merged_meta_stats_overall<-community_cov_merged_meta %>%
  group_by(group)%>% 
  dplyr::summarise(across(
    .cols = prokaryotic_cov,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))

# per sampletype per group
community_cov_merged_meta_stats_sampletype<-community_cov_merged_meta %>%
  group_by(mfd_sampletype,group)%>% 
  dplyr::summarise(across(
    .cols = prokaryotic_cov,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))

# per sampletype_areatype per group
community_cov_merged_meta_stats_sample_areatype<-community_cov_merged_meta %>%
  mutate(sample_area=paste(mfd_sampletype,mfd_areatype,sep="_"))%>%
  group_by(sample_area,group)%>% 
  dplyr::summarise(across(
    .cols = prokaryotic_cov,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))

# per Habitat1 per group
community_cov_merged_meta_stats_sample_area_Habitat<-community_cov_merged_meta %>%
  group_by(complex,group)%>% 
  dplyr::summarise(across(
    .cols = prokaryotic_cov,
    .fns =  list(N=~dplyr::n(),mean = ~round(mean(., na.rm = T),2), 
                 MD = ~round(median(., na.rm = T),2),sd=~round(sd(.,na=T),2),
                 quantile_10=~round(quantile(.,0.1, na.rm = TRUE),2),
                 quantile_25=~round(quantile(.,0.25, na.rm = TRUE),2),
                 quantile_75=~round(quantile(.,0.75, na.rm = TRUE),2),
                 quantile_90=~round(quantile(.,0.9, na.rm = TRUE),2))))%>%
  tidyr::separate(complex, into = c("mfd_sampletype","mfd_areatype" ,"area_hab1"), 
                  sep = ",", remove = FALSE,extra = "merge")

############################################################
################### Plot jitter plot ####################### 
############################################################
# different sampletypes
plot_jitter_shortreadMAGcov_merged_sample<-ggplot2::ggplot(community_cov_merged_meta ,
                                                           aes(x = reorder(mfd_sampletype, prokaryotic_cov,median), 
                                                               y=prokaryotic_cov,color=mfd_sampletype)) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2)+
  geom_jitter(shape=16,size=0.2, position=position_jitter(0.2))+  
  theme_bw()+
  scale_color_manual(values=sampletype.palette)+
  labs(y="Short-read MAG prokaryotic coverage (%)", x="",color="Sample types")+
  geom_text(data = community_cov_merged_meta_stats_sampletype, 
            aes(mfd_sampletype, y = max(prokaryotic_cov_MD)+65 ,
                label = paste("MD =", round(prokaryotic_cov_MD, 1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+ 
  geom_text(data = community_cov_merged_meta_stats_sampletype, 
            aes(mfd_sampletype, y = max(prokaryotic_cov_MD)+60 ,
                label = paste("N =",prokaryotic_cov_N), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+ 
  geom_text(data = community_cov_merged_meta_stats_sampletype, 
            aes(mfd_sampletype, y = max(prokaryotic_cov_MD)+55 ,
                label = paste("sd =",round(prokaryotic_cov_sd,1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+
  theme(axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        legend.text = element_text(size=9),
        legend.title = element_text(size=12),
        axis.title.x  = element_text(size=12),
        axis.title.y =  element_text(size=8),
        plot.title = element_text(),
        axis.text.y= element_text(size =9))

# different sample_areatypes
# reorder x axis by the median of each group first
community_cov_merged_meta_new<-community_cov_merged_meta %>%
  mutate(sample_area=paste(mfd_sampletype,mfd_areatype,sep="_"))
community_cov_merged_meta_new_medians <- aggregate(prokaryotic_cov ~ sample_area, 
                                                   community_cov_merged_meta_new, median)
ordered_sample_area <- community_cov_merged_meta_new_medians$sample_area[order(community_cov_merged_meta_new_medians$prokaryotic_cov)]
community_cov_merged_meta_new$sample_area <- factor(community_cov_merged_meta_new$sample_area, levels = ordered_sample_area)

plot_jitter_shortreadMAGcov_merged_sample_area<-ggplot2::ggplot(community_cov_merged_meta_new,
                                                                aes(x = sample_area, 
                                                                    y=prokaryotic_cov,color=sample_area,
                                                                )) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2)+
  geom_jitter(shape=16,size=0.2, position=position_jitter(0.2))+  
  theme_bw()+
  scale_color_manual(values=sample_area_palette)+
  labs(y="Short-read MAG prokaryotic coverage (%)", x="")+
  geom_text(data = community_cov_merged_meta_stats_sample_areatype, 
            aes(sample_area, y = max(prokaryotic_cov_MD)+63 ,
                label = paste("MD =", round(prokaryotic_cov_MD, 1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+ 
  geom_text(data = community_cov_merged_meta_stats_sample_areatype, 
            aes(sample_area, y = max(prokaryotic_cov_MD)+58 ,
                label = paste("N =",prokaryotic_cov_N), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+ 
  geom_text(data = community_cov_merged_meta_stats_sample_areatype, 
            aes(sample_area, y = max(prokaryotic_cov_MD)+53 ,
                label = paste("sd =",round(prokaryotic_cov_sd,1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+
  theme(axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        legend.text = element_text(size=9),
        legend.title = element_text(size=12),
        axis.title  = element_text(size=8),
        plot.title = element_text(),
        axis.text.y= element_text(size =9))

# different habitat
level_complex<-c("Soil, Natural, Bogs, mires and fens","Soil, Natural, Coastal","Sediment, Urban, Freshwater",
                 "Sediment, Natural, Freshwater","Soil, Urban, Greenspaces","Soil, Agriculture, Fields",
                 "Soil, Natural, Grassland formations","Water, Urban, Drinking water","Soil, Natural, Dunes",
                 "Sediment, Natural, Saltwater","Soil, Natural, Forests","Sediment, Urban, Saltwater",
                 "Soil, Natural, Rocky habitats and caves","Soil, Natural, Temperate heath and scrub",
                 "Sediment, Subterranean, Saltwater","Water, Natural, Saltwater","Water, Urban, Wastewater",
                 "Water, Urban, Biogas")
plot_jitter_shortreadMAGcov_merged<-ggplot2::ggplot(community_cov_merged_meta %>%
                                                      # order complex
                                                      mutate(complex = factor(complex,levels = level_complex, ordered = TRUE)),
                                                    aes(# order by value
                                                      forcats::fct_reorder(complex, 
                                                                           prokaryotic_cov,median), 
                                                      y=prokaryotic_cov,color=complex)) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2)+
  geom_jitter(shape=16,size=0.2, position=position_jitter(0.2))+  
  theme_bw()+
  scale_color_manual(values=mfdo1.palette)+
  scale_y_continuous(limits = c(0,100))+
  labs(y="Short-read MAG prokaryotic coverage (%)", x="Habitats",color="Sample types")+
  geom_text(data = community_cov_merged_meta_stats_sample_area_Habitat, 
            aes(complex, y = max(prokaryotic_cov_MD)+29 ,
                label = paste("MD =", round(prokaryotic_cov_MD, 1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.4, vjust = 0)+ 
  geom_text(data = community_cov_merged_meta_stats_sample_area_Habitat, 
            aes(complex, y = max(prokaryotic_cov_MD)+19 ,
                label = paste("N =",prokaryotic_cov_N), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.4, vjust = 0)+ 
  geom_text(data = community_cov_merged_meta_stats_sample_area_Habitat,
            aes(complex, y = max(prokaryotic_cov_MD)+24 ,
                label = paste("sd =",round(prokaryotic_cov_sd,1)), fontface = "bold"),
            position = position_dodge(width = 0.8), size = 1.4, vjust = 0)+
  theme(axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none",
    legend.text = element_text(size=9),
    legend.title = element_text(size=12),
    axis.title  = element_text(size=10),
    plot.title = element_text(),
    axis.text.y= element_text(size =9))#+
coord_flip()


######################
#### merge figure ####
######################

plot_jitter_recovery_merged_prokaryotes<-ggdraw() +
  draw_plot(plot_jitter_shortreadMAGcov_merged_sample, x = 0, y = 0.62, width = .3, height = 0.37) +
  draw_plot(plot_jitter_shortreadMAGcov_merged_sample_area, x = .35, y = 0.55, width = .6, height = 0.45) +
  draw_plot(plot_jitter_shortreadMAGcov_merged, x = 0, y = 0, width = 1, height = 0.55) +
  draw_plot_label(label = c("a", "b", "c"), size = 15,
                  x = c(0, 0.30, 0), y = c(1, 1, 0.6))
ggsave(plot_jitter_recovery_merged_prokaryotes,
       filename = "02-output/SIFigS7.png",
       width = 7,height = 8)

save.image("02-SIFigS7.Rdata")
