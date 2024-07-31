# This script recreates Supplementary Figure 6 in the big MfD paper 
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

### Load microbial fraction results from SingleM
read_frac<-fread("00-inputs/SIFig6_merged_singlem_microbial_frac.txt",
                 header=T, sep="\t")
# Remove "%" and convert to numeric
read_frac$microbial_fraction <- as.numeric(gsub("%", "", read_frac$microbial_fraction))
# Keep samples without microbial fraction warnings and merge with metadata information
read_frac_select<-read_frac %>%
  dplyr::select(c(1,4,5))%>%
  filter(WARNING=="")
read_frac_select_metadata<-merge(read_frac_select, mfd_metadata, 
                                 by.x="sample",by.y="flat_name",all=F)
read_frac_metadata<-read_frac_select_metadata

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

################################################################
###### plot boxplot per habitat: sampletype_areatype_hab1 ###### 
################################################################
# Calculate stats of bac_arch read_frac in different samples
read_frac_metadata_Median_full <- read_frac_metadata %>%
  # if remove samples with WARNING in microbial_fraction
  filter(WARNING=="")%>%
  #dplyr::filter(mfd_sampletype!="soil")%>% # remove 10 samples currently labelled as "soil", instead of "Soil" in metadata
  dplyr::select(sample,microbial_fraction,complex,mfd_sampletype,mfd_hab1,mfd_hab2,mfd_hab3,mfd_areatype,WARNING,complex)%>%
  dplyr::summarise(dplyr::group_by(., mfd_sampletype),
                   total_MD = median(microbial_fraction), 
                   na.rm = T)%>%
  mutate(Habitat = paste(mfd_sampletype,mfd_areatype,mfd_hab1, sep = "_"))%>%
  mutate(sample_area = paste(mfd_sampletype,mfd_areatype, sep = "_"))

read_frac_metadata_Median_habitat<-read_frac_metadata_Median_full%>% 
  dplyr::group_by(complex,mfd_sampletype,mfd_areatype,mfd_hab1) %>% 
  dplyr::summarise(across(
    .cols = microbial_fraction,
    .fns =  list(n=~dplyr::n(), mean = ~mean(.,na.rm = TRUE),MD = ~median(.,na.rm = TRUE),
                 sd=~sd(.,na.rm=T))))%>%
  mutate(Habitat = paste(mfd_sampletype,mfd_areatype,mfd_hab1, sep = "_"))

plot_jitter_readfrac_complex<-ggplot2::ggplot(read_frac_metadata_Median_full %>%
                                                #if remove those samples with WARNING in microbial_fraction
                                                filter(WARNING=="")%>%
                                                dplyr::select(sample,microbial_fraction,complex,mfd_sampletype,
                                                              mfd_areatype,mfd_hab1,mfd_hab2,mfd_hab3,Habitat),
                                              aes(x = reorder(complex, microbial_fraction,median), # mfd_sampletype or Habitat
                                                  y=microbial_fraction,color=complex)) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2)+
  geom_jitter(shape=16,size=0.2, position=position_jitter(0.2))+  
  theme_bw()+
  scale_color_manual(values=mfdo1.palette)+
  labs(y="Microbial read fraction (%)", x="",color="Sample types")+
  geom_text(data = read_frac_metadata_Median_habitat, 
            aes(complex, y = max(microbial_fraction_MD)+26 ,# mfd_sampletype or Habitat
                label = paste("MD =", round(microbial_fraction_MD, 1)), 
                fontface = "bold",), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+ 
  geom_text(data = read_frac_metadata_Median_habitat, 
            aes(complex, y = max(microbial_fraction_MD)+21 ,# mfd_sampletype or Habitat
                label = paste("N =",microbial_fraction_n), fontface = "bold",
                ), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+ 
  geom_text(data = read_frac_metadata_Median_habitat, 
            aes(complex, y = max(microbial_fraction_MD)+16 ,# mfd_sampletype or Habitat
                label = paste("sd =",round(microbial_fraction_sd,1)), fontface = "bold",
                ), 
            position = position_dodge(width = 0.8), size = 1.7, vjust = 0)+
  theme(axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        legend.text = element_text(size=9),
        legend.title = element_text(size=12),
        axis.title  = element_text(size=12),
        plot.title = element_text(),
        axis.text.y= element_text(size =9))

###########################################################
###### plot boxplot per habitat: sampletype_areatype ###### 
###########################################################
# Calculate stats of bac_arch read_frac in different sample_areas
read_frac_metadata_Median_sample_area<-read_frac_metadata_Median_full%>% 
  dplyr::group_by(sample_area) %>% 
  dplyr::summarise(across(
    .cols = microbial_fraction,
    .fns =  list(n=~dplyr::n(), mean = ~mean(.,na.rm = TRUE),MD = ~median(.,na.rm = TRUE),
                 sd=~sd(.,na.rm=T))))
# reorder x axis by the median of each group first
read_frac_metadata_Median_sample_area_medians <- aggregate(microbial_fraction ~ sample_area, 
                                                           read_frac_metadata_Median_full, median)
ordered_sample_area <- read_frac_metadata_Median_sample_area_medians$sample_area[order(read_frac_metadata_Median_sample_area_medians$microbial_fraction)]
read_frac_metadata_Median_full$sample_area <- factor(read_frac_metadata_Median_full$sample_area, levels = ordered_sample_area)

plot_jitter_readfrac_sample_area<-ggplot2::ggplot(read_frac_metadata_Median_full %>%
                                                    #if remove those samples with WARNING in microbial_fraction
                                                    filter(WARNING=="")%>%
                                                    dplyr::select(sample,microbial_fraction,mfd_sampletype,mfd_hab1,mfd_hab2,mfd_hab3,Habitat,sample_area),
                                                  aes(x = sample_area, # mfd_sampletype or Habitat
                                                      y=microbial_fraction,color=sample_area
                                                  )) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2)+
  geom_jitter(shape=16,size=0.2, position=position_jitter(0.2))+  
  theme_bw()+
  scale_color_manual(values=sample_area_palette)+
  #scale_color_manual(values=sampletype.palette)+
  labs(y="Microbial read fraction (%)",x="",#,color="Sample types"
  )+
  geom_text(data = read_frac_metadata_Median_sample_area, 
            aes(sample_area, y = max(microbial_fraction_MD)+31 ,# mfd_sampletype or Habitat
                label = paste("MD =", round(microbial_fraction_MD, 1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.5, vjust = 0)+ 
  geom_text(data = read_frac_metadata_Median_sample_area, 
            aes(sample_area, y = max(microbial_fraction_MD)+26 ,# mfd_sampletype or Habitat
                label = paste("N =",microbial_fraction_n), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.5, vjust = 0)+ 
  geom_text(data = read_frac_metadata_Median_sample_area, 
            aes(sample_area, y = max(microbial_fraction_MD)+21 ,# mfd_sampletype or Habitat
                label = paste("sd =",round(microbial_fraction_sd,1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.5, vjust = 0)+
  theme(#axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=9),
    legend.title = element_text(size=12),
    axis.title.y  = element_text(size=12),
    axis.title.x  = element_blank(),
    plot.title = element_text(),
    axis.text.y= element_text(size =9),
    axis.text.x= element_text(size =9,angle = 90, vjust = 0.5, hjust=1))


#############################################
###### plot boxplot per mfd_sampletype ###### 
#############################################
read_frac_metadata_Median_sampletype<-read_frac_metadata_Median_full%>% 
  dplyr::group_by(mfd_sampletype) %>% 
  dplyr::summarise(across(
    .cols = microbial_fraction,
    .fns =  list(n=~dplyr::n(), mean = ~mean(.,na.rm = TRUE),MD = ~median(.,na.rm = TRUE),
                 sd=~sd(.,na.rm=T))))

plot_jitter_readfrac_sample<-ggplot2::ggplot(read_frac_metadata_Median_full %>%
                                               #if remove those samples with WARNING in microbial_fraction
                                               filter(WARNING=="")%>%
                                               #filter(mfd_sampletype!="soil")%>% # remove 10 samples currently labelled as "soil", instead of "Soil" in metadata
                                               dplyr::select(sample,microbial_fraction,mfd_sampletype,mfd_hab1,mfd_hab2,mfd_hab3,Habitat),
                                             aes(x = mfd_sampletype, # mfd_sampletype or Habitat
                                                 y=microbial_fraction,color=mfd_sampletype)) + 
  geom_boxplot(outlier.shape =16,outlier.size=0.2,lwd=0.2)+
  geom_jitter(shape=16,size=0.2, position=position_jitter(0.2))+  
  theme_bw()+
  scale_color_manual(values=sampletype.palette)+
  labs(y="Microbial read fraction (%)",color="Sample types",x="")+
  geom_text(data = read_frac_metadata_Median_sampletype, 
            aes(mfd_sampletype, y = max(microbial_fraction_MD)+31 ,# mfd_sampletype or Habitat
                label = paste("MD =", round(microbial_fraction_MD, 1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.5, vjust = 0)+ 
  geom_text(data = read_frac_metadata_Median_sampletype, 
            aes(mfd_sampletype, y = max(microbial_fraction_MD)+26 ,# mfd_sampletype or Habitat
                label = paste("N =",microbial_fraction_n), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.5, vjust = 0)+ 
  geom_text(data = read_frac_metadata_Median_sampletype, 
            aes(mfd_sampletype, y = max(microbial_fraction_MD)+21 ,# mfd_sampletype or Habitat
                label = paste("sd =",round(microbial_fraction_sd,1)), fontface = "bold"), 
            position = position_dodge(width = 0.8), size = 1.5, vjust = 0)+
  theme(#axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=9),
    legend.title = element_text(size=12),
    axis.title.y  = element_text(size=12),
    axis.title.x  = element_blank(),
    plot.title = element_text(),
    axis.text.y= element_text(size =9),
    axis.text.x= element_text(size =9,angle = 90, vjust = 0.5, hjust=1))


######################
#### merge figure ####
######################
top_row <- plot_grid(plot_jitter_readfrac_sample, plot_jitter_readfrac_sample_area, 
                     labels = c('a', 'b'), label_size = 12,
                     align = 'h',
                     rel_widths = c(0.8,2))
plot_jitter_readfrac_merged_grid<-plot_grid(top_row, 
                                            plot_jitter_readfrac_complex,
                                            rel_heights = c(1.2,1.7), 
                                            labels = c('', 'c'), label_size = 12, ncol = 1)

ggsave(plot_jitter_readfrac_merged_grid,
       filename = "02-output/SIFigS6.png",
       width = 8,height = 8)

save.image(file="03-SIFigS6.Rdata")

