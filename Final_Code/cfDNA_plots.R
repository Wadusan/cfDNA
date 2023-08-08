###Libraries
library(readr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)
library(kableExtra)
library(janitor)
library(readxl)
library(ggrepel)
library(viridis)
library(magick)
library(gridExtra)
###Load data
metadata <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/final_metadata.csv")
fragment_lengths <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/fragment_lengths/fragment_lengths.csv")
ichorCNA_summary_finaledb <- read_table("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/ichorCNA/ichorCNA_summary_finaledb_13.tsv")
cristiano_id_sample_metadata <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/cristiano_id_sample_metadata.csv")
cristiano_sup_table <- read_delim("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/cristiano_sup_table.csv", 
                                  delim = ";", escape_double = FALSE, trim_ws = TRUE, 
                                  skip = 1)
ichorCNA_unfiltered <- read_table("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/ichorCNA/ichorCNA_summary_finaledb_04_unfiltered.tsv")
ichorCNA_filtered <- read_table("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/ichorCNA/ichorCNA_summary_finaledb_05_filtered.tsv")
cristiano_sup_table_6 <- read_excel("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/ichorCNA/41586_2019_1272_MOESM2_ESM.xlsx", 
                                          sheet = "6", skip = 1)
liq_unfiltered <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/LIQUORICE/LIQUORICE_537.csv")
liq_150 <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/LIQUORICE/summary_across_samples_and_ROIS.csv")
liq_top50 <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/LIQUORICE/liq_top50.csv")
SAEC_hg38 <- read_delim("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/LIQUORICE/SAEC_hg38.bed", 
                        delim = "\t", escape_double = FALSE, 
                        col_names = FALSE, trim_ws = TRUE)
FrEIA_mouliere_tnc_and_panel <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/FrEIA/Cristiano_FrEIA_score_mouliere_tnc_and_panel.csv")
FrEIA_my_panel_and_mouliere_tnc <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/FrEIA/Cristiano_FrEIA_score_my_panel_and_mouliere_tnc.csv")
FrEIA_my_tnc_mouliere_panel <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/FrEIA/Cristiano_FrEIA_score_my_tnc_mouliere_panel.csv")
FrEIA_my_panel_and_tnc <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/FrEIA/Cristiano_FrEIA_score_my_panel_and_tnc.csv")
Cristiano_logFC_p_my_panel_and_tnc <- read_csv("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/FrEIA/Cristiano_logFC_p_my_panel_and_tnc.csv")
###FRAGMENT LENGTH###

###Fragment length Data wrangling
#Distribution of healthy vx cancer
frag_count <- fragment_lengths[,-1]
frag_ratio <- frag_count %>% mutate_at(vars(-length), funs(./sum(.)))
names(metadata)[1] <- "sample"
names(metadata)[2] <- "group"

frag_ratiolong <- melt(frag_ratio, 
                       id.vars=c("length"),
                       variable.name="sample",
                       value.name="ratio")
frag_ratiolong <- merge(frag_ratiolong,metadata,by="sample")
frag_ratiolong <- frag_ratiolong[- grep("EE87920", frag_ratiolong$sample),]
frag_ratiolong <- frag_ratiolong[- grep("EE87921", frag_ratiolong$sample),]

frag_ratiolong2groups <- frag_ratiolong

frag_ratiolong2groups$group[frag_ratiolong2groups$group != "Healthy"] <- "Cancer"
#Ratio of short fragments overall(90-150)
frag_ratio_short <- frag_ratio[1:61,]
frag_ratio_short_sum <- as.data.frame(colSums(frag_ratio_short[,-1]))
frag_ratio_short_sum$sample <- row.names(frag_ratio_short_sum)
rownames(frag_ratio_short_sum) <- NULL
names(frag_ratio_short_sum)[1] <- "ratio_short_sum"
frag_ratio_short_sum <- frag_ratio_short_sum[c("sample", "ratio_short_sum")]
frag_ratio_short_sum <- merge(frag_ratio_short_sum,metadata,by="sample")
frag_ratio_short_sum$group[frag_ratio_short_sum$group != "Healthy"] <- "Cancer"
frag_ratio_short_sum <- frag_ratio_short_sum[- grep("EE87920", frag_ratio_short_sum$sample),]
frag_ratio_short_sum <- frag_ratio_short_sum[- grep("EE87921", frag_ratio_short_sum$sample),]
c_frag_ratio_short_sum <- frag_ratio_short_sum[- grep("Healthy", frag_ratio_short_sum$group),]
h_frag_ratio_short_sum <- frag_ratio_short_sum[- grep("Cancer", frag_ratio_short_sum$group),]
median(c_frag_ratio_short_sum$ratio_short_sum)
median(h_frag_ratio_short_sum$ratio_short_sum)

#average length
t_frag_count <- as.data.frame(t(frag_count))
t_frag_count <- t_frag_count %>% row_to_names(row_number = 1)
t_frag_count$sample <- row.names(t_frag_count)
t_frag_count <- merge(t_frag_count,metadata,by="sample")
t_frag_count$group[t_frag_count$group != "Healthy"] <- "Cancer"
c_frag_count <- t_frag_count[- grep("Healthy", t_frag_count$group),]
c_frag_count <- select(c_frag_count, -c(group))
c_frag_count <- select(c_frag_count, -c(sample))
c_frag_count <- as.data.frame(colSums(c_frag_count))
c_frag_count$length <- row.names(c_frag_count)
names(c_frag_count)[1] <- "sum"
c_frag_count$length <- as.numeric(c_frag_count$length)
c_frag_count$prod <- c_frag_count$sum * c_frag_count$length
c_avg_bp <- sum(c_frag_count$prod)/sum(c_frag_count$sum)

h_frag_count <- t_frag_count[- grep("Cancer", t_frag_count$group),]
h_frag_count <- select(h_frag_count, -c(group))
h_frag_count <- select(h_frag_count, -c(sample))
h_frag_count <- as.data.frame(colSums(h_frag_count))
h_frag_count$length <- row.names(h_frag_count)
names(h_frag_count)[1] <- "sum"
h_frag_count$length <- as.numeric(h_frag_count$length)
h_frag_count$prod <- h_frag_count$sum * h_frag_count$length
h_avg_bp <- sum(h_frag_count$prod)/sum(h_frag_count$sum)
###Fragmentlength Plots
#Distribution of healthy vs cancer
frag_len_p1<-ggplot(frag_ratiolong2groups, 
                aes(x=length, y=ratio, group=sample)) +
  geom_line(aes(color=group),alpha=0.15) + xlim(80,250) + theme_classic() +
  scale_color_manual(values = c("red","blue")) +
  labs(y= "Frequency", x = "Fragment length(bp)") +
  guides(col = guide_legend(override.aes = list(alpha = 1))) + labs(tag = "A")

frag_len_p1
#Ratio of short fragments overall(90-150)
frag_len_p2 <- ggplot(frag_ratio_short_sum, aes(x=group, y=ratio_short_sum,color=group,))+ theme_classic() + 
  scale_color_manual(values = c("red","blue")) +
  labs(y= "Frequency of short fragments", x = "") +
  geom_violin() + labs(tag = "B")

frag_len_p2 <- frag_len_p2 +geom_boxplot(width=0.1)+ stat_compare_means(method = "t.test",label.x = 1.3) + theme(legend.position = "none")
#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

###ichorCNA###

###ichorCNA Data wrangling
#Tumor fraction overview of all cohorts
ichorCNA_summary_finaledb$sample <- gsub(".sortByCoord","",as.character(ichorCNA_summary_finaledb$sample))
ichordf <-  merge(x=ichorCNA_summary_finaledb,y=metadata,by="sample",all.x = TRUE)
ichordf <- ichordf[- grep("EE87920", ichordf$sample),]
ichordf <- ichordf[- grep("EE87921", ichordf$sample),]
ichordf1 <- ichordf
group_ordered <- with(ichordf1,
                      reorder(group,
                              tfx,
                              median))
#Comparison cancer stage and tumor fraction
ichordf <- merge(x=ichordf,y=cristiano_id_sample_metadata,by="sample",all.x = TRUE)
ichordf <- merge(x=ichordf,y=cristiano_sup_table,by="Patient",all.x = TRUE)
ichordf2<- subset(ichordf, Stage %in% c("I", "II", "III", "IV"))
#Longitudinal data
ichordf3 <- ichordf %>% select("Patient","sample","tfx","Patient Type","Timepoint")
ichordf3 <- ichordf3[grep("Day", ichordf3$Timepoint), ]
ichordf3$Timepoint <- gsub("Pre-treatment, Day ","",as.character(ichordf3$Timepoint))
ichordf3$Timepoint <- gsub("Post-treatment, Day ","",as.character(ichordf3$Timepoint))

ichordf3 <- transform(ichordf3, Timepoint = as.numeric(Timepoint))
cristiano_sup_table_6 <- cristiano_sup_table_6 %>% select("Patient","Maximum Mutant Allele Fractionᶥ")
ichordf3 <- merge(x=ichordf3,y=cristiano_sup_table_6,by="Patient",all.x = TRUE)
freia_time <- FrEIA_my_panel_and_tnc %>% select("sample","FrEIA_score")
ichordf3 <- merge(x=ichordf3,y=freia_time,by="sample",all.x = TRUE)
size_time <- frag_ratio_short_sum  %>% select("sample","ratio_short_sum")
ichordf3 <- merge(x=ichordf3,y=size_time,by="sample",all.x = TRUE)
liq_150 <- select(liq_150,"sample","region-set","Dip depth: z-score vs controls in same region set")
names(liq_150)[names(liq_150) == "Dip depth: z-score vs controls in same region set"] <- "Dip_depth"
names(liq_150)[names(liq_150) == "region-set"] <- "region"
liq_150 <- liq_150 %>% pivot_wider(names_from = region, values_from = c(Dip_depth))
liq_150$sample <- gsub(".sortByCoord","",as.character(liq_150$sample))
liq_time <- liq_150 %>% select("sample","SAEC_hg38")
ichordf3 <- merge(x=ichordf3,y=liq_time,by="sample",all.x = TRUE)
names(ichordf3)[names(ichordf3) == "Maximum Mutant Allele Fractionᶥ"] <- "maf"
names(ichordf3)[names(ichordf3) == "FrEIA_score"] <- "freia"
names(ichordf3)[names(ichordf3) == "ratio_short_sum"] <- "shortratio"
names(ichordf3)[names(ichordf3) == "SAEC_hg38"] <- "dipsaec"


CGPLLU13 <-ichordf3[grep("CGPLLU13", ichordf3$Patient), ]
CGPLLU14 <-ichordf3[grep("CGPLLU14", ichordf3$Patient), ]
CGPLLU244 <-ichordf3[grep("CGPLLU244", ichordf3$Patient), ]
CGPLLU245 <-ichordf3[grep("CGPLLU245", ichordf3$Patient), ]
CGPLLU246 <-ichordf3[grep("CGPLLU246", ichordf3$Patient), ]
CGPLLU264 <-ichordf3[grep("CGPLLU264", ichordf3$Patient), ]
CGPLLU265 <-ichordf3[grep("CGPLLU265", ichordf3$Patient), ]
CGPLLU266 <-ichordf3[grep("CGPLLU266", ichordf3$Patient), ]
CGPLLU267 <-ichordf3[grep("CGPLLU267", ichordf3$Patient), ]
CGPLLU269 <-ichordf3[grep("CGPLLU269", ichordf3$Patient), ]
CGPLLU271 <-ichordf3[grep("CGPLLU271", ichordf3$Patient), ]
CGPLLU43 <-ichordf3[grep("CGPLLU43", ichordf3$Patient), ]
CGPLLU86 <-ichordf3[grep("CGPLLU86", ichordf3$Patient), ]
CGPLLU88 <-ichordf3[grep("CGPLLU88", ichordf3$Patient), ]
CGPLLU89 <-ichordf3[grep("CGPLLU89", ichordf3$Patient), ]
library(pandas)
CGPLLU13 <- CGPLLU13[order(CGPLLU13$Timepoint),]
CGPLLU14 <- CGPLLU14[order(CGPLLU14$Timepoint),]
CGPLLU244 <- CGPLLU244[order(CGPLLU244$Timepoint),]
CGPLLU245 <- CGPLLU245[order(CGPLLU245$Timepoint),]
CGPLLU246 <- CGPLLU246[order(CGPLLU246$Timepoint),]
CGPLLU264 <- CGPLLU264[order(CGPLLU264$Timepoint),]
CGPLLU265 <- CGPLLU265[order(CGPLLU265$Timepoint),]
CGPLLU266 <- CGPLLU266[order(CGPLLU266$Timepoint),]
CGPLLU267 <- CGPLLU267[order(CGPLLU267$Timepoint),]
CGPLLU269 <- CGPLLU269[order(CGPLLU269$Timepoint),]
CGPLLU271 <- CGPLLU271[order(CGPLLU271$Timepoint),]
CGPLLU43 <- CGPLLU43[order(CGPLLU43$Timepoint),]
CGPLLU86 <- CGPLLU86[order(CGPLLU86$Timepoint),]
CGPLLU88 <- CGPLLU88[order(CGPLLU88$Timepoint),]
CGPLLU89 <- CGPLLU89[order(CGPLLU89$Timepoint),]
#In silico size selection effect
ichorCNA_unfiltered <- ichorCNA_unfiltered[,-1]
names(ichorCNA_unfiltered)[names(ichorCNA_unfiltered) == "tfx"] <- "sample"
names(ichorCNA_unfiltered)[names(ichorCNA_unfiltered) == "gender"] <- "tfx"
ichor_unfiltered <-  merge(x=ichorCNA_unfiltered,y=metadata,by="sample",all.x = TRUE)
ichor_unfiltered <- ichor_unfiltered[- grep("EE87920", ichor_unfiltered$sample),]
ichor_unfiltered <- ichor_unfiltered[- grep("EE87921", ichor_unfiltered$sample),]
unfiltered_group_ordered <- with(ichor_unfiltered,
                      reorder(group,
                              tfx,
                              median))

ichorCNA_filtered <- ichorCNA_filtered[,-1]
names(ichorCNA_filtered)[names(ichorCNA_filtered) == "tfx"] <- "sample"
names(ichorCNA_filtered)[names(ichorCNA_filtered) == "gender"] <- "tfx"
ichorCNA_filtered$sample <- gsub("_filtered.sortByCoord","",as.character(ichorCNA_filtered$sample))
ichor_filtered <-  merge(x=ichorCNA_filtered,y=metadata,by="sample",all.x = TRUE)
ichor_filtered <- ichor_filtered[- grep("EE87920", ichor_filtered$sample),]
ichor_filtered <- ichor_filtered[- grep("EE87921", ichor_filtered$sample),]
filtered_group_ordered <- with(ichor_filtered,
                                 reorder(group,
                                         tfx,
                                         median))
###ichorCNA plots
#Tumor fraction overview of all cohorts
my_comparisons <- list( c("Healthy", "Breast"), c("Healthy", "Bile duct"), c("Healthy", "Pancreas"),c("Healthy", "Ovarian"), c("Healthy", "Lung"),c("Healthy", "Colorectal"), c("Healthy", "Stomach") )

ichor_p1 <- ggplot(ichordf1, aes(x=group_ordered, y=tfx,color=group))+ theme_classic() + 
  labs(y= "Tumor fraction", x = "") +
  geom_boxplot()+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",label.y = 0.35,step.increase = 0.065,tip.length = 0.01)

ichor_p1 <- ichor_p1 + theme(legend.position = "none") + labs(tag="A")
#Comparison cancer stage and tumor fraction
ichor_p2 <- ggplot(ichordf2, aes(x=Stage, y=tfx,color=Stage))+ theme_classic() + 
  labs(y= "Tumor fraction", x = "Cancer Stage") + 
  geom_boxplot()

ichor_p2 <- ichor_p2  + theme(legend.position = "none") + labs(tag="B")
#Longitudinal data
q_tfx <- ggplot(CGPLLU13,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=20,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU13") +
  geom_point(data = CGPLLU13[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() + 
  ylim(0,0.2) 

w_tfx <- ggplot(CGPLLU14,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU14") +
  geom_point(data = CGPLLU14[4, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.05) 

e_tfx <- ggplot(CGPLLU244,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU244") +
  geom_point(data = CGPLLU244[2, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() + 
  ylim(0,0.2)

r_tfx <- ggplot(CGPLLU245,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU245") +
  geom_point(data = CGPLLU245[2, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() + 
  ylim(0,0.2)
                                                                                                  
t_tfx <- ggplot(CGPLLU246,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU246") +
  geom_point(data = CGPLLU246[2, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.15)

z_tfx <- ggplot(CGPLLU264,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU264") +
  geom_point(data = CGPLLU264[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.006)

u_tfx <- ggplot(CGPLLU265,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU265") + 
  geom_point(data = CGPLLU265[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.01)

i_tfx <- ggplot(CGPLLU266,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU266") +  
  geom_point(data = CGPLLU266[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.005)

o_tfx <- ggplot(CGPLLU267,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU267") +  
  geom_point(data = CGPLLU267[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.1)
  
p_tfx <- ggplot(CGPLLU269,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU269") +
  geom_point(data = CGPLLU269[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.02)

a_tfx <- ggplot(CGPLLU271,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU271") +  
  geom_point(data = CGPLLU271[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.07)

s_tfx <- ggplot(CGPLLU43,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU43") +
  geom_point(data = CGPLLU43[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.008)

d_tfx <- ggplot(CGPLLU86,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU86") +
  geom_point(data = CGPLLU86[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.015)
  
f_tfx <- ggplot(CGPLLU88,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU88") +
  geom_point(data = CGPLLU88[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.2)

g_tfx <- ggplot(CGPLLU89,aes(x = Timepoint, y= tfx)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = tfx), color = "blue") + ggtitle("CGPLLU89") +
  geom_point(data = CGPLLU89[1, ], aes(x = Timepoint, y= tfx), fill = "red", size = 3,shape=23) + 
  labs(y = "Tumor fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.015)

q_maf <- ggplot(CGPLLU13,aes(x = Timepoint, y= maf)) +
  geom_point(shape=20,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU13") +
  geom_point(data = CGPLLU13[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() + 
  ylim(0,0.2)#+ labs(tag="A")

w_maf <- ggplot(CGPLLU14,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU14") +
  geom_point(data = CGPLLU14[4, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.05)#+ labs(tag="B")

e_maf <- ggplot(CGPLLU244,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU244") +
  geom_point(data = CGPLLU244[2, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() + 
  ylim(0,0.2)#+ labs(tag="C")

r_maf <- ggplot(CGPLLU245,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU245") +
  geom_point(data = CGPLLU245[2, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() + 
  ylim(0,0.2)#+ labs(tag="D")

t_maf <- ggplot(CGPLLU246,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU246") +
  geom_point(data = CGPLLU246[2, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.15)#+ labs(tag="E")

z_maf <- ggplot(CGPLLU264,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU264") +
  geom_point(data = CGPLLU264[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.006)#+ labs(tag="F")

u_maf <- ggplot(CGPLLU265,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU265") + 
  geom_point(data = CGPLLU265[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.01)#+ labs(tag="G")

i_maf <- ggplot(CGPLLU266,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU266") +  
  geom_point(data = CGPLLU266[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.005)#+ labs(tag="H")

o_maf <- ggplot(CGPLLU267,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU267") +  
  geom_point(data = CGPLLU267[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.1)#+ labs(tag="I")

p_maf <- ggplot(CGPLLU269,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU269") +
  geom_point(data = CGPLLU269[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.02)#+ labs(tag="J")

a_maf <- ggplot(CGPLLU271,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU271") +  
  geom_point(data = CGPLLU271[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.07)#+ labs(tag="K")

s_maf <- ggplot(CGPLLU43,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU43") +
  geom_point(data = CGPLLU43[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.008)#+ labs(tag="L")

d_maf <- ggplot(CGPLLU86,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU86") +
  geom_point(data = CGPLLU86[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.015)#+ labs(tag="M")

f_maf <- ggplot(CGPLLU88,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU88") +
  geom_point(data = CGPLLU88[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.2)#+ labs(tag="N")

g_maf <- ggplot(CGPLLU89,aes(x = Timepoint, y= maf)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = maf), color = "blue") + ggtitle("CGPLLU89") +
  geom_point(data = CGPLLU89[1, ], aes(x = Timepoint, y= maf), fill = "red", size = 3,shape=23) + 
  labs(y = "Mutant allele fraction",x="Timepoint(day)") +
  theme_classic() +
  ylim(0,0.015)#+ labs(tag="O")

q_freia <- ggplot(CGPLLU13,aes(x = Timepoint, y= freia)) +
  geom_point(shape=20,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU13") +
  geom_point(data = CGPLLU13[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

w_freia <- ggplot(CGPLLU14,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU14") +
  geom_point(data = CGPLLU14[4, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

e_freia <- ggplot(CGPLLU244,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU244") +
  geom_point(data = CGPLLU244[2, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

r_freia <- ggplot(CGPLLU245,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU245") +
  geom_point(data = CGPLLU245[2, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

t_freia <- ggplot(CGPLLU246,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU246") +
  geom_point(data = CGPLLU246[2, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

z_freia <- ggplot(CGPLLU264,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU264") +
  geom_point(data = CGPLLU264[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

u_freia <- ggplot(CGPLLU265,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU265") + 
  geom_point(data = CGPLLU265[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

i_freia <- ggplot(CGPLLU266,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU266") +  
  geom_point(data = CGPLLU266[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

o_freia <- ggplot(CGPLLU267,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU267") +  
  geom_point(data = CGPLLU267[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

p_freia <- ggplot(CGPLLU269,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU269") +
  geom_point(data = CGPLLU269[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

a_freia <- ggplot(CGPLLU271,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU271") +  
  geom_point(data = CGPLLU271[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

s_freia <- ggplot(CGPLLU43,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU43") +
  geom_point(data = CGPLLU43[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

d_freia <- ggplot(CGPLLU86,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU86") +
  geom_point(data = CGPLLU86[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic() 

f_freia <- ggplot(CGPLLU88,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU88") +
  geom_point(data = CGPLLU88[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

g_freia <- ggplot(CGPLLU89,aes(x = Timepoint, y= freia)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = freia), color = "blue") + ggtitle("CGPLLU89") +
  geom_point(data = CGPLLU89[1, ], aes(x = Timepoint, y= freia), fill = "red", size = 3,shape=23) + 
  labs(y = "FrEIA score",x="Timepoint(day)") +
  theme_classic()

q_shortratio <- ggplot(CGPLLU13,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=20,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU13") +
  geom_point(data = CGPLLU13[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

w_shortratio <- ggplot(CGPLLU14,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU14") +
  geom_point(data = CGPLLU14[4, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

e_shortratio <- ggplot(CGPLLU244,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU244") +
  geom_point(data = CGPLLU244[2, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

r_shortratio <- ggplot(CGPLLU245,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU245") +
  geom_point(data = CGPLLU245[2, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

t_shortratio <- ggplot(CGPLLU246,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU246") +
  geom_point(data = CGPLLU246[2, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

z_shortratio <- ggplot(CGPLLU264,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU264") +
  geom_point(data = CGPLLU264[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

u_shortratio <- ggplot(CGPLLU265,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU265") + 
  geom_point(data = CGPLLU265[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

i_shortratio <- ggplot(CGPLLU266,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU266") +  
  geom_point(data = CGPLLU266[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

o_shortratio <- ggplot(CGPLLU267,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU267") +  
  geom_point(data = CGPLLU267[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

p_shortratio <- ggplot(CGPLLU269,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU269") +
  geom_point(data = CGPLLU269[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

a_shortratio <- ggplot(CGPLLU271,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU271") +  
  geom_point(data = CGPLLU271[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

s_shortratio <- ggplot(CGPLLU43,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU43") +
  geom_point(data = CGPLLU43[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

d_shortratio <- ggplot(CGPLLU86,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU86") +
  geom_point(data = CGPLLU86[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic() 

f_shortratio <- ggplot(CGPLLU88,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU88") +
  geom_point(data = CGPLLU88[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

g_shortratio <- ggplot(CGPLLU89,aes(x = Timepoint, y= shortratio)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = shortratio), color = "blue") + ggtitle("CGPLLU89") +
  geom_point(data = CGPLLU89[1, ], aes(x = Timepoint, y= shortratio), fill = "red", size = 3,shape=23) + 
  labs(y = "Short fragment ratio",x="Timepoint(day)") +
  theme_classic()

q_dipsaec <- ggplot(CGPLLU13,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=20,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU13") +
  geom_point(data = CGPLLU13[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

w_dipsaec <- ggplot(CGPLLU14,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU14") +
  geom_point(data = CGPLLU14[4, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

e_dipsaec <- ggplot(CGPLLU244,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU244") +
  geom_point(data = CGPLLU244[2, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

r_dipsaec <- ggplot(CGPLLU245,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU245") +
  geom_point(data = CGPLLU245[2, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

t_dipsaec <- ggplot(CGPLLU246,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU246") +
  geom_point(data = CGPLLU246[2, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

z_dipsaec <- ggplot(CGPLLU264,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU264") +
  geom_point(data = CGPLLU264[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

u_dipsaec <- ggplot(CGPLLU265,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU265") + 
  geom_point(data = CGPLLU265[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

i_dipsaec <- ggplot(CGPLLU266,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU266") +  
  geom_point(data = CGPLLU266[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

o_dipsaec <- ggplot(CGPLLU267,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU267") +  
  geom_point(data = CGPLLU267[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

p_dipsaec <- ggplot(CGPLLU269,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) +
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU269") +
  geom_point(data = CGPLLU269[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

a_dipsaec <- ggplot(CGPLLU271,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU271") +  
  geom_point(data = CGPLLU271[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

s_dipsaec <- ggplot(CGPLLU43,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU43") +
  geom_point(data = CGPLLU43[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

d_dipsaec <- ggplot(CGPLLU86,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU86") +
  geom_point(data = CGPLLU86[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic() 

f_dipsaec <- ggplot(CGPLLU88,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU88") +
  geom_point(data = CGPLLU88[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

g_dipsaec <- ggplot(CGPLLU89,aes(x = Timepoint, y= dipsaec)) +
  geom_point(shape=18,color = "blue",fill="blue",size=3) + 
  geom_line(aes(y = dipsaec), color = "blue") + ggtitle("CGPLLU89") +
  geom_point(data = CGPLLU89[1, ], aes(x = Timepoint, y= dipsaec), fill = "red", size = 3,shape=23) + 
  labs(y = "Dip depth SAEC",x="Timepoint(day)") +
  theme_classic()

time_all <-ggarrange(q_maf,q_tfx,q_freia,q_shortratio,q_dipsaec,
                     w_maf,w_tfx,w_freia,w_shortratio,w_dipsaec,
                     e_maf,e_tfx,e_freia,e_shortratio,e_dipsaec,
                     r_maf,r_tfx,r_freia,r_shortratio,r_dipsaec,
                     t_maf,t_tfx,t_freia,t_shortratio,t_dipsaec,
                     z_maf,z_tfx,z_freia,z_shortratio,z_dipsaec,
                     u_maf,u_tfx,u_freia,u_shortratio,u_dipsaec,
                     i_maf,i_tfx,i_freia,i_shortratio,i_dipsaec,
                     o_maf,o_tfx,o_freia,o_shortratio,o_dipsaec,
                     p_maf,p_tfx,p_freia,p_shortratio,p_dipsaec,
                     a_maf,a_tfx,a_freia,a_shortratio,a_dipsaec,
                     s_maf,s_tfx,s_freia,s_shortratio,s_dipsaec,
                     d_maf,d_tfx,d_freia,d_shortratio,d_dipsaec,
                     f_maf,f_tfx,f_freia,f_shortratio,f_dipsaec,
                     g_maf,g_tfx,g_freia,g_shortratio,g_dipsaec,
                     ncol = 5, nrow = 15)


annotate_figure(time_qwert, top = text_grob("Longitudinal lung cancer data", 
                                      color = "black", size = 14))
 
#in silico size selection
ichor_p4 <- ggplot(ichor_unfiltered, aes(x=unfiltered_group_ordered, y=tfx,color=group))+ theme_classic() + 
  labs(y= "Tumor fraction", x = "Unfiltered") +
  geom_boxplot() + theme(legend.position = "none")+ labs(tag="C") + ylim(0,0.4)

ichor_p5 <- ggplot(ichor_filtered, aes(x=filtered_group_ordered, y=tfx,color=group))+ theme_classic() + 
  labs(y= "", x = "Filtered") + labs(tag="") +
  geom_boxplot() + theme(legend.position = "none")+ ylim(0,0.4)

ichor_p6 <- ggarrange(ichor_p4,ichor_p5,ncol = 2, nrow = 1) 

###LIQUORICE###

###LIQUORICE Data wrangling
#liq_unfiltered
liq_unfiltered <- select(liq_unfiltered,"sample","region-set","Dip depth: z-score vs controls in same region set")
names(liq_unfiltered)[names(liq_unfiltered) == "Dip depth: z-score vs controls in same region set"] <- "Dip_depth"
names(liq_unfiltered)[names(liq_unfiltered) == "region-set"] <- "region"
liq_unfiltered <- liq_unfiltered %>% 
  pivot_wider(names_from = region, values_from = c(Dip_depth))

liq_unfiltered$sample <- gsub(".sortByCoord","",as.character(liq_unfiltered$sample))
liq_unfiltered <-  merge(x=liq_unfiltered,y=metadata,by="sample",all.x = TRUE)

liq_unfiltered_hl <- liq_unfiltered[liq_unfiltered$group =="Healthy" | liq_unfiltered$group =="Lung",]
#liq_150
liq_150 <- select(liq_150,"sample","region-set","Dip depth: z-score vs controls in same region set")
names(liq_150)[names(liq_150) == "Dip depth: z-score vs controls in same region set"] <- "Dip_depth"
names(liq_150)[names(liq_150) == "region-set"] <- "region"
liq_150 <- liq_150 %>% 
  pivot_wider(names_from = region, values_from = c(Dip_depth))

liq_150$sample <- gsub(".sortByCoord","",as.character(liq_150$sample))
liq_150 <-  merge(x=liq_150,y=metadata,by="sample",all.x = TRUE)

liq_150_hl <- liq_150[liq_150$group =="Healthy" | liq_150$group =="Lung",]
#liq_top50
liq_top50 <- select(liq_top50,"sample","region-set","Dip depth: z-score vs controls in same region set")
names(liq_top50)[names(liq_top50) == "Dip depth: z-score vs controls in same region set"] <- "Dip_depth"
names(liq_top50)[names(liq_top50) == "region-set"] <- "region"
liq_top50 <- liq_top50 %>% 
  pivot_wider(names_from = region, values_from = c(Dip_depth))

liq_top50$sample <- gsub(".sortByCoord","",as.character(liq_top50$sample))
liq_top50 <-  merge(x=liq_top50,y=metadata,by="sample",all.x = TRUE)

liq_top50_hl <- liq_top50[liq_top50$group =="Healthy" | liq_top50$group =="Lung",]
###LIQUORICE plots
#show different regions
liq_p1 <- ggplot(liq_150_hl, aes(x=group, y=hematopoietic_all_hg38,color=group))+ theme_classic() + 
  labs(y= "Coverage dip depth", x = "hematopoietic_all_hg38") + theme(legend.position = "none") +
  geom_boxplot() + stat_compare_means(method = "t.test",label.x = 1.3)+ labs(tag="B") +
  scale_color_manual(values=c("blue", "red"))
  

liq_p2 <- ggplot(liq_150_hl, aes(x=group, y=SAEC_hg38,color=group))+ theme_classic() + 
  labs(y= "", x = "SAEC_hg38") + theme(legend.position = "none")  +
  geom_boxplot() + stat_compare_means(method = "t.test",label.x = 1.3)+ labs(tag="") +
  scale_color_manual(values=c("blue", "red"))

liq_p3 <- ggplot(liq_150_hl, aes(x=group, y=A549_hg38,color=group))+ theme_classic() + 
  labs(y= "", x = "A549_hg38") + theme(legend.position = "none") +
  geom_boxplot() + stat_compare_means(method = "t.test",label.x = 1.3) + labs(tag="") +
  scale_color_manual(values=c("blue", "red"))

liq_p7 <- ggarrange(liq_p1,liq_p2,liq_p3,nrow=1,ncol=3) 
#show different filtering
liq_p4 <- ggplot(liq_unfiltered_hl, aes(x=group, y=SAEC_hg38,color=group))+ theme_classic() + 
  labs(y= "Coverage dip depth SAEC_hg38", x = "Unfiltered") + theme(legend.position = "none") +
  geom_boxplot() + stat_compare_means(method = "t.test",label.x = 1.3) + ylim(-4,6)+ labs(tag="D") +
  scale_color_manual(values=c("blue", "red"))

liq_p5 <- ggplot(liq_150_hl, aes(x=group, y=SAEC_hg38,color=group))+ theme_classic() + 
  labs(y= "", x = "Above 150 Quality Score") + theme(legend.position = "none") +
  geom_boxplot() + stat_compare_means(method = "t.test",label.x = 1.3) + ylim(-4,6)+ labs(tag="") +
  scale_color_manual(values=c("blue", "red"))

liq_p6 <- ggplot(liq_top50_hl, aes(x=group, y=SAEC_hg38,color=group))+ theme_classic() + 
  labs(y= "", x = "Top 50 Regions") + theme(legend.position = "none") +
  geom_boxplot() + stat_compare_means(method = "t.test",label.x = 1.3) + ylim(-4,6)+ labs(tag="") +
  scale_color_manual(values=c("blue", "red"))

liq_p8 <- ggarrange(liq_p4,liq_p5,liq_p6,nrow=1,ncol=3) 
#quality scores saec
saec_qual <- ggplot(SAEC_hg38, aes(x=X4)) +
  geom_histogram(color="darkgreen", fill="lightgreen") +
  labs(y= "Count", x = "Qualitiy scores of SAEC region-set") +
  theme_classic() + labs(tag="C")

###FrEIA###

###FrEIA Data wrangling
#Overview of cohorts
names(FrEIA_mouliere_tnc_and_panel)[names(FrEIA_mouliere_tnc_and_panel) == "sample_name"] <- "sample"
names(FrEIA_my_panel_and_mouliere_tnc)[names(FrEIA_my_panel_and_mouliere_tnc) == "sample_name"] <- "sample"
names(FrEIA_my_tnc_mouliere_panel)[names(FrEIA_my_tnc_mouliere_panel) == "sample_name"] <- "sample"
names(FrEIA_my_panel_and_tnc)[names(FrEIA_my_panel_and_tnc) == "sample_name"] <- "sample"

FrEIA_mouliere_tnc_and_panel$sample<-gsub(".sortByCoord","",as.character(FrEIA_mouliere_tnc_and_panel$sample))
FrEIA_my_panel_and_mouliere_tnc$sample<-gsub(".sortByCoord","",as.character(FrEIA_my_panel_and_mouliere_tnc$sample))
FrEIA_my_tnc_mouliere_panel$sample<-gsub(".sortByCoord","",as.character(FrEIA_my_tnc_mouliere_panel$sample))
FrEIA_my_panel_and_tnc$sample<-gsub(".sortByCoord","",as.character(FrEIA_my_panel_and_tnc$sample))

FrEIA_mouliere_tnc_and_panel <- merge(x=FrEIA_mouliere_tnc_and_panel,y=metadata,by="sample",all.x = TRUE)
FrEIA_my_panel_and_mouliere_tnc <- merge(x=FrEIA_my_panel_and_mouliere_tnc,y=metadata,by="sample",all.x = TRUE)
FrEIA_my_tnc_mouliere_panel <- merge(x=FrEIA_my_tnc_mouliere_panel,y=metadata,by="sample",all.x = TRUE)
FrEIA_my_panel_and_tnc <- merge(x=FrEIA_my_panel_and_tnc,y=metadata,by="sample",all.x = TRUE)

FrEIA_mouliere_tnc_and_panel_group_ordered <- with(FrEIA_mouliere_tnc_and_panel,reorder(group,FrEIA_score,median))
FrEIA_my_panel_and_mouliere_tnc_group_ordered <- with(FrEIA_my_panel_and_mouliere_tnc,reorder(group,FrEIA_score,median))
FrEIA_my_tnc_mouliere_panel_group_ordered <- with(FrEIA_my_tnc_mouliere_panel,reorder(group,FrEIA_score,median))
FrEIA_my_panel_and_tnc_ordered <- with(FrEIA_my_panel_and_tnc,reorder(group,FrEIA_score,median))
#Correlation FrEIA and ichorCNA
ichor_corr <- ichordf %>% select('sample','tfx')
freia_corr <- FrEIA_my_panel_and_tnc %>% select('sample','FrEIA_score')
ichor_freia_corr <-  merge(x=ichor_corr,y=freia_corr,by="sample",all.x = TRUE)
ichor_freia_corr <-  merge(x=ichor_freia_corr,y=metadata,by="sample",all.x = TRUE)

my_comparisons_1 <- list( c("Healthy", "Colorectal"),c("Healthy", "Pancreas"),c("Healthy", "Lung"),c("Healthy", "Ovarian"),c("Healthy", "Stomach"),c("Healthy", "Breast"), c("Healthy", "Bile duct") )
my_comparisons_2 <- list( c("Healthy", "Stomach"),c("Healthy", "Lung"),c("Healthy", "Breast"),c("Healthy", "Pancreas"),c("Healthy", "Ovarian"),c("Healthy", "Colorectal"), c("Healthy", "Bile duct") )
my_comparisons_3 <- list( c("Healthy", "Pancreas"),c("Healthy", "Colorectal"),c("Healthy", "Breast"),c("Healthy", "Ovarian"), c("Healthy", "Bile duct"),c("Healthy", "Stomach"),c("Healthy", "Lung") )
my_comparisons_4 <- list( c("Healthy", "Pancreas"),c("Healthy", "Colorectal"),c("Healthy", "Breast"),c("Healthy", "Lung"),c("Healthy", "Ovarian"),c("Healthy", "Stomach"), c("Healthy", "Bile duct") )

###FrEIA plots
#Overview of cohorts
FrEIA_p1 <- ggplot(FrEIA_mouliere_tnc_and_panel, aes(x=FrEIA_mouliere_tnc_and_panel_group_ordered, y=FrEIA_score,color=group))+ theme_classic() + 
  labs(y= "FrEIA_score", x = "") + theme(legend.position = "none") + 
  #ggtitle("Original trinucleotide selection and original panel of median") +
  geom_violin() + geom_boxplot(width=0.1) + ylim(0,11.5) + stat_compare_means(comparisons = my_comparisons_1,label = "p.signif",label.y = 6,step.increase = 0.1,tip.length = 0.01) + 
  labs(tag="A")

FrEIA_p2 <- ggplot(FrEIA_my_panel_and_mouliere_tnc, aes(x=FrEIA_my_panel_and_mouliere_tnc_group_ordered, y=FrEIA_score,color=group))+ theme_classic() + 
  labs(y= "FrEIA_score", x = "") + theme(legend.position = "none") + 
  #ggtitle("Adjusted panel of median and original trinucleotide selection") +
  geom_violin() + geom_boxplot(width=0.1) + ylim(0,11.5) + stat_compare_means(comparisons = my_comparisons_2,label = "p.signif",label.y = 6,step.increase = 0.08,tip.length = 0.01)+
  labs(tag="B")
FrEIA_p3 <- ggplot(FrEIA_my_tnc_mouliere_panel, aes(x=FrEIA_my_tnc_mouliere_panel_group_ordered, y=FrEIA_score,color=group))+ theme_classic() + 
  labs(y= "FrEIA_score", x = "") + theme(legend.position = "none") + 
  #ggtitle("Adjusted trinucleotide selection and original panel of median") +
  geom_violin() + geom_boxplot(width=0.1) + ylim(0,11.5) + stat_compare_means(comparisons = my_comparisons_3,label = "p.signif",label.y = 6,step.increase = 0.1,tip.length = 0.01)+
  labs(tag="C")
FrEIA_p4 <- ggplot(FrEIA_my_panel_and_tnc, aes(x=FrEIA_my_panel_and_tnc_ordered, y=FrEIA_score,color=group))+ theme_classic() + 
  labs(y= "FrEIA_score", x = "") + theme(legend.position = "none") + 
  #ggtitle("Adjusted trinucleotide selection and adjusted panel of median") +
  geom_violin() + geom_boxplot(width=0.1) + ylim(0,11.5) + stat_compare_means(comparisons = my_comparisons_4,label = "p.signif",label.y = 6,step.increase = 0.1,tip.length = 0.01)+
  labs(tag="D")
FrEIA_p5 <- ggarrange(FrEIA_p1,FrEIA_p2,FrEIA_p3,FrEIA_p4,ncol = 2, nrow = 2)

FrEIA_mouliere_tnc_and_panel.aov <- aov(FrEIA_score ~ group, data = FrEIA_mouliere_tnc_and_panel)
FrEIA_my_panel_and_mouliere_tnc.aov <- aov(FrEIA_score ~ group, data = FrEIA_my_panel_and_mouliere_tnc)
FrEIA_my_tnc_mouliere_panel.aov <- aov(FrEIA_score ~ group, data = FrEIA_my_tnc_mouliere_panel)
FrEIA_my_panel_and_tnc.aov <- aov(FrEIA_score ~ group, data = FrEIA_my_panel_and_tnc)

summary(FrEIA_mouliere_tnc_and_panel.aov)
summary(FrEIA_my_panel_and_mouliere_tnc.aov)
summary(FrEIA_my_tnc_mouliere_panel.aov)
summary(FrEIA_my_panel_and_tnc.aov)

TukeyHSD(FrEIA_mouliere_tnc_and_panel.aov)
TukeyHSD(FrEIA_my_panel_and_mouliere_tnc.aov)
TukeyHSD(FrEIA_my_tnc_mouliere_panel.aov)
TukeyHSD(FrEIA_my_panel_and_tnc.aov)
#ichorCNA and FrEIA correlation
FrEIA_p6 <- ggplot(ichor_freia_corr, aes(x=tfx, y=FrEIA_score,color=group))+ theme_classic() + 
  labs(y= "FrEIA_score", x = "Tumor fraction") + 
  #ggtitle("Correlation FrEIA score and Tumor fraction") +
  geom_point() + labs(tag="A")

shapiro.test(ichor_freia_corr$tfx)
ggqqplot(ichor_freia_corr$FrEIA_score)
cor.test(ichor_freia_corr$tfx,ichor_freia_corr$FrEIA_score,method = "spearman")

#FrEIA and cancer stage
freia_stage <- merge(x=FrEIA_my_panel_and_tnc,y=cristiano_id_sample_metadata,by="sample",all.x = TRUE)
freia_stage <- merge(x=freia_stage,y=cristiano_sup_table,by="Patient",all.x = TRUE)
freia_stage<- subset(freia_stage, Stage %in% c("I", "II", "III", "IV"))
freia_stage_p <- ggplot(freia_stage, aes(x=Stage, y=FrEIA_score,color=Stage))+ theme_classic() + 
  labs(y= "FrEIA score", x = "Cancer Stage") + 
  #ggtitle("Correlation between FrEIA score and cancer stage") +
  geom_boxplot() + labs(tag="B")

FrEIA_p7 <- freia_stage_p  + theme(legend.position = "none")
#freia volcano
freia_volcano <- Cristiano_logFC_p_my_panel_and_tnc
#freia_volcano <- freia_volcano %>% remove_rownames %>% column_to_rownames(var="base")
names(freia_volcano)[names(freia_volcano) == "p-val"] <- "pval"
# add a column of NAs
freia_volcano$diffexpressed <- "No change"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
freia_volcano$diffexpressed[freia_volcano$Log10FC > 0.0132836968 & freia_volcano$pval < 0.01] <- "Increased"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
freia_volcano$diffexpressed[freia_volcano$Log10FC < -0.0073361298 & freia_volcano$pval < 0.01] <- "Decreased"
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
freia_volcano$delabel <- NA
freia_volcano$delabel[freia_volcano$diffexpressed != "No change"] <- freia_volcano$base[freia_volcano$diffexpressed != "No change"]

freia_volcano_p <- ggplot(data=freia_volcano, aes(x=Log10FC, y=-log10(pval), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "red", "black")) +
  geom_vline(xintercept=c(-0.0073361298, 0.0132836968), col="black",linetype="dotted") +
  geom_hline(yintercept=-log10(0.01), col="black",linetype="dotted") + labs(tag="C",col = "Proportions in cancer")

###Heatmap###

###Heatmap Data wrangling
#ratio
ratio <- read_delim("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/Heatmap/5mb_binned_SLratio.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
ratio$binid <- paste0(ratio$c1,"_", ratio$start, "_", ratio$end)

ratio <- subset(ratio, select=-c(c1,start, end))
ratio[is.na(ratio)] = 0
ratio <- ratio %>% remove_rownames %>% column_to_rownames(var="binid")
ratio <- as.data.frame(t(ratio))
ratio <- tibble::rownames_to_column(ratio, "sample")
ratio <-  merge(x=ratio,y=metadata,by="sample",all.x = TRUE)
healthy_ratio <- ratio[ratio$group =="Healthy",]
healthy_ratio <- healthy_ratio[-c(1,2),]
ref_ratio <- healthy_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
ref_ratio <- ref_ratio[,-566]
ref_ratio <- colMeans(ref_ratio)

healthy_ratio <- healthy_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
healthy_ratio <- healthy_ratio[,-566]
t_healthy_ratio <- t(healthy_ratio)
n_healthy_ratio <- t_healthy_ratio/ref_ratio
healthy_ratio <- scale(n_healthy_ratio)
m_healthy_ratio <- melt(healthy_ratio)
names(m_healthy_ratio)[names(m_healthy_ratio) == "Var1"] <- "Genome"
names(m_healthy_ratio)[names(m_healthy_ratio) == "Var2"] <- "Sample"

lung_ratio <- ratio[ratio$group =="Lung",]
lung_ratio <- lung_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
lung_ratio <- lung_ratio[,-566]
t_lung_ratio <- as.data.frame(t(lung_ratio))
n_lung_ratio <- as.matrix(t_lung_ratio/ref_ratio)
lung_ratio <- scale(n_lung_ratio)
m_lung_ratio <- melt(lung_ratio)
names(m_lung_ratio)[names(m_lung_ratio) == "Var1"] <- "Genome"
names(m_lung_ratio)[names(m_lung_ratio) == "Var2"] <- "Sample"

breast_ratio <- ratio[ratio$group =="Breast",]
breast_ratio <- breast_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
breast_ratio <- breast_ratio[,-566]
t_breast_ratio <- as.data.frame(t(breast_ratio))
n_breast_ratio <- as.matrix(t_breast_ratio/ref_ratio)
breast_ratio <- scale(n_breast_ratio)
m_breast_ratio <- melt(breast_ratio)
names(m_breast_ratio)[names(m_breast_ratio) == "Var1"] <- "Genome"
names(m_breast_ratio)[names(m_breast_ratio) == "Var2"] <- "Sample"

colorectal_ratio <- ratio[ratio$group =="Colorectal",]
colorectal_ratio <- colorectal_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
colorectal_ratio <- colorectal_ratio[,-566]
t_colorectal_ratio <- as.data.frame(t(colorectal_ratio))
n_colorectal_ratio <- as.matrix(t_colorectal_ratio/ref_ratio)
colorectal_ratio <- scale(n_colorectal_ratio)
m_colorectal_ratio <- melt(colorectal_ratio)
names(m_colorectal_ratio)[names(m_colorectal_ratio) == "Var1"] <- "Genome"
names(m_colorectal_ratio)[names(m_colorectal_ratio) == "Var2"] <- "Sample"

stomach_ratio <- ratio[ratio$group =="Stomach",]
stomach_ratio <- stomach_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
stomach_ratio <- stomach_ratio[,-566]
t_stomach_ratio <- as.data.frame(t(stomach_ratio))
n_stomach_ratio <- as.matrix(t_stomach_ratio/ref_ratio)
stomach_ratio <- scale(n_stomach_ratio)
m_stomach_ratio <- melt(stomach_ratio)
names(m_stomach_ratio)[names(m_stomach_ratio) == "Var1"] <- "Genome"
names(m_stomach_ratio)[names(m_stomach_ratio) == "Var2"] <- "Sample"

bile_duct_ratio <- ratio[ratio$group =="Bile duct",]
bile_duct_ratio <- bile_duct_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
bile_duct_ratio <- bile_duct_ratio[,-566]
t_bile_duct_ratio <- as.data.frame(t(bile_duct_ratio))
n_bile_duct_ratio <- as.matrix(t_bile_duct_ratio/ref_ratio)
bile_duct_ratio <- scale(n_bile_duct_ratio)
m_bile_duct_ratio <- melt(bile_duct_ratio)
names(m_bile_duct_ratio)[names(m_bile_duct_ratio) == "Var1"] <- "Genome"
names(m_bile_duct_ratio)[names(m_bile_duct_ratio) == "Var2"] <- "Sample"

pancreas_ratio <- ratio[ratio$group =="Pancreas",]
pancreas_ratio <- pancreas_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
pancreas_ratio <- pancreas_ratio[,-566]
t_pancreas_ratio <- as.data.frame(t(pancreas_ratio))
n_pancreas_ratio <- as.matrix(t_pancreas_ratio/ref_ratio)
pancreas_ratio <- scale(n_pancreas_ratio)
m_pancreas_ratio <- melt(pancreas_ratio)
names(m_pancreas_ratio)[names(m_pancreas_ratio) == "Var1"] <- "Genome"
names(m_pancreas_ratio)[names(m_pancreas_ratio) == "Var2"] <- "Sample"

ovarian_ratio <- ratio[ratio$group =="Ovarian",]
ovarian_ratio <- ovarian_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
ovarian_ratio <- ovarian_ratio[,-566]
t_ovarian_ratio <- as.data.frame(t(ovarian_ratio))
n_ovarian_ratio <- as.matrix(t_ovarian_ratio/ref_ratio)
ovarian_ratio <- scale(n_ovarian_ratio)
m_ovarian_ratio <- melt(ovarian_ratio)
names(m_ovarian_ratio)[names(m_ovarian_ratio) == "Var1"] <- "Genome"
names(m_ovarian_ratio)[names(m_ovarian_ratio) == "Var2"] <- "Sample"

duodenum_ratio <- ratio[ratio$group =="Duodenum",]
duodenum_ratio <- duodenum_ratio %>% remove_rownames %>% column_to_rownames(var="sample")
duodenum_ratio <- duodenum_ratio[,-566]
t_duodenum_ratio <- as.data.frame(t(duodenum_ratio))
n_duodenum_ratio <- as.matrix(t_duodenum_ratio/ref_ratio)
duodenum_ratio <- scale(n_duodenum_ratio)
m_duodenum_ratio <- melt(duodenum_ratio)
names(m_duodenum_ratio)[names(m_duodenum_ratio) == "Var1"] <- "Genome"
names(m_duodenum_ratio)[names(m_duodenum_ratio) == "Var2"] <- "Sample"
#coverage
cov <- read_delim("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/Heatmap/5mb_binned_coverage.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
cov$binid <- paste0(cov$c1,"_", cov$start, "_", cov$end)

cov <- subset(cov, select=-c(c1,start, end))
cov[is.na(cov)] = 0
cov <- cov %>% remove_rownames %>% column_to_rownames(var="binid")
cov <- as.data.frame(t(cov))
cov <- tibble::rownames_to_column(cov, "sample")
cov <-  merge(x=cov,y=metadata,by="sample",all.x = TRUE)
healthy_cov <- cov[cov$group =="Healthy",]
ref_cov <- healthy_cov %>% remove_rownames %>% column_to_rownames(var="sample")
ref_cov <- ref_cov[,-566]
ref_cov <- colMeans(ref_cov)

healthy_cov <- healthy_cov %>% remove_rownames %>% column_to_rownames(var="sample")
healthy_cov <- healthy_cov[,-566]
#t_healthy_cov <- t(healthy_cov)
#n_healthy_cov <- t_healthy_cov/ref_cov

lung_cov <- cov[cov$group =="Lung",]
lung_cov <- lung_cov %>% remove_rownames %>% column_to_rownames(var="sample")
lung_cov <- lung_cov[,-566]
#t_lung_cov <- as.data.frame(t(lung_cov))
#n_lung_cov <- as.matrix(t_lung_cov/ref_cov)

breast_cov <- cov[cov$group =="Breast",]
breast_cov <- breast_cov %>% remove_rownames %>% column_to_rownames(var="sample")
breast_cov <- breast_cov[,-566]
#t_breast_cov <- as.data.frame(t(breast_cov))
#n_breast_cov <- as.matrix(t_breast_cov/ref_cov)

colorectal_cov <- cov[cov$group =="Colorectal",]
colorectal_cov <- colorectal_cov %>% remove_rownames %>% column_to_rownames(var="sample")
colorectal_cov <- colorectal_cov[,-566]
#t_colorectal_cov <- as.data.frame(t(colorectal_cov))
#n_colorectal_cov <- as.matrix(t_colorectal_cov/ref_cov)

stomach_cov <- cov[cov$group =="Stomach",]
stomach_cov <- stomach_cov %>% remove_rownames %>% column_to_rownames(var="sample")
stomach_cov <- stomach_cov[,-566]
#t_stomach_cov <- as.data.frame(t(stomach_cov))
#n_stomach_cov <- as.matrix(t_stomach_cov/ref_cov)

bile_duct_cov <- cov[cov$group =="Bile duct",]
bile_duct_cov <- bile_duct_cov %>% remove_rownames %>% column_to_rownames(var="sample")
bile_duct_cov <- bile_duct_cov[,-566]
#t_bile_duct_cov <- as.data.frame(t(bile_duct_cov))
#n_bile_duct_cov <- as.matrix(t_bile_duct_cov/ref_cov)

pancreas_cov <- cov[cov$group =="Pancreas",]
pancreas_cov <- pancreas_cov %>% remove_rownames %>% column_to_rownames(var="sample")
pancreas_cov <- pancreas_cov[,-566]
#t_pancreas_cov <- as.data.frame(t(pancreas_cov))
#n_pancreas_cov <- as.matrix(t_pancreas_cov/ref_cov)

ovarian_cov <- cov[cov$group =="Ovarian",]
ovarian_cov <- ovarian_cov %>% remove_rownames %>% column_to_rownames(var="sample")
ovarian_cov <- ovarian_cov[,-566]
#t_ovarian_cov <- as.data.frame(t(ovarian_cov))
#n_ovarian_cov <- as.matrix(t_ovarian_cov/ref_cov)

duodenum_cov <- cov[cov$group =="Duodenum",]
duodenum_cov <- duodenum_cov %>% remove_rownames %>% column_to_rownames(var="sample")
duodenum_cov <- duodenum_cov[,-566]
#t_duodenum_cov <- as.data.frame(t(duodenum_cov))
#n_duodenum_cov <- as.matrix(t_duodenum_cov/ref_cov)
###Heatmap plots 
#coverage
cov_healthy <- heatmap(n_healthy_cov,Colv = NA, Rowv = NA,labRow="", main="Healthy")
cov_lung <-heatmap(n_lung_cov,Colv = NA, Rowv = NA,labRow="", main="Lung")
cov_breast <-heatmap(n_breast_cov,Colv = NA, Rowv = NA,labRow="", main="Breast")

cov_colorectal <-heatmap(n_colorectal_cov,Colv = NA, Rowv = NA,labRow="", main="Colorectal")
cov_stomach <-heatmap(n_stomach_cov,Colv = NA, Rowv = NA,labRow="", main="Stomach")
cov_bile_duct <-heatmap(n_bile_duct_cov,Colv = NA, Rowv = NA,labRow="", main="Bile duct")

cov_pancreas <-heatmap(n_pancreas_cov,Colv = NA, Rowv = NA,labRow="", main="Pancreas")
cov_ovarian <-heatmap(n_ovarian_cov,Colv = NA, Rowv = NA,labRow="", main="Ovarian")
cov_duodenum <-heatmap(n_duodenum_cov,Colv = NA, Rowv = NA,labRow="", main="Duodenum")
#ratio
hm_healthy <- ggplot(m_healthy_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) +  scale_fill_viridis(discrete=FALSE) + labs(title="Healthy",tag="C")
hm_lung <- ggplot(m_lung_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) + scale_fill_viridis(discrete=FALSE) + labs(title="Lung")
hm_breast <- ggplot(m_breast_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) + scale_fill_viridis(discrete=FALSE) + labs(title="Breast")

hm_colorectal <- ggplot(m_colorectal_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) +  scale_fill_viridis(discrete=FALSE) + labs(title="Colorectal")
hm_stomach <- ggplot(m_stomach_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) + scale_fill_viridis(discrete=FALSE) + labs(title="Stomach")
hm_bile_duct <- ggplot(m_bile_duct_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) + scale_fill_viridis(discrete=FALSE) + labs(title="Bile duct")

hm_pancreas <- ggplot(m_pancreas_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) +  scale_fill_viridis(discrete=FALSE) + labs(title="Pancreas",tag="D")
hm_ovarian <- ggplot(m_ovarian_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) + scale_fill_viridis(discrete=FALSE) + labs(title="Ovarian")
#hm_duodenum <- ggplot(m_duodenum_ratio, aes(Genome,Sample, fill= value)) + geom_tile()+ theme(axis.ticks = element_blank(),axis.text = element_blank()) + scale_fill_viridis(discrete=FALSE) + labs(title="Lung")
###Figures###
##Fragment length
frag_len_panel <- ggarrange(frag_len_p1,frag_len_p2,hm_healthy,hm_pancreas, ncol=2,nrow=2)
#frag_len_panel <- annotate_figure(frag_len_panel,fig.lab= "Figure X: cfDNA fragment lengths in healthy individuals and cancer patients")
##ichorCNA
#CNV
EE87865 <- image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/ichorCNA/EE87865.jpg")
EE87923 <- image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Data/ichorCNA/EE87923.jpg")

EE87865 <- image_ggplot(EE87865)
EE87923 <- image_ggplot(EE87923)
EE87923 <- EE87923 + ggtitle("A")
EE87865 <- EE87865 + ggtitle("B")
ichor_panel_1 <- ggarrange(EE87923,EE87865,ncol=1,nrow=2)
#ichor_panel_1 <- annotate_figure(ichor_panel_1,fig.lab= "Figure X: Copy number variants in a healthy individuals and a cancer patient")
#TFX
ichor_panel_2 <- ggarrange(ichor_p1,ggarrange(ichor_p2,ichor_p6,ncol=2,nrow=1),ncol=1,nrow=2)
#ichor_panel_2 <- annotate_figure(ichor_panel_2,fig.lab= "Figure X: Tumor fraction estimates in cfDNA fragments")
##Time series
time_panel <- annotate_figure(time_qwert,fig.lab= "Figure X: Time series analysis of lung cancer samples")
##Freia
freia_panel_1 <- annotate_figure(FrEIA_p5,fig.lab= "Figure X: Effects of parameter tuning on FrEIA score")
freia_panel_2 <- ggarrange(ggarrange(FrEIA_p6,FrEIA_p7,ncol=2,nrow=1),freia_volcano_p,ncol=1,nrow=2)
#freia_panel_2 <- annotate_figure(freia_panel_2,fig.lab= "Figure X: Nucleotide endmotif evaluation")
##LIQUORICE
SAEC <- image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/LIQUORICE/SAEC_overlay.jpg")
SAEC <- image_ggplot(SAEC)
SAEC <- SAEC + labs(tag="A")
liq_panel <- ggarrange(SAEC,liq_p7,saec_qual,liq_p8,ncol=2,nrow=2)
liq_panel <-  annotate_figure(liq_panel,fig.lab= "Figure X: Coverage loss in healthy individuals and cancer patients")
##ML
hvc_roc <- image_ggplot(image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/ML model/hvc_roc.svg")) + ggtitle("B")
hvc_fi <- image_ggplot(image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/ML model/feature_importance.svg")) + ggtitle("C")
hvc_table <- image_ggplot(image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/ML model/hvc.png"))
hvc_table <- hvc_table + labs(tag="A")
hvc_panel <- ggarrange(hvc_table,ggarrange(hvc_roc,hvc_fi,nrow=1,ncol=2),ncol=1,nrow=2)
#hvc_panel <-  annotate_figure(hvc_panel,fig.lab= "Figure X: Evaluation of cancer detection")

gnb_too_roc <- image_ggplot(image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/ML model/gnb_too_roc.svg")) + ggtitle("B")
gnb_too_cm <- image_ggplot(image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/ML model/gnb_too_cm.svg")) + ggtitle("C")
too_table <- image_ggplot(image_read("C:/Users/wadus/OneDrive/Desktop/Studium/Masterarbeit/Projects/Final_Plots/ML model/too.png"))
too_table <- too_table + labs(tag="A")
too_panel <- ggarrange(too_table,ggarrange(gnb_too_roc,gnb_too_cm,nrow=1,ncol=2),ncol=1,nrow=2)

#too_panel <-  annotate_figure(too_panel,fig.lab= "Figure X: Evaluation of cancer identification")

###ML model

#Preparing df
frag_df <- data.frame (p76_90 = apply(frag_count[(frag_count$length %in% seq(76, 90)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p91_120 = apply(frag_count[(frag_count$length %in% seq(91, 120)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p121_150 = apply(frag_count[(frag_count$length %in% seq(121, 150)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p151_180 = apply(frag_count[(frag_count$length %in% seq(151, 180)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p181_200 = apply(frag_count[(frag_count$length %in% seq(181, 200)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p201_250 = apply(frag_count[(frag_count$length %in% seq(201, 250)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p251_300 = apply(frag_count[(frag_count$length %in% seq(251, 300)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p301_350 = apply(frag_count[(frag_count$length %in% seq(301, 350)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p351_400 = apply(frag_count[(frag_count$length %in% seq(351, 400)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p401_500 = apply(frag_count[(frag_count$length %in% seq(401, 500)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p501_700 = apply(frag_count[(frag_count$length %in% seq(501, 700)), colnames(frag_count) !="length"], 2, sum, na.rm=T),
                       p701_999 = apply(frag_count[(frag_count$length %in% seq(701, 999)), colnames(frag_count) !="length"], 2, sum, na.rm=T))
frag_ratiodf <- t(apply(frag_df,1, function(x) x/sum(x)))

frag_df <- merge(x = frag_df, y = metadata, by.x = 0, by.y = "sample", all.x = TRUE)
names(frag_df)[names(frag_df) == 'Row.names'] <- 'sample'

frag_a <- apply(frag_ratiodf[rownames(frag_ratiodf) %in% 
                               metadata$sample[metadata$group=="Healthy"],], 
                2, median, na.rm=T)

frag_ratiodf <- merge(x=frag_ratiodf, y=metadata, by.x=0, by.y="sample", all.x=TRUE)
names(frag_ratiodf)[names(frag_ratiodf) == 'Row.names'] <- 'sample'

frag_ratiodf[nrow(frag_ratiodf)+1, ] <- c("median", as.vector(frag_a), "Healthy")

setDT(frag_ratiodf)
set(frag_ratiodf, j = "p76_90", value = as.numeric(frag_ratiodf[["p76_90"]]))
set(frag_ratiodf, j = "p91_120", value = as.numeric(frag_ratiodf[["p91_120"]]))
set(frag_ratiodf, j = "p121_150", value = as.numeric(frag_ratiodf[["p121_150"]]))
set(frag_ratiodf, j = "p151_180", value = as.numeric(frag_ratiodf[["p151_180"]]))
set(frag_ratiodf, j = "p181_200", value = as.numeric(frag_ratiodf[["p181_200"]]))
set(frag_ratiodf, j = "p201_250", value = as.numeric(frag_ratiodf[["p201_250"]]))
set(frag_ratiodf, j = "p251_300", value = as.numeric(frag_ratiodf[["p251_300"]]))
set(frag_ratiodf, j = "p301_350", value = as.numeric(frag_ratiodf[["p301_350"]]))
set(frag_ratiodf, j = "p351_400", value = as.numeric(frag_ratiodf[["p351_400"]]))
set(frag_ratiodf, j = "p401_500", value = as.numeric(frag_ratiodf[["p401_500"]]))
set(frag_ratiodf, j = "p501_700", value = as.numeric(frag_ratiodf[["p501_700"]]))
set(frag_ratiodf, j = "p701_999", value = as.numeric(frag_ratiodf[["p701_999"]]))

frag_ratiodf[, `:=`(p76_90 = p76_90/p76_90[sample == "median"],
                    p91_120 = p91_120/p91_120[sample == "median"],
                    p121_150 = p121_150/p121_150[sample == "median"],
                    p151_180 = p151_180/p151_180[sample == "median"],
                    p181_200 = p181_200/p181_200[sample == "median"],
                    p201_250 = p201_250/p201_250[sample == "median"],
                    p251_300 = p251_300/p251_300[sample == "median"],
                    p301_350 = p301_350/p301_350[sample == "median"],
                    p351_400 = p351_400/p351_400[sample == "median"],
                    p401_500 = p401_500/p401_500[sample == "median"],
                    p501_700 = p501_700/p501_700[sample == "median"],
                    p701_999 = p701_999/p701_999[sample == "median"])]
frag_ratiodf[is.na(frag_ratiodf)] <- 0

ml_df <- frag_ratiodf[,-c(8,9,9,10,11,12,13)]
ml_df <- ml_df[-c(539),]
ml_df <- ml_df[- grep("EE87920", ml_df$sample),]
ml_df <- ml_df[- grep("EE87921", ml_df$sample),]

ml_ichor <- ichordf %>% select("sample","tfx")

ml_df <- merge(x = ml_df, y =ml_ichor , by = "sample", all.x = TRUE)

ml_liq <- liq_150[,-15]
ml_liq <- ml_liq[- grep("EE87920", ml_liq$sample),]
ml_liq <- ml_liq[- grep("EE87921", ml_liq$sample),]

ml_df <- merge(x = ml_df, y =ml_liq , by = "sample", all.x = TRUE)

ml_freia <- FrEIA_my_panel_and_tnc %>% select("sample","FrEIA_score")
ml_freia <- ml_freia[- grep("EE87920", ml_freia$sample),]
ml_freia <- ml_freia[- grep("EE87921", ml_freia$sample),]

ml_df <- merge(x = ml_df, y =ml_freia , by = "sample", all.x = TRUE)
ml_df <- ml_df %>% relocate(group, .after = FrEIA_score)
ml_df <- ml_df %>% relocate(tfx, .before = FrEIA_score)


ml_hvc <- ml_df
ml_too <- ml_df

ml_hvc$group[which(ml_hvc$group!="Healthy")] <- "1"
ml_hvc$group[which(ml_hvc$group=="Healthy")] <- "0"

ml_too$group[which(ml_too$group=="Healthy")] <- "0"
ml_too$group[which(ml_too$group=="Lung")] <- "1"
ml_too$group[which(ml_too$group=="Breast")] <- "2"

ml_too$group[which(ml_too$group=="Colorectal")] <- "3"
ml_too$group[which(ml_too$group=="Pancreas")] <- "4"
ml_too$group[which(ml_too$group=="Bile duct")] <- "5"

ml_too$group[which(ml_too$group=="Stomach")] <- "6"
ml_too$group[which(ml_too$group=="Ovarian")] <- "7"
ml_too$group[which(ml_too$group=="Duodenum")] <- "8"

write.csv(ml_hvc, "C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Data\\ML_model\\ml_hvc.csv",row.names=FALSE)

ml_too_2 <- ml_too[- grep("8", ml_too$group),]
write.csv(ml_too, "C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Data\\ML_model\\ml_too.csv",row.names=FALSE)
write.csv(ml_too_2, "C:\\Users\\wadus\\OneDrive\\Desktop\\Studium\\Masterarbeit\\Projects\\Final_Data\\ML_model\\ml_too_2.csv",row.names=FALSE)
#ML results
log_hvc <- c(0.76,0.78, 0.73, 0.75)
knn_hvc <- c(0.63,0.71, 0.49, 0.57)
gnb_hvc <- c(0.79,0.94, 0.64,0.76)
dtc_hvc <- c(0.73,0.75, 0.73, 0.74)
gbc_hvc <- c(0.81,0.83,0.78, 0.81)
ml_hvc_results<- data.frame(log_hvc,knn_hvc,gnb_hvc,dtc_hvc,gbc_hvc)
colnames(ml_hvc_results) <- c('Logistic Regression','K-Nearest-Neighbour','Gaussian Naive Bayes','Decision Tree','Gradient Boosting')
rownames(ml_hvc_results) <- c('Mean Validation Accuracy', 'Mean Validation Precision', 'Mean Validation Recall', 'Mean Validation F1 Score')
hvc_table <-kable(ml_hvc_results) %>% kable_styling(latex_options = 'striped')

log_too <- c(0.58,0.38, 0.30, 0.32)
knn_too <- c(0.51,0.29, 0.17, 0.17)
gnb_too <- c(0.58,0.40, 0.34,0.35)
dtc_too <- c(0.47,0.25, 0.26, 0.25)
gbc_too <- c(0.55,0.26,0.25, 0.24)
ml_too_results<- data.frame(log_too,knn_too,gnb_too,dtc_too,gbc_too)
colnames(ml_too_results) <- c('Logistic Regression','K-Nearest-Neighbour','Gaussian Naive Bayes','Decision Tree','Gradient Boosting')
rownames(ml_too_results) <- c('Mean Validation Accuracy', 'Mean Validation Precision', 'Mean Validation Recall', 'Mean Validation F1 Score')
too_table <-kable(ml_too_results) %>% kable_styling(latex_options = 'striped')
