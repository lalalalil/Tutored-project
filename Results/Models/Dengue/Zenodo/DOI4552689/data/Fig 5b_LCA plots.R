library(tidyverse)
library(ggrepel)

Day0v8lca <-as.data.frame(read.csv("~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Manuscript/ExcelFiles_4_Sean_June10_2020/LCA_Day_8_v_0.csv",header = TRUE, sep = ","))
head(Day0v8lca)
Day0v8lca2<- Day0v8lca[-grep('MSTRG|AC|AL|Or', Day0v8lca$Gene),]
Day0v8lca3 <-mutate(Day0v8lca2, Cutoff = ifelse(Log2CPM >= 1, "Threshold", "Non-Threshold"))
write.csv(Day0v8lca3,"~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Day0v8lca3.csv")


Day8v28lca <-as.data.frame(read.csv("~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Manuscript/ExcelFiles_4_Sean_June10_2020/LCA_Day_28_v_8.csv",header = TRUE, sep = ","))
head(Day8v28lca)
Day8v28lca2<- Day8v28lca[-grep('MSTRG|AC|AL|Or', Day8v28lca$Gene),]
Day8v28lca3 <-mutate(Day8v28lca2, Cutoff = ifelse(Log2CPM >= 1, "Threshold", "Non-Threshold"))
write.csv(Day8v28lca3,"~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Day8v28lca3.csv")

Day0v28lca <-as.data.frame(read.csv("~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Manuscript/ExcelFiles_4_Sean_June10_2020/LCA_Day_28_v_0.csv",header = TRUE, sep = ","))
head(Day8v28lca)
Day0v28lca2<- Day0v28lca[-grep('MSTRG|AC|AL|Or', Day0v28lca$Gene),]
Day0v28lca3 <-mutate(Day0v28lca2, Cutoff = ifelse(Log2CPM >= 1, "Threshold", "Non-Threshold"))
write.csv(Day0v28lca3,"~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Day0v28lca3.csv")

ggplot(Day0v8a, aes(x=Log2CPM,y=Log2FC, label=Gene))+
  geom_hline(yintercept = 1, linetype="solid")+
  geom_hline(yintercept = -1, linetype="solid")+
  scale_colour_manual(values = c("Non-Threshold"= "grey", "Threshold"="black"))+
  scale_size_manual(values=c(1, 3))+
  scale_fill_manual(values=c("Non-Threshold"= "grey", "Threshold"="#00BFC4"))+
   geom_point(aes(size=Group,fill=Group,color=Group), pch=21)+
  coord_cartesian(xlim = c(-3,14),ylim = c(-3,3))+
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))+
  geom_text_repel(
    data = subset(Day0v8a, 
                  Group =="Threshold" & abs(Log2FC)>1),
    aes(label = Gene),
    size=5,box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
    force=5)+
  labs(title = "Day 8 vs Day 0")+
  theme(axis.text.y = element_text(size=24, color="black"),
        axis.text.x   = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        plot.title = element_text(size=24, color="black"),
        panel.border = element_rect(color = "black",fill=NA, size = 1),
        axis.ticks.length=unit(3, "mm"),
        legend.key = element_rect(fill = "grey90", colour = "black"))

setwd("~/OneDrive - UVM Larner College of Medicine/DENV2-RNASeq/Manuscript/Figures")
ggsave("i0v8v2.png", width=8, height=6.45, units="in", dpi=300,device='png')


#plot day 8 vs 28 individual subjects
Day8v28 <-as.data.frame(read.csv("/Volumes/sadiehl/MedMyDocs/DENV2-RNASeq/VIGR/For_Sean_Nov13_2019/Indiv_Day8_28_Plot_df.csv",header = TRUE, sep = ","))
Day8v28a<- Day8v28[-grep('MSTRG|AC|AL|Or', Day8v28$Gene),]

ggplot(Day8v28a, aes(x=Log2CPM,y=Log2FC, label=Gene))+
  geom_hline(yintercept = 1, linetype="solid")+
  geom_hline(yintercept = -1, linetype="solid")+
  scale_colour_manual(values = c("Non-Threshold"= "grey", "Threshold"="black"))+
  scale_size_manual(values=c(1, 3))+
  scale_fill_manual(values=c("Non-Threshold"= "grey", "Threshold"="#00BFC4"))+
  geom_point(aes(size=Group,fill=Group,color=Group), pch=21)+
  coord_cartesian(xlim = c(-3,14),ylim = c(-3,3))+
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))+
  geom_text_repel(
    data = subset(Day8v28a, 
                  Group =="Threshold" & abs(Log2FC)>1),
    aes(label = Gene),
    size=5,box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
    force=5)+
  labs(title = "Day 28 vs Day 8")+
  theme(axis.text.y = element_text(size=24, color="black"),
        axis.text.x   = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        plot.title = element_text(size=24, color="black"),
        panel.border = element_rect(color = "black",fill=NA, size = 1),
        axis.ticks.length=unit(3, "mm"),
        legend.key = element_rect(fill = "grey90", colour = "black"))

ggsave("i8v28v2.png", width=8, height=6.45, units="in", dpi=300,device='png')


#plot day 0 vs 28 individual subjects
Day0v28 <-as.data.frame(read.csv("/Volumes/sadiehl/MedMyDocs/DENV2-RNASeq/VIGR/For_Sean_Nov13_2019/Indiv_Day0_28_Plot_df.csv",header = TRUE, sep = ","))
Day0v28a<- Day0v28[-grep('MSTRG|AC|AL', Day0v28$Gene),]

ggplot(Day0v28a, aes(x=Log2CPM,y=Log2FC, label=Gene))+
  geom_hline(yintercept = 1, linetype="solid")+
  geom_hline(yintercept = -1, linetype="solid")+
  scale_colour_manual(values = c("Non-Threshold"= "grey", "Threshold"="black"))+
  scale_size_manual(values=c(1, 3))+
  scale_fill_manual(values=c("Non-Threshold"= "grey", "Threshold"="#00BFC4"))+
  geom_point(aes(size=Group,fill=Group,color=Group), pch=21)+
  coord_cartesian(xlim = c(-3,14),ylim = c(-3,3))+
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))+
  geom_text_repel(
    data = subset(Day0v28a, 
                  Group =="Threshold"),
    aes(label = Gene),
    size=5,box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
    force=5)+
  labs(title = "Day 28 vs Day 8")+
  theme(axis.text.y = element_text(size=24, color="black"),
        axis.text.x   = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        plot.title = element_text(size=24, color="black"),
        panel.border = element_rect(color = "black",fill=NA, size = 1),
        axis.ticks.length=unit(3, "mm"),
        legend.key = element_rect(fill = "grey90", colour = "black"))
ggsave("i0v28v2.png", width=8, height=6.45, units="in", dpi=300,device='png')
