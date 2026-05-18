library(tidyverse)
library(ggrepel)
library(magrittr)
Day0v8g <-as.data.frame(read.csv("/Volumes/Home/sadiehl/MedMyDocs/DENV2-RNASeq/VIGR/For_Sean_Nov13_2019/Group_Day0_8_Plot_df.csv",header = TRUE, sep = ","))
head(Day0v8g)

Day0v8ga<- Day0v8g[-grep('MSTRG|Or|AC|AL', Day0v8g$Gene),]
Day0v8ga2 <-mutate(Day0v8ga, Cutoff = ifelse(Log2CPM >= 1, "Threshold", "Non-Threshold"))
write.csv(Day0v8ga2,"~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Day0v8ga2.csv")
Day0v8ga3 <- subset(Day0v8ga2,Cutoff =="Threshold")                            

##Fig 3c
ggplot(Day0v8ga3, aes(x=Log2CPM,y=Log2FC, label=Gene))+
  scale_colour_manual(values = c("Non-Threshold"= "grey", "Threshold"="red"))+
    geom_point(aes(color=Group),alpha=1/2)+
  scale_x_continuous(breaks = seq(0, 14, 2),limits=c(0, 14))+
  scale_y_continuous(breaks = seq(-7, 7, 2),limits=c(-7, 7))+
        geom_text_repel(
                  data = subset(Day0v8ga, 
                  Group =="Threshold" & Log2CPM >=2 & abs(Log2FC)>3),
                  aes(label = Gene),
                  size=5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  force=2)+
  geom_hline(yintercept = 0.6, linetype="solid")+
  geom_hline(yintercept = -0.6, linetype="solid")+
  geom_hline(yintercept = 0, linetype="dashed")+
   labs(title = "Day 0 vs Day 8")+
  ylab(bquote(Fold-Change~(log[2])))+
  theme(axis.text.y = element_text(size=24, color="black"),
        axis.text.x   = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        plot.title = element_text(size=24, color="black"),
        panel.border = element_rect(color = "black",fill=NA, size = 1),
        axis.ticks.length=unit(3, "mm"),
        legend.key = element_rect(fill = "grey90", colour = "black"))

setwd("~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq")
ggsave("g0v8v3.png", width=8, height=6.45, units="in", dpi=300,device='png')

#Original version Fig 2c

Day8v28g <-as.data.frame(read.csv("/Volumes/Home/sadiehl/MedMyDocs/DENV2-RNASeq/VIGR/For_Sean_Nov13_2019/Group_Day8_28_Plot_df.csv",header = TRUE, sep = ","))
Day8v28ga<- Day8v28g[-grep('MSTRG|Or|AC|AL', Day8v28g$Gene),]
Day8v28ga2 <-mutate(Day8v28ga, Cutoff = ifelse(Log2CPM >= 1, "Threshold", "Non-Threshold"))
write.csv(Day8v28ga2,"~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Day8v28ga2.csv")
Day8v28ga3 <- subset(Day8v28ga2,Cutoff =="Threshold")      

ggplot(Day8v28ga3, aes(x=Log2CPM,y=Log2FC, label=Gene))+
  scale_colour_manual(values = c("Non-Threshold"= "grey", "Threshold"="red"))+
  geom_point(aes(color=Group),alpha=1/2)+
  scale_x_continuous(breaks = seq(0, 14, 2),limits=c(0, 14))+
  scale_y_continuous(breaks = seq(-7, 7, 2),limits=c(-7, 7))+
  geom_text_repel(
    data = subset(Day8v28ga3, 
                  Group =="Threshold" & Log2CPM >=2 & Log2FC>1.3 | Log2FC < -3 ),
    aes(label = Gene),size=5, 
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    force=2)+
  geom_hline(yintercept = 0.6, linetype="solid")+
  geom_hline(yintercept = -0.6, linetype="solid")+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(title = "Day 8 vs Day 28")+
  ylab(bquote(Fold-Change~(log[2])))+
  theme(axis.text.y = element_text(size=24, color="black"),
        axis.text.x   = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        plot.title = element_text(size=24, color="black"),
        panel.border = element_rect(color = "black",fill=NA, size = 1),
        axis.ticks.length=unit(3, "mm"),
        legend.key = element_rect(fill = "grey90", colour = "black"))
ggsave("g8v28v3.png", width=8, height=6.45, units="in", dpi=300,device='png')


Day0v28g <-as.data.frame(read.csv("/Volumes/Home/sadiehl/MedMyDocs/DENV2-RNASeq/VIGR/For_Sean_Nov13_2019/Group_Day0_28_Plot_df.csv",header = TRUE, sep = ","))
Day0v28ga<- Day0v28g[-grep('MSTRG|Or|AC|AL', Day0v28g$Gene),]
Day0v28ga2 <-mutate(Day0v28ga, Cutoff = ifelse(Log2CPM >= 1, "Threshold", "Non-Threshold"))
write.csv(Day0v28ga2,"~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Day0v28ga2.csv")

Day0v28ga3 <- subset(Day0v28ga2,Cutoff =="Threshold")    

ggplot(Day0v28ga3, aes(x=Log2CPM,y=Log2FC, label=Gene))+
  scale_colour_manual(values = c("Non-Threshold"= "grey", "Threshold"="red"))+
  geom_point(aes(color=Group),alpha=1/2)+
  scale_x_continuous(breaks = seq(0, 14, 2),limits=c(0, 14))+
  scale_y_continuous(breaks = seq(-7, 7, 2),limits=c(-7, 7))+
  geom_text_repel(
    data = subset(Day0v28ga3, 
                  Group =="Threshold" & Log2CPM >=2 & Log2FC>1 | Log2FC < -2.5 ),
    aes(label = Gene),size=5, 
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    force=2)+
  geom_hline(yintercept = 0.6, linetype="solid")+
  geom_hline(yintercept = -0.6, linetype="solid")+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(title = "Day 0 vs Day 28")+
  ylab(bquote(Fold-Change~(log[2])))+
  theme(axis.text.y = element_text(size=24, color="black"),
        axis.text.x   = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        plot.title = element_text(size=24, color="black"),
        panel.border = element_rect(color = "black",fill=NA, size = 1),
        axis.ticks.length=unit(3, "mm"),
        legend.key = element_rect(fill = "grey90", colour = "black"))
  ggsave("g0v28v3.png", width=8, height=6.45, units="in", dpi=300,device='png')