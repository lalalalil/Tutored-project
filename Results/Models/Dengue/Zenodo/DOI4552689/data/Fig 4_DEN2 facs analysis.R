library(readxl)
library(tidyverse)
library(reshape2)
library(lemon)

DEN2facs <-read_excel("~/OneDrive - UVM Larner College of Medicine/Manuscripts/DENV2-RNASeq/Results/FACS/DEN2_pbmc_facs112420.xlsx",sheet = 2)
                      
DEN2facs.m <- melt(DEN2facs, id.vars = c("day","DeID"),  measure.vars = c("bcells","pb","s_mbc","naive_b","c_mono","i_mono","nc_mono",
                                                                          "nkt","NK1","NK2","NK3","NK4","NK5","DC","DC_act","tcells","t_regs","cd4t_naive","cd4t_cm","cd4t_em","cd4t_emra","CD4_homeo","CD4_cyto","CD4_act",
                                                                          "cd8t_naive","cd8t_cm","cd8t_em","cd8t_emra","CD8_homeo","CD8_cyto","CD8_act","CD4_AIM","CD4_AIM_old"), 
                   variable.name = "testname", value.name = "testresult")

DEN2facs.m.s <-DEN2facs.m %>%
  select(DeID, day, testname, testresult) %>%
  filter(DeID %in% c("B","C","D","E","F","G","H","I") & day %in% c("0","8","28"))

DEN2facs.m.s$day<-as.factor(DEN2facs.m.s$day)
DEN2facs.m.s$testname_f = factor(DEN2facs.m.s$testname, levels=c("bcells","pb","naive_b","s_mbc","c_mono","i_mono","nc_mono","nkt","NK1","NK2","NK3","NK4","NK5","DC","DC_act", "tcells", 
                                                             "cd4t_naive","cd4t_cm","cd4t_em","cd4t_emra","CD4_homeo","CD4_cyto","CD4_act","t_regs", "CD4_AIM","CD4_AIM_old", 
                                                             "cd8t_naive","cd8t_cm","cd8t_em","cd8t_emra","CD8_homeo","CD8_cyto","CD8_act"))
ggplot(DEN2facs.m.s,aes(day, testresult,fill=day)) + geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("red", "green", "blue"))+
  facet_rep_wrap(~testname_f, repeat.tick.labels = 'bottom',ncol=5,scales="free_y")

#statistics

B<- DEN2facs %>% filter(day %in% c("0","8","28") & DeID %in% c("B","C","D","E","F","G","H","I"))
B$day<-as.factor(B$day)

bcells <- kruskal.test(bcells ~ day, data = B)
bcells
pb <- kruskal.test(pb ~ day, data = B)
pb
naive_b <-kruskal.test(naive_b ~ day, data =B)
naive_b
s_mbc <-kruskal.test(s_mbc~day,data=B)
s_mbc

DC_act <- kruskal.test(DC_act ~ day, data = B)
DC_act
DCs <- kruskal.test(DC ~ day, data = B)
DCs

c_mono <- kruskal.test(c_mono ~ day, data = B)
c_mono
i_mono <- kruskal.test(i_mono ~ day, data = B)
i_mono
nc_mono <- kruskal.test(nc_mono ~ day, data = B)
nc_mono

nkt <- kruskal.test(nkt ~ day, data = B)
nkt
nk1 <- kruskal.test(NK1 ~ day, data = B)
nk1
nk2 <- kruskal.test(NK2 ~ day, data = B)
nk2
nk3 <- kruskal.test(NK3 ~ day, data = B)
nk3
nk4 <- kruskal.test(NK4 ~ day, data = B)
nk4
nk5 <- kruskal.test(NK5 ~ day, data = B)
nk5

tcells <-kruskal.test(tcells ~ day, data = B)
tcells
CD4 <- kruskal.test(cd4 ~ day, data = B)
CD4
CD4_act <- kruskal.test(CD4_act ~ day, data = B)
CD4_act
CD4_cyto <- kruskal.test(CD4_cyto ~ day, data = B)
CD4_cyto
t_regs <- kruskal.test(t_regs ~ day, data = B)
t_regs
CD4_AIM_old <-kruskal.test(CD4_AIM_old ~ day, data = B)
CD4_AIM_old
CD4_AIM <-kruskal.test(CD4_AIM ~ day, data = B)
CD4_AIM
cd4_homeo <-kruskal.test(CD4_homeo ~ day, data = B)
cd4_homeo
CD4t_naive <- kruskal.test(cd4t_naive ~ day, data = B)
CD4t_naive
CD4t_cm <- kruskal.test(cd4t_cm ~ day, data = B)
CD4t_cm
CD4t_em <- kruskal.test(cd4t_em ~ day, data = B)
CD4t_em
CD4t_emra <- kruskal.test(cd4t_emra ~ day, data = B)
CD4t_emra

CD8 <- kruskal.test(cd8 ~ day, data = B)
CD8
CD8_act <- kruskal.test(CD8_act ~ day, data = B)
CD8_act
CD8_cyto <- kruskal.test(CD8_cyto ~ day, data = B)
CD8_cyto
CD8t_naive <- kruskal.test(cd8t_naive ~ day, data = B)
CD8t_naive
CD8t_cm <- kruskal.test(cd8t_cm ~ day, data = B)
CD8t_cm
CD8t_em <- kruskal.test(cd8t_em ~ day, data = B)
CD8t_em
CD8t_emra <- kruskal.test(cd8t_emra ~ day, data = B)
CD8t_emra
cd8_homeo <-kruskal.test(CD8_homeo ~ day, data = B)
cd8_homeo
