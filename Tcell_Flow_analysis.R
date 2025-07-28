

# Code Written by Laila M Rad
# Lonnie Shea Lab
# Updated July 2025


#Flow Gating
#CD4 and CD8 T cells
#Treg: CD4+CD25+FoxP3+
#Gut-homing markers: CCR9 and LPAM-1

# File path ----
path = "~/Documents/UMich Research/Allergy/Cohort 17/"

setwd(path)

# Load Libraries ----
library(ggplot2)
library(ggprism)
library(stringr)
library(stats)

library(dplyr)
library(reshape2) 

library(tidyverse)
library(rstatix)
library(ggpubr)
library(rlang)
library(glue)
library(multcomp)
library(emmeans)
library(car)

#Load Data----

#find dataset in location
file = list.files(path=path, pattern = ".csv" )
file


data <-read.csv("TcellflowjoCounts.csv", check.names=FALSE,header=T, sep=",")
data = data[,1:15]
flowlabels = read.csv( "flowsamplelabels.csv" , check.names=FALSE,header=T, sep=",")
colnames(data)[1] = "Sample"
data$Name = sub("[[:alpha:]]", "", data$Sample)
data$Name = sub("[[:digit:]]", "", data$Name)
data$Name = sub("^[[:digit:]]", "", data$Name)
data$Name = sub(" ", "", data$Name)
data$Name = sub("\\.fcs*", "", data$Name)

data$Tissue = substr(data$Name, start = 1, stop = 2)




counts = data[which(data$Tissue == "SI"),]
counts$"Flow ID" = gsub("[^0-9]","",counts$Name)
counts = merge(counts, flowlabels, by = "Flow ID")
counts$Depletion = counts$Condition
counts[which(counts$Depletion == "Isotype"),"Depletion"] = "Ctrl Ab"
counts$Depletion = factor(counts$Depletion, levels = c("Ctrl Ab", "anti-CD20"))

counts[which(counts$Tx == "PLG(OVA)"),"Tx"] = "OVA NPs"
counts$Tx = factor(counts$Tx, levels = c("PBS", "OVA NPs"))

#absolute counts
HCcounts <-read.csv("HemocytometerCounts.csv", check.names=FALSE,header=T, sep=",")

counts = merge(counts,HCcounts[,c(2,8)], by = 'Name')

counts$Condition = paste0(counts$Depletion, "_", counts$Tx)
#counts = counts[-which(counts$Name == "SILP 6"),]
Counts = counts #save original data

#Analyze counts data----
counts$"% Live Cells" = (counts$`Live Cell Count`/counts$`Single Cell Count`)*100

# make percents columns 
#colnum = 23 #column number of live cells %
colnum = 24
colname_list = list()
#live cell count to last count 
#for (i in (5:16)){
for (i in (6:17)){
  colnum = colnum + 1
  colname = paste0("% ", colnames(counts)[i])
  colname = sub(" Count", "", colname)
  colname_list[[length(colname_list) + 1]] = colname
  counts[,colnum] = (counts[,i]/counts$`Live Cell Count`)*100
}
colnames(counts)[(colnum-length(colname_list)+1):colnum] = colname_list




colnum = 36 # last column
colname_list = list()
for (i in c(25:36)){
  colnum = colnum + 1
  colname = paste0("Total ", colnames(counts)[i], " Counts")
  colname = gsub('% ', '', colname)
  colname_list[[length(colname_list) + 1]] = colname
  #counts[,colnum] = ((counts[,i]/100)*(counts$`% Live Cells`/100)*counts$`Total Cell Count`)
  counts[,colnum] = (counts[,i]/100)*counts$"Total Cell Count"
}
colnames(counts)[(colnum-length(colname_list)+1):colnum] = colname_list

colnum = 48
colname_list = list()
#first  Tcell counts
for (i in c(7:10,15:17)){
  colnum = colnum + 1
  colname = gsub(' Count', '', colnames(counts)[i])
  colname = paste0("% ", colname, "s (of CD4+ T Cells)")
  colname_list[[length(colname_list) + 1]] = colname
  counts[,colnum] = ((counts[,i])/(counts$`CD4+ T Cell Count`))*100
}
colnames(counts)[(colnum-length(colname_list)+1):colnum] = colname_list

colnum = 55
colname_list = list()
#first  Tcell counts
for (i in c(12:14)){
  colnum = colnum + 1
  colname = gsub(' Count', '', colnames(counts)[i])
  colname = paste0("% ", colname, "s (of CD8+ T Cells)")
  colname_list[[length(colname_list) + 1]] = colname
  counts[,colnum] = ((counts[,i])/(counts$`CD8+ T Cell Count`))*100
}
colnames(counts)[(colnum-length(colname_list)+1):colnum] = colname_list


#Two way Anova Analysis ----


#for(i in c(19:45)){
#for (i in c(4:17,23:58)){
for (i in c(26:29,38:41,49:52)){
  #for (i in c(33)){
  colnum = i
  counts_data = counts
  counts_data = counts_data[-which(counts$`Mouse ID` == "C3M5" | counts$`Mouse ID` == "C4M1"),]
  colname = colnames(counts_data)[colnum]
  print(paste(i, colname))
  
  data_subset <- counts_data[c("Tx", "Depletion" ,paste(colname))]
  colnames(data_subset) = c("Tx","Depletion", "value")
  #data_subset$Depletion = factor(data_subset$Depletion)
  #data_subset$Tx = factor(data_subset$Tx)
  
  type3 <- list(Depletion = contr.sum, Tx = contr.sum)
  m2 <- lm(value ~ Tx * Depletion , data = data_subset, contrasts = type3)
  print(car::Anova(m2, type = 3 )) # correct type 3
  
  
  #check assumptions
  #check the equality of variances
  print(car::leveneTest(value ~ Tx * Depletion, data = data_subset))
  
  #check normality 
  print(m2$residuals %>% 
          shapiro_test())
  
  
  with(data_subset, interaction.plot(Tx, Depletion, value, fun = mean,
                                     main = "Interaction Plot"))
  
  #res.aov3 <- aov(value ~ Tx * Condition , data = data_subset, type = "III")
  #print(summary(res.aov3))
  
  # Post-hoc tests using Tukey's HSD
  library(multcomp)
  
  #print(tukey_hsd(m2),        which = "Tx:Depletion")
  
  EMM <- emmeans(m2,  ~ Depletion * Tx)   # where treat has 2 levels
  pairs(EMM, adjust = "sidak")   # adjustment is ignored - only 1 test per group
  print(summary(pairs(EMM)[c(1:2,5:6)], by = NULL, adjust = "sidak") )  # all are in one group now
  
  #posthoc <- glht(m2, linfct = mcp( Tx = "Tukey"))
  #summary(posthoc)
  
}
dev.off()

#plot data 2-way ANOVA stats ----

plotFlowdataA2 = function(colnum){
  
  
  data_subset = counts
  data_subset = counts[-which(counts$`Mouse ID` == "C3M5" | counts$`Mouse ID` == "C4M1"),]
  
  colnum = colnum
  colname = colnames(data_subset)[colnum]
  colname
  colnames(data_subset)[colnum] = "value"
  
  tab_cs = data_subset %>%
    group_by(Tx, Depletion)%>%
    dplyr::summarise(n = n(), value = mean(value)) 
  
  
  max.value = max(data_subset$value)
  tab_cs$value = max.value/25
  tab_cs$print = paste0("n = ", tab_cs$n)
  tab_cs
  
  type3 <- list(Depletion = contr.sum, Tx = contr.sum)
  m2_cs <- lm(value ~ Tx * Depletion , data = data_subset, contrasts = type3)
  a2_cs = Anova(m2_cs, type = 3 )
  a2_cs = a2_cs %>% as.data.frame() 
  
  
  EMMcs <- emmeans(m2_cs,  ~ Depletion * Tx)   # where treat has 2 levels
  stat.test_cs = summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak") 
  stat.test_cs
  
  stat.test_cs_pwc = stat.test_cs %>% tibble(".y." = "value", 
                                             "group1" = "Isotype", 
                                             "group2" = "anti-CD20",
                                             "p.adj" = round(p.value,4),
                                             "method" = "sidak") %>%
    add_x_position() %>% mutate(xmin = c(0.85,0.8,1.8,1.15), xmax = c(1.85,1.2,2.2,2.15)) %>%
    add_significance() 
  
  ggplot(data_subset, aes(x = Depletion, y = value)) +
    stat_summary(aes(color = Tx,group = Tx),
                 fun.min = function(x) mean(x), 
                 fun.max = function(x) mean(x) + sd(x), 
                 geom = "errorbar", width = 0.2,
                 position = position_dodge(0.6),
                 color = "black"
    ) +
    stat_summary(
      fun = mean, geom="col", 
      position = "dodge",
      aes(fill = Tx),
      color = "black", 
      size = 0.5, 
      width = 0.6
    ) +
    scale_fill_manual(values = c("grey", "#89CFF0")) +
    geom_point(aes(fill = Tx),
               position = position_jitterdodge(
                 jitter.width = 0.12,
                 jitter.height = 0,
                 dodge.width = 0.6),
               color = "black",show.legend = F, size = 0.5
    ) + 
    geom_text(data = tab_cs, aes(label = print, group = Tx), 
              position = position_dodge(width = 0.6),
              color = "white", size = 1.5) +
    ggtitle(str_wrap(colname, 20)) + 
    labs(y = "",
         x = "",
         caption =  paste0("Two-Way Anova, post-hoc: Sidak Test \n[",rownames(a2_cs)[2],"] F(", a2_cs["Tx","Df"],",",a2_cs["Residuals","Df"],") = ",round(a2_cs["Tx","F value"], 4) ,", p = ",round(a2_cs["Tx","Pr(>F)"], 4) ," \n[",rownames(a2_cs)[3],"] F(", a2_cs["Depletion","Df"],",",a2_cs["Residuals","Df"],") = "
                           ,round(a2_cs["Depletion","F value"], 4),", p = ",  round(a2_cs["Depletion","Pr(>F)"], 4), " \n[",rownames(a2_cs)[4],
                           "] F(", a2_cs["Tx:Depletion","Df"],",",a2_cs["Residuals","Df"],") = ", round(a2_cs["Tx:Depletion","F value"], 4),", p = ",  round(a2_cs["Tx:Depletion","Pr(>F)"], 4), " "
         ))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    coord_cartesian(clip = "off")+
    theme_prism(base_family = "Arial", base_fontface = "bold", base_size = 8) + 
    scale_colour_prism(palette = "black_and_white")+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          legend.position = "top", legend.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0,face = "plain", size = 4),
          plot.caption.position = "plot",
          plot.title.position = "plot"
          #plot.subtitle = element_text(size = 12),
    ) + 
    stat_pvalue_manual(
      stat.test_cs_pwc[1,],  y.position = max.value, bracket.nudge.y = max.value*.5,
      step.group.by = "contrast",
      label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, 
      lineheight = 0.1, size = 2
    ) + 
    stat_pvalue_manual(
      stat.test_cs_pwc[2,],  y.position = max.value, bracket.nudge.y = max.value*.1,
      step.group.by = "contrast",
      label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
      lineheight = 0, size = 2
    ) + 
    stat_pvalue_manual(
      stat.test_cs_pwc[3,],  y.position = max.value, bracket.nudge.y = max.value*.1,
      step.group.by = "contrast",
      label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
      lineheight = 0.1, size = 2
    ) + 
    stat_pvalue_manual(
      stat.test_cs_pwc[4,],  y.position = max.value, bracket.nudge.y = max.value*.3,
      step.group.by = "contrast",
      label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
      lineheight = 0.1, size = 2
    ) 
  
}

plotFlowdataA2(38)
length( c(26:29,38:41,49:52))
myplots = lapply( c(26:29,38:41,49:52), FUN = plotFlowdataA2)

cairo_pdf("TregFlowdata_twoANOVA.pdf", width = 8, height =9)
png("cohort17TregFlowdata_2wayAnova_Sidak_removelowAb_byhand.png",    units = "in", width = 8 , height = 9, res = 400)
ggarrange(plotlist = myplots, ncol=4, nrow=3, 
          common.legend = T, legend = "top")
dev.off()


#Plot specific columns: ----

#Plot %tregs ----
data_subset = counts
data_subset = counts[-which(counts$`Mouse ID` == "C3M5" | counts$`Mouse ID` == "C4M1"),]
data_subset$`Mouse ID`

type3 <- list(Depletion = contr.sum, Tx = contr.sum)
colnum = 26
colname = colnames(data_subset)[colnum]
colname
colnames(data_subset)[colnum] = "value"
m2_cs <- lm(value ~ Tx * Depletion , data = data_subset, contrasts = type3)
print(car::Anova(m2_cs, type = 3 )) # correct type 3

a2_cs = car::Anova(m2_cs, type = 3 )
#check assumptions
#check the equality of variances
print(car::leveneTest(value ~ Tx * Depletion, data = data_subset))

#check normality 
print(m2_cs$residuals %>% 
        shapiro_test())


with(data_subset, interaction.plot(Tx, Depletion, value, fun = mean,
                                   main = "Interaction Plot"))


# Post-hoc tests using Tukey's HSD

print(tukey_hsd(m2_cs, which = "Tx:Depletion"))

#posthoc <- glht(m2, linfct = mcp( Tx = "Tukey"))
#summary(posthoc)

EMMcs <- emmeans(m2_cs,  ~ Depletion * Tx)   # where treat has 2 levels
pairs(EMMcs, adjust = "sidak")   # adjustment is ignored - only 1 test per group
summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak")   # all are in one group now
stat.test_cs = summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak") 
stat.test_cs

stat.test_cs_pwc = stat.test_cs %>% tibble(".y." = "value", 
                                           "group1" = "Isotype", 
                                           "group2" = "anti-CD20",
                                           "p.adj" = round(p.value,4),
                                           "method" = "sidak") %>%
  add_x_position() %>% mutate(xmin = c(0.85,0.8,1.8,1.15), xmax = c(1.85,1.2,2.2,2.15)) %>%
  add_significance()


stat.test_cs_pwc



a2_cs = a2_cs %>% as.data.frame() 


ggboxplot(data_subset, x = "Depletion", y = "value", 
          color = "Tx",group = "Tx") +
  stat_anova_test(aes(group = Tx), method = "two_way",
                  type = 3)


tab_cs = data_subset %>%
  group_by(Tx, Depletion)%>%
  dplyr::summarise(n = n(), value = mean(value)) 

tab_cs$value = max(data_subset$value)/25
tab_cs$print = paste0("n = ", tab_cs$n)
tab_cs
max(data_subset$value)


cairo_pdf("cohort17TregPer_2wayAnova_Sidak_removelowAb_byhand.pdf",          width = 2 , height = 3)
#png("cohort17TregPer_2wayAnova_Sidak_removelowAb_byhand.png",    units = "in", width = 4.2 , height = 6, res = 400)
ggplot(data_subset, aes(x = Depletion, y = value)) +
  stat_summary(aes(color = Tx,group = Tx),
               fun.min = function(x) mean(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2,
               position = position_dodge(0.6),
               color = "black"
  ) +
  stat_summary(
    fun = mean, geom="col", 
    position = "dodge",
    aes(fill = Tx),
    color = "black", 
    size = 0.5, 
    width = 0.6
  ) +
  scale_fill_manual(values = c("grey", "#89CFF0")) +
  geom_point(aes(fill = Tx),
             position = position_jitterdodge(
               jitter.width = 0.12,
               jitter.height = 0,
               dodge.width = 0.6),
             color = "black",show.legend = F, size = 0.5
  ) + 
  geom_text(data = tab_cs, aes(label = print, group = Tx), 
            position = position_dodge(width = 0.6),
            color = "white", size = 1.5) +
  labs(#subtitle ="Two-way Anova, post-hoc: Sidak Test", 
    #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001\n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
    caption =  paste0("Two-Way Anova, post-hoc: Sidak Test \n[",rownames(a2_cs)[2],"] F(", a2_cs["Tx","Df"],",",a2_cs["Residuals","Df"],") = ",round(a2_cs["Tx","F value"], 4) ,", p = ",round(a2_cs["Tx","Pr(>F)"], 4) ," .\n[",rownames(a2_cs)[3],"] F(", a2_cs["Depletion","Df"],",",a2_cs["Residuals","Df"],") = "
                      ,round(a2_cs["Depletion","F value"], 4),", p = ",  round(a2_cs["Depletion","Pr(>F)"], 4), " *\n[",rownames(a2_cs)[4],
                      "] F(", a2_cs["Tx:Depletion","Df"],",",a2_cs["Residuals","Df"],") = ", round(a2_cs["Tx:Depletion","F value"], 4),", p = ",  round(a2_cs["Tx:Depletion","Pr(>F)"], 4), " ."
    ),
    y = "Percent") +
  ggtitle(paste(colname)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),  breaks = seq(0, 1.2, 0.4))+
  coord_cartesian(clip = "off")+
  theme_prism(base_family = "Arial", base_fontface = "bold", base_size = 8) + 
  scale_colour_prism(palette = "black_and_white")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0,face = "plain", size = 4),
        plot.caption.position = "plot",
        plot.title.position = "plot"
        #plot.subtitle = element_text(size = 12),
  ) +
  #stat_anova_test(aes(group = Tx), method = "two_way",type = 3) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[1,],  y.position =1.8, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 700000,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[2,],  y.position = 1.2, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,#bracket.nudge.y = 200000,
    lineheight = 0, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[3,],  y.position = 1.2, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[4,],  y.position = 1.5, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 400000,
    lineheight = 0.1, size = 2
  ) 
dev.off()
#Plot %tregs of CD4 ----
data_subset = counts
data_subset = counts[-which(counts$`Mouse ID` == "C3M5" | counts$`Mouse ID` == "C4M1"),]
data_subset$`Mouse ID`

type3 <- list(Depletion = contr.sum, Tx = contr.sum)
colnum = 49
colname = colnames(data_subset)[colnum]
colname
colnames(data_subset)[colnum] = "value"
m2_cs <- lm(value ~ Tx * Depletion , data = data_subset, contrasts = type3)
print(car::Anova(m2_cs, type = 3 )) # correct type 3

a2_cs = car::Anova(m2_cs, type = 3 )
#check assumptions
#check the equality of variances
print(car::leveneTest(value ~ Tx * Depletion, data = data_subset))

#check normality 
print(m2_cs$residuals %>% 
        shapiro_test())


with(data_subset, interaction.plot(Tx, Depletion, value, fun = mean,
                                   main = "Interaction Plot"))


# Post-hoc tests using Tukey's HSD

print(tukey_hsd(m2_cs, which = "Tx:Depletion"))

#posthoc <- glht(m2, linfct = mcp( Tx = "Tukey"))
#summary(posthoc)

EMMcs <- emmeans(m2_cs,  ~ Depletion * Tx)   # where treat has 2 levels
pairs(EMMcs, adjust = "sidak")   # adjustment is ignored - only 1 test per group
summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak")   # all are in one group now
stat.test_cs = summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak") 
stat.test_cs

stat.test_cs_pwc = stat.test_cs %>% tibble(".y." = "value", 
                                           "group1" = "Isotype", 
                                           "group2" = "anti-CD20",
                                           "p.adj" = round(p.value,4),
                                           "method" = "sidak") %>%
  add_x_position() %>% mutate(xmin = c(0.85,0.8,1.8,1.15), xmax = c(1.85,1.2,2.2,2.15)) %>%
  add_significance()


stat.test_cs_pwc



a2_cs = a2_cs %>% as.data.frame() 


ggboxplot(data_subset, x = "Depletion", y = "value", 
          color = "Tx",group = "Tx") +
  stat_anova_test(aes(group = Tx), method = "two_way",
                  type = 3)


tab_cs = data_subset %>%
  group_by(Tx, Depletion)%>%
  dplyr::summarise(n = n(), value = mean(value)) 

tab_cs$value = max(data_subset$value)/25
tab_cs$print = paste0("n = ", tab_cs$n)
tab_cs
max(data_subset$value)


cairo_pdf("cohort17TregofCD4_2wayAnova_Sidak_removelowAb_byhand.pdf",          width = 2 , height = 3)
#png("cohort17TregPer_2wayAnova_Sidak_removelowAb_byhand.png",    units = "in", width = 4.2 , height = 6, res = 400)
ggplot(data_subset, aes(x = Depletion, y = value)) +
  stat_summary(aes(color = Tx,group = Tx),
               fun.min = function(x) mean(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2,
               position = position_dodge(0.6),
               color = "black"
  ) +
  stat_summary(
    fun = mean, geom="col", 
    position = "dodge",
    aes(fill = Tx),
    color = "black", 
    size = 0.5, 
    width = 0.6
  ) +
  scale_fill_manual(values = c("grey", "#89CFF0")) +
  geom_point(aes(fill = Tx),
             position = position_jitterdodge(
               jitter.width = 0.12,
               jitter.height = 0,
               dodge.width = 0.6),
             color = "black",show.legend = F, size = 0.5
  ) + 
  geom_text(data = tab_cs, aes(label = print, group = Tx), 
            position = position_dodge(width = 0.6),
            color = "white", size = 1.5) +
  labs(#subtitle ="Two-way Anova, post-hoc: Sidak Test", 
    #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001\n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
    caption =  paste0("Two-Way Anova, post-hoc: Sidak Test \n[",rownames(a2_cs)[2],"] F(", a2_cs["Tx","Df"],",",a2_cs["Residuals","Df"],") = ",round(a2_cs["Tx","F value"], 4) ,", p = ",round(a2_cs["Tx","Pr(>F)"], 4) ," \n[",rownames(a2_cs)[3],"] F(", a2_cs["Depletion","Df"],",",a2_cs["Residuals","Df"],") = "
                      ,round(a2_cs["Depletion","F value"], 4),", p = ",  round(a2_cs["Depletion","Pr(>F)"], 4), " \n[",rownames(a2_cs)[4],
                      "] F(", a2_cs["Tx:Depletion","Df"],",",a2_cs["Residuals","Df"],") = ", round(a2_cs["Tx:Depletion","F value"], 4),", p = ",  round(a2_cs["Tx:Depletion","Pr(>F)"], 4), " "
    ),
    y = "Percent") +
  ggtitle(paste(colname)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),  breaks = seq(0, 6, 1))+
  coord_cartesian(clip = "off")+
  theme_prism(base_family = "Arial", base_fontface = "bold", base_size = 8) + 
  scale_colour_prism(palette = "black_and_white")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0,face = "plain", size = 4),
        plot.caption.position = "plot",
        plot.title.position = "plot"
        #plot.subtitle = element_text(size = 12),
  ) +
  #stat_anova_test(aes(group = Tx), method = "two_way",type = 3) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[1,],  y.position =8, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 700000,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[2,],  y.position = 5.6, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,#bracket.nudge.y = 200000,
    lineheight = 0, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[3,],  y.position = 5.6, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[4,],  y.position = 6.8, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 400000,
    lineheight = 0.1, size = 2
  ) 
dev.off()
# Plot "51 % CCR9+ Tregs (of CD4+ T Cells)"----
data_subset = counts
data_subset = counts[-which(counts$`Mouse ID` == "C3M5" | counts$`Mouse ID` == "C4M1"),]
data_subset$`Mouse ID`

type3 <- list(Depletion = contr.sum, Tx = contr.sum)
colnum = 51
colname = colnames(data_subset)[colnum]
colname
colnames(data_subset)[colnum] = "value"
m2_cs <- lm(value ~ Tx * Depletion , data = data_subset, contrasts = type3)
print(car::Anova(m2_cs, type = 3 )) # correct type 3

a2_cs = car::Anova(m2_cs, type = 3 )
#check assumptions
#check the equality of variances
print(car::leveneTest(value ~ Tx * Depletion, data = data_subset))

#check normality 
print(m2_cs$residuals %>% 
        shapiro_test())


with(data_subset, interaction.plot(Tx, Depletion, value, fun = mean,
                                   main = "Interaction Plot"))


# Post-hoc tests using Tukey's HSD

print(tukey_hsd(m2_cs, which = "Tx:Depletion"))

#posthoc <- glht(m2, linfct = mcp( Tx = "Tukey"))
#summary(posthoc)

EMMcs <- emmeans(m2_cs,  ~ Depletion * Tx)   # where treat has 2 levels
pairs(EMMcs, adjust = "sidak")   # adjustment is ignored - only 1 test per group
summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak")   # all are in one group now
stat.test_cs = summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak") 
stat.test_cs

stat.test_cs_pwc = stat.test_cs %>% tibble(".y." = "value", 
                                           "group1" = "Isotype", 
                                           "group2" = "anti-CD20",
                                           "p.adj" = round(p.value,4),
                                           "method" = "sidak") %>%
  add_x_position() %>% mutate(xmin = c(0.85,0.8,1.8,1.15), xmax = c(1.85,1.2,2.2,2.15)) %>%
  add_significance()


stat.test_cs_pwc



a2_cs = a2_cs %>% as.data.frame() 


ggboxplot(data_subset, x = "Depletion", y = "value", 
          color = "Tx",group = "Tx") +
  stat_anova_test(aes(group = Tx), method = "two_way",
                  type = 3)


tab_cs = data_subset %>%
  group_by(Tx, Depletion)%>%
  dplyr::summarise(n = n(), value = mean(value)) 

tab_cs$value = max(data_subset$value)/25
tab_cs$print = paste0("n = ", tab_cs$n)
tab_cs
max(data_subset$value)

colname
cairo_pdf("cohort17CCR9TregPerofCD4_2wayAnova_Sidak_removelowAb_byhand.pdf",          width = 2 , height = 3)
#png("cohort17TregPer_2wayAnova_Sidak_removelowAb_byhand.png",    units = "in", width = 4.2 , height = 6, res = 400)
ggplot(data_subset, aes(x = Depletion, y = value)) +
  stat_summary(aes(color = Tx,group = Tx),
               fun.min = function(x) mean(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2,
               position = position_dodge(0.6),
               color = "black"
  ) +
  stat_summary(
    fun = mean, geom="col", 
    position = "dodge",
    aes(fill = Tx),
    color = "black", 
    size = 0.5, 
    width = 0.6
  ) +
  scale_fill_manual(values = c("grey", "#89CFF0")) +
  geom_point(aes(fill = Tx),
             position = position_jitterdodge(
               jitter.width = 0.12,
               jitter.height = 0,
               dodge.width = 0.6),
             color = "black",show.legend = F, size = 0.5
  ) + 
  geom_text(data = tab_cs, aes(label = print, group = Tx), 
            position = position_dodge(width = 0.6),
            color = "white", size = 1.5) +
  labs(#subtitle ="Two-way Anova, post-hoc: Sidak Test", 
    #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001\n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
    caption =  paste0("Two-Way Anova, post-hoc: Sidak Test \n[",rownames(a2_cs)[2],"] F(", a2_cs["Tx","Df"],",",a2_cs["Residuals","Df"],") = ",round(a2_cs["Tx","F value"], 4) ,", p = ",round(a2_cs["Tx","Pr(>F)"], 4) ," \n[",rownames(a2_cs)[3],"] F(", a2_cs["Depletion","Df"],",",a2_cs["Residuals","Df"],") = "
                      ,round(a2_cs["Depletion","F value"], 4),", p = ",  round(a2_cs["Depletion","Pr(>F)"], 4), " \n[",rownames(a2_cs)[4],
                      "] F(", a2_cs["Tx:Depletion","Df"],",",a2_cs["Residuals","Df"],") = ", round(a2_cs["Tx:Depletion","F value"], 4),", p = ",  round(a2_cs["Tx:Depletion","Pr(>F)"], 4), " *"
    ),
    y = "Percent") +
  ggtitle(gsub("[:(:]", "\n(", colname)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),  breaks = seq(0, 2.4, 0.4))+
  coord_cartesian(clip = "off")+
  theme_prism(base_family = "Arial", base_fontface = "bold", base_size = 8) + 
  scale_colour_prism(palette = "black_and_white")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0,face = "plain", size = 4),
        plot.caption.position = "plot",
        plot.title.position = "panel"
        #plot.subtitle = element_text(size = 12),
  ) +
  #stat_anova_test(aes(group = Tx), method = "two_way",type = 3) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[1,],  y.position =3, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 700000,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[2,],  y.position = 2.2, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,#bracket.nudge.y = 200000,
    lineheight = 0, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[3,],  y.position = 1.6, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[4,],  y.position = 2.5, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 400000,
    lineheight = 0.1, size = 2
  ) 
dev.off()
#Plot tregs counts ----
names(counts)
data_subset = counts
data_subset = counts[-which(counts$`Mouse ID` == "C3M5" | counts$`Mouse ID` == "C4M1"),]
data_subset$`Mouse ID`

type3 <- list(Depletion = contr.sum, Tx = contr.sum)
colnum = 38
colname = colnames(data_subset)[colnum]
colname
colnames(data_subset)[colnum] = "value"
m2_cs <- lm(value ~ Tx * Depletion , data = data_subset, contrasts = type3)
print(Anova(m2_cs, type = 3 )) # correct type 3

a2_cs = Anova(m2_cs, type = 3 )
#check assumptions
#check the equality of variances
print(car::leveneTest(value ~ Tx * Depletion, data = data_subset))

#check normality 
print(m2_cs$residuals %>% 
        shapiro_test())


with(data_subset, interaction.plot(Tx, Depletion, value, fun = mean,
                                   main = "Interaction Plot"))


# Post-hoc tests using Tukey's HSD

print(tukey_hsd(m2_cs, which = "Tx:Depletion"))

#posthoc <- glht(m2, linfct = mcp( Tx = "Tukey"))
#summary(posthoc)

EMMcs <- emmeans(m2_cs,  ~ Depletion * Tx)   # where treat has 2 levels
pairs(EMMcs, adjust = "sidak")   # adjustment is ignored - only 1 test per group
summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak")   # all are in one group now
stat.test_cs = summary(pairs(EMMcs)[c(1:2,5:6)], by = NULL, adjust = "sidak") 
stat.test_cs

stat.test_cs_pwc = stat.test_cs %>% tibble(".y." = "value", 
                                           "group1" = "Isotype", 
                                           "group2" = "anti-CD20",
                                           "p.adj" = round(p.value,4),
                                           "method" = "sidak") %>%
  add_x_position() %>% mutate(xmin = c(0.85,0.8,1.8,1.15), xmax = c(1.85,1.2,2.2,2.15)) %>%
  add_significance()


stat.test_cs_pwc



a2_cs = a2_cs %>% as.data.frame() 


ggboxplot(data_subset, x = "Depletion", y = "value", 
          color = "Tx",group = "Tx") +
  stat_anova_test(aes(group = Tx), method = "two_way",
                  type = 3)


tab_cs = data_subset %>%
  group_by(Tx, Depletion)%>%
  dplyr::summarise(n = n(), value = mean(value)) 

tab_cs$value = max(data_subset$value)/25
tab_cs$print = paste0("n = ", tab_cs$n)
tab_cs
max(data_subset$value)


cairo_pdf("cohort17TregCounts_2wayAnova_Sidak_removelowAb_byhand.pdf",          width = 2 , height = 3)
#png("cohort17TregPer_2wayAnova_Sidak_removelowAb_byhand.png",    units = "in", width = 4.2 , height = 6, res = 400)
ggplot(data_subset, aes(x = Depletion, y = value)) +
  stat_summary(aes(color = Tx,group = Tx),
               fun.min = function(x) mean(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2,
               position = position_dodge(0.6),
               color = "black"
  ) +
  stat_summary(
    fun = mean, geom="col", 
    position = "dodge",
    aes(fill = Tx),
    color = "black", 
    size = 0.5, 
    width = 0.6
  ) +
  scale_fill_manual(values = c("grey", "#89CFF0")) +
  geom_point(aes(fill = Tx),
             position = position_jitterdodge(
               jitter.width = 0.12,
               jitter.height = 0,
               dodge.width = 0.6),
             color = "black",show.legend = F, size = 0.5
  ) + 
  geom_text(data = tab_cs, aes(label = print, group = Tx), 
            position = position_dodge(width = 0.6),
            color = "white", size = 1.5) +
  labs(#subtitle ="Two-way Anova, post-hoc: Sidak Test", 
    #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001\n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
    caption =  paste0("Two-Way Anova, post-hoc: Sidak Test \n[",rownames(a2_cs)[2],"] F(", a2_cs["Tx","Df"],",",a2_cs["Residuals","Df"],") = ",round(a2_cs["Tx","F value"], 4) ,", p = ",round(a2_cs["Tx","Pr(>F)"], 4) ," *\n[",rownames(a2_cs)[3],"] F(", a2_cs["Depletion","Df"],",",a2_cs["Residuals","Df"],") = "
                      ,round(a2_cs["Depletion","F value"], 4),", p = ",  round(a2_cs["Depletion","Pr(>F)"], 4), " \n[",rownames(a2_cs)[4],
                      "] F(", a2_cs["Tx:Depletion","Df"],",",a2_cs["Residuals","Df"],") = ", round(a2_cs["Tx:Depletion","F value"], 4),", p = ",  round(a2_cs["Tx:Depletion","Pr(>F)"], 4), " ."
    ),
    y = "Counts") +
  ggtitle(paste(colname)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),  breaks = seq(0, 6e4, 2e4))+
  coord_cartesian(clip = "off")+
  theme_prism(base_family = "Arial", base_fontface = "bold", base_size = 8) + 
  scale_colour_prism(palette = "black_and_white")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0,face = "plain", size = 4),
        plot.caption.position = "plot",
        plot.title.position = "panel"
        #plot.subtitle = element_text(size = 12),
  ) +
  #stat_anova_test(aes(group = Tx), method = "two_way",type = 3) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[1,],  y.position =5.5e4, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 700000,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[2,],  y.position = 4.2e4, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,#bracket.nudge.y = 200000,
    lineheight = 0, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[3,],  y.position = 6e4, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3,
    lineheight = 0.1, size = 2
  ) + 
  stat_pvalue_manual(
    stat.test_cs_pwc[4,],  y.position = 7.2e4, step.group.by = "contrast",
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 400000,
    lineheight = 0.1, size = 2
  ) 
dev.off()

Max_67_data = read_csv("Severity_ClinicalOutcomes_Cohort17.csv")
colnames(Max_67_data)[1] ="Mouse ID"
Max_67_data$Severity = factor(Max_67_data$Severity, levels = c("No Reactivity", "Mild", "Severe"))
Max_67_data = merge(Max_67_data[,c(1:4,9)],counts,by="Mouse ID")

data_subset = Max_67_data
colnum = 42
colname = colnames(data_subset)[colnum]
colname
colnames(data_subset)[colnum] = "value"

stat.test <- as.data.frame(data_subset) %>% 
  wilcox_test( value ~ Condition) %>%
  add_significance() %>% add_xy_position(x = "Condition")
stat.test

ggplot(data_subset, aes(x = Severity, y = value)) +
  facet_grid(~Tx) + 
  stat_summary(aes(color = Depletion,group = Depletion),
               fun.min = function(x) mean(x) , 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2,
               position = position_dodge(0.6, preserve = "single"),
               color = "black"
  ) +
  stat_summary(
    fun = mean, geom="col", 
    #position = "dodge",
    position = position_dodge(preserve = "single"),
    aes(fill = Depletion),
    color = "black", 
    size = 0.5, 
    width = 0.6
  ) +
  scale_fill_manual(values = c("grey", "darkgreen")) +
  geom_point(aes(fill = Depletion, group = Depletion),
             position = position_jitterdodge(
               jitter.width = 0.1,
               jitter.height = 0,
               dodge.width = 0.6),
             color = "black",show.legend = F
  ) +
  #geom_hline(yintercept = 0, size = 1,linetype = 1, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  coord_cartesian(clip = "off")+
  labs(subtitle =" ", 
       #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001 * \n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
       caption = " ",
       y = colname) +
  ggtitle(colname) + 
  theme_prism(base_family = "Arial", base_fontface = "bold") + 
  scale_colour_prism(palette = "black_and_white")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0,face = "plain", size = 8),
        plot.subtitle = element_text(size = 12)) 
dev.off()

myplots_list =  list()
j = 1
for (i in c(29:62)){
  data_subset = Max_67_data
  colnum = i
  colname = colnames(data_subset)[colnum]
  colname
  colnames(data_subset)[colnum] = "value"
  
  
  plot = ggplot(data_subset, aes(x = Severity, y = value)) +
    facet_grid(~Tx) + 
    stat_summary(aes(color = Depletion,group = Depletion),
                 fun.min = function(x) mean(x) , 
                 fun.max = function(x) mean(x) + sd(x), 
                 geom = "errorbar", width = 0.2,
                 position = position_dodge(0.6, preserve = "single"),
                 color = "black"
    ) +
    stat_summary(
      fun = mean, geom="col", 
      #position = "dodge",
      position = position_dodge(preserve = "single"),
      aes(fill = Depletion),
      color = "black", 
      size = 0.5, 
      width = 0.6
    ) +
    scale_fill_manual(values = c("grey", "darkgreen")) +
    geom_point(aes(fill = Depletion),
               position = position_jitterdodge(
                 jitter.width = 0.1,
                 jitter.height = 0,
                 dodge.width = 0.5),
               color = "black",show.legend = F
    ) +
    #geom_hline(yintercept = 0, size = 1,linetype = 1, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    coord_cartesian(clip = "off")+
    labs(subtitle =" ", 
         #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001 * \n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
         caption = " ",
         y = colname) +
    ggtitle(colname) + 
    theme_prism(base_family = "Arial", base_fontface = "bold") + 
    scale_colour_prism(palette = "black_and_white")+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          legend.position = "top", legend.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0,face = "plain", size = 8),
          plot.subtitle = element_text(size = 12)) 
  
  myplots_list[[j]] = plot
  j = j + 1
}

getwd()
png(file = "Severity_analysis_TregPanel_byhand.png",
    units = "in",width = 40, height = 30, res = 400)
ggarrange(plotlist = myplots_list, ncol=5, nrow=7, 
          common.legend = T, legend = "top")
dev.off()

png(file = "Severity_analysis_TregPanel_byhand_Treg.png",
    units = "in",width = 20, height = 4, res = 400)
ggarrange(plotlist = myplots_list[c(2,14,25)], ncol=3, nrow=1, 
          common.legend = T, legend = "top")
dev.off()

myplots_list =  list()
j = 1
for (i in c(29:50)){
  data_subset = Max_67_data
  colnum = i
  colname = colnames(data_subset)[colnum]
  colname
  colnames(data_subset)[colnum] = "value"
  
  stat.test <- data_subset %>%
    wilcox_test(value ~ Severity, paired = F) %>%
    add_significance()  %>% add_xy_position(x = "Condition")
  stat.test
  
  
  plot =  ggplot(data_subset, aes(x = Severity, y = value)) +
    #facet_grid(~Tx) + 
    stat_summary(aes(color = Tx, group = Tx),
                 fun.min = function(x) mean(x) , 
                 fun.max = function(x) mean(x) + sd(x), 
                 geom = "errorbar", width = 0.2,
                 position = position_dodge(0.6),
                 color = "black"
    ) +
    stat_summary(
      fun = mean, geom="col", 
      position = "dodge",
      aes(fill = Tx),
      color = "black", 
      size = 0.5, 
      width = 0.6
    ) +
    scale_fill_manual(values = c("grey", "#89CFF0")) +
    geom_point(aes(fill = Tx, group = Tx),
               position = position_jitterdodge(
                 jitter.width = 0.12,
                 jitter.height = 0,
                 dodge.width = 0.6),
               color = "black",show.legend = F
    ) + 
    #geom_hline(yintercept = 0, size = 1,linetype = 1, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    coord_cartesian(clip = "off")+
    labs(subtitle =" ", 
         #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001 * \n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
         caption = " ",
         y = colname) +
    ggtitle(colname) + 
    theme_prism(base_family = "Arial", base_fontface = "bold") + 
    scale_colour_prism(palette = "black_and_white")+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          legend.position = "top", legend.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0,face = "plain", size = 8),
          plot.subtitle = element_text(size = 12)) 
  
  myplots_list[[j]] = plot
  j = j + 1
}

png(file = "Severity_analysis_DCPanel_byTx_noDep_byhand.png",
    units = "in",width = 40, height = 20, res = 400)
ggarrange(plotlist = myplots_list, ncol=5, nrow=5, 
          common.legend = T, legend = "top")
dev.off()

ggplot(data_subset, aes(x = Severity, y = value)) +
  facet_grid(~Tx) + 
  stat_summary(aes(color = Severity,group = Severity),
               fun.min = function(x) mean(x) , 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2,
               position = position_dodge(0.6),
               color = "black"
  ) +
  stat_summary(
    fun = mean, geom="col", 
    position = "dodge",
    aes(fill = Severity),
    color = "black", 
    size = 0.5, 
    width = 0.6
  ) +
  #scale_fill_manual(values = c("grey", "#89CFF0")) +
  geom_point(aes(fill = Severity),
             position = position_jitterdodge(
               jitter.width = 0.12,
               jitter.height = 0,
               dodge.width = 0.6),
             color = "black",show.legend = F
  ) + 
  
  #geom_hline(yintercept = 0, size = 1,linetype = 1, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  coord_cartesian(clip = "off")+
  labs(subtitle =" ", 
       #caption = "Two-Way Anova \n [Tx] F(1,32) = 25.08, p<0.0001 * \n [Depletion] F(1,32) = 0.8316, p = 0.3686 ns \n [Tx:Depletion] F(1,32) = 0.4076, p = 0.5277 ns",
       caption = " ",
       y = colname) +
  ggtitle(colname) + 
  theme_prism(base_family = "Arial", base_fontface = "bold") + 
  scale_colour_prism(palette = "black_and_white")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0,face = "plain", size = 8),
        plot.subtitle = element_text(size = 12)) 
dev.off()

cor(Max_67_data[,c(2:4,20:49)][,unlist(lapply(Max_67_data[,c(2:4,20:49)], is.numeric))])

#load psych package
library(psych)

#create pairs plot
#png(file = "Severity_analysis_DCPanel_Correlation.png",units = "in",width = 40, height = 20, res = 400)
pairs.panels(Max_67_data[,c(2:4,20:49)][,unlist(lapply(Max_67_data[,c(2:4,20:49)], is.numeric))], lm = T)
dev.off()
x = Max_67_data[,c(2:4,53:62)]
c()
x = Max_67_data[,c(2:4,30:33,42:45,54:56)]
colnames(x)
png(file = "Severity_analysis_TregPanel_Tregs_Correlation_byhand.png",
    units = "in",width = 30, height = 25, res = 400)
corPlot(r = x[,unlist(lapply(x, is.numeric))], 
        las=2,
        main = "Correlation of SILP T cells and clinical outcomes",
        colors = T, show.legend = F, scale = F,
        cex.axis = 1.5, cex = 1,diag = T,
        #cex.main = 20,
        min.length = 40,
        MAR = 20,
        gr = colorRampPalette(c( "#2171B5","white", "#B52127")) )
dev.off()

