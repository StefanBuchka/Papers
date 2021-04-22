# On the optimistic performance evaluationof newly introduced bioinformatic methods
#The code to generate the quantitative results from the Commentary by Buchka et al.
#This code: Generating the figure. 
#This code must be run after running analysis_1_level_of_paper.R and analysis_2_level_of_sutbstudy.R
#Latest Date: 31.03.2021

rm(list=ls())

library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)
library(cowplot)

load(file="comparisons_papers.RData")
load(file="comparisons_substudies.RData")

######################### Creation of figure

### Adjustment of both datasets (substudy level and paper level) to be able to create the figure
# 1.) Paper level
comparisons_papers$introducing$study_type <- "introducing"
comparisons_papers$neutral$study_type     <- "neutral"

comparisons_papers$introducing$unit <- "(paper as unit)"
comparisons_papers$neutral$unit     <- "(paper as unit)"

figure_paper <- rbind(comparisons_papers$introducing,comparisons_papers$neutral)
names(figure_paper)[which(names(figure_paper) == "n_papers")] <- "n_comparisons"

# 2.) Substudy level
comparisons_substudies$introducing$study_type <- "introducing"
comparisons_substudies$neutral$study_type     <- "neutral"

comparisons_substudies$introducing$unit <- "(substudy as unit)"
comparisons_substudies$neutral$unit     <- "(substudy as unit)"

figure_substudies <- rbind(comparisons_substudies$introducing,comparisons_substudies$neutral)
names(figure_substudies)[which(names(figure_substudies) == "n_substudies")] <- "n_comparisons"

# 3.) Combination of the two data sets
figure        <- rbind(figure_paper,figure_substudies)
figure$colour <- paste0(figure$study_type," ",figure$unit)

figure$perc_new_better <- 100*figure$perc_new_better

#Creation of the plot itself

plot_base <- figure %>%
  group_by(desc(study_type), desc(unit)) %>%
  arrange(desc(perc_new_better),.by_group =T) %>%
  ungroup() %>% 
  mutate(order = 1:dim(figure)[1]) %>% 
  mutate(comparisons = fct_reorder(comparisons,order,min)) %>%
  ggplot(.,aes(comparisons,perc_new_better)) +
  geom_line(aes(group = comparisons))


position <- position_dodge(.5)

plot_final <- plot_base +
  geom_point(aes(comparisons,perc_new_better,
                 colour = colour,
                 size = n_comparisons),
             position = position,
             alpha = 0.6) + 
  ggtitle("Percentage calculated over studies") +
  scale_color_manual(values = c("#F8766D","darkred","#619CFF","blue4")) +
  labs(size ="Number of comparisons", colour="Type of comparison") +
  theme(legend.text = element_text(size=15, face="bold"), legend.title = element_text(size=18, face="bold")) +
  theme(legend.text = element_text(size = 15)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) + 
  guides(colour = guide_legend(override.aes = list(size=4))) + 
  xlab("Pairs") + 
  ylab("New method better (frequency, %)" ) +
  theme(legend.position="right") +
  theme(plot.title = element_text(color = "salmon", size = 20, face = "bold"),
        axis.title.x = element_text(color = "salmon", size = 20, face = "bold", vjust = -1),
        axis.title.y = element_text(color = "salmon", size = 20, face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.2, face = "bold", size = 15)) +
  theme(axis.text.y = element_text(size = 15, face = "bold")) +
  scale_size(breaks = c(1,5,10,15,seq(20, 70, 10)),range = c(3, 16.1))



jpeg("figure_final.jpeg", res = 1000, units = "in", width = 18, height = 10)

ggdraw() +
  draw_plot(plot_final, x = 0, y = 0, width = 0.99, height = 0.98) 

dev.off()


  