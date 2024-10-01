#locate the data files collection
setwd("E:/206/ESKAPE/Actual_data")
#Find out the MAle Female Age difference
library(dplyr)
normality_test_M_vs_F <- file1 %>%
  group_by(Sex) %>%
  summarize(p_value = shapiro.test(Age)$p.value)
  #non-normal data
library(rstatix)
wilcox_test_M_vs_F <- file1 %>% 
     wilcox_test(Age ~ Sex)
  #plot the data
ggplot(file1, aes(Sex,Age))+
  geom_boxplot()+
  theme_minimal()+
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c("M", "F")), label.y = 0.75, label.x = 1.05, label.size = 2)  # Add t-test and display p-value

population <- file1%>% 
  group_by(Sex)%>%
  summarise(n = n())
  
  
  
#Find out the biofilm percentage
      #load the metadata
file1 <-read.csv("metadata.csv")
biofilmogram <- file1 %>%
  select("Bacteria", "Biofilm", "Class") %>%
  group_by(Bacteria, Biofilm) %>%  
  summarise(n = n()) %>%          
  mutate(percentage = n/sum(n)*100) 

# Calculate the overall percentage and add it as a new column  
overall_summary <- file1 %>%
  group_by(Biofilm) %>%           # Group only by Biofilm for overall summary
  summarise(n = n()) %>%          
  mutate(percentage = n/sum(n)*100)%>%
  mutate(Bacteria = "Overall") 
#calculate 
class_biofilm <- file1 %>%
  group_by(Class, Biofilm) %>%
  summarise(n = n()) %>%          
  mutate(percentage = n/sum(n)*100)%>%
  rename(Bacteria = Class)
  
#Join both tables
combined_data_1 <- bind_rows(biofilmogram, overall_summary, class_biofilm)

#Arranging the data
combined_data_1$Biofilm <- factor(combined_data_1$Biofilm, 
                                  levels = c("no-biofilm", "weak", "medium", "Strong"))
combined_data_1$Bacteria <- factor(combined_data_1$Bacteria, 
                                  levels = c("A. baumannii", "K. pneumoniae", "P. aeruginosa", "S. aureus", "E. feacium", "gram-negative", "gram-positive", "Overall"))
# plot the figure: Biofilm Percentage
Figure_4 <- ggplot(combined_data_1, aes(Bacteria, percentage, fill = Biofilm))+
  geom_bar(stat = "identity", color = "black", position = "dodge", linewidth = 0.8) +
  scale_fill_brewer()+
  labs(x = "Bacteria", y = "Percentage") +
  theme_test()+
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 13),
    strip.placement = "outside", strip.clip = "off",
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", margin = margin(b = 5), size = 13),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 13), 
    axis.title = element_text(size = 15),
    panel.border = element_rect(linewidth = 1.2),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),  # Corrected font setting
    legend.background = element_rect(color = "black", linewidth = 0.8),
    legend.position = c(0, 1),     # Position: top-right corner (1, 1)
    legend.justification = c(0, 1),  # Justification: top-right corner (1, 1)
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Corrected border setting
  )


#p value annotatio
custom_annotation <- function(pvalue) {
  ifelse(pvalue < 0.001, "p < 0.001", paste0("p = ", formatC(pvalue, format = "f", digits = 3)))
}
# Biofilm_variation in gram_negative and gram_positive isolates
         #Need to check two assumption for one way ANOVa
        #first assumption: Normality
normality_test_plot1 <- file1 %>%
  group_by(Bacteria) %>%
  filter(Class == "gram-negative")%>%
  summarize(p_value = shapiro.test(meanOD.ODcut)$p.value)
       #p values less thann 0.05 in each group, deviate significantly from a normal distribution.
       # normality no assumed, so we will need to perform Kruskall-wallis test
#Create a file only for gram negative Bacteria
library(dplyr)
file2 <- file1 %>%
  filter(Class == "gram-negative")
library(ggpubr)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(rstatix)
library(RColorBrewer)

# Perform Kruskal-Wallis test
kruskal_test_result_1 <- kruskal.test(meanOD.ODcut ~ Bacteria, data = file2)

# Perform Dunn's test (with Bonferroni correction as an example)
dunn_test_results_1 <- file2 %>% 
  dunn_test(meanOD.ODcut ~ Bacteria, p.adjust.method = "bonferroni") %>%
  mutate(p.adj.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Bacteria")

# Print Dunn's test results to inspect
print(dunn_test_results_1)
# Convert p.adj from scientific notation to decimal
dunn_test_results_1$p.adj <- as.numeric(dunn_test_results_1$p.adj) 

# Format p.adj values (now in decimal) and create a new column
dunn_test_results_1$p.adj.format <- formatC(dunn_test_results_1$p.adj, format = "f", digits = 6)
pvalue_p1 <- dunn_test_results_1$p.adj.format
# Create the plot
p1 <- ggplot(file2, aes(x = Bacteria, y = meanOD.ODcut)) +
      geom_point(color = "#9B9B9B",shape = 95, size = 30, alpha = .60)+
      theme_test() +
      #stat_pvalue_manual(dunn_test_results_1, label = "p.adj.format", hide.ns = FALSE) + # Display Dunn's test results
  geom_signif(comparisons = list(c("A. baumannii","K. pneumoniae"),c("A. baumannii", "P. aeruginosa"), c("P. aeruginosa", "K. pneumoniae")), 
              
              step_increase = 0.2,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_p1, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_p1),
              y_position = max(file2$meanOD.ODcut) * 1.08) +
      labs(y = "RBF")+
  theme(axis.text.x = element_text(angle = 12, size = 13, face = "italic"), legend.position = "none",
        panel.border = element_rect(linewidth = 1),
        axis.text.y = element_text(size = 13), 
        axis.title = element_text(size = 15))





#Gram-positive file
file3 <- file1 %>%
  filter(Class == "gram-positive")
#Gram-positive Biofilm Difference
       #Only two bacteria so t-test will work
       #Need to check normality
normality_test_plot2 <- file3 %>%
  group_by(Bacteria) %>%
  summarize(p_value = shapiro.test(meanOD.ODcut)$p.value)
       #non-normal data
        #need to do wilcox.t-test
       # Create the plot
wilcox_test_1 <- file3 %>% 
  wilcox_test(meanOD.ODcut ~ Bacteria) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Bacteria")
pvalue_p2 <- wilcox_test_1$p

p2 <- ggplot(file3, aes(x = Bacteria, y = meanOD.ODcut))+
  geom_point(color = "#9B9B9B",shape = 95, size = 30, alpha = .60)+
      theme_test()+
  geom_signif(comparisons = list(c("E. feacium", "S. aureus")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_p2, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_p2),
              y_position = max(file3$meanOD.ODcut) * 1.05) +
      #stat_pvalue_manual(wilcox_test_1, label = "p.signif", hide.ns = FALSE) + 
      theme(axis.title.y = element_blank(), axis.text.x = element_text(size = 13, face = "italic"), legend.position = "none",
            panel.border = element_rect(linewidth = 1),
            axis.text = element_text(size = 13),
            axis.title.x = element_text(size = 15))


# Gram-posistive vs Gram- negative
      #normality-check
 normality_test_plot3 <- file1 %>%
   group_by(Class) %>%
   summarize(p_value = shapiro.test(meanOD.ODcut)$p.value)  
 #non-normal data
 #need to do wilcox.t-test
 # Create the plot
 wilcox_test_2 <- file1 %>% 
   wilcox_test(meanOD.ODcut ~ Class) %>%
   mutate(p.signif = case_when(
     p <= 0.001 ~ "***",
     p <= 0.01  ~ "**",
     p <= 0.05  ~ "*",
     TRUE           ~ "NS"
   ))%>%
   add_xy_position(x = "Class")
 pvalue_p3 <-  wilcox_test_2$p
p3 <- ggplot(file1, aes(x = Class, y = meanOD.ODcut)) +
  geom_point(color = "#9B9B9B",shape = 95, size = 30, alpha = .60)+
  theme_test()+
  geom_signif(comparisons = list(c("gram-negative", "gram-positive")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_p3, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_p3),
              y_position = max(file1$meanOD.ODcut) * 1.05) +
      #stat_pvalue_manual(wilcox_test_2,, vjust = 1.5, label = "p.signif", hide.ns = FALSE)+
  theme(axis.title.y = element_blank(), legend.position = "none",
        panel.border = element_rect(linewidth = 1),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 15))
      #labs(y = "meanOD/ODcut") 

 

#Carbapenemase Producer vs non-producer meaOD/ ODcut
 
 normality_test_plot4 <- file2 %>%
   group_by(Carbapenemase) %>%
   summarize(p_value = shapiro.test(meanOD.ODcut)$p.value)  
 #non-normal data
 #need to do wilcox.t-test
 # Create the plot
 wilcox_test_3 <- file2 %>% 
   wilcox_test(meanOD.ODcut ~ Carbapenemase) %>%
   mutate(p.signif = case_when(
     p <= 0.001 ~ "***",
     p <= 0.01  ~ "**",
     p <= 0.05  ~ "*",
     TRUE           ~ "NS"
   ))%>%
   add_xy_position(x = "Carbapenemase")
 pvalue_p4 <- wilcox_test_3$p
p4 <- ggplot(file2, aes(x = Carbapenemase, y = meanOD.ODcut)) +
  geom_point(color = "#9B9B9B",shape = 95, size = 30, alpha = .60)+
  theme_test()+
  geom_signif(comparisons = list(c("NP", "P")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_p4, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_p4),
              y_position = max(file1$meanOD.ODcut) * 1.05) +
      #stat_pvalue_manual(wilcox_test_3, label = "p.signif", hide.ns = FALSE) + 
      labs(y = "RBF")+
  theme(legend.position = "none",
        panel.border = element_rect(linewidth = 1),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))
 
 #Find out the sample source percentage
source <- file1 %>%
  group_by(Source) %>%           
  summarise(n = n()) %>%          
  mutate(percentage = n/sum(n)*100)
# Filter the predominant source
file4 <- file1 %>%
  filter(Source %in% c("W/S", "U", "T/A", "P", "S"))
   #normality test
normality_test_plot5 <- file4 %>%
  group_by(Source) %>%
  summarize(p_value = shapiro.test(meanOD.ODcut)$p.value)  
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(rstatix)
library(RColorBrewer)

# Perform Kruskal-Wallis test
kruskal_test_result_2 <- kruskal.test(meanOD.ODcut ~ Source, data = file4)

# Perform Dunn's test (with Bonferroni correction as an example)
dunn_test_results_2 <- file4 %>% 
  dunn_test(meanOD.ODcut ~ Source, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Source")

# Print Dunn's test results to inspect
print(dunn_test_results_2)
# significant difference only between T/A and U sample source
# So again evaluate by wilcox_test
#Source only T/A and U
file4 <- file1 %>%
  filter(Source %in% c("U", "T/A"))
wilcox_test_4 <- file4 %>% 
  wilcox_test(meanOD.ODcut ~ Source) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Source")
pvalue_p5 <- wilcox_test_4$p
p5 <- ggplot(file4, aes(x = Source, y = meanOD.ODcut)) +
  geom_point(color = "#9B9B9B",shape = 95, size = 30, alpha = .60)+
  theme_test()+
  geom_signif(comparisons = list(c("T/A", "U")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_p5, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_p5),
              y_position = max(file4$meanOD.ODcut) * 1.05) +
      #stat_pvalue_manual(wilcox_test_4, label = "p.signif", hide.ns = FALSE) + 
  theme(axis.title.y = element_blank(), legend.position = "none",
        panel.border = element_rect(linewidth = 1),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 15))
# MDR and non-MDR meanOD differece
normality_test_plot6 <- file1 %>%
  group_by(MDR) %>%
  summarize(p_value = shapiro.test(meanOD.ODcut)$p.value)  
wilcox_test_5 <- file1 %>% 
  wilcox_test(meanOD.ODcut ~ MDR) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "MDR")
pvalue_p6 <- wilcox_test_5$p
p6 <- ggplot(file1, aes(x = MDR, y = meanOD.ODcut)) +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.3) + 
  geom_point(color = "#9B9B9B",shape = 95, size = 30, alpha = .60)+
  theme_test()+
  geom_signif(comparisons = list(c("MDR", "non-MDR")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_p6, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_p6),
              y_position = max(file1$meanOD.ODcut) * 1.05) +
  #stat_pvalue_manual(wilcox_test_5, label = "p.signif", vjust = 1.5, hide.ns = FALSE) + 
  theme(axis.title.y = element_blank(), legend.position = "none",
        panel.border = element_rect(linewidth = 1),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 15))
library(patchwork)
Figure_5 <- (p1 | p2 | p3 ) / (p4| p5 | p6) 
Figure_5 <- Figure_5 + plot_annotation(tag_levels = 'a')+
  theme(plot.tag = element_text(face = "bold")) 
Figure_5




#Antibiogram: Figure 1
file5 <- read.csv("resistance_percentage.csv")
file5$Bacteria <- factor(file5$Bacteria, levels = c("A. baumannii", "P. aeruginosa", "K. pneumoniae", "S. aureus", "E. feacium"))
file5$Resistance <- factor(file5$Resistance, levels = c("S", "I", "R"))
file5$Antibiotics <- factor(file5$Antibiotics, levels = c("IPM", "MEM", "CN","AMK","CAZ", "FEP", "TZP", "CT", "CIP", "LEV",
                                                          "TE", "DO","ATM", "SXT", "NFT", "VA", "C", "LIN"))
library(ggplot2)
Figure_1 <-ggplot(file5, aes(x = Bacteria, y = Percentage, fill = Resistance)) +
  geom_bar(stat = "identity", position= "stack", alpha=0.75, width = 0.7, color = "black") +
  theme_minimal() +
  facet_wrap(~Antibiotics, nrow = 6, ncol = 6, labeller = label_wrap_gen(width = 2)) +
  xlab("Bacteria") +
  ylab("Percentage") +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 13),
    strip.placement = "outside", strip.clip = "off",
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", margin = margin(b = 5), size = 13),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 13), 
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),  # Corrected font setting
    legend.background = element_rect(color = "black", linewidth = 1),
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Corrected border setting
  )+ scale_fill_brewer()
        #scale_fill_manual(values = c("white", "grey", "black"))






#cMARI in Gram-positive isolates
   #normaility test
ormality_test_plot7 <- file3 %>%
  group_by(Bacteria) %>%
  summarize(p_value = shapiro.test(cMARI)$p.value)
    #non-normal data
wilcox_test_6 <- file3 %>% 
  wilcox_test(cMARI ~ Bacteria) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Bacteria")
pvalue_q1 <- wilcox_test_6$p
q1 <- ggplot(file3, aes(x = Bacteria, y = cMARI)) +
  geom_violin(aes(fill = Bacteria)) +
  scale_fill_brewer() +
  theme_test() +
  geom_signif(comparisons = list(c("E. feacium", "S. aureus")), 
              step_increase = 0.1,
              vjust = 1.5,
              annotations =  custom_annotation(pvalue_q1),
              y_position = max(file3$cMARI) * 1.06) +
  #stat_pvalue_manual(wilcox_test_6, label = "p.signif", hide.ns = FALSE) + 
  labs(y = "cMARI", x = "gram-positive")+
  theme(axis.text.x = element_text(face = "italic", hjust = 1, margin = margin(b = 5), size = 12),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 15),
        panel.border = element_rect(linewidth = 1.2),
        legend.position = "none")

#cMARI in Gram-negative isolates
normality_test_plot8 <- file2 %>%
  group_by(Bacteria) %>%
  summarize(p_value = shapiro.test(cMARI)$p.value)  
#non-normal data
# Perform Kruskal-Wallis test
kruskal_test_result_3 <- kruskal.test(cMARI ~ Bacteria, data = file2)
# Perform Dunn's test (with Bonferroni correction as an example)
dunn_test_results_3 <- file2 %>% 
  dunn_test(cMARI ~ Bacteria, p.adjust.method = "bonferroni") %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Bacteria")
  #Plot the figure
pvalue_q2 <- dunn_test_results_3$p.adj

q2 <- ggplot(file2, aes(x = Bacteria, y = cMARI)) +
  geom_violin(aes(fill = Bacteria))+
  scale_fill_brewer() +
  theme_test() +
  geom_signif(comparisons = list(c("A. baumannii","K. pneumoniae"),c("A. baumannii", "P. aeruginosa"), c("P. aeruginosa", "K. pneumoniae")), 
              step_increase = 0.1,
              vjust = 1.5,
              annotations =  custom_annotation(pvalue_q2),
              y_position = max(file2$cMARI) * 1.06) +
  #stat_pvalue_manual(dunn_test_results_3, label = "p.signif", hide.ns = FALSE,  y.position = c(1.05, 1.10, 1.15)) + 
  labs(y = "cMARI", x = "gram-negative")+
  theme(axis.text.x = element_text(face = "italic", hjust = 1, margin = margin(b = 5), size = 12),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 15),
        panel.border = element_rect(linewidth = 1.2),
        legend.position = "none")
  #combine the plot
library(patchwork)
Figure_2 <- q1 | q2 
Figure_2 <- Figure_2 + plot_annotation(tag_levels = 'a')+
  theme(plot.tag = element_text(face = "bold")) 
Figure_2
  
 




#gram-positive specific resistance

library(tidyr)
library(ggplot2)
library(ggpubr)
file6 <- read.csv("specific_resistance_gp.csv")%>%
  gather(Antibiotics,Susceptibility,4:5)

ggplot(file6, aes(Susceptibility,RBF))+
  geom_boxplot()+
  facet_wrap(~Antibiotics)+
  theme_minimal()+
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c("S", "NS")), label.y = 0.75, label.x = 1.05, label.size = 2) +  # Add t-test and display p-value
  labs(title = "Relation between Susceptibility and RBF by Antibiotics") 

#gram-negative specific resistance
  
file7 <- read.csv("specific_resistance_gn.csv")%>%
  na.omit(file7)%>%
  gather(Antibiotics,Susceptibility,7:16)%>%
  mutate(Antibiotics = factor(Antibiotics, levels = c("IPM", "MEM", "CN", "AK", "CAZ", "FEP", "CIP", "LEV", "TZP", "CT")))%>%
  mutate(Antibiotic_Class = case_when(
    Antibiotics %in% c("IPM", "MEM") ~ "Carbapenems",
    Antibiotics %in% c("CN", "AK")   ~ "Aminoglycosides",
    Antibiotics %in% c("CAZ", "FEP") ~ "Cephalosporins",
    Antibiotics %in% c("CIP", "LEV") ~ "Fluoroquinolones",
    Antibiotics == "TZP"             ~ "BLI",
    Antibiotics == "CT"              ~ "Colistin"
  ))

wilcox_results_r <- file7 %>%
  group_by(Antibiotics) %>%
  do(tidy(wilcox.test(RBF ~ Susceptibility, data = .))) 
pvalue_r <- wilcox_results_r$p.value  
# Reshape wilcox_results_r for geom_signif
wilcox_results_r_signif <- wilcox_results_r %>%
  select(Antibiotics, p.value) %>%
  mutate(
    group1 = "S",
    group2 = "NS"
  )

# Create the plot
Figure_6 <- ggplot(file7, aes(Susceptibility, RBF)) +
  geom_boxplot(aes(fill = Antibiotic_Class)) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Antibiotics) +
  theme_test() +
  geom_signif(
    data = wilcox_results_r_signif, # Use modified wilcox_results_r
    aes(xmin = group1, xmax = group2, annotations = custom_annotation(p.value), y_position = max(file1$meanOD.ODcut) * 1.06), 
    manual = TRUE, # Indicate manual input for comparisons
    step_increase = 0.1,
    vjust = 1.5,
    textsize = 4,     # Adjust text size as needed
    tip_length = 0.02 # Adjust tip length as needed
  ) +
  #stat_compare_means(method = "wilcox.test",  aes(label = after_stat(custom_annotation(p))), comparisons = list(c("S", "NS")), label.y = 6, label.x = 1.05, label.size = 5) + 
  #geom_text(data = wilcox_results_r, aes(x = 1.5, y = 5.75, label = custom_annotation(p.value)), size = 7)+ #use p value directly
  #stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c("S", "NS")), label.y = 0.75, label.x = 1.05, label.size = 2) + #only p value
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 13),
    strip.placement = "outside", strip.clip = "off",
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(hjust = 1, margin = margin(b = 5), size = 15),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 15), 
    axis.title = element_text(size = 15),
    panel.border = element_rect(linewidth = 1.2),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),  # Corrected font setting
    legend.background = element_rect(color = "black", linewidth = 1),
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Corrected border setting
  )


#MARI by Source
#MARI by Source
#MARI by Source
#MARI by Source
library("ggdist")
library("rstatix")
file4 <- file1 %>%
  filter(Source %in% c("W/S", "U", "T/A", "P", "S"))
#normality test
normality_test_plot9 <- file4 %>%
  group_by(Source) %>%
  summarize(p_value = shapiro.test(MARI)$p.value)  

# Perform Kruskal-Wallis test
kruskal_test_result_4 <- kruskal.test(MARI ~ Source, data = file4)


# Perform Dunn test and add p-values to dunn_test_results_4
dunn_test_results_4 <- file4 %>% 
  dunn_test(MARI ~ Source, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Source") %>%
  mutate(custom_label = custom_annotation(p.adj))  # Create custom labels using custom_annotation2

# Now create the ggplot with stat_pvalue_manual using the custom labels
MARI1 <- ggplot(file4, aes(x = Source, y = MARI)) +
  ggdist::stat_gradientinterval(aes(fill = Source), width = 0.5, color = "black", linewidth = 1) +
  scale_fill_brewer(palette = "BuPu") +
  theme_test() +
  stat_pvalue_manual(
    dunn_test_results_4,
    label = "custom_label", # Specify the custom labels here
    hide.ns = TRUE,
    size = 5,
    bracket.size = 0.7,
    step.increase = 0.1,
    y_position = max(file4$MARI) * 1.05
  ) +
  theme(
    axis.text.x = element_text(hjust = 1, margin = margin(b = 5), size = 12),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 15),
    panel.border = element_rect(linewidth = 1.2),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13),
    legend.background = element_rect(color = "black", linewidth = 1),
    legend.margin = margin(t = 30, r = 5, b = 5, l = 5, unit = "pt")
  )+
  ylim(c(min(file4$MARI), max(file4$MARI) * 1.7))
  #theme(axis.title.y = element_blank())

#Carbapenemase P NP and MARI
#Carbapenemase P NP and MARI
#Carbapenemase P NP and MARI
#Carbapenemase P NP and MARI
normality_test_plot10 <- file2 %>%
  group_by(Carbapenemase) %>%
  summarize(p_value = shapiro.test(MARI)$p.value)  
#non-normal data
#need to do wilcox.t-test
# Create the plot
wilcox_test_7 <- file2 %>% 
  wilcox_test(MARI ~ Carbapenemase) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Carbapenemase")

pvalue_MARI2 <- wilcox_test_7$p
MARI2 <- ggplot(file2, aes(x = Carbapenemase, y = MARI)) +
  ggdist::stat_gradientinterval(aes(fill =Carbapenemase),
                                width = .5, color = "black", linewidth = 1
  ) +
  scale_fill_brewer(palette = "BuPu") +
  theme_test() +
  geom_signif(comparisons = list(c("P", "NP")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_MARI2, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_MARI2),
              y_position = max(file2$MARI) * 1.06)+
  #stat_pvalue_manual(wilcox_test_7, label = "p.signif", hide.ns = FALSE,size = 5, bracket.size = 0.7, v.just = 1.2)+
  theme(axis.text.x = element_text(hjust = 1, margin = margin(b = 5), size = 12),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 15),
        panel.border = element_rect(linewidth = 1.2),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),  # Corrected font setting
        legend.background = element_rect(color = "black", linewidth = 1),
        legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Corrected border setting
  )
#Biofilm group vs MARI
#Biofilm group vs MARI
#Biofilm group vs MARI
#Biofilm group vs MARI
normality_test_plot11 <- file1 %>%
  group_by(Formation) %>%
  summarize(p_value = shapiro.test(MARI)$p.value)
wilcox_test_8 <- file1 %>% 
  wilcox_test(MARI ~ Formation) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE           ~ "NS"
  ))%>%
  add_xy_position(x = "Formation")
pvalue_MARI3 <- wilcox_test_8$p
MARI3 <- ggplot(file1, aes(x = Formation, y = MARI)) +
  ggdist::stat_gradientinterval(aes(fill = Formation),
                                width = .5, color = "black", linewidth = 1
  ) +
  scale_fill_brewer(palette = "BuPu") +
  theme_test() +
  geom_signif(comparisons = list(c("Considerable", "Negligible")), 
              step_increase = 0.1,
              vjust = 1.3,
              #annotations = paste0("p = ", format(pvalue_MARI1, scientific = FALSE, digits = 3)),
              annotations = custom_annotation(pvalue_MARI3),
              y_position = max(file1$MARI) * 1.06)+
  #stat_pvalue_manual(wilcox_test_8, label = "p.signif", hide.ns = FALSE, size = 5, bracket.size = 0.7, vjust= 1.2)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, margin = margin(b = 5), size = 12),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 15),
        panel.border = element_rect(linewidth = 1.2),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),  # Corrected font setting
        legend.background = element_rect(color = "black", linewidth = 1),
        legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Corrected border setting
  )
#combine the plot
library(patchwork)
Figure_3 <- MARI1 / (MARI2|MARI3)
Figure_3 <- Figure_3 + plot_annotation(tag_levels = "a")
Figure_3



file8$Percentage <- as.numeric(file8$Percentage)
#Figure_7
file8 <- read.csv("file8.csv")%>%
  na.omit()
ggplot(file8, aes(Biofilm_Gene, Percentage, fill = Source))+
  geom_bar(stat = "identity", width =0.7, just= 1, color = "black", position = "dodge", linewidth = 0.8) +
  scale_fill_brewer(palette = "BuPu")+
  labs(x = "Biofilm Gene", y = "Percentage") +
  theme_test()+
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100))+
  facet_wrap(~Bacteria, scales = "free_x" ,nrow = 2,ncol = 3, labeller = label_wrap_gen(width = 2))+
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(face = "italic",size = 12),
    strip.placement = "inside", strip.clip = "off",
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", margin = margin(b = 5), size = 12),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 14),
    panel.border = element_rect(linewidth = 1.2),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),  # Corrected font setting
    legend.background = element_rect(color = "black", linewidth = 0.8),
    legend.position = c(0.9, 0),     # Position: top-right corner (1, 1)
    legend.justification = c(0.9, 0),  # Justification: top-right corner (1, 1)
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt") # Corrected border setting
  )
  




#correlation between biofilm formation and MARI
#how to view correlation between two variable
Correlation <- file1 %>%
  group_by(Bacteria) %>%
  summarize(
    correlation = cor(meanOD.ODcut, MARI, method = "pearson"),
    p_value = cor.test(meanOD.ODcut, MARI, method = "pearson")$p.value
  )
Correlation <- Correlation %>%
  mutate(Significance = ifelse(p_value < 0.05, "Significant", "Not Significant")
#mostly non significant

ggplot(file1, aes(meanOD.ODcut,MARI, color = Bacteria))+
  geom_point( size = 3, alpha = 0.5  )+
  scale_fill_brewer()+
  #geom_smooth(se = TRUE)+
  #facet_wrap(~Bacteria)+
  theme_test()





#MDR percentage in each bacterial group
MDR_percentage <- file1%>%
  group_by(Bacteria, MDR)%>%
  summarise(n = n()) %>%          
  mutate(percentage = n/sum(n)*100)
  















### Rough
#gram-negative specific resistance
facet_labels <- c(
  "a", "b", "c", "d", "e",
  "f", "g", "h", "i", "j"
)

library(tidyr)
library(ggplot2)
library(ggpubr)
file7 <- read.csv("specific_resistance_gn.csv")%>%
  na.omit(file7)%>%
  gather(Antibiotics,Susceptibility,7:16)%>%
  mutate(Antibiotics = factor(Antibiotics, levels = c("IPM", "MEM", "CN", "AK", "CAZ", "FEP", "CIP", "LEV", "TZP", "CT")))%>%
  mutate(Antibiotic_Class = case_when(
    Antibiotics %in% c("IPM", "MEM") ~ "Carbapenems",
    Antibiotics %in% c("CN", "AK")   ~ "Aminoglycosides",
    Antibiotics %in% c("CAZ", "FEP") ~ "Cephalosporins",
    Antibiotics %in% c("CIP", "LEV") ~ "Fluoroquinolones",
    Antibiotics == "TZP"             ~ "BLI",
    Antibiotics == "CT"              ~ "Colistin"
  ))
wilcox_results_mmm <- file7 %>%
  group_by(Antibiotics) %>%
  do(tidy(wilcox.test(RBF ~ Susceptibility, data = .))) %>%
  ungroup() %>%
  mutate(
    p.signif = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      TRUE             ~ "NS"
    ),
    start = 1,  # Assuming group 1 is "NS"
    end = 2     # Assuming group 2 is "S"
  )
ggplot(file7, aes(x = Susceptibility, y = RBF)) +
  geom_boxplot(aes(fill = Antibiotic_Class)) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Antibiotics + labeller = labeller(Antibiotics = facet_labels)) +) +
  theme_test() +
  geom_signif(
    data = wilcox_results_mmm,
    aes(xmin = start, xmax = end, annotations = p.signif),
    manual = TRUE,
    y_position = max(file7$RBF),  # Adjust y_position as needed
    tip_length = 2,
    textsize = 4
  ) +
  labs(y = "RBF", x = "Susceptibility")


stat_compare_means(
  method = "kruskal.test", 
  label = "p.format", 
  comparisons = list(c("A. baumannii", "P. aeruginosa"), 
                     c("P. aeruginosa", "K. pneumoniae"),
                     c("A. baumannii","K. pneumoniae")),
  label.y = c(6.2, 5.2, 4.2)   # Adjust y positions for each comparison
) +
  