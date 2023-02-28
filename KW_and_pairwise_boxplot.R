library (ggpubr)

# Kruskal Wallis + pairwise comparisons added automatically, but need to define which pairwise comparisons needed first
my_comparisons <- list(c("group_1", "group_2"), c("group_1", "group_3"), c("group_2","group_3")) #need to define all groups you want to compare with each other

ggboxplot(dataframe, x = "group", y = "marker_expression", color = "group",
          add = "jitter", legend = "none") +
  stat_compare_means(label.y = 5)+        # Add global annova p-value, default is Kruskal-Wallis (method="anova" if anova required)
  stat_compare_means(comparisons=my_comparisons, label = "p.signif", method = "t.test") #adds pairwise comparisons. Change label to "p.adj" to display exact adjusted p values 


#Kruskal Wallis + Dunns post test. Will do all copmarisons but need to adjust position of p values on plot
library(FSA)
DT <- dunnTest(group~marker_expression, data = dataframe, method = "bh")

ggboxplot(dataframe, x = "group", y = "marker_expression", color = "group",
          add = "jitter", legend = "none") +
  stat_compare_means(label.y = 5)+        # Add global annova p-value, default is Kruskal-Wallis (method="anova" if anova required)
  stat_pvalue_manual(DT, label = "p.adj", y.position = c(3,4,2)) #adds pairwise comparisons. 
