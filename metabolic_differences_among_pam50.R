library(tidyverse)
library(hrbrthemes)
library(viridis)


#### from https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
# Load myeloma data from survminer package
library("survminer")


GSE25065 = transform(GSE25065,
                    SLC1A4 = as.numeric(SLC1A4))

# Perform the test
SLC1A4_p_adj_GSE25065 = compare_means(SLC1A4~ pam50_class,  data = GSE25065,
              ref.group = ".all.",
              method = "t.test", 
              p.adjust.method = "bonferroni")

SLC1A4_p_adj_GSE25065_anabolic_path


ggboxplot(GSE25065, x = "pam50_class", y = "SLC1A4", color = "pam50_class", 
          add = "jitter", 
          legend = "none", 
          palette = "jco",
          xlab = "Subtipos moleculares",
          ylab = expression(~bolditalic("SLC1A4"))) +
  stat_compare_means(method = "anova", 
                     label.y = 1.60) +
  theme(axis.title.y = element_text(size = 14, face = "italic", 
                                    margin = margin(t = 0, r = 18, b = 0, l = 0)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 18, r = 0, b = 0, l = 0))) +# Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all

GSE= SLC1A4_GSE25065_anabolic_path
##### GET two columns: a column with pam50 and a 
#### column with gene expression by bimodal gene


## 1. transpose, change colnames and convert to dataframe

GSE25065 = t(GSE25065)



## select metabolic bimodal genes
# note: ELOVL2, SLC1A40, FADS1, DPYSL4, UPP1, SLC16A3 doesnt exist
GSE25065_anabolic_path  = GSE1456 %>%
  select(pam50_class,
         SLC1A4
  )
  # myscela

#GSE25065_anabolic_path = GSE25065_anabolic_path %>% 
 # mutate(pam50_class = ifelse(as.character(pam50_class) == "Normal", "Normal-like", as.character(pam50_class)))

## make the first row my colnames
colnames(GSE1456) <- GSE1456[1,]
GSE1456 <- GSE1456[-1, ] 
GSE1456 = as.data.frame(GSE1456)
## write table
write.table(GSE25065_anabolic_path, file = "GSE25065_anabolic_path.txt", sep = "\t",
            row.names = FALSE)
