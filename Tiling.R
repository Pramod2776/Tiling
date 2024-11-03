setwd("~/path/to/workdir/")
library("Biostrings")
library("tidyverse")
library("seqinr")
library(ggplot2)

## read fasta or read sequence file and counts

reference1 = as.character(read.table("./CDS.txt"))
reference1_rev = reverseComplement(DNAString(reference1))
reference1_rev = as.character(reference1_rev)

sgRNA_5 = readxl::read_excel("./read_counts.xls", sheet = 5) %>%
  as.data.frame()

sgRNA_sequences_5 = sgRNA_5 %>%
  dplyr::filter(Tiling == "Tiling")

list_forward = lapply(sgRNA_sequences_5$sgRNA_sequence, function(x){
  print(x)
  m1 = matchPattern(x, reference1)
  if (length(m1) == 0){
    S1 = 0
    E1 = 0
  } else{
    S1 = start(m1)
    E1 = end(m1)
  }
  
  dat = data.frame(sgRNA_sequence = x,
                   S1 = as.character(print(S1)),
                   E1 =  as.character(print(E1)),
                   S1_S = as.character(print(paste0(S1,"_F"))),
                   E1_S =  as.character(print(paste0(E1,"_F"))))
  
}) %>%
  bind_rows()

list_rev = lapply(sgRNA_sequences_5$sgRNA_sequence, function(x){
  print(x)
  m1 = matchPattern(x, reference1_rev)
  if (length(m1) == 0){
    S1 = 0
    E1 = 0
  } else{
    S1 = start(m1)
    E1 = end(m1)
  }

  dat = data.frame(sgRNA_sequence = x,
                   S1 = as.character(print(S1)),
                   E1 =  as.character(print(E1)),
                   S1_S = as.character(print(paste0(S1,"_R"))),
                   E1_S =  as.character(print(paste0(E1,"_R"))))
  
}) %>%
  bind_rows()


df_for = list_forward %>%
  left_join(sgRNA_sequences_5, by = "sgRNA_sequence") %>%
  dplyr::select("sgRNA_sequence", "S1", "E1",  "S1_S", "E1_S", "Tiling", "RATIO")

df_for = df_for[
  with(df_for, order(S1, decreasing = T)),
]

df1_for = df_for %>%
  dplyr::filter(S1 > 0)

df_rev = list_rev %>%
  left_join(sgRNA_sequences_5, by = "sgRNA_sequence") %>%
  dplyr::select("sgRNA_sequence", "S1", "E1",  "S1_S", "E1_S", "Tiling", "RATIO")

df_rev = df_rev[
  with(df_for, order(S1, decreasing = T)),
]

df1_rev = df_rev %>%
  dplyr::filter(S1 > 0)
dim(df1_rev)

df1 = rbind(df1_for, df1_rev)
dim(df1)


df2 = read.csv("df5.csv")

df3 = df2 %>%
  dplyr::select("sgRNA_sequence", "S1", "RATIO")



df4 = melt(df3, id = c("sgRNA_sequence", "S1"))
colnames(df4) = c("sgRNA_sequence", "S1", "case","Ratio")


# rolling median function
rolling_median <- function(formula, data, n_roll = 10, ...) {
  formula = y ~ x
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}

predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}


ggplot(df4, aes(x=S1, y= Ratio, colour = case)) +
  geom_point(size = 2)+
  geom_smooth(formula = y ~ x, method = "rolling_median", se = FALSE,  n = 70)+
  ##expand_limits(x = 0, y = 0)+
  scale_x_continuous(breaks=seq(0, 2200, 100))+
  scale_y_continuous(breaks=seq(0, 4, 1))+
  coord_cartesian(expand = FALSE)+
  ylim(0,4)+
  xlim (0, 2200)+
  theme_classic() 


df2 = read.csv("df2.csv")
df3 = read.csv("df3.csv")
df4 = read.csv("df4.csv")
df5 = read.csv("df5.csv")


library(reshape)
plot_data = df2 %>%
  left_join(df3, by = "sgRNA_sequence") %>%
  left_join(df4, by = "sgRNA_sequence") %>%
  left_join(df5, by = "sgRNA_sequence")
colnames(plot_data)

plot_data = plot_data %>%
  dplyr::select("sgRNA_sequence", "S1.x", "E1.x", "S1_S.x", "E1_S.x", "Tiling.x", "Ratio1",
               "RATIO2", "RATIO3",
                "RATIO4")

write.csv(plot_data, "plot_data.csv")


library(dplyr)

99_data = plot_data %>% 
  rowwise() %>% 
  do(data.frame (Ratio_1= .$Ratio1, value = .$S1.x:.$ E1.x))

21_data = plot_data %>% 
  rowwise() %>% 
  do(data.frame (RATIO2 = .$RATIO2, value = .$S1.x:.$ E1.x))

46_data = plot_data %>% 
  rowwise() %>% 
  do(data.frame (RATIO3= .$RATIO3 value = .$S1.x:.$ E1.x))

DO_data = plot_data %>% 
  rowwise() %>% 
  do(data.frame (RATIO4= .$RATIO4 value = .$S1.x:.$ E1.x))


sgRNA_sequence_data = plot_data %>% 
  rowwise() %>% 
  do(data.frame (sgRNA_sequence = .$sgRNA_sequence, value = .$S1.x:.$ E1.x))


combined_data = cbind(sgRNA_sequence_data, DO_data, 46_data, 21_data, 99_data)

write.csv(combined_data, "plot_data_combined.csv")


DO_sort = combined_data %>%
  dplyr::select("sgRNA_sequence", "value", "RATIO4") %>%
  group_by(value) %>% slice(which.max(RATIO4))


46_sort = combined_data %>%
  dplyr::select("sgRNA_sequence", "value", "RATIO3) %>%
  group_by(value) %>% slice(which.max(RATIO3))



21_sort = combined_data %>%
  dplyr::select("sgRNA_sequence", "value", "RATIO2") %>%
  group_by(value) %>% slice(which.max(RATIO2))


99_sort = combined_data %>%
  dplyr::select("sgRNA_sequence", "value", "Ratio1") %>%
  group_by(value) %>% slice(which.max(Ratio1)

combined_data_sort = cbind(DO_sort, 46_sort, 21_sort, 99_sort)
write.csv(combined_data_sort, "pl_data_sort_nuc.csv")


colnames(plot_data)

plot_data_all = plot_data %>%
  dplyr::select("sgRNA_sequence", "S1.x",  "Ratio1", "RATIO2", 
                "RATIO3", "RATIO4")

plot_data_all_log = plot_data_all[, c(3:6)]
plot_data_all_log = log10(plot_data_all_log)

plot_data_all_log = cbind(plot_data[, c(1:2)], plot_data_all_log)


df_plot = melt(plot_data_all, id = c("sgRNA_sequence", "S1.x"))
head(df_plot)
colnames(df_plot) = c("sgRNA_sequence", "S1", "case","Ratio")
levels(df_plot$case)
head(df_plot)


df_plot_log = melt(plot_data_all_log, id = c("sgRNA_sequence", "S1.x"))
colnames(df_plot_log) = c("sgRNA_sequence", "S1", "case","Ratio")


plot_median_S = ggplot(df6, aes(x=S1, y= Ratio, colour = case)) +
  geom_point(size = 2)+
  geom_smooth(formula = y ~ x, method = "rolling_median", se = FALSE,  n = 70)+
  ##expand_limits(x = 0, y = 0)+
  scale_x_continuous(breaks=seq(0, 2000, 200))+
  ##scale_y_continuous(breaks=seq(0, 4, 1))+
  #coord_cartesian(expand = FALSE)+
  ##ylim(0,20)+
  xlim (0, 2000)+
  xlab("nucleotide_sequence")+
  ylab("Treatment_ratios")+
  theme_classic() 
ggsave('ratios_complete.png', plot_median_S, height = 8, width = 8)

plot_median = ggplot(df_plot, aes(x=S1, y= Ratio, colour = case)) +
  geom_point(size = 2)+
  geom_smooth(formula = y ~ x, method = "rolling_median", se = FALSE,  n = 70)+
  ##expand_limits(x = 0, y = 0)+
  scale_x_continuous(breaks=seq(0, 1800, 200))+
  scale_y_continuous(breaks=seq(0, 4, 1))+
  #coord_cartesian(expand = FALSE)+
  ylim(0,4)+
  xlim (0, 1800)+
  xlab("nucleotide_sequence")+
  ylab("Treatment_ratios")+
  theme_classic() 
ggsave('ratios_complement.png', plot_median, height = 8, width = 8)

plot_median_log = ggplot(df_plot_log, aes(x=S1, y= Ratio, colour = case)) +
  geom_point(size = 2)+
  geom_smooth(formula = y ~ x, method = "rolling_median", se = FALSE,  n = 70)+
  ##expand_limits(x = 0, y = 0)+
  #scale_x_continuous(breaks=seq(0, 1800, 200))+
  #scale_y_continuous(breaks=seq(0, 4, 1))+
  #coord_cartesian(expand = FALSE)+
  #ylim(-10,4)+
  #xlim (0, 1800)+
  xlab("nucleotide_sequence")+
  ylab("Treatment_log_ratios")+
  theme_classic() 
ggsave('log_ratios_complement.png', plot_median_log, height = 8, width = 8)




