#############################################
###### Import packages and bandwidths #######
#############################################

library(kableExtra)
library(dplyr)

Ac = c("A10","A20","A40")
bopt = c(0.5, 400^(-1/6), 1600^(-1/6))
tb.list.M1 = tb.list.M2 = tb.list.M3 = list()
for (i in seq_along(Ac)){
  tb.list.M1[[i]] = readRDS(paste0("../../estimate/M1/", Ac[i], "/M1_", Ac[i], "_optbandwidth.rds"))
  tb.list.M2[[i]] = readRDS(paste0("../../estimate/M2/", Ac[i], "/M2_", Ac[i], "_optbandwidth.rds"))
  tb.list.M3[[i]] = readRDS(paste0("../../estimate/M3/", Ac[i], "/M3_", Ac[i], "_optbandwidth.rds"))
}
names(tb.list.M1) = names(tb.list.M2) = names(tb.list.M3) = Ac




###############################
###### Generate table 1 #######
###############################

df.band = cbind(lapply(tb.list.M1, summary)[[1]], lapply(tb.list.M2, summary)[[1]], lapply(tb.list.M3, summary)[[1]],
                lapply(tb.list.M1, summary)[[2]], lapply(tb.list.M2, summary)[[2]], lapply(tb.list.M3, summary)[[2]],
                lapply(tb.list.M1, summary)[[3]], lapply(tb.list.M2, summary)[[3]], lapply(tb.list.M3, summary)[[3]])
df.band = data.frame(df.band)
colnames(df.band) = c(paste0(rep(Ac, each = 3) ,paste0("M", 1:3)))
df.band = rbind(df.band, rep(bopt, each = 3))
df.band$val = c("min", "Q1", "Q2", "Mean", "Q3", "max", "$b_\\text{opt}$")
df.band = df.band[df.band$val %in% c("Q1","Q2","Q3","Mean","$b_\\text{opt}$"),]
df.band = df.band %>% slice(match(c("$b_\\text{opt}$","Q1","Q2","Q3","Mean"), val))
df.band = df.band[, c("val", colnames(df.band)[-ncol(df.band)])]
row.names(df.band) = NULL

tb.tex = df.band %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  kbl(format = "latex",
      booktabs = T,
      escape = F,
      longtable = F,
      align = "c",
      col.names = c("", rep(c("M1","M2","M3"), 3)),
      caption = "Table 1: The three quartiles (second to fourth row) and the average (fifth row) of $b_\\text{CV}$ from 500
replications for each window and model. The first row indicates the optimal bandwidth.") %>%
  add_header_above(c(" " = 1,
                     r"($D = [-5,5]^2$)" = 3,
                     r"($D = [-10,10]^2$)" = 3,
                     r"($D = [-20,20]^2$)" = 3),
                   escape = F) 
tb.tex = gsub("\\caption", "\\caption*", x = as.character(tb.tex), ignore.case = F, fixed = T)

# Copy below result to LaTeX editor to generate table 1
# Remember to add `\usepackage{amsmath, booktabs, caption}` in the preamble to make it works
cat(tb.tex)

# \begin{table}
# 
# \caption*{Table 1: The three quartiles (second to fourth row) and the average (fifth row) of $b_\text{CV}$ from 500
#   replications for each window and model. The first row indicates the optimal bandwidth.}
# \centering
# \begin{tabular}[t]{cccccccccc}
# \toprule
# \multicolumn{1}{c}{ } & \multicolumn{3}{c}{$D = [-5,5]^2$} & \multicolumn{3}{c}{$D = [-10,10]^2$} & \multicolumn{3}{c}{$D = [-20,20]^2$} \\
# \cmidrule(l{3pt}r{3pt}){2-4} \cmidrule(l{3pt}r{3pt}){5-7} \cmidrule(l{3pt}r{3pt}){8-10}
# & M1 & M2 & M3 & M1 & M2 & M3 & M1 & M2 & M3\\
# \midrule
# $b_\text{opt}$ & 0.50 & 0.50 & 0.50 & 0.37 & 0.37 & 0.37 & 0.29 & 0.29 & 0.29\\
# Q1 & 1.13 & 1.06 & 1.08 & 0.59 & 0.54 & 0.55 & 0.37 & 0.28 & 0.28\\
# Q2 & 1.51 & 1.18 & 1.18 & 0.76 & 0.60 & 0.60 & 0.40 & 0.30 & 0.30\\
# Q3 & 1.88 & 1.55 & 1.53 & 0.84 & 0.77 & 0.75 & 0.40 & 0.32 & 0.30\\
# Mean & 1.61 & 1.26 & 1.29 & 0.80 & 0.66 & 0.65 & 0.39 & 0.31 & 0.30\\
# \bottomrule
# \end{tabular}
# \end{table}
