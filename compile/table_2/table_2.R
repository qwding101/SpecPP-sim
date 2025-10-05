##################################################################
###### Import packages, IBIAS and IMSE of the three models #######
##################################################################

library(tidyr)
library(dplyr)
library(stringr)
library(kableExtra)

Mc = paste0("M", 1:3)
tb.list.A10 = tb.list.A20 = tb.list.A40 = list()
for (i in seq_along(Mc)){
  tb.list.A10[[i]] = t(readRDS(paste0(Mc[i],"_summary.rds"))$A10)
  tb.list.A20[[i]] = t(readRDS(paste0(Mc[i],"_summary.rds"))$A20)
  tb.list.A40[[i]] = t(readRDS(paste0(Mc[i],"_summary.rds"))$A40)
}



###############################
###### Generate table 2 #######
###############################

tb.df = data.frame()
for (m in 1:3){
  df.m = as.data.frame(rbind(tb.list.A10[[m]], tb.list.A20[[m]], tb.list.A40[[m]]))
  df.m$ij = rep(row.names(tb.list.A10[[m]]),3)
  df.m$model = Mc[m]
  df.m$window = rep(c("$[-5,5]^2$", "$[-10,10]^2$", "$[-20,20]^2$"), each = nrow(tb.list.A10[[m]]))
  
  tb.df = rbind(tb.df, df.m)
  row.names(tb.df) = NULL
}
tb.df$ij = factor(tb.df$ij, levels = c("11","22","12"),
                  labels = c("$F^{(1,1)}_h$","$F^{(2,2)}_h$","$F^{(1,2)}_h$"))
tb.df = tb.df[, c(colnames(tb.df)[c(8,7,9)], colnames(tb.df)[1:6])] %>%
  arrange(model, ij)

tb.tex = tb.df %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  kbl(format = "latex",
      booktabs = T,
      escape = F,
      longtable=F,
      align = "c",
      col.names = c("Model", "Pseudo-spectrum", "Window", rep(c("IBIAS","IMSE"),3)),
      caption = "Table 2: IBIAS$^2$ and IMSE for the three pseudo-spectrum estimators.") %>%
  add_header_above(c(" " = 3,
                     r"($\\hat{I}_{h,n}$)" = 2,
                     r"($\\hat{F}_\\text{opt}$)" = 2,
                     r"($\\hat{F}_\\text{CV}$)" = 2),
                   escape = F) %>%
  collapse_rows(columns = 1:2, latex_hline = "custom", custom_latex_hline=1:2)

tb.tex = gsub("\\centering\\arraybackslash ", "", x = as.character(tb.tex), ignore.case = F, fixed = T)
tb.tex = gsub("\\caption", "\\caption*", x = as.character(tb.tex), ignore.case = F, fixed = T)


# Copy below result to LaTeX editor to generate table 2
# Remember to add `\usepackage{amsmath, multirow, booktabs, caption}` in the preamble to make it works

cat(tb.tex)

# \begin{table}
# 
# \caption*{Table 2: IBIAS$^2$ and IMSE for the three pseudo-spectrum estimators.}
# \centering
# \begin{tabular}[t]{ccccccccc}
# \toprule
# \multicolumn{3}{c}{ } & \multicolumn{2}{c}{$\hat{I}_{h,n}$} & \multicolumn{2}{c}{$\hat{F}_\text{opt}$} & \multicolumn{2}{c}{$\hat{F}_\text{CV}$} \\
# \cmidrule(l{3pt}r{3pt}){4-5} \cmidrule(l{3pt}r{3pt}){6-7} \cmidrule(l{3pt}r{3pt}){8-9}
# Model & Pseudo-spectrum & Window & IBIAS & IMSE & IBIAS & IMSE & IBIAS & IMSE\\
# \midrule
# &  & $[-5,5]^2$ & 0.00 & 1.14 & 0.00 & 0.83 & 0.00 & 0.24\\
# 
# &  & $[-10,10]^2$ & 0.00 & 1.07 & 0.00 & 0.29 & 0.00 & 0.12\\
# 
# & \multirow{-3}{*}{$F^{(1,1)}_h$} & $[-20,20]^2$ & 0.00 & 1.01 & 0.00 & 0.12 & 0.00 & 0.12\\
# \cmidrule{2-9}
# &  & $[-5,5]^2$ & 0.00 & 1.04 & 0.00 & 0.75 & 0.00 & 0.18\\
# 
# &  & $[-10,10]^2$ & 0.00 & 1.02 & 0.00 & 0.26 & 0.00 & 0.10\\
# 
# & \multirow{-3}{*}{$F^{(2,2)}_h$} & $[-20,20]^2$ & 0.00 & 1.00 & 0.00 & 0.11 & 0.00 & 0.11\\
# \cmidrule{2-9}
# &  & $[-5,5]^2$ & 0.39 & 145.85 & 0.27 & 101.55 & 0.04 & 16.72\\
# 
# &  & $[-10,10]^2$ & 0.18 & 117.77 & 0.05 & 28.29 & 0.02 & 9.95\\
# 
# \multirow{-9}{*}[1\dimexpr\aboverulesep+\belowrulesep+\cmidrulewidth]{M1} & \multirow{-3}{*}{$F^{(1,2)}_h$} & $[-20,20]^2$ & 0.21 & 102.52 & 0.02 & 10.97 & 0.02 & 10.97\\
# \cmidrule{1-9}
# &  & $[-5,5]^2$ & 0.01 & 1.19 & 0.01 & 0.90 & 0.01 & 0.38\\
# 
# &  & $[-10,10]^2$ & 0.00 & 1.03 & 0.00 & 0.31 & 0.00 & 0.16\\
# 
# & \multirow{-3}{*}{$F^{(1,1)}_h$} & $[-20,20]^2$ & 0.00 & 1.01 & 0.00 & 0.13 & 0.00 & 0.13\\
# \cmidrule{2-9}
# &  & $[-5,5]^2$ & 0.01 & 1.02 & 0.01 & 0.73 & 0.01 & 0.24\\
# 
# &  & $[-10,10]^2$ & 0.00 & 1.00 & 0.00 & 0.25 & 0.00 & 0.12\\
# 
# & \multirow{-3}{*}{$F^{(2,2)}_h$} & $[-20,20]^2$ & 0.00 & 0.98 & 0.00 & 0.11 & 0.00 & 0.11\\
# \cmidrule{2-9}
# &  & $[-5,5]^2$ & 0.06 & 27.80 & 0.04 & 19.62 & 0.02 & 5.48\\
# 
# &  & $[-10,10]^2$ & 0.04 & 22.52 & 0.01 & 5.69 & 0.01 & 2.55\\
# 
# \multirow{-9}{*}[1\dimexpr\aboverulesep+\belowrulesep+\cmidrulewidth]{M2} & \multirow{-3}{*}{$F^{(1,2)}_h$} & $[-20,20]^2$ & 0.04 & 20.16 & 0.01 & 2.29 & 0.01 & 2.29\\
# \cmidrule{1-9}
# &  & $[-5,5]^2$ & 0.01 & 1.16 & 0.01 & 0.87 & 0.01 & 0.31\\
# 
# &  & $[-10,10]^2$ & 0.00 & 1.05 & 0.00 & 0.31 & 0.00 & 0.15\\
# 
# & \multirow{-3}{*}{$F^{(1,1)}_h$} & $[-20,20]^2$ & 0.00 & 1.00 & 0.00 & 0.13 & 0.00 & 0.13\\
# \cmidrule{2-9}
# &  & $[-5,5]^2$ & 0.01 & 1.06 & 0.01 & 0.75 & 0.01 & 0.19\\
# 
# &  & $[-10,10]^2$ & 0.00 & 1.00 & 0.00 & 0.25 & 0.00 & 0.12\\
# 
# & \multirow{-3}{*}{$F^{(2,2)}_h$} & $[-20,20]^2$ & 0.00 & 0.99 & 0.00 & 0.11 & 0.00 & 0.11\\
# \cmidrule{2-9}
# &  & $[-5,5]^2$ & 0.40 & 288.41 & 0.28 & 198.14 & 0.10 & 34.83\\
# 
# &  & $[-10,10]^2$ & 0.48 & 208.98 & 0.11 & 49.69 & 0.05 & 20.46\\
# 
# \multirow{-9}{*}[1\dimexpr\aboverulesep+\belowrulesep+\cmidrulewidth]{M3} & \multirow{-3}{*}{$F^{(1,2)}_h$} & $[-20,20]^2$ & 0.34 & 174.90 & 0.03 & 19.10 & 0.03 & 19.10\\
# \bottomrule
# \end{tabular}
# \end{table}