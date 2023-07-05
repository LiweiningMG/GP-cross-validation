#!/work/apps/tools/conda/minconda3/20230202/bin/Rscript

## 检查需要的R包是否都已安装  ##
## liwn 2023-07-05 ##

# 要检查的 R 包列表
packages <- c("data.table", "LAVA", "getopt", "MASS", "RcppEigen", "Matrix", "ggplot2", "reshape2",
  "pedigree", "RcppArmadillo", "dplyr", "coda", "lattice", "stringr", "Cairo", "Rcpp", "reshape")

# 检查并安装缺失的包
missing_packages <- packages[!packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  message("以下包未安装：")
  message(missing_packages)
  message("请在R语言中运行以下命令安装：")
  message(paste0("install.packages(\"", missing_packages, "\")"))
}
