#!/BIGDATA1/app/R/4.1.0/bin/Rscript
## phenotype simulation  ##
## liwn 2021-09-22 ##

## 获取命令行参数
spec <- matrix(
  c(
    "gt",       "A", 1, "character", "[Required] qtl plink file prefix of A",
    "h2",       "1", 1, "character", "[Optional] heritability of populations [0.5]",
    "rg",       "r", 1, "character", "[Optional] genetic correlation between breeds [0.2]",
    "win",      "w", 1, "integer",   "[Optional] number of snp in each bins [100]",
    "binf",     "b", 1, "character", "[Optional] ld block file",
    "binc",     "C", 1, "integer",   "[Optional] column number of nSNP in block file [5]",
    "miss",     "m", 1, "character", "[Optional] miss value in phenotype [-99]",
    "dist_cor", "c", 1, "character", "[Optional] distribution of genetic correlation of qtl [identical/uniform]",
    "qtlbins",  "q", 1, "integer",   "[Optional] [0.01]",
    "qtlnsnp",  "Q", 1, "integer",   "[Optional] [10]",
    "qtlpct",   "p", 1, "integer",   "[Optional] Proportion of QTL to total SNPs [0.01]",
    "seed",     "s", 1, "integer",   "[Optional] seed of random effect sampling",
    "mean",     "e", 1, "character", "[Optional] number of poopulations [1.0]",
    "overlap",  "v", 2, "character", "[Optional] allow QTLs and markers to overlap",
    "out",      "o", 1, "character", "[Optional] output phenotype file name",
    "qtlf",     "t", 1, "character", "[Optional] output QTL effects file name",
    "fid",      "f", 0, "logical",   "[Optional] whether output family id",
    "help",     "h", 0, "logical",   "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt::getopt(spec = spec)

## 检查参数
n_need <- length(c(opt$gt))
if (!is.null(opt$help) || n_need < 1) {
  cat(paste(getopt::getopt(spec = spec, usage = TRUE), "\n"))
  quit()
}

## 自定义函数
argp_parse <- function(char, num, sep = " ", label = "parameters") {
  argv <- as.numeric(unlist(strsplit(char, sep)))

  if (length(argv) == 1) {
    return(rep(argv, num))
  } else if (length(argv) != num) {
    cat("The number of", label, "should be:", num, "!\n")
    quit()
  } else {
    return(argv)
  }
}

## 加载需要的程序包（未安装的话请提前安装）
cat("Loading required packages... \n")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("dplyr"))

## 默认参数
if (is.null(opt$seed)) opt$seed <- Sys.time()
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$h2)) opt$h2 <- "0.5"
if (is.null(opt$rg)) opt$rg <- "0.2"
if (is.null(opt$qtlbins)) opt$qtlbins <- 0.1 ## 大于1为具体存在相关的区间数，小于1为相关的区间占总区间的比例
if (is.null(opt$qtlpct)) opt$qtlpct <- 0.01 ## QTL数占总SNP数的比例
if (is.null(opt$qtlnsnp)) opt$qtlnsnp <- 10
if (is.null(opt$dist_cor)) opt$dist_cor <- "identical"
if (is.null(opt$mean)) opt$mean <- "1.0"
if (is.null(opt$binc)) opt$binc <- 5
if (is.null(opt$win)) opt$win <- 100

## 随机数种子
set.seed(as.integer(opt$seed))

## 检查基因型文件是否存在
bfiles <- trimws(opt$gt, which = c("both"), whitespace = "[ \t\r\n]")  ## 去除首尾空格
bfiles <- unlist(strsplit(bfiles, " "))
for (prefix in bfiles) {
  for (suffix in c(".map", ".ped")) {
    f <- paste0(prefix, suffix)
    if (!file.exists(f)) {
      cat("file needed not found:", f, "\n")
      quit()
    }
  }
}

## 品种数
np <- length(bfiles)

## 相关系数矩阵上三角元素个数
nrg <- choose(np, 2)

## 参数解析
rgs <- argp_parse(opt$rg, nrg, label = "rg")
h2s <- argp_parse(opt$h2, np, label = "h2")
means <- argp_parse(opt$mean, np, label = "mean")

## 读取plink格式基因型
cat("Loading plink genotype files... \n")
maps <- list()
peds <- list()
nsnps <- c()
for (i in 1:np) {
  maps[[i]] <- fread(paste0(bfiles[i], ".map"))
  peds[[i]] <- fread(paste0(bfiles[i], ".ped"))
  names(maps[[i]]) <- c("CHR", "SNP", "cM", "POS")
  nsnps <- c(nsnps, nrow(maps[[i]]))
}

## 共有标记位点
if (length(unique(nsnps)) > 1) {
  cat("The number of markers in the genotype file is inconsistent with that in the ldblock file!\n")
  quit()
}

## snp所属区间
map <- maps[[1]]
if (!is.null(opt$binf)) {
  bin <- fread(paste0(opt$binf))
  names(bin)[opt$binc] <- "nsnp"
} else {
  bin <- nsnp_bin <- NULL
  chrs <- unique(map$CHR)
  for (i in seq_along(chrs)) {
    nsnp_i <- sum(map$CHR == chrs[i])
    win_num <- floor(nsnp_i / opt$win)
    if (win_num < 1) win_num <- 1
    nsnp_bini <- rep(opt$win, win_num)
    nsnp_bini[win_num] <- nsnp_bini[win_num] + nsnp_i - sum(nsnp_bini)
    nsnp_bin <- c(nsnp_bin, nsnp_bini)
  }
  bin <- data.frame(bin = seq_len(length(nsnp_bin)), nsnp = nsnp_bin)
}
map$bin <- rep(seq_len(nrow(bin)), times = bin$nsnp)

### 挑选QTL ###
## 每个区间挑选一个标记位点
qmap <- map %>%
  group_by(bin) %>%
  slice_sample()

## 随机挑选qtlbins个区间，再挑选qtlnsnp-1个标记位点
bins_sel <- c()
if (opt$qtlbins > 0) {
  if (opt$qtlbins < 1) {
    opt$qtlbins <- floor(opt$qtlbins * nrow(bin))
    cat("The number of bins with non-zero correaltion is:", opt$qtlbins, "\n")
  }
  bins_sel <- sample(qmap$bin, opt$qtlbins)
  qmap2 <- subset(map, bin %in% bins_sel & !SNP %in% qmap$SNP) %>%
    group_by(bin) %>%
    slice_sample(n = opt$qtlnsnp - 1)
  qmap <- bind_rows(qmap, qmap2)
}

## 保证qtl数符合要求
nqtl <- round(opt$qtlpct * nrow(map))
over <- nrow(qmap) - nqtl
if (over > 0) {
  keep_single <- nqtl - sum(qmap$bin %in% bins_sel)
  keep_qtl <- sample(qmap$SNP[!qmap$bin %in% bins_sel], keep_single)
  qmap <- subset(qmap, SNP %in% keep_qtl | bin %in% bins_sel)
}
cat("total number of qtl: ", nrow(qmap), "\n")

## 提取基因型信息
qtl <- ebvs <- list()
fids <- c()
for (i in 1:np) {
  qtl_index <- which(maps[[i]]$SNP %in% qmap$SNP)
  col_index <- c(1:6, 2 * qtl_index + 6, 2 * qtl_index - 1 + 6)

  ## QTL基因含量矩阵
  qtl[[i]] <- subset(peds[[i]], select = col_index)
  ebvs[[i]] <- qtl[[i]][, 2:1]
  names(ebvs[[i]]) <- c("iid", "fid")

  ## 品种ID
  fids <- c(fids, unique(qtl[[i]]$V1))
}

## 判断fid合法性
if (length(fids) != np) {
  cat("The 'fid' in the ped file is not a population identifier!\n")
  quit()
}

## 加性效应协方差阵
cormat <- diag(np)
index <- 1
for (i in 1:(np - 1)) {
  for (j in (i + 1):np) {
    cormat[i, j] <- cormat[j, i] <- rgs[index]
    index <- index + 1
  }
}

## 存在遗传相关的标记
if (opt$qtlbins > 0) { # nolint
  ## 遗传相关的列
  rg_cols <- c()
  for (i in 1:(np - 1)) {
    for (j in (i + 1):np) {
      rg_cols <- c(rg_cols, paste(fids[i], fids[j], sep = "_"))
    }
  }

  if (opt$dist_cor == "identical") {
    qtl_cor <- which(qmap$bin %in% bins_sel)
    nsample <- length(qtl_cor)
    qmap[qtl_cor, fids] <- mvrnorm(nsample, rep(0, np), cormat)
    qmap[qtl_cor, rg_cols] <- matrix(rep(rgs, times = nsample), nsample, nrg, byrow = TRUE)
  } else {
    ## 各区间遗传相关
    if (opt$dist_cor == "normal") {
      rg_bins <- rnorm(opt$qtlbins, mean = 0, sd = 0.5)
    } else if (opt$dist_cor == "uniform") {
      rg_bins <- runif(opt$qtlbins, -1, 1)
    } else {
      cat("dist_cor can noly be identical, normal or uniform.\n")
      quit()
    }

    ## 抽样
    qmap[, fids] <- NA
    rgsi <- rep(0, nrg)
    for (i in seq_along(bins_sel)) {
      ## 标记效应协方差阵
      cormati <- cormat
      ## 赋值
      index <- 1
      for (j in 1:(np - 1)) {
        for (k in (j + 1):np) {
          if (cormat[j, k] == 0) next
          rgij <- rg_bins[i] - mean(rg_bins) + cormat[j, k]

          ## 保证合法性
          if (rgij > 1) rgij <- 1.0
          if (rgij < -1) rgij <- -1.0

          rgsi[index] <- cormati[k, j] <- cormati[j, k] <- rgij
          index <- index + 1
        }
      }

      ## 保证矩阵正定性
      cormat_pd <- nearPD(cormati, keepDiag = TRUE) # default
      if (!cormat_pd$converged) {
        stop("warning! The positive definiteness of genetic covariance matrix cannot be guaranteed.\n")
      } else if (cormat_pd$iterations > 1) {
        cat("add a small value to genetic effect covariance matrix\n")
        cormati <- cormat_pd$mat
      }

      ## 标记效应抽样
      qtl_cor <- which(qmap$bin == bins_sel[i])
      nsample <- length(qtl_cor)
      qmap[qtl_cor, fids] <- mvrnorm(nsample, rep(0, np), cormati)
      
      ## 保存相关系数
      qmap[qtl_cor, rg_cols] <- matrix(rep(rgsi, times = nsample), nsample, nrg, byrow = TRUE)
    }
  }

  ## 不存在遗传相关的标记
  qtl_no_cor <- which(is.na(qmap[[fids[1]]]))
  qmap[qtl_no_cor, fids] <- mvrnorm(length(qtl_no_cor), rep(0, np), diag(np))
  qmap[qtl_no_cor, rg_cols] <- 0
}
cat("global genetical correlation is:\n")
cor(subset(qmap, select = fids))

## 生成表型信息
cat("Creating phenotypes files... \n")
pheno <- data.frame()
for (i in seq_len(np)) {
  gene <- qtl[[i]][, 7:ncol(qtl[[i]])]

  ## 如果基因型为字母表示，转换为1、2表示
  if (gene[1, 1] %in% LETTERS) {
    cat("Converting genotype data to 012 format...\n")

    ## 参考碱基
    ref_allele <- unlist(subset(gene, c(TRUE, rep(FALSE, nrow(gene) - 1)),
      select = seq(1, ncol(gene), 2)
    ))

    ## 根据参考碱基转换为1、2标识符形式
    gene12 <- as.matrix(gene)
    for (j in 1:(ncol(gene) / 2)) {
      gene12[, 2 * j - 1] <- as.integer(subset(gene12, select = 2 * j - 1) == ref_allele[j]) + 1
      gene12[, 2 * j] <- as.integer(subset(gene12, select = 2 * j) == ref_allele[j]) + 1
    }

    ## 字符转数值
    gene <- apply(gene12, 2, as.numeric)
  }

  ## 基因含量转成012格式
  # cat("Converting genotype data to 012 format...\n")
  gene <- sapply(
    seq(1, ncol(gene) - 1, by = 2),
    function(i) rowSums(gene[, i:(i + 1)]) - 2
  )

  ## 计算加性效应
  ebvs[[i]]$tbv_raw <- apply(gene, 1, function(x) sum(x * qmap[[fids[i]]]))

  ## 校正为标准正态
  ebv_mean <- mean(ebvs[[i]]$tbv_raw)
  ebv_sd <- sd(ebvs[[i]]$tbv_raw)
  ebvs[[i]]$tbv_nor <- (ebvs[[i]]$tbv_raw - ebv_mean) / ebv_sd

  ## 残差效应(不相关)
  ve <- ((1 / h2s[i]) - 1) * var(ebvs[[i]]$tbv_nor)
  ebvs[[i]]$envir <- rnorm(nrow(ebvs[[i]]), 0, sqrt(ve))

  ## 加性 + 残差 + 群体均值 = 表型
  ebvs[[i]]$phe <- ebvs[[i]]$tbv_nor + ebvs[[i]]$envir + means[i]

  ## 合并各个品种的表型
  pheno <- rbind(pheno, ebvs[[i]])
}

## 群(品种)效应标识符
pheno$breed <- as.integer(as.factor(pheno$fid))

## 如不需要，则去掉fid列
if (is.null(opt$fid)) {
  pheno <- subset(pheno, select = colnames(pheno) != "fid")
}

## 列排序
pheno <- pheno[, c("iid", "fid", "breed", "tbv_raw", "tbv_nor", "envir", "phe")]

## 输出表型文件
if (is.null(opt$out)) opt$out <- paste(c(fids, "pheno.txt"), collapse = "_")
write.table(pheno, opt$out, row.names = FALSE, quote = FALSE)
cat("phenotypes output to file:", opt$out, "\n")

## 加性效应真值
if (is.null(opt$qtlf)) opt$qtlf <- paste(c(fids, "qtl.txt"), collapse = "_")
qmap <- qmap[order(qmap$bin), ]
write.table(qmap, opt$qtlf, row.names = FALSE, quote = FALSE)
cat("true qtl effects output to file:", opt$qtlf, "\n")

## 标记与QTL是否重叠
if (is.null(opt$overlap) && !is.null(opt$binf)) {
  ## 更新每个区间内的标记数
  qtl_num <- as.data.frame(table(as.factor(qmap$bin)), stringsAsFactors = FALSE)
  qtl_num$Var1 <- as.integer(qtl_num$Var1)
  if (nrow(qtl_num) < nrow(bin)) {
    bins <- seq_len(nrow(bin))
    add_qtl_num <- data.frame(Var1 = bins[!bins %in% qmap$bin], Freq = 0)
    qtl_num <- bind_rows(qtl_num, add_qtl_num)
  }
  qtl_num <- qtl_num[order(qtl_num$Var1), ]
  bin$nsnp <- bin$nsnp - qtl_num$Freq
  newbinf <- paste0(basename(opt$binf), ".noqtl")
  write.table(bin$nsnp, newbinf, row.names = FALSE, quote = FALSE, col.names = FALSE)
  cat("new bins file output to:", newbinf, "\n")
}

quit()

## debug
# setwd('/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/identical/cor0.2')
opt <- list()
opt$h2 <- "0.3 0.3 0.1 0.1"
opt$mean <- "1.0 1.0 0.5 0.5"
opt$rg <- "0.2"
opt$gt <- " /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/Am /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/Bm /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/Cm /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/Dm" # nolint
# opt$binf <- "/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Three/M_bin.txt"
opt$fid <- "true"
opt$win <- 100
opt$dist_cor <- "identical"
opt$seed <- 168181
opt$out <- "pheno_sim.txt"
opt$qtlf <- "qtl_info.txt"
