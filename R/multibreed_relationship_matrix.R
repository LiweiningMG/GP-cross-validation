#!/apps/local/software/program/R-4.0.2/bin/Rscript
##!/usr/bin/Rscript
# 李伟宁 2021-12-08
## 计算多群体关系矩阵 ###
## 输入两个群体合并基因型文件(plink的012格式，*.raw) ###
## 输出不同方法计算的多群体关系矩阵[(idA+idB) x (idA+idB)] ###
## 输出格式为三列，id1 id2 value，tab分隔符
## 使用前需要安装相应的包，以及调用c++时需要的组件
## debug
# opt = list(rawf='popA_popB.raw', idAf='popA.ids', idBf='popB.ids')

# Load packages
cat('Loading packages needed...\n')
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("RcppArmadillo"))
suppressPackageStartupMessages(library("RcppEigen"))

## command parameters
spec <- matrix(
  c("rawf",    "R", 2, "character", "[Required] Raw file after merging two populations",
    "idAf",    "A", 2, "character", "[Required] ID file of population A with only 1 column",
    "idBf",    "B", 2, "character", "[Required] ID file of population B with only 1 column",
    "phef",    "P", 2, "character", "[Optional] If columns 1 and 2 in G matrix are row/column index, the ID in the phenotype file is renumbered",
    "method",  "M", 2, "character", "[Optional] Method of calculating G matrix/[Yvonne(default)/No/Chen/Mean]",
    "workdir", "D", 2, "character", "[Optional] working directory",
    "out",     "O", 2, "character", "[Optional] output file name(Can include path)",
    "tol",     "T", 0, "double",    "[Optional] tol of make matrix positive define",
    "inv",     "V", 0, "logical",   "[Optional] Output inverse matrix",
    "logdet",  "L", 0, "logical",   "[Optional] Add a row of '0 0 log(det)' (for DMU)",
    "header",  "H", 0, "logical",   "[Optional] include header in output files/NULL",
    "index",   "I", 0, "logical",   "[Optional] The first two columns of the output file are line index instead of ID\n",
    "help",    "h", 0, "logical",   "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt::getopt(spec=spec)

## 检查必要参数
if(!is.null(opt$help) || is.null(opt$rawf) || is.null(opt$idAf) || is.null(opt$idBf)){
   ## print help message
  cat(paste(getopt::getopt(spec=spec, usage = TRUE), "\n"))
  quit(status = -1)
}

## 默认参数设置
if (is.null(opt$header)) header = FALSE else header = TRUE
if (is.null(opt$tol)) opt$tol = 1e-6
if (is.null(opt$out)) opt$out = 'Yvonne2017'
if (is.null(opt$method)) opt$method = 'Yvonne'

## 检查method参数
if (!opt$method %in% c('Yvonne', 'Chen', 'Mean', 'No')) {
  cat('Method can only have one of these options: Yvonne/No/Chen/Mean\n')
  quit(status = -1)
}

## 更换工作路径
if (!is.null(opt$workdir)) setwd(opt$workdir)

## 准备加速矩阵乘法的c++脚本
cpp <- '// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]'
cpp <- c(cpp, '#include <RcppArmadillo.h>')
cpp <- c(cpp, '#include <RcppEigen.h>')
cpp <- c(cpp, '// [[Rcpp::export]]')
cpp <- c(cpp, 'SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){')
cpp <- c(cpp, '    Eigen::MatrixXd C = A * B;')
cpp <- c(cpp, '    return Rcpp::wrap(C);')
cpp <- c(cpp, '}')
cpp <- c(cpp, '// [[Rcpp::export]]')
cpp <- c(cpp, 'arma::mat getInv(arma::mat M) {')
cpp <- c(cpp, '    return arma::inv(M);')
cpp <- c(cpp, '}')

write.table(cpp, 'Matrix_multiplication.cpp', sep = '\n', col.names = FALSE,
            row.names = FALSE, quote = FALSE)

## 加载c++函数
cat('Loading C++ functions...\n')
sourceCpp('Matrix_multiplication.cpp')

## 读取基因型文件
cat('Loading genotype...\n')
raw <- fread(opt$rawf)

## 获取不同群体id
idA <- read.table(opt$idAf)
idB <- read.table(opt$idBf)

## 群体个体数
nIndA <- nrow(idA)
nIndB <- nrow(idB)
nIndT <- nIndA + nIndB

## 不同群体基因含量矩阵
geno_A <- as.matrix(raw[raw$IID %in% idA$V1, -c(1:6)])
geno_B <- as.matrix(raw[raw$IID %in% idB$V1, -c(1:6)])

## 标记数
nSNPA <- ncol(geno_A)
nSNPB <- ncol(geno_B)
nSNPT <- ncol(geno_A)

## 基因频率
cat('Calculating gene frequency...\n')
AF_A <- apply(geno_A, 2, sum)/nIndA/2
Sum_2pq_A <- sum(2 * AF_A * (1 - AF_A))

AF_B <- apply(geno_B, 2, sum)/nIndB/2
Sum_2pq_B <- sum(2 * AF_B * (1 - AF_B))

Sum_2pq_AB <- 2 * sum(sqrt(AF_A * (1 - AF_A) * AF_B * (1 - AF_B)))

AF_T <- apply(raw[, -c(1:6)], 2, sum)/nIndT/2
Sum_2pq_T <- sum(2 * AF_T * (1 - AF_T))

## Z-matrices
PA <- matrix(2 * AF_A, byrow = TRUE, nrow = nIndA, ncol = nSNPA)
ZA <- geno_A - PA

PB <- matrix(2 * AF_B, byrow = TRUE, nrow = nIndB, ncol = nSNPB)
ZB <- geno_B - PB

PT <- matrix(2 * AF_T, byrow = TRUE, nrow = nIndT, ncol = nSNPT)
ZT <- as.matrix(raw[, -c(1:6)]) - PT

## 删除变量，释放内存
rm(geno_A)
rm(geno_B)
rm(raw)
gclog <- gc(verbose = FALSE)

## 关系矩阵容器
Gmat <- matrix(0, nrow = nIndT, ncol = nIndT)

## 群体在Gmat中的位置索引
idA_index <- c(rep(TRUE, nIndA), rep(FALSE, nIndB))

## 基因组关系矩阵校正因子
scaleA <- Sum_2pq_A
scaleB <- Sum_2pq_B
if (opt$method == 'Yvonne') {
  scaleAB <- sqrt(Sum_2pq_A) * sqrt(Sum_2pq_B)
} else if (opt$method == 'Chen') {    
  scaleAB <- Sum_2pq_AB
} else if (opt$method == 'Mean') {  
  scaleA <- scaleB <- scaleAB <- Sum_2pq_T
} else if (opt$method == 'No') {
  scaleA <- scaleB <- scaleAB <- 1
}

## 生成G矩阵
cat('Generating G matrix...\n')
if (opt$method != 'Mean') {
  Gmat[idA_index,  idA_index ] <- eigenMapMatMult(ZA, t(ZA)) / scaleA
  Gmat[!idA_index, !idA_index] <- eigenMapMatMult(ZB, t(ZB)) / scaleB
  Gmat[!idA_index, idA_index ] <- eigenMapMatMult(ZB, t(ZA)) / scaleAB
  Gmat[idA_index,  !idA_index] <- eigenMapMatMult(ZA, t(ZB)) / scaleAB
} else {
  Gmat[, ] <- eigenMapMatMult(ZT, t(ZT)) / Sum_2pq_T
}

## 保证矩阵正定
make.positive.definite <- function(M, tol=1e-6) {
  eig <- eigen(M, symmetric = TRUE)
  rtol <- tol * eig$values[1]
  if(min(eig$values) < rtol) {
    cat('Make the matrix positive definite...\n')
    vals <- eig$values
    vals[vals < rtol] <- rtol
    srev <- eigenMapMatMult(eig$vectors, vals * t(eig$vectors))
    dimnames(srev) <- dimnames(M)
    return(srev)
  } else {
    return(M)
  }
}

## 保证矩阵正定
Gmat_pd <- make.positive.definite(Gmat, tol = opt$tol)

## 是否输出逆矩阵
if (!is.null(opt$inv)) {
  cat('Inverting matrix...\n')
  opt$out <- paste0(opt$out, '.grm.inv')
  Gmat_out <- getInv(Gmat_pd)
} else {
  opt$out <- paste0(opt$out, '.grm')
  Gmat_out <- Gmat_pd
}

## 是否以id命名
if (is.null(opt$index)) {
  ## 输出文件中1、2列为个体基因型文件中id
  rownames(Gmat_out) <- colnames(Gmat_out) <- c(idA$V1, idB$V1)
} else {
  if (!is.null(opt$phef)) {
    ## 输出的基因型矩阵中1、2列为行/列号，对表型文件中id进行重编
    phe <- fread(opt$phef)  ## 默认第一列为id
    names(phe)[1] <- 'org'
    mtab <- data.frame(org=c(idA$V1, idB$V1), new=1:nIndT)
    pcol <- ncol(phe)
    phe2 <- merge(phe, mtab, by = 'org', sort = FALSE)
    phe <- subset(phe2, select = c(pcol+1, 2:pcol, 1))
    fwrite(phe, paste0(basename(opt$phef), '.new'), col.names = FALSE, sep=' ')
    fwrite(mtab, 'grm_id_match.txt', sep=' ')
    cat('The ID in phenotype file has been replaced\n')
  }
}

## 将矩阵转换成三列形式(下三角)
Gmat_out[upper.tri(Gmat_out)] <- NA
Gmat3col <- melt(Gmat_out, na.rm = TRUE)

## 是否输出矩阵的log(determinant)
if (!is.null(opt$logdet)) {
  cat('Calculating determinant of matrix...\n')
  logdet <- determinant(Gmat_pd)$modulus
  row1add <- matrix(c(0, 0, logdet), nrow = 1, ncol = 3)
  if (is.na(row1add[,3])) {
    cat('The determinant of matrix is NA, and log(det) is not output\n')
  } else {
    colnames(row1add) <- colnames(Gmat3col)
    Gmat3col <- rbind(row1add, Gmat3col)
  } 
}

## 写出矩阵
cat('Writing out matrix...\n')
fwrite(Gmat3col, opt$out, sep=' ', col.names = header)

cat('G matrix has been exported to: ', opt$out, '\n')
