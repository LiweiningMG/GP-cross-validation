#!/bin/usr/bash

## liwn 2023-05-29
## 联系方式：liwn@cau.edu.cn
## 功能：基于公开（模拟）数据集验证多品种评估模型：

## 主脚本路径
code=/home/liujf/WORKSPACE/liwn/mbGS/code
GP_cross=/home/liujf/WORKSPACE/liwn/mbGS/code/GP_cross_validation.sh

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

######################## 1. 猪数据集 ########################

## Xie, L., J. Qin, L. Rao, X. Tang and D. Cui et al., 2021 Accurate prediction and genome-wide association analysis of digital intramuscular fat content in longissimus muscle of pigs. Animal Genetics 52: 633-644. https://doi.org/10.1111/age.13121

## 路径
pro=/home/liujf/WORKSPACE/liwn/mbGS/data/Xie2021
bfile=${pro}/Genotype.id.qc
phef=${pro}/phenotypes_dmu.txt
breeds=(YY LL)
traits=(PFAI MS)

## 计算校正表型
sbatch -c2 --mem=2G $GP_cross \
  --proj ${pro} \
  --breeds "${breeds[*]}" \
  --traits "${traits[*]}" \
  --bfile ${bfile} \
  --phef ${phef} \
  --code ${code} \
  --type adj

## 计算品种内评估准确性
sbatch -c26 --mem=25G $GP_cross \
  --proj ${pro} \
  --breeds "${breeds[*]}" \
  --traits "${traits[*]}" \
  --bfile ${bfile} \
  --phef ${phef} \
  --code ${code} \
  --rep 5 \
  --fold 5 \
  --type within

## 计算多品种GBLUP评估准确性
for type in blend union; do
  for trait in "${traits[@]}"; do
    sbatch -c26 --mem=25G $GP_cross \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${trait}" \
      --trait "${trait}" \
      --bfile ${bfile} \
      --phef ${phef} \
      --code ${code} \
      --suffix \
      --type ${type}
  done
done

## 计算多品种多性状Bayes评估准确性
for bin in lava cubic; do
  for trait in "${traits[@]}"; do
    sbatch -c26 --mem=26G $GP_cross \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${traits[*]}" \
      --trait "${trait}" \
      --bfile ${bfile} \
      --phef ${phef} \
      --code ${code} \
      --bin ${bin} \
      --suffix \
      --type multi
    sleep 1m
  done
done

## 准确性结果统计
$GP_cross \
  --proj ${pro} \
  --breeds "${breeds[*]}" \
  --traits "${traits[*]}" \
  --bin "fix cubic lava" \
  --code ${code} \
  --type accur

######################## 模拟数据 ########################

pro=/home/liujf/WORKSPACE/liwn/mbGS/data/Two
bfile=/home/liujf/WORKSPACE/liwn/mbGS/data/Two
breeds="A B C"
means="0.5 1.0 1.5"
h2s="0.5 0.4 0.3"

for r in 1; do # r=1;dist=identical;cor=0.2
  ## 生成模拟基因型
  sbatch -c20 --mem=100G $GP_cross \
    --type gsim \
    --code "${code}" \
    --proj "${pro}" \
    --breeds "${breeds}" \
    --nchr "5" \
    --thread "20" \
    --sim_dir "rep${r}" \
    --founder_sel "rnd,rnd,rnd" \
    --seg_sel "phen,phen /l,rnd" \
    --last_sel "rnd,rnd,rnd" \
    --seg_gens "40 40 40" \
    --extentLDs "10 10 10" \
    --last_males "40 40 40" \
    --last_females "200 200 200" \
    --geno_gen "8-10" \
    --nmloc "100000" \
    --nginds "600 600 600"

  ## 从模拟群体中筛选个体和SNP标记
  sbatch -c1 --mem=10G $GP_cross \
    --type geno \
    --code ${code} \
    --proj "${pro}/rep${r}" \
    --breeds "${breeds}" \
    --last_females "200 200 200" \
    --nginds "600 600 600" \
    --geno_gen "8-10" \
    --nsnp "14000" \
    --out merge

  ## 随机数种子
  seed=$(cat ${pro}/rep${r}/random.seed)

  for dist in identical uniform; do # dist=identical;cor=0.2;bin=fix
    for cor in 0.2 0.4; do
      proi=${pro}/rep${r}/${dist}/cor${cor}
      mkdir -p ${proi}

      ## 生成模拟表型
      sbatch -c1 --mem=10G $GP_cross \
        --type psim \
        --proj ${proi} \
        --bfile ${pro}/rep${r}/merge \
        --breeds "${breeds}" \
        --code ${code} \
        --means "${means}" \
        --h2s "${h2s}" \
        --rg_sim "${cor}" \
        --rg_dist ${dist} \
        --seed ${seed}

      ## 计算品种内评估准确性
      sbatch -c26 --mem=50G $GP_cross \
        --type within \
        --proj ${proi} \
        --breeds "${breeds}" \
        --bfile ${pro}/rep${r}/merge \
        --phef ${proi}/pheno_sim.txt \
        --code ${code} \
        --tbv_col 6 \
        --seed ${seed} \
        --rep 5 \
        --fold 5

      ## 多品种评估准确性(GBLUP)
      for type in blend union; do
        sbatch -c26 --mem=25G $GP_cross \
          --type ${type} \
          --proj ${proi} \
          --breeds "${breeds}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --code ${code} \
          --seed ${seed} \
          --tbv_col 6
        sleep 30
      done

      ## 计算多品种多性状Bayes评估准确性
      for bin in fix lava cubic; do
        sbatch -c40 --mem=80G $GP_cross \
          --type multi \
          --proj ${proi} \
          --breeds "${breeds}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --thread 20 \
          --code ${code} \
          --seed ${seed} \
          --bin ${bin}
        sleep 30
      done
    done
  done
done
