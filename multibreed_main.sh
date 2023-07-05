#!/bin/usr/bash

########################################################################################################################
## 版本: 1.1.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 功能：
## 基于公开（模拟）数据集验证多品种评估模型：
## 
## 使用: multibreed_main.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################



## 主脚本路径
code=/work/home/ljfgroup01/WORKSPACE/liwn/code/GitHub/GP-cross-validation
GP_cross=${code}/shell/GP_cross_validation.sh
combination=${code}/R/array_combination.R

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH


###########################################################################################
## 
##  公开数据集
## 
###########################################################################################


######################## 数据集1 - 猪 ########################

## Xie, L., J. Qin, L. Rao, X. Tang and D. Cui et al., 2021 Accurate prediction and genome-wide association analysis of digital intramuscular fat content in longissimus muscle of pigs. Animal Genetics 52: 633-644. https://doi.org/10.1111/age.13121

## 路径(需根据实际数据文件路径修改参数)
pro=/work/home/ljfgroup01/WORKSPACE/liwn/mbGS/Real/Xie2021
bfile=${pro}/data/Genotype.id.qc
phef=${pro}/data/phenotypes_dmu.txt
breeds=(YY LL)
traits_all=(PFAI MS)
traits_cal=(PFAI MS)


######################## 数据集2 - 菜豆 ########################

## KELLER B, ARIZA-SUAREZ D, PORTILLA-BENAVIDES A E, et al. Improving Association Studies and Genomic Predictions for Climbing Beans With Data From Bush Bean Populations[J]. Frontiers in Plant Science, 2022,13

## 路径(需根据实际数据文件路径修改参数)
pro=/work/home/ljfgroup01/WORKSPACE/liwn/mbGS/Real/Keller2022
bfile=${pro}/data/merge
phef=${pro}/data/adj_dmu.txt
breeds=(VEC ADP AxM MIP VEF)
traits_all=(100SdW DF DPM SdFe Yield)
traits_cal=(100SdW DF Yield)


######################## 数据分析 ########################

## 随机数种子
if [[ ! -s "${pro}/random.seed" ]]; then
  seed=$RANDOM
  echo ${seed} >"${pro}/random.seed"
else
  seed=$(cat "${pro}/random.seed")
fi

## 计算校正表型
for trait in "${traits_cal[@]}"; do
  ## 需根据Linux服务器特性和需求调整参数
  sbatch -c2 -p XiaoYueHe --mem=4G $GP_cross \
    --proj ${pro} \
    --breeds "${breeds[*]}" \
    --traits "${traits_all[*]}" \
    --trait "${trait}" \
    --bfile ${bfile} \
    --phef ${phef} \
    --code ${code} \
    --type adj
  sleep 5
done

## 计算品种内评估准确性GBLUP
# tbv_col=" --tbv_col same "
for trait in "${traits_cal[@]}"; do
  sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
    --proj ${pro} \
    --breeds "${breeds[*]}" \
    --traits "${traits_all[*]}" \
    --trait "${trait}" \
    --bfile ${bfile} \
    --phef ${phef} \
    --seed ${seed} \
    ${tbv_col} \
    --code ${code} \
    --rep 5 \
    --fold 5 \
    --type within
  sleep 10
done

## 品种的各种可能组合
$combination --array "${breeds[*]}" --label VEC --max 3 --out ${pro}/breeds_combination.txt
echo "VEC ADP AxM MIP VEF" >>${pro}/breeds_combination.txt

## 计算多品种GBLUP评估准确性
for trait in "${traits_cal[@]}"; do
  for type in blend union; do
    while IFS= read -r breed_comb; do
      # 判断是否已完成分析
      [[ -d ${pro}/${trait}/${type}_${breed_comb// /_} ]] && \
        continue

      # echo "${pro} ${trait} ${type} ${breed_comb}"
      sbatch -c25 -p XiaoYueHe --mem=80G $GP_cross \
        --type ${type} \
        --proj ${pro} \
        --breeds "${breed_comb}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --code ${code} \
        --phef ${phef} \
        ${tbv_col} \
        --seed ${seed} \
        --thread 26 \
        --suffix
      sleep 20
    done <"${pro}/breeds_combination.txt"
  done
done

## 计算多品种多性状Bayes评估准确性
for trait in "${traits_cal[@]}"; do
  while IFS= read -r breed_comb; do
    for bin in fix lava cubic; do
      ## 判断是否已完成分析
      [[ -s ${pro}/${trait}/multi_${breed_comb// /_}/val1/rep1/EBV_${bin}_y1.txt ]] && \
        continue

      sbatch -c25 -p XiaoYueHe --mem=80G $GP_cross \
        --type multi \
        --proj ${pro} \
        --breeds "${breed_comb}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --code ${code} \
        ${tbv_col} \
        --seed ${seed} \
        --phef ${phef} \
        --thread 26 \
        --bin ${bin} \
        --suffix
      sleep 10
    done
  done <"${pro}/breeds_combination.txt"
done

## 计算品种内评估准确性BayesAS
tbv_col=""
## 菜豆数据集不用计算校正表型
[[ ${pro} =~ "Keller" ]] && tbv_col=" --tbv_col same "
for trait in "${traits_cal[@]}"; do
  ## 各个参考群组合中的fix类型文件应该一致，找出第一个符合的即可
  binf=$(find ${pro} -name fix_100.txt -type f -print -quit)
  [[ ! -f ${binf} ]] && continue
  for b in "${breeds[@]}"; do
    sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
      --proj ${pro} \
      --breeds "${b}" \
      --traits "${traits_all[*]}" \
      --trait "${trait}" \
      --bfile ${bfile} \
      --phef ${phef} \
      --seed ${seed} \
      --method BayesAS \
      --binf ${binf} \
      ${tbv_col} \
      --code ${code} \
      --rep 5 \
      --fold 5 \
      --type within \
      --out accur_BayesAS.txt
    sleep 10
  done
done

## 统计准确性和方差组分结果
for type in accur var; do # type=var
  $GP_cross \
    --type ${type} \
    --proj ${pro} \
    --breeds "${breeds[*]}" \
    --bin "fix lava cubic" \
    --traits "${traits_cal[*]}" \
    --code ${code}
done


###########################################################################################
## 
##  模拟数据集
## 
###########################################################################################


######################## QMSim模拟(大+小) ########################
## 路径需根据实际需要修改
pro=/work/home/ljfgroup01/WORKSPACE/liwn/mbGS/QMSim/Two
breed_sim="A B"
means="1.0 0.5"
h2s="0.5 0.3"

## 基因型数据模拟
sbatch -c40 -p XiaoYueHe --mem=100G $GP_cross \
  --type gsim \
  --code "${code}" \
  --proj "${pro}" \
  --breeds "${breed_sim}" \
  --nchr "18" \
  --thread "40" \
  --sim_dir "$(seq -s ' ' -f "rep%.0f" 1 10)" \
  --founder_sel "rnd,rnd" \
  --seg_sel "phen /h,phen /l" \
  --seg_gens "40 10" \
  --last_sel "rnd,rnd" \
  --extentLDs "10 10" \
  --last_males "100 20" \
  --last_females "500 100" \
  --geno_gen "8-10" \
  --nmloc "200000"

for r in {1..10}; do # r=6;dist=uniform;cor=0.8;bin=lava
  ## 从模拟群体中筛选个体和SNP标记
  [[ ! -s ${pro}/rep${r}/merge.fam ]] && \
    $GP_cross \
      --type geno \
      --code ${code} \
      --proj "${pro}/rep${r}" \
      --breeds "${breed_sim}" \
      --last_females "500 100" \
      --nginds "3000 600" \
      --binDiv "pos" \
      --binThr "10" \
      --geno_gen "8-10" \
      --nsnp "50000" \
      --out "merge"

  ## 随机数种子
  seed=$(cat ${pro}/rep${r}/random.seed)

  for dist in identical uniform; do # dist=identical;cor=0.2;bin=cubic
    for cor in 0.2 0.5 0.8; do
      proi=${pro}/rep${r}/${dist}/cor${cor}
      mkdir -p ${proi}

      ## 生成模拟表型
      # sbatch -c1 --mem=10G \
      $GP_cross \
        --type psim \
        --proj ${proi} \
        --bfile ${pro}/rep${r}/merge \
        --breeds "${breed_sim}" \
        --code ${code} \
        --means "${means}" \
        --h2s "${h2s}" \
        --rg_sim "${cor}" \
        --rg_dist ${dist} \
        --nqtl "400" \
        --nsnp_cor "10" \
        --nbin_cor "10" \
        --nsnp_sim "50" \
        --seed ${seed}
      sleep 10

      ## 计算品种内评估准确性GBLUP
      for b in ${breed_sim}; do
        [[ -s "${proi}/${b}/accur_GBLUP.txt" ]] && continue
        sbatch -c25 -p XiaoYueHe --mem=100G \
          $GP_cross \
            --type within \
            --proj ${proi} \
            --breeds "${b}" \
            --bfile ${pro}/rep${r}/merge \
            --phef ${proi}/pheno_sim.txt \
            --code ${code} \
            --thread 26 \
            --debug \
            --tbv_col 6 \
            --seed ${seed} \
            --rep 5 \
            --fold 5
      done

      ## 等待品种内评估执行完毕
      while [[ ! -f "${proi}/${breed_sim##* }/val5/rep5/pheno.txt" ]]; do
        sleep 3
      done

      ## 多品种评估准确性(GBLUP)
      for type in blend union; do
        [[  -d ${proi}/${type}_${breed_sim// /_} ]] && continue
        sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
          --type ${type} \
          --proj ${proi} \
          --breeds "${breed_sim}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --code ${code} \
          --suffix \
          --thread 26 \
          --seed ${seed} \
          --tbv_col 6
        sleep 20
      done

      ## 计算多品种多性状Bayes评估准确性
      for bin in fix lava cubic; do
        [[ -s ${proi}/multi_${breed_sim// /_}/val1/rep1/EBV_${bin}_y1.txt ]] && continue
        sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
          --type multi \
          --proj ${proi} \
          --breeds "${breed_sim}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --thread 26 \
          --code ${code} \
          --seed ${seed} \
          --debug \
          --suffix \
          --tbv_col 6 \
          --bin ${bin}
        sleep 20
      done

      ## 计算品种内评估准确性BayesAS
      binf=${proi}/multi_${breed_sim// /_}/fix_100.txt
      [[ ! -f ${binf} ]] && continue
      for b in ${breed_sim}; do
        sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
          --type within \
          --proj ${proi} \
          --breeds "${b}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --code ${code} \
          --binf ${binf} \
          --method BayesAS \
          --thread 26 \
          --tbv_col 6 \
          --seed ${seed} \
          --rep 5 \
          --fold 5
        sleep 10
      done
    done
  done
done

## 统计准确性结果
for type in accur var; do # type=accur
  $GP_cross \
    --type ${type} \
    --proj ${pro} \
    --rep "$(seq -s " " 1 10)" \
    --rg_dist "identical uniform" \
    --rg_sim "0.2 0.5 0.8" \
    --bin "fix lava cubic" \
    --breeds "${breed_sim}" \
    --code ${code}
done


######################## QMSim模拟(三个小数量品种) ########################
## 路径需根据实际需要修改
pro=/work/home/ljfgroup01/WORKSPACE/liwn/mbGS/QMSim/Three
breed_sim="A B C"
means="1.5 1.0 0.5"
h2s="0.5 0.4 0.3"

## 生成模拟基因型
sbatch -c40 -p XiaoYueHe --mem=100G $GP_cross \
  --type gsim \
  --code "${code}" \
  --proj "${pro}" \
  --breeds "${breed_sim}" \
  --nchr "18" \
  --thread "40" \
  --sim_dir "$(seq -s ' ' -f "rep%.0f" 1 10)" \
  --founder_sel "rnd,rnd,rnd" \
  --seg_sel "rnd,phen /h,phen /l" \
  --seg_gens "70 40 10" \
  --last_sel "rnd,rnd,rnd" \
  --extentLDs "10 10 10" \
  --last_males "20 20 20" \
  --last_females "100 100 100" \
  --geno_gen "8-10" \
  --nmloc "200000"

# breed_sim="A B C"
for r in {1..10}; do # r=6;dist=identical;cor=0.5;bin=cubic
  ## 从模拟群体中筛选个体和SNP标记
  # sbatch -c1 --mem=10G \
  [[ ! -s ${pro}/rep${r}/merge.fam ]] && \
    $GP_cross \
      --type geno \
      --code ${code} \
      --proj "${pro}/rep${r}" \
      --breeds "${breed_sim}" \
      --last_females "100 100 100" \
      --nginds "600 600 600" \
      --binDiv "pos" \
      --binThr "10" \
      --geno_gen "8-10" \
      --nsnp "50000" \
      --out "merge"

  ## 随机数种子
  seed=$(cat ${pro}/rep${r}/random.seed)

  for dist in identical uniform; do # dist=identical;cor=0.5;bin=cubic
    for cor in 0.5; do
      proi=${pro}/rep${r}/${dist}/cor${cor}
      mkdir -p ${proi}

      ## 生成模拟表型
      # sbatch -c1 --mem=10G \
      [[ ! -s ${proi}/pheno_sim.txt ]] && \
      $GP_cross \
        --type psim \
        --proj ${proi} \
        --bfile ${pro}/rep${r}/merge \
        --breeds "${breed_sim}" \
        --code ${code} \
        --means "${means}" \
        --h2s "${h2s}" \
        --rg_sim "${cor}" \
        --rg_dist ${dist} \
        --nqtl "400" \
        --nsnp_cor "10" \
        --nbin_cor "10" \
        --nsnp_sim "50" \
        --seed ${seed}
      sleep 10

      ## 计算品种内评估准确性
      for b in ${breed_sim}; do
        [[ -s "${proi}/${b}/accur_GBLUP.txt" ]] && continue
        sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
          --type within \
          --proj ${proi} \
          --breeds "${b}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --code ${code} \
          --thread 26 \
          --tbv_col 6 \
          --seed ${seed} \
          --rep 5 \
          --fold 5
      done

      ## 等待最后一个品种内评估执行完毕
      while [[ ! -f "${proi}/${breed_sim##* }/val1/rep1/pheno.txt" ]]; do
        sleep 3
      done

      ## 包含A的所有可能参考群组合
      $combination --array "${breed_sim}" --label A --out ${pro}/breeds_combination.txt

      ## 多品种评估准确性(GBLUP)
      while IFS= read -r breed_comb; do
        for type in blend union; do
          [[ -d ${proi}/${type}_${breed_comb// /_} ]] && continue
          sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
            --type ${type} \
            --proj ${proi} \
            --breeds "${breed_comb}" \
            --bfile ${pro}/rep${r}/merge \
            --phef ${proi}/pheno_sim.txt \
            --code ${code} \
            --dense \
            --suffix \
            --thread 26 \
            --seed ${seed} \
            --tbv_col 6
          sleep 20
        done

        ## 计算多品种多性状Bayes评估准确性
        for bin in fix lava cubic; do
          [[ -d ${proi}/multi_${breed_comb// /_} ]] && continue
          sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
            --type multi \
            --proj ${proi} \
            --breeds "${breed_comb}" \
            --bfile ${pro}/rep${r}/merge \
            --phef ${proi}/pheno_sim.txt \
            --thread 26 \
            --code ${code} \
            --seed ${seed} \
            --suffix \
            --tbv_col 6 \
            --bin ${bin}
          sleep 20
        done
      done <"${pro}/breeds_combination.txt"

      ## 计算品种内评估准确性BayesAS
      binf=${proi}/multi_${breed_sim// /_}/fix_100.txt
      [[ ! -f ${binf} ]] && continue
      for b in $breed_sim; do
        sbatch -c25 -p XiaoYueHe --mem=100G $GP_cross \
          --type within \
          --proj ${proi} \
          --breeds "${b}" \
          --bfile ${pro}/rep${r}/merge \
          --phef ${proi}/pheno_sim.txt \
          --code ${code} \
          --binf ${binf} \
          --method BayesAS \
          --thread 26 \
          --tbv_col 6 \
          --seed ${seed} \
          --rep 5 \
          --fold 5 \
          --out accur_BayesAS.txt
        sleep 20
      done
    done
  done
done

## 统计准确性结果
for type in accur var; do # type=accur
  $GP_cross \
    --type ${type} \
    --proj ${pro} \
    --rep "$(seq -s " " 1 10)" \
    --rg_dist "identical uniform" \
    --rg_sim "0.5" \
    --bin "fix lava cubic" \
    --breeds "${breed_sim}" \
    --code ${code}
done
