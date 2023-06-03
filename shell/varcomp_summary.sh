#!/bin/bash
#SBATCH --job-name=accur_GBLUP

########################################################################################################################
## 版本: 1.0.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-05-30
## 
## 统计各种情形下的方差组分结果
## 
## 使用: ./varcomp_summary.sh --proj ...
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long proj:,breeds:,traits:,bin:,dirPre:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必要参数]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'YY DD' [必要参数]
    --traits )   traits="$2";   shift 2 ;; ## 性状名称，如"DF DPM" [NULL]
    --dirPre )   dirPre="$2";   shift 2 ;; ## JWAS输出文件夹增加的前缀 [NULL]
    --bin )      bin="$2";      shift 2 ;; ## 多品种评估时区间划分方法，fix/frq/ld/lava/cubic ["fix lava cubic"]
    --out )      out="$2";      shift 2 ;; ## 准确性输出文件名 [accuracy_$date.txt]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 检查必要参数是否提供
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
  exit 1
elif [[ ${breeds} ]]; then
  echo "${breeds} is reduired! "
  exit 1
fi

## 时间
today=$(date +%Y%m%d)

## 默认参数
out=${out:=VarComp_${today}.txt}
bin=${bin:="fix lava cubic"}

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## DMU参数卡名称
DIR_full=phe_adj_PBLUP
DIR_within=within
DIR_blend=blend
DIR_union=union

##############  方差组分整理  ##############
###########################################
echo "模型 品种 性状 交叉重复 交叉折数 类型 值" >${proj}/VarComp_${today}.txt
for traiti in "${traits[@]}"; do # traiti=${traits[0]};b=${breeds[0]};r=1;f=1;bin=fix
  path=${proj}/${traiti}
  [[ ! -d ${path} ]] && continue

  ## 品种内遗传力
  for b in "${breeds[@]}"; do
    ## 交叉验证参数
    rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
    fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

    ## 每个子集方差组分
    for r in $(seq 1 ${rep}); do # r=1;f=1
      for f in $(seq 1 ${fold}); do
          ## GBLUP
          lst=${path}/${b}/val${f}/rep${r}/${DIR_within}.lst
          [[ ! -s ${lst} ]] && continue
          h2_within=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
          echo "w-GBLUP ${b} ${traiti} ${r} ${f} h2 ${h2_within}" >>${proj}/VarComp_${today}.txt

          ## Bayes
          ebvf=${path}/${b}/val${f}/rep${r}/EBV_fix_y1.txt
          [[ ! -s ${ebvf} ]] && continue
          varf=${path}/${b}/val${f}/rep${r}/var_fix.txt
          varg=$(sed '1d' ${ebvf} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
          vare=$(tail -n 1 ${varf} | awk -v col=$(((i + 1) * (i + 1))) '{print $col}')
          h2_Bayes=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
          echo "w-fix ${b} ${traiti} ${r} ${f} h2 ${h2_Bayes}" >>${proj}/VarComp_${today}.txt
      done
    done
  done

  ## 品种组合
  mapfile -t blends < <(find ${path} -name "blend_*" -type d)
  finished=""

  for t in "${blends[@]}"; do # t=${blends[0]}
    paths=$(basename ${t})
    types=${paths/blend_/}
    IFS='_' read -r -a breeds_sub <<<"$types"
    n=${#breeds_sub[@]}

    ## 组合中每个品种
    for i in "${!breeds_sub[@]}"; do # i=0
      ## 子集
      for r in $(seq 1 ${rep}); do # r=1;f=1
        for f in $(seq 1 ${fold}); do
          ## 单性状GBLUP联合评估
          lst=${path}/blend_${types}/val${f}/rep${r}/${DIR_blend}.lst
          if [[ -s ${lst} ]]; then
            h2_blend=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
          # else
          #   h2_blend=""
          #   echo "${lst} not found! "
          fi

          ## 双性状联合评估
          lst=${path}/union_${types}/val${f}/rep${r}/${DIR_union}.lst
          if [[ -s ${lst} ]]; then
            message="Correlation matrix for random"
            h2_union=$(grep -A $((i + 2)) 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')

            for j in $(seq $((i + 1)) $((n - 1))); do
              rg_union=$(grep -A $((j + 2)) "${message}" ${lst} | tail -n 1 | awk -v col=$((i + 2)) '{print $col}')
              echo "union ${breeds_sub[i]}_${breeds_sub[j]} ${traiti} ${r} ${f} rg ${rg_union}" >>${proj}/VarComp_${today}.txt
            done
          # else
          #   h2_union=""
          #   echo "${lst} not found! "
          fi

          for bin in ${bin}; do
            for dirPre in ${dirPre}; do
              multi_path=${path}/multi_${types}/val${f}/rep${r}
              ebvfi=$(find ${multi_path} -name "EBV_${bin}_y$((i + 1)).txt")
              varfi=$(find ${multi_path} -name "var_${bin}*.txt")
              [[ ! -s ${ebvfi} ]] && continue
              ## 遗传力
              varg=$(sed '1d' ${ebvfi} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
              vare=$(tail -n 1 ${varfi} | awk -v col=$(((i + 1) * (i + 1))) '{print $col}')
              h2_multi=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
              echo "${bin} ${breeds_sub[i]} ${traiti} ${r} ${f} h2 ${h2_multi}" >>${proj}/VarComp_${today}.txt

              ## 遗传相关
              for j in $(seq $((i + 1)) $((n - 1))); do
                ebvfj=$(find ${multi_path} -name "EBV_${bin}_y$((j + 1)).txt")
                [[ ! -s ${ebvfj} ]] && continue

                # 计算A文件中第2列减去均值的值
                meanA=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfi})
                awk -v meanA="$meanA" 'NR>1{print $2-meanA}' ${ebvfi} >A.tmp
                # 计算B文件中第2列减去均值的值
                meanB=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfj})
                awk -v meanB="$meanB" 'NR>1{print $2-meanB}' ${ebvfj} >B.tmp
                # 计算协方差
                cov=$(paste A.tmp B.tmp | awk '{sum+=($1*$2)}END{print sum/(NR - 1)}')
                var2=$(sed '1d' ${ebvfj} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
                rg_multi=$(echo "scale=4; $cov / sqrt($varg * $var2)" | bc | xargs printf "%.4f")

                echo "${bin} ${breeds_sub[i]}_${breeds_sub[j]} ${traiti} ${r} ${f} rg ${rg_multi}" >>${proj}/VarComp_${today}.txt
              done
            done
          done

          ## 写出到文件
          {
            echo "blend ${breeds_sub[i]} ${traiti} ${r} ${f} h2 ${h2_blend}"
            echo "union ${breeds_sub[i]} ${traiti} ${r} ${f} h2 ${h2_union}"
          } >>${proj}/VarComp_${today}.txt
        done
      done
    done
  done
done

## 去除na值行
sed -i '/ $/d' ${proj}/VarComp_${today}.txt
sed -i 's/ /\t/g' ${proj}/VarComp_${today}.txt

rm A.tmp B.tmp
