#!/bin/bash
#SBATCH --job-name=accuracy

########################################################################################################################
## 版本: 1.0.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-05-30
## 
## 统计各种情形下的交叉验证准确性结果
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
TEMP=$(getopt -o h --long code:,proj:,breeds:,traits:,bin:,soft:,dirPre:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必要参数]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'YY DD' [必要参数]
    --traits )   traits="$2";   shift 2 ;; ## 性状名称，如"DF DPM" [NULL]
    --dirPre )   dirPre="$2";   shift 2 ;; ## JWAS输出文件夹增加的前缀 [""]
    --bin )      bins="$2";     shift 2 ;; ## 多品种评估时区间划分方法，fix/frq/ld/lava/cubic ["fix lava cubic"]
    --soft)      software="$2"; shift 2 ;; ## 育种值估计程序，可为JWAS/C [JWAS]
    --code )     code="$2";     shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
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
elif [[ ! ${breeds} ]]; then
  echo "para --breeds is reduired! "
  exit 1
fi

## 时间
today=$(date +%Y%m%d)

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code="${script_path%%code*}code"
fi

## 路径/脚本
accur_cal=${code}/R/accuracy_bias_calculation.R
[[ ! -s ${accur_cal} ]] && echo "Error: ${accur_cal} not found! " && exit 5

## 默认参数
out=${out:=${proj}/accuracy_${today}.txt}
bins=${bins:="fix lava cubic"}
dirPre=${dirPre:=""}
software=${software:="C"}

## 解析参数
read -ra breeds_array <<<"$breeds"
read -ra traits_array <<<"$traits"
read -ra bins_array <<<"$bins"

##############  准确性结果统计(各种组合)  ##########
################################################
echo "模型 参考群 品种 性状 重复 交叉折数 准确性 无偏性 秩相关 验证群大小" >${out}
for traiti in "${traits_array[@]}"; do # traiti=${traits_array[0]};b=${breeds_array[0]}
  for b in "${breeds_array[@]}"; do
    path=${proj}/${traiti}
    ## 交叉验证参数
    rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
    fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

    ## 群体内
    accf_within=${path}/${b}/accur_GBLUP.txt
    accf_BayesAS=${path}/${b}/accur_BayesAS.txt
    # if [[ -s ${accf_within} ]]; then
    #   sed '1d' ${accf_within} | awk '{print "within","'${b}'","'${b}'","'${traiti}'", $0}' >>${out}
    # # else
    #   # echo "${accf} not found! "
    # fi

    ## 群体合并单性状
    accf=$(find ${path}/blen* -name "accur_*${b}*txt")
    # accf=${path}/blend/accur_GBLUP_${b}.txt
    for accfi in ${accf}; do
      if [[ -s ${accfi} ]]; then
        ## 文件夹名
        type=$(dirname ${accfi})
        type=$(basename ${type})
        type=${type/blend_/}
        {
          awk '{print "ST-GBLUP","'${type}'","'${b}'","'${traiti}'",$0}' ${accfi}
          [[ -s ${accf_within} ]] && \
            awk '{print "GBLUP","'${type}'","'${b}'","'${traiti}'", $0}' ${accf_within}
          [[ -s ${accf_BayesAS} ]] && \
            awk '{print "BayesAS","'${type}'","'${b}'","'${traiti}'", $0}' ${accf_BayesAS}
        } >>${out}
      # else
        # echo "${accf} not found! "
      fi
    done

    ## 群体合并多性状GBLUP
    accf=$(find ${path}/unio* -name "accur_*${b}*txt")
    # accf=${path}/union/accur_GBLUP_${b}.txt
    for accfi in ${accf}; do
      if [[ -s ${accfi} ]]; then
        ## 文件夹名
        type=$(dirname ${accfi})
        type=$(basename ${type})
        type=${type/union_/}
        {
          awk '{print "MT-GBLUP","'${type}'","'${b}'","'${traiti}'",$0}' ${accfi}
          [[ -s ${accf_within} ]] && \
            awk '{print "GBLUP","'${type}'","'${b}'","'${traiti}'", $0}' ${accf_within}
          [[ -s ${accf_BayesAS} ]] && \
            awk '{print "BayesAS","'${type}'","'${b}'","'${traiti}'", $0}' ${accf_BayesAS}
        } >>${out}
      else
        echo "${accf} not found! "
      fi
    done

    ## MT-bayesAS
    for soft in ${software}; do # soft=C;bin=fix;dirPre=""
      for bin in "${bins_array[@]}"; do
        # accf=${path}/multi/accur_bayes_${dirPre}${soft}_${bin}_${b}.txt
        accf=$(find ${path}/mult* -name "accur*${bin}*${b}*txt")
        for accfi in ${accf}; do # accfi=${accf[0]}
          type=$(dirname ${accfi})
          type=$(basename ${type})
          type=${type/multi_/}
          type="${dirPre}${type}"
          # type="${dirPre}${soft}_${bin}"
          if [[ -s ${accfi} ]]; then
            {
              awk '{print "'${bin}'","'${type}'","'${b}'","'${traiti}'",$0}' ${accfi}
              [[ -s ${accf_within} ]] && \
                awk '{print "GBLUP","'${type}'","'${b}'","'${traiti}'", $0}' ${accf_within}
              [[ -s ${accf_BayesAS} ]] && \
                awk '{print "BayesAS","'${type}'","'${b}'","'${traiti}'", $0}' ${accf_BayesAS}
            } >>${out}
          else
            echo "${accfi} not found! "
          fi
        done
      done
    done
  done
done
## 去掉首位行
sed -i '/rep/d' ${out}
sed -i '/mean/d' ${out}
## 替换空格
sed -i 's/ /\t/g' ${out}

## debug
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Xie2021
breeds="YY LL"
traits="PFAI MS"
