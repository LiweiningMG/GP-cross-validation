#!/usr/bin/bash

## 对提供的任意多个群体进行平均遗传距离计算
## 需要用到的软件：R、Plink1.9
## 需要提供的信息，各个群体的plink文件前缀（plink文件需转为fam、bim、bed的二进制格式）
## 存在不在当前工作路径的文件时，提供全路径，但不用带后缀(.fam)，如/home/user/pop1

## 用法：./distance_multi_pop.sh.sh -P "prefix1 prefix2 ..." -N "N1 N2 ..."
# -P 需要作图的群体的plink文件路径（需为二进制文件，且路径中文件名不用加.fam .bim .bed等后缀）
# -N 每个群体在PCA作图中的群体标识符 (两个参数均需加双引号)
# -O 输出文件名
# -h 帮助

###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long bfile:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --bfile )   bfile="$2";  shift 2 ;;  ## plink文件前缀，如"/public/home/popA /public/home/popB"
    --out )     out="$2" ;   shift 2 ;;  ## 输出文件名
    -h | --help )     grep " shift " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

## 默认参数
out=${out:=distance}

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code="${script_path%%code*}code"
fi

## 脚本
dist=${code}/R/mean_distance.R
func=${code}/shell/function.sh

## 载入自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查plink二进制文件是否存在
check_plink ${bfile}

## 工作文件夹
dir=$(dirname ${bfile})
cd ${dir} || exit

## 个体按照FID排序
plink --bfile ${bfile} --indiv-sort natural --make-bed --out distance_tmp

## 计算距离矩阵
plink --bfile distance_tmp --distance square0 1-ibs --out distance_tmp

## 计算矩阵中品种配对对应的块的均值
$dist --prefix distance_tmp --out ${out}

## 删除中间文件
rm distance_tmp.*

## debug
# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Xie2021
bfile="/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Xie2021/Genotype.id.qc.fam"
