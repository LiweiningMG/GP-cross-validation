#!/bin/bash
#SBATCH --job-name=accuracy

########################################################################################################################
## 版本: 1.0.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-06-03
## 
## 根据指定的参数，生成QMSim软件所需的参数文件*.prm
## 
## 使用: QMSim_genome_process.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long proj:,breeds:,rep:,nsnp:,geno_gen:,geno_sel:,binDiv:,binThr:,maf:,nginds:,litters:,females:,code:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必要参数]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'A B C' ["A B"]
    --rep )      rep="$2";      shift 2 ;; ## 需要处理的是QMSim输出的第几次重复模拟，如lm_mrk_001.txt中的1 [1]
    --nsnp )     nsnp="$2";     shift 2 ;; ## 需要选择的标记个数 [50000]
    --geno_gen ) geno_gen="$2"; shift 2 ;; ## 输出基因型个体的世代 [8-10]
    --geno_sel ) geno_sel="$2"; shift 2 ;; ## 从哪几个世代中挑选基因型个体 [geno_gen]
    --binDiv )   binDiv="$2";   shift 2 ;; ## 抽样SNP时区间划分的依据，pos/frq [pos]
    --binThr )   binThr="$2";   shift 2 ;; ## 抽样SNP时区间划分长度，物理位置cM或基因频率步长 [10]
    --maf )      maf="$2";      shift 2 ;; ## 抽样SNP时允许的最小等位基因频率 [10]
    --nginds )   nginds="$2";   shift 2 ;; ## 每个品种选择的基因型个体数 ["600 600 ..."]
    --last_litters )  litters="$2";  shift 2 ;; ## 各品种在LD稳定阶段的每窝个体数 ["10 10 ..."]
    --last_females )  females="$2";  shift 2 ;; ## 各品种最后一个阶段的群体中的雌性个体数 ["500 500 ..."]
    --code )     code="$2";     shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --out )      out="$2";      shift 2 ;; ## 最终输出基因型文件前缀 [merge]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 检查必要参数是否提供
if [[ ! -d ${proj} ]]; then
  echo "project path ${proj} not found! "
  exit 1
fi

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
func=${code}/shell/function.sh
mrk_sel=${code}/R/QMSim_mrk_select.R
gind_sel=${code}/R/geno_individuals_select.R
PCA_cal=${code}/shell/pca_multi_pop.sh
geno_dist=${code}/shell/distance_multi_pop.sh
PCA_plot=${code}/R/PCA_plot.R
LD_cor=${code}/R/LD_decay_plot.R
# block_LD=${code}/R/block_LD_cor.R
corr_cal=${code}/R/columns_correlation.R

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 加载自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink QMSim_selected

## 检查需要的脚本文件是否存在且具有执行权限
check_command $gind_sel $PCA_cal $geno_dist $PCA_plot $LD_cor $mrk_sel $corr_cal

## 默认参数
out=${out:=merge}
windows=${windows:="10000"}
r2=${r2:="0"}
inter=${inter:="99999"}
nsnp=${nsnp:="50000"}
geno_gen=${geno_gen:="8-10"}
maf=${maf:="0.01"}
binDiv=${binDiv:="pos"}
binThr=${binThr:="10"}
rep=${rep:="1"}
breeds=${breeds:="A B"}
geno_sel=${geno_sel:="${geno_gen}"}

## 获取品种个数
read -ra breeds <<<"$breeds"
np=${#breeds[@]}

## 由品种数确定的默认参数
nginds=${nginds:=$(printf "%${np}s" | sed "s/ /600 /g" | sed 's/ *$//')}
litters=${litters:=$(printf "%${np}s" | sed "s/ /10 /g" | sed 's/ *$//')}
females=${females:=$(printf "%${np}s" | sed "s/ /200 /g" | sed 's/ *$//')}

## 解析参数
read -ra nginds <<<"$nginds"
read -ra litters <<<"$litters"
read -ra females <<<"$females"
qrep=$(printf "%03d" "$rep")
mapfile -t gid_gen < <(seq ${geno_gen:0:1} ${geno_gen:2:3})

## 工作路径
cd ${proj} || exit

## 从QMSim中获取第一个随机数为后面随机过程的种子
seed=$(sed -n '2p' seed | awk '{print $2}')
seed=${seed:0:6}
echo ${seed} >random.seed

## 选择标记
$mrk_sel \
  --popN "${breeds[*]}" \
  --freqf "%pop%_freq_mrk_${qrep}.txt" \
  --mapf "lm_mrk_${qrep}.txt" \
  --nsel ${nsnp} \
  --binDiv ${binDiv} \
  --binThr ${binThr} \
  --maf ${maf} \
  --out mrk_sel_index.txt \
  --outmapf ${breeds[0]}.map \
  --seed ${seed}

## 选择个体和提取标记信息
for j in $(seq 0 $((np - 1))); do
  ## 系谱文件
  sed '1d' ${breeds[${j}]}_data_${qrep}.txt | awk '{print $1,$2,$3}' >${breeds[${j}]}_pedi.txt

  ## 选择基因型个体
  $gind_sel \
    --dataf ${breeds[${j}]}_data_${qrep}.txt \
    --gen_all ${geno_gen} \
    --gen_sel ${geno_sel} \
    --nsel ${nginds[${j}]} \
    --outIndex \
    --out ${breeds[${j}]}_Ind_sel_index.txt

  ## map文件
  [[ ${j} != "0" ]] && cp ${breeds[0]}.map ${breeds[${j}]}.map

  ## 基因型个体数
  ngid=$((${#gid_gen[@]} * ${litters[${j}]} * ${females[${j}]}))

  ## ped文件
  QMSim_selected \
    --indexf mrk_sel_index.txt \
    --indIndexf ${breeds[${j}]}_Ind_sel_index.txt \
    --mrkf ${breeds[${j}]}_mrk_${qrep}.txt \
    --nInd ${ngid} \
    --fid ${breeds[${j}]} \
    --out ${breeds[${j}]}.ped &
done
wait

## 质控
for b in "${breeds[@]}"; do
  plink --file ${b} --maf 0.01 --make-bed --out ${b}q
done

## 筛选出在所有品种中均通过质控的标记位点
awk '{print $2}' "${breeds[@]/%/q.bim}" | sort | uniq -c | awk -v n=${np} '$1==n {print $2}' >common.snp
for b in "${breeds[@]}"; do
  plink --file ${b} --extract common.snp --make-bed --out ${b}q
done

## pca计算
qc_files=("${breeds[@]/%/q}")
$PCA_cal --pre_list "${qc_files[*]}" --fids "${breeds[*]}" --fid

## pca作图
$PCA_plot \
  --eigv "$(printf '%s_' "${qc_files[@]}")pca.txt" \
  --out "$(printf '%s_' "${breeds[@]}")pca"

## LD计算
for b in "${breeds[@]}"; do
  plink \
    --bfile ${b}q \
    --freq \
    --r2 \
    --ld-window-kb ${windows} \
    --ld-window ${inter} \
    --ld-window-r2 ${r2} \
    --out ${b}q
done

## LD结果统计、作图
$LD_cor \
  --files "$(printf '%sq ' "${breeds[@]}")" \
  --popN "${breeds[*]}" \
  --bin1 50 \
  --breaks 1000 \
  --bin2 100 \
  --max 5000 \
  --out g_"$(printf '%s_' "${breeds[@]}")"_5Mb

## 品种对之间的基因频率和LD相关性大小
for t in frq ld; do
  :>${t}_cor.txt
  for bi in $(seq 0 $((np - 1))); do
    for bj in $(seq ${bi} $((np - 1))); do
      [[ ${bi} == "${bj}" ]] && continue
      cor=$($corr_cal --file1 ${breeds[${bi}]}q.${t} --file2 ${breeds[${bj}]}q.${t})
      echo "${breeds[${bi}]} ${breeds[${bj}]} ${cor}" >>${t}_cor.txt
    done
  done
  echo "correlation of ${t}:"
  cat ${t}_cor.txt
done

## 合并所有品种的基因型
: >plink_merge_list.txt
for b in "${breeds[@]}"; do
  [[ ${b} == "${breeds[0]}" ]] && continue
  echo "${b}q" >>plink_merge_list.txt
done
plink \
  --bfile ${breeds[0]}q \
  --merge-list plink_merge_list.txt \
  --maf 0.05 \
  --make-bed \
  --out ${out}

## 群体间的遗传距离
$geno_dist --bfile ${out} --out ${out}.dist.summ

## debug
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Two/rep3
breeds=(A B C)
geno_gen=8-10
maf=0.01
