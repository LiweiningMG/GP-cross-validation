#!/bin/bash
#SBATCH --job-name=accuracy

########################################################################################################################
## 版本: 1.0.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-06-02
## 
## 根据指定的参数，生成QMSim软件所需的参数文件*.prm，并运行QMSim软件进行群体模拟
## 
## 使用: run_QMSim.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long proj:,breeds:,thread:,nginds:,seg_gens:,extentLDs:,last_males:,last_females:,founder_sel:,seg_sel:,last_sel:,last_litters:,geno_gen:,nchr:,nmloc:,nqloci:,QMSim_h2:,QMSim_qtlh2:,QMSim_rep:,bottleneck:,code:,prmpath:,sim_dir:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )         proj="$2";         shift 2 ;; ## 项目目录 [必要参数]
    --breeds )       breeds="$2";       shift 2 ;; ## 群体/品种标识符，如'YY DD' ["A B"]
    --thread )       thread="$2";       shift 2 ;; ## QMSim运行时的线程数 [10]
    --nginds )       nginds="$2";       shift 2 ;; ## 各品种在输出基因型的群体中选择的基因型个体数 ["3000 300"]
    --seg_gens )     seg_gens="$2";     shift 2 ;; ## 各品种从历史群体中分离后，分离的世代数 ["40 10"]
    --extentLDs )    extentLDs="$2";    shift 2 ;; ## 各品种最后稳定LD经历的世代数 ["10 10"]
    --last_males )   last_males="$2";   shift 2 ;; ## 各品种最后一个阶段的群体中的雄性个体数 ["100 10"]
    --last_females ) last_females="$2"; shift 2 ;; ## 各品种最后一个阶段的群体中的雌性个体数 ["500 50"]
    --founder_sel )  founder_sel="$2";  shift 2 ;; ## 各品种从历史群体中选择个体的依据 ["tbv /h,tbv /l"]
    --seg_sel )      seg_sel="$2";      shift 2 ;; ## 各品种在世代选择阶段的个体选留依据 ["phen /h,phen /l"]
    --last_sel )     last_sel="$2";     shift 2 ;; ## 各品种在LD稳定阶段的个体选留依据 ["rnd,rnd"]
    --last_litters ) last_litters="$2"; shift 2 ;; ## 各品种在LD稳定阶段的每窝个体数 ["10 10"]
    --geno_gen )     geno_gen="$2";     shift 2 ;; ## 输出基因型个体的世代 [8-10]
    --nchr )         nchr="$2";         shift 2 ;; ## 染色体数目 [18]
    --nmloc )        nmloc="$2";        shift 2 ;; ## 每条染色体上的标记数 [300000]
    --nqloci )       nqloci="$2";       shift 2 ;; ## 每条染色体上的QTL数 [100]
    --QMSim_h2 )     QMSim_h2="$2";     shift 2 ;; ## 性状的广义遗传力大小 [0.3]
    --QMSim_qtlh2 )  QMSim_qtlh2="$2";  shift 2 ;; ## 性状的狭义遗传力大小 [0.3]
    --QMSim_rep )    QMSim_rep="$2";    shift 2 ;; ## QMSim模拟重复次数 [1]
    --bottleneck )   bottleneck="$2";   shift 2 ;; ## 历史群体模拟中的瓶颈阶段个体数 [250]
    --code )         code="$2";         shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [自动获取]
    --prmpath )      prmpath="$2";      shift 2 ;; ## 参数卡模板文件所在路径 [${code}/prm]
    --sim_dir )      sim_dir="$2";      shift 2 ;; ## 模拟结果文件输出文件夹名，注意该文件夹不能已存在 [rep1]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 检查必要参数是否提供
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
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

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 路径/脚本
func=${code}/shell/function.sh
gind_sel=${code}/R/geno_individuals_select.R
PCA_cal=${code}/shell/pca_multi_pop.sh
geno_dist=${code}/shell/distance_multi_pop.sh
PCA_plot=${code}/R/PCA_plot.R
LD_cor=${code}/R/LD_decay_plot.R
mrk_sel=${code}/R/QMSim_mrk_select.R
block_LD=${code}/R/block_LD_cor.R

## 参数卡模板文件
prm_hist=templete_gloabal_history.prm
prm_sub=templete_subpopulation.prm
prm_geno=templete_genome_output.prm

## 加载自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink QMSim

## 检查需要的脚本文件是否存在且具有执行权限
check_command $gind_sel $PCA_cal $geno_dist $PCA_plot $LD_cor $mrk_sel $block_LD

## 默认参数
breeds=${breeds:="A B"}
sim_dir=${sim_dir:=rep1}
QMSim_rep=${QMSim_rep:="1"}
nchr=${nchr:="18"}
nmloc=${nmloc:="300000"}
geno_gen=${geno_gen:="8-10"}
nqloci=${nqloci:="100"}
QMSim_h2=${QMSim_h2:="0.3"}
QMSim_qtlh2=${QMSim_qtlh2:="${QMSim_h2}"}
bottleneck=${bottleneck:="250"}
thread=${thread:="10"}
prmpath=${prmpath:="${code}/prm"}

## 获取品种个数
read -ra breeds <<<"$breeds"
np=${#breeds[@]}

## 由品种数确定的默认参数
founder_sel=${founder_sel:=$(printf "%${np}s" | sed "s/ /rnd,/g" | sed 's/,$//')}
seg_sel=${seg_sel:=${founder_sel}}
last_sel=${last_sel:=${founder_sel}}
last_litters=${last_litters:=$(printf "%${np}s" | sed "s/ /10 /g" | sed 's/ *$//')}
extentLDs=${extentLDs:=$(printf "%${np}s" | sed "s/ /10 /g" | sed 's/ *$//')}
last_males=${last_males:=$(printf "%${np}s" | sed "s/ /40 /g" | sed 's/ *$//')}
last_females=${last_females:=$(printf "%${np}s" | sed "s/ /200 /g" | sed 's/ *$//')}
nginds=${nginds:=$(printf "%${np}s" | sed "s/ /600 /g" | sed 's/ *$//')}
seg_gens=${seg_gens:=$(printf "%${np}s" | sed "s/ /40 /g" | sed 's/ *$//')}

## 解析参数
read -ra nginds <<<"$nginds"
read -ra seg_gens <<<"$seg_gens"
read -ra last_males <<<"$last_males"
read -ra last_females <<<"$last_females"
read -ra extentLDs <<<"$extentLDs"
read -ra last_litters <<<"$last_litters"
IFS=, read -ra founder_sel <<<"$founder_sel"
IFS=, read -ra seg_sel <<<"$seg_sel"
IFS=, read -ra last_sel <<<"$last_sel"

mapfile -t gid_gen < <(seq ${geno_gen:0:1} ${geno_gen:2:3})
gen_all=${geno_gen:2:3}
rand=$RANDOM

## 查看参数卡模板文件是否存在
if [[ ! -d ${prmpath} ]]; then
  echo "${prmpath} not exists! "
  exit 5
else
  for f in $prm_hist $prm_sub $prm_geno; do
    [[ ! -s ${prmpath}/$f ]] && echo "$f not found! " && exit 4
  done
fi

###################### 模拟群体参数卡 ######################
cd ${prmpath} || exit

## 历史群体
sed "s/%nthread%/${thread}/" $prm_hist >temp_${rand}.prm
sed -i "s/%rep%/${QMSim_rep}/" temp_${rand}.prm
sed -i "s/%h2%/${QMSim_h2}/" temp_${rand}.prm
sed -i "s/%qtlh2%/${QMSim_qtlh2}/" temp_${rand}.prm
sed -i "s/%bottleneck%/${bottleneck}/" temp_${rand}.prm

## 子群体
for i in $(seq 0 $((np - 1))); do
  sed "s/%pop%/${breeds[${i}]}/" $prm_sub >${breeds[${i}]}.prm
  sed -i "s#%founder_select%#${founder_sel[${i}]}#" ${breeds[${i}]}.prm
  sed -i "s#%seg_select%#${seg_sel[${i}]}#" ${breeds[${i}]}.prm
  sed -i "s#%last_select%#${last_sel[${i}]}#" ${breeds[${i}]}.prm
  sed -i "s/%md%/rnd/" ${breeds[${i}]}.prm
  sed -i "s/%accur%/${QMSim_qtlh2}/" ${breeds[${i}]}.prm
  sed -i "s/%seg_ng%/${seg_gens[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%last_male%/${last_males[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%last_female%/${last_females[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%extentLD%/${extentLDs[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%last_litter_size%/${last_litters[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%geno_gen%/${gid_gen[*]}/g" ${breeds[${i}]}.prm
  sed -i "s/%freq_gen%/${gid_gen[*]}/g" ${breeds[${i}]}.prm

  ## 合并到主参数文件
  cat ${breeds[${i}]}.prm >>temp_${rand}.prm
  rm ${breeds[${i}]}.prm
done

## 基因组和输出参数
sed "s/%nmloci%/${nmloc}/" $prm_geno >genome.prm
sed -i "s/%nchr%/${nchr}/" genome.prm
sed -i "s/%nqloci%/${nqloci}/" genome.prm

## 最终模拟参数卡
cat genome.prm >>temp_${rand}.prm
rm genome.prm

## 文件路径
cd ${proj} || exit 5

## 输出文件夹
sed "s#%out_dir%#${sim_dir}#" ${prmpath}/temp_${rand}.prm >QMSim.prm
rm ${prmpath}/temp_${rand}.prm

## QMSim模拟程序运行
QMSim QMSim.prm

## 更换工作路径
cd ${proj}/${sim_dir} || exit

## 保存QMSim模拟种子
cp seed QMSim_${sim_dir}_seed.prv

## debug
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Two
