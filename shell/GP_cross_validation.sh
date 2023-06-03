#!/bin/bash
#SBATCH --job-name=GP_cross_validation
#SBATCH --output=GP_cross_%j.log


########################################################################################################################
## 版本: 1.0.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-06-02
## 
## 功能：
## 1.根据提供的plink格式的基因型文件进行表型模拟
## 2.根据提供的表型和基因型(系谱)信息进行品种内遗传评估并计算评估准确性
## 3.根据提供的表型和基因型(系谱)信息进行多品种遗传评估并计算评估准确性
## 
## 使用: GP_cross_validation.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################



###################  参数处理  #####################
###################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long code:,proj:,type:,breeds:,thread:,traits:,trait:,h2s:,rg_sim:,rg_pri:,rg_dist:,means:,method:,phef:,pedf:,binf:,nqtl:,nbin_cor:,nsnp_cor:,nsnp_win:,prior:,bin:,bin_sim:,tbv_col:,soft:,all_eff:,ran_eff:,iter:,burnin:,ref:,dirPre:,nbin:,bfile:,seed:,fold:,rep:,gen:,nsnp:,nsnp_sim:,out:,sim_dir:,nginds:,seg_gens:,extentLDs:,last_males:,last_females:,founder_sel:,seg_sel:,last_sel:,last_litters:,geno_gen:,maf:,binDiv:,binThr:,nchr:,nmloc:,nqloci:,QMSim_h2:,debug,suffix,noCov,append,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必需]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'YY DD' [必需]
    --traits )   traits="$2";   shift 2 ;; ## 表型文件中所有性状名，如'AGE BF' [必需]
    --type )     type="$2";     shift 2 ;; ## 分析类型，bin/adj/gsim/geno/psim/within/blend/union/multi/accur/var [必需]
    --trait )    trait="$2";    shift 2 ;; ## 需要进行分析的性状名，如'AGE' [traits]
    --bfile )    bfile="$2";    shift 2 ;; ## plink文件前缀，如"/public/home/merge" [NULL]
    --phef )     phef="$2";     shift 2 ;; ## 表型文件 [NULL]
    --pedf )     phef="$2";     shift 2 ;; ## 系谱文件 [NULL]
    --rg_dist )  rg_dist="$2";  shift 2 ;; ## 加性遗传相关大小服从的分布，uniform/normal/identical [identical]
    --rg_pri )   rg_pri="$2";   shift 2 ;; ## 多性状模型中群体间的加性遗传相关大小先验 [NULL]
    --h2s )      h2s="$2";      shift 2 ;; ## 性状模拟时各品种性状的遗传力大小 [0.2]
    --means )    means="$2";    shift 2 ;; ## 性状模拟时各品种群体均值 [NULL]
    --iter )     iter="$2";     shift 2 ;; ## MCMC 总循环次数 [30000]
    --burnin )   burnin="$2";   shift 2 ;; ## MCMC burn-in循环数 [20000]
    --ref )      ref="$2";      shift 2 ;; ## 用于区间划分的SNP面板，M/1/2/... [M]
    --dirPre )   dirPre="$2";   shift 2 ;; ## JWAS输出文件夹增加的前缀 [NULL]
    --method)    method="$2";   shift 2 ;; ## 育种值估计方法，可为PBLUP/GBLUP/ssGBLUP/BayesAS [GBLUP]
    --nsnp_cor ) nsnp_cor="$2"; shift 2 ;; ## 当设定指定数目区间内品种间存在遗传相关时，区间内的QTL数目 [10]
    --nbin_cor ) nbin_cor="$2"; shift 2 ;; ## 设定指定数目区间内品种间存在遗传相关 [10]
    --rg_sim )   rg_sim="$2";   shift 2 ;; ## 性状模拟时群体间的加性遗传相关大小 [0.2]
    --nsnp_sim ) nsnp_sim="$2"; shift 2 ;; ## 性状模拟时，以固定数目进行区间划分时，每个区间内所含的标记数 [60]
    --bin_sim )  bin_sim="$2";  shift 2 ;; ## 性状模拟时区间划分方法，binf路径/win/chr [win]
    --nsnp_win ) nsnp_win="$2"; shift 2 ;; ## 多品种评估以固定数目进行区间划分时，每个区间内所含的标记数 [100]
    --nqtl )     nqtl="$2";     shift 2 ;; ## 性状模拟时的QTL数目 [300]
    --nbin )     nbin="$2";     shift 2 ;; ## 多品种评估以固定数目SNP为区间划分依据时，尽可能接近这个设置的区间数 [NULL]
    --tbv_col )  tbv_col="$2";  shift 2 ;; ## 不进行校正表型计算，真实育种值在表型文件中tbv_col列 [NULL]
    --bin )      bin="$2";      shift 2 ;; ## 多品种评估时区间划分方法，fix/frq/ld/lava/cubic [NULL]
    --prior )    prior="$2";    shift 2 ;; ## 多品种评估时方差组分先验文件 [NULL]
    --all_eff )  all_eff="$2";  shift 2 ;; ## 固定效应和随机效应所在列，如"2 1" ["2 1"]
    --ran_eff )  ran_eff="$2";  shift 2 ;; ## 随机效应所在列，如"1" [1]
    --binf )     binf="$2";     shift 2 ;; ## 区间划分文件 [NULL]
    --soft)      software="$2"; shift 2 ;; ## 育种值估计程序，可为JWAS/C [JWAS]
    --rep)       rep="$2";      shift 2 ;; ## 交叉验证重复数 [1]
    --fold)      fold="$2";     shift 2 ;; ## 交叉验证折数 [1]
    --gen)       gen="$2";      shift 2 ;; ## 有表型个体中倒数gen个世代个体为验证群 [ALL]
    --seed )     seed="$2";     shift 2 ;; ## MCMC抽样及验证群划分时的随机种子 [40296]
    --code )     code="$2";     shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --thread )   thread="$2";   shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --sim_dir )      sim_dir="$2";      shift 2 ;; ## 模拟结果文件输出文件夹名，注意该文件夹不能已存在 [rep1]
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
    --maf )          maf="$2";          shift 2 ;; ## 品种内模拟基因型标记选择标准 [0.01]
    --binDiv )       binDiv="$2";       shift 2 ;; ## 抽样SNP时区间划分的依据，pos/frq [pos]
    --binThr )       binThr="$2";       shift 2 ;; ## 抽样SNP时区间划分长度，物理位置cM或基因频率步长 [10]
    --nsnp )         nsnp="$2";         shift 2 ;; ## 模拟基因型中需要选择的标记个数 [50000]
    --rg_local )  rg_local=true;shift   ;; ## 每个基因组分区的协方差不同，为局部ld、frq相关系数
    --noCov )    noCov=true;    shift   ;; ## 性状间的残差效应约束为0
    --suffix )   suffix=true;   shift   ;; ## 在union/blend/multi等文件夹后添加品种名称后缀，如blend_YY_LL
    --append )   append=true;   shift   ;; ## 计算校正表型时，把各个群体的校正表型进行合并
    --debug )    debug=true;    shift   ;; ## 不跑DMU、gmatrix、bayes等时间长的步骤
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
  echo "parameter --breeds is reduired! "
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

## 脚本
pheno_sim=${code}/R/pheno_simulation.R
phe_adj=${code}/shell/dmu_get_pheno_adj.sh
GP_single=${code}/shell/GP_single_breed.sh
GP_multi=${code}/shell/GP_multi_breed.sh
run_QMSim=${code}/shell/run_QMSim.sh
genome_process=${code}/shell/QMSim_genome_process.sh
accuracy_summ=${code}/shell/accuracy_summary.sh
varComp_summ=${code}/shell/varcomp_summary.sh
block_define=${code}/shell/lava_cubic_bolck.sh
func=${code}/shell/function.sh

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 加载自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink gmatrix mbBayesAS LD_mean_r2 run_dmu4 run_dmuai

## 检查需要的脚本文件是否存在且具有执行权限
check_command $pheno_sim $GP_single $GP_multi $block_define $phe_adj

## 默认参数
method=${method:=GBLUP}
all_eff=${all_eff:="2 1"}
maf=${maf:="0.01"}
ran_eff=${ran_eff:=1}
seed=${seed:=8123}
rg_sim=${rg_sim:=0.2}
h2s=${h2s:=0.2}
iter=${iter:=30000}
burnin=${burnin:=20000}
ref=${ref:=M}
nsnp_cor=${nsnp_cor:=10}
nbin_cor=${nbin_cor:=10}
nsnp_win=${nsnp_win:=100}
nsnp_sim=${nsnp_sim:=60}
nsnp=${nsnp:=50000}
nqtl=${nqtl:=300}
bin_sim=${bin_sim:=win}
type=${type:=sim}
trait=${trait:=${traits}}
rg_dist=${rg_dist:=identical}
QMSim_rep=${QMSim_rep:="1"}
nchr=${nchr:="18"}
nmloc=${nmloc:="300000"}
nqloci=${nqloci:="100"}
QMSim_h2=${QMSim_h2:="0.3"}
QMSim_qtlh2=${QMSim_qtlh2:="${QMSim_h2}"}
bottleneck=${bottleneck:="250"}
prmpath=${prmpath:="${code}/prm"}
sim_dir=${sim_dir:=rep1}
geno_gen=${geno_gen:="8-10"}
SLURM_JOB_ID=${SLURM_JOB_ID:=$RANDOM}

## 从参数获取信息
read -ra breeds_array <<<"$breeds"
read -ra traits_array <<<"$traits"
read -ra trait_array <<<"$trait"
nbreed=${#breeds_array[@]}
ntrait=${#trait_array[@]}
[[ ${out} ]] && out=" --out ${out} "
[[ ${debug} ]] && debug=" --debug "
[[ ${suffix} ]] && suffix=" --suffix "
[[ ${noCov} ]] && noCov=" --res_const "
[[ ${rg_local} ]] && rg_local=" --rg_local "
[[ ${gen} ]] && gen=" --gen ${gen} "
[[ ${binf} ]] && binf=" --binf ${binf} "
[[ ${pedf} ]] && pedf=" --pedf ${pedf} "
[[ ${tbv_col} ]] && tbv_col=" --tbv_col ${tbv_col} "
[[ ${nbin} ]] && nbin=" --nbin ${nbin} "
[[ ${rg} ]] && rg=" --rg ${rg} "
[[ ${dirPre} ]] && dirPre=" --dirPre ${dirPre} "
[[ ${prior} ]] && prior=" --priorVar ${prior} "
[[ ${bin} ]] && bin=" --bin ${bin} "
[[ ${soft} ]] && software=" --software ${software} "
[[ ${founder_sel} ]] && founder_sel=" --founder_sel ${founder_sel} "
[[ ${seg_sel} ]] && seg_sel=" --seg_sel ${seg_sel} "
[[ ${last_sel} ]] && last_sel=" --seg_sel ${last_sel} "
[[ ${nginds} ]] && nginds=" --nginds ${nginds} "
[[ ${last_males} ]] && last_males=" --last_males ${last_males} "
[[ ${last_females} ]] && last_females=" --last_females ${last_females} "
[[ ${last_litters} ]] && last_litters=" --last_females ${last_litters} "
[[ ${extentLDs} ]] && extentLDs=" --extentLDs ${extentLDs} "
[[ ${seg_gens} ]] && seg_gens=" --seg_gens ${seg_gens} "

## 模拟表型时，性状数为1，无性状名
if [[ ${ntrait} == "0" ]]; then
  ntrait=1
  trait_array=('/')
  traits=/
fi

## 修改作业名
if [[ ${SLURM_CPUS_ON_NODE} ]]; then 
  thread=${thread:=${SLURM_CPUS_ON_NODE}}

  ## 修改作业名
  scontrol update jobid=${SLURM_JOB_ID} name=${type}
fi

## 并行作业数
thread=${thread:=$((rep * fold))}

## 保存当前的工作路径
workdir=$(pwd)

## 修改工作文件夹
cd ${proj} || exit 5

## 从合并群中提取各个品种基因型信息(需要map/ped文件用于生成模拟表型)
if [[ ${bfile} ]]; then
  check_plink ${bfile}
  unset bfiles
  for b in "${breeds_array[@]}"; do
    echo ${b} >fid.txt
    [[ ! -s ${b}m.map ]] &&
      plink --bfile ${bfile} --keep-fam fid.txt --freq --recode --out ${b}m &>/dev/null
    bfiles="${bfiles} ${proj}/${b}m"
  done
  rm fid.txt
fi

## 文件夹
# mkdir -p ${proj}
# cd ${proj} || exit

## 日志文件
logf=${proj}/${type}_${SLURM_JOB_ID}.log

if [[ ${type} == "bin" ]]; then
  ## 区间划分
  $block_define --bfile ${bfile} --minSize ${nsnp_sim} --out M_bin.txt
elif [[ ${type} == "adj" ]]; then
  ## 计算校正表型（ebv+re）
  for pi in $(seq 1 ${ntrait}); do # pi=1;b=${breeds_array[0]}
    phedir=${proj}/${trait_array[$((pi - 1))]}

    ## 分析指定表型的文件夹
    mkdir -p ${phedir}
    cd "${phedir}" || exit

    [[ -s ${phedir}/phe_adj_PBLUP.txt ]] && \
      echo "Warning: phe_adj_PBLUP.txt existe! "

    for b in "${breeds_array[@]}"; do
      ## 分析指定表型中指定品种的文件夹
      mkdir -p ${phedir}/${b}
      cd ${phedir}/${b} || exit

      if [[ ! -s ${b}_dmu.txt ]]; then
        ## 检查基因型文件是否存在及是否为二进制文件
        check_plink ${bfile}

        ## 检查plink文件中有无该品种标识的家系id
        if [[ $(grep -c "${b}" ${bfile}.fam) -eq 0 ]]; then
          echo "no family id ${b} in ${bfile}.fam file! "
          exit 1
        fi

        ## 提取指定品种基因型(可能需要对基因型文件进行修改，所以也是复制基因型)
        echo ${b} >tmp_fid.txt
        plink --bfile ${bfile} --keep-fam tmp_fid.txt --chr-set 30 --make-bed --out ${b} &>>${logf}
        rm tmp_fid.txt

        ## 品种相应表型文件
        awk 'FNR==NR{a[$2];next} $1 in a' ${b}.fam ${phef} >${b}_dmu.txt
        # grep ${b} ${phef} | awk '{$NF="";print}' > ${b}_dmu.txt # 表型文件中最后一列为品种

        ## 计算校正表型
        $phe_adj \
          --phereal ${pi} \
          --bfile ${b} \
          --DIR phe_adj_PBLUP \
          --phef ${b}_dmu.txt \
          --all_eff "2 1" \
          --ran_eff "1" \
          --out ${phedir}/phe_adj_PBLUP.txt \
          --append &>>${logf}
      fi
    done
  done
elif [[ ${type} == "gsim" ]]; then
  srun -c10 --mem=50G  $run_QMSim \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --thread "${thread}" \
    --sim_dir "${sim_dir}" \
    --geno_gen "${geno_gen}" \
    --nmloc "${nmloc}" \
    --nchr "${nchr}" \
    "${extentLDs}" \
    "${seg_gens}" \
    "${last_males}" \
    "${last_females}" \
    "${last_litters}" \
    "${nginds}" \
    "${founder_sel}" \
    "${seg_sel}" \
    "${last_sel}"
elif [[ ${type} == "geno" ]]; then
  $genome_process \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --geno_gen "${geno_gen}" \
    --maf "${maf}" \
    --nsnp "${nsnp}" \
    "${last_females}" \
    "${last_litters}" \
    "${nginds}" \
    "${binDiv}" \
    "${binThr}"
elif [[ ${type} == "psim" ]]; then
  ## 模拟表型
  if [[ ! -s pheno_sim.txt ]]; then
    $pheno_sim \
      --h2 "${h2s}" \
      --mean "${means}" \
      --rg "${rg_sim}" \
      --gt "${bfiles}" \
      --bin ${bin_sim} \
      --win ${nsnp_sim} \
      --nqtl ${nqtl} \
      --nbin_cor ${nbin_cor} \
      --nsnp_cor ${nsnp_cor} \
      --dist_cor ${rg_dist} \
      --seed ${seed} \
      --out pheno_sim.txt \
      --fid \
      ${binf} \
      --qtlf qtl_info.txt &>>${logf}
    [[ $? -ne 0 ]] && echo "phenotypes simulation error! " && exit 1

    ## 去除qtl
    awk '{print $2}' qtl_info.txt >qtl_snpid.txt
  fi
elif [[ ${type} == "within" ]]; then
  ## 群体内评估
  for pi in $(seq 1 ${ntrait}); do # pi=1;b=${breeds_array[0]}
    phedir=${proj}/${trait_array[$((pi - 1))]}

    ## 真实育种值所在文件
    if [[ -s ${phedir}/phe_adj_PBLUP.txt ]]; then
      tbvf=${phedir}/phe_adj_PBLUP.txt
    else
      tbvf=${phef}
    fi

    for b in "${breeds_array[@]}"; do
      ## 工作文件夹
      mkdir -p ${phedir}/${b}
      cd ${phedir}/${b} || exit

      if [[ -s ${phedir}/qtl_snpid.txt ]]; then
        ## 在基因型文件中去除qtl
        plink \
          --file ${phedir}/${b}m \
          --exclude ${phedir}/qtl_snpid.txt \
          --make-bed \
          --out ${phedir}/${b}mq &>>${logf}

        ## 提取品种b的表型
        grep "${b}" ../pheno_sim.txt | awk '{print $2="", $0}' | awk '$2="1"' >${b}_dmu_pheno.txt
        bfile=${phedir}/${b}mq
        phef=${b}_dmu_pheno.txt
      else
        bfile=${phedir}/${b}/${b}
      fi

      ## dmu计算准确性
      $GP_single \
        --label ${b} \
        --phef ${phef} \
        --bfile ${bfile} \
        --method ${method} \
        --seed ${seed} \
        --all_eff "${all_eff}" \
        --ran_eff "${ran_eff}" \
        --rep ${rep} \
        --fold ${fold} \
        --phereal ${pi} \
        --thread ${thread} \
        --tbvf ${tbvf} \
        ${software} \
        ${pedf} \
        ${debug} \
        ${gen} \
        ${binf} \
        ${tbv_col} \
        --out accur_GBLUP.txt &>>${logf}
    done
  done
elif [[ ${type} == "blend" || ${type} == "union" || ${type} == "multi" ]]; then
  for ti in "${trait_array[@]}"; do # ti=MS;pi=1
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)

    phedir=${proj}/${trait_array[$((pi - 1))]}
    cd ${phedir} || exit 5

    ## 真实育种值所在文件
    if [[ -s ${phedir}/phe_adj_PBLUP.txt ]]; then
      ## 校正表型
      tbvf=${phedir}/phe_adj_PBLUP.txt
    else
      ## 表型文件中含"真实"育种值
      tbvf=${phef}
    fi

    ## 多品种评估
    # $GP_multi \
    #   --pops "${breeds}" \
    #   --type ${type} \
    #   --tbvf ${tbvf} \
    #   --phereal ${pi} \
    #   --thread ${thread} \
    #   ${suffix} \
    #   ${tbv_col} \
    #   ${debug} \
    #   --out accur_GBLUP &>>${logf}
    $GP_multi \
      --pops "${breeds}" \
      --type ${type} \
      --tbvf ${tbvf} \
      --phereal ${pi} \
      --nsnp_win ${nsnp_win} \
      --thread ${thread} \
      --iter ${iter} \
      --burnin ${burnin} \
      --ref ${ref} \
      --seed ${seed} \
      ${software} \
      ${bin} \
      ${prior} \
      ${dirPre} \
      ${rg_local} \
      ${binf} \
      ${rg} \
      ${suffix} \
      ${nbin} \
      ${tbv_col} \
      ${noCov} \
      ${debug} &>>${logf}
  done
elif [[ ${type} == "multi" ]]; then
  ## 真实育种值所在文件
  if [[ -s ${phedir}/phe_adj_PBLUP.txt ]]; then
    ## 校正表型
    tbvf=${phedir}/phe_adj_PBLUP.txt
  else
    ## 表型文件中含"真实"育种值
    tbvf=${phef}
  fi

  ## 多品种评估 Bayes
  $GP_multi \
    --pops "${breeds}" \
    --type ${type} \
    --tbvf ${tbvf} \
    --phereal ${pi} \
    --nsnp_win ${nsnp_win} \
    --thread ${thread} \
    --iter ${iter} \
    --burnin ${burnin} \
    --ref ${ref} \
    --seed ${seed} \
    ${software} \
    ${bin} \
    ${prior} \
    ${dirPre} \
    ${rg_local} \
    ${binf} \
    ${rg} \
    ${suffix} \
    ${nbin} \
    ${tbv_col} \
    ${noCov} \
    ${debug} \
    --out accur_bayes &>>${logf}
elif [[ ${type} == "accur" ]]; then
  ## 准确性统计
  $accuracy_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --bin "${bin}" \
    --code "${code}" \
    ${out} \
    ${software}
elif [[ ${type} == "var" ]]; then
  ## 遗传参数统计
  $varComp_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --bin "${bin}" \
    --code "${code}" \
    ${out} \
    ${software}
fi

## 删除没有信息的日志文件
# rm ${workdir}/GP_cross_${SLURM_JOB_ID}.log

## within
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Xie2021
breeds="YY LL"
traits="PFAI MS"
bfile=${proj}/Genotype.id.qc
phef=${proj}/phenotypes_dmu.txt
type=adj

code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Xie2021
# phef=${proj}/phenotypes_dmu.txt
breeds="YY LL"
traits="PFAI MS"
# trait="MS"
type=accur
bin="fix lava cubic"

## gsim
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Two/rep3
breeds="A B C"
thread="10"
sim_dir="rep3"
founder_sel="rnd,rnd,rnd"
seg_sel="phen,phen /l,rnd"
last_sel="rnd,rnd,rnd"
seg_gens="40 40 40"
extentLDs="10 10 10"
last_males="40 40 40"
last_females="200 200 200"
nmloc="100000"
nginds="600 600 600"
type=gsim

## psim
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Two/rep1/identical/cor0.2
bfile=/home/liujf/WORKSPACE/liwn/mbGS/data/Two/rep1/merge
breeds="A B"
means="0.5 1.0"
h2s="0.5 0.4"
rg_sim="0.2"
rg_dist=identical
seed=149994
type=sim

## sim within 
code=/home/liujf/WORKSPACE/liwn/mbGS/code
proj=/home/liujf/WORKSPACE/liwn/mbGS/data/Two/rep1/identical/cor0.2
bfile=/home/liujf/WORKSPACE/liwn/mbGS/data/Two/rep1/merge
phef=${proj}/pheno_sim.txt
breeds="A B"
tbv_col=6
# rep=5
# fold=5
type=blend
