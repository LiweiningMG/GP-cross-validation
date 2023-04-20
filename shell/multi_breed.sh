#!/bin/bash
#SBATCH --job-name=multi
#SBATCH --output=/BIGDATA2/cau_jfliu_2/liwn/mbGS/logs/multi_%j.log

## date
## liwn: 2021-12-05
## 用于估计联合评估ssGBLUP的预测准确性，需要先完成群体内准确性计算(即不同群体目录下有val*文件夹)再运行此脚本
## 需要的文件：
##   1.表型文件(id在第一列)
##   2.系谱文件
##   3.plink二进制文件
##   4.参数卡模板
##   5.R脚本文件(可以复制到你指定文件夹，并修改脚本路径)
## 以上文件均可以不在当前文件夹，不在时提供全路径即可
## 需要用到的软件(需保证在控制台中直接键入以下命令能正常运行)
## Rscript; plink; gmatrix BayesAS_LD

###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
TEMP=$(getopt -o h4 --long ref:,rg_local,keep_all,nbin:,binf:,summs:,dirPre:,res_const,all_samps,debug,traceplot,bfileW:,bfileM:,pops:,pedf:,priorVar:,DIR:,vidf:,method:,rg:,re:,type:,GmatM:,gmat:,gidf:,invA:,all_eff:,ran_eff:,tbvf:,tbv_col:,phereal:,add_rf:,fold:,rep:,miss:,bincol:,num_int:,win:,nsnp_win:,r2_merge:,bin:,LD_maxD:,r2:,inter:,fix_snp:,seed:,iter:,burnin:,thin:,report_sep:,code:,thread:,alpha:,out:,software:,updatepri:,VarA,dmu4,overwri,suffix,dense \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"

## 解析参数
while true; do
    case "$1" in
    --pops )      pops="$2";       shift 2 ;; ## 群体A标识符，如'YY DD'
    --bfileW )    bfileW="$2";     shift 2 ;; ## 每个品种的基因型文件，如'./data/YY ./data/DD'
    --bfileM )    bfileM="$2";     shift 2 ;; ## plink二进制文件前缀，群体合并后的
    --binf )      binf="$2";       shift 2 ;; ## 区间文件
    --bincol )    bincol="$2";     shift 2 ;; ## 区间文件中指示区间内标记数的列 [1]
    --pedf )      pedf="$2";       shift 2 ;; ## 系谱文件，可为合并的单个系谱，也可为以空格分隔的每个群体的系谱文件
    --summs )     summs="$2";      shift 2 ;; ## GEMMA结果文件前缀
    --priorVar )  priorVar="$2";   shift 2 ;; ## 方差组分初始值来源 null/predict/pheno/A/B/... [pheno]
    --h2 )        h2="$2";         shift 2 ;; ## 遗传力，用于生成初始方差组分，两个性状不同时空格隔开，如 '0.3 0.1' [0.5]
    --priorRg )   rg="$2";         shift 2 ;; ## 群体A和B的加性相关大小先验(两性状模型) [0.001]
    --priorRe )   re="$2";         shift 2 ;; ## 群体A和B的残差相关大小先验(两性状模型) [0.001]
    --DIR )       DIR="$2";        shift 2 ;; ## 参数卡前缀 [type]
    --method )    method="$2";     shift 2 ;; ## dmu模型
    --type )      type="$2";       shift 2 ;; ## 参考群合并方式 blend/union/multi(两性状模型) [blend]
    --GmatM )     GmatM="$2";      shift 2 ;; ## 基因型矩阵构建方式[multi/single]
    --gmat )      gmat="$2";       shift 2 ;; ## 用户提供关系矩阵或逆矩阵(id id value)
    --gidf )      gidf="$2";       shift 2 ;; ## 基因型个体id，与用户指定的G阵文件中个体id一致
    --invA )      invA="$2";       shift 2 ;; ## A逆构建方式(1/2/3/4/6)，1为考虑近交，2为不考虑近交，其他见DMU说明书 [1]
    --all_eff  )  all_eff="$2";    shift 2 ;; ## DIR中$MODEL第3行(所有效应)，前3位不用，如"2 1"，品种间不同时用-分隔
    --ran_eff  )  ran_eff="$2";    shift 2 ;; ## DIR中$MODEL第4行(随机效应分)，第1位不用，如"1"，品种间不同时用-分隔
    --tbvf )      tbvf="$2";       shift 2 ;; ## 包含真实育种值的文件
    --tbv_col )   tbv_col="$2";    shift 2 ;; ## 真实育种值在育种值文件的第几列
    --phereal )   phereal="$2";    shift 2 ;; ## 表型在表型文件中实数列的位置
    --add_rf )    add_rf="$2";     shift 2 ;; ## 加性效应所在分组
    --fold )      fold="$2";       shift 2 ;; ## 交叉验证倍数
    --rep )       rep="$2";        shift 2 ;; ## 重复计算次数
    --miss )      miss="$2";       shift 2 ;; ## 缺失表型表示符
    --num_int )   num_int="$2";    shift 2 ;; ## 表型文件中整型变量的列数 [从A品种表型文件中获取]
    --nsnp_win )  nsnp_win="$2";   shift 2 ;; ## 合并前每个窗口SNP数
    --win )       win="$2";        shift 2 ;; ## 计算指定SNP侧翼R2均值时的窗口大小 [50]
    --bin )       bin="$2";        shift 2 ;; ## 是否合并临近窗口，fix/frq/ld/lava/cubic [lava]
    --maf )       maf="$2";        shift 2 ;; ## 当bin为cubic时，根据指定次等位基因频率阈值对SNP进行过滤 [-0.01]
    --nbin )      nbin="$2";       shift 2 ;; ## 当bin为fix时，划分区间数约束为与nbin相近
    --ref )       ref="$2";        shift 2 ;; ## Reference panel for dividing blocks M/1[index of breed]/2/... [M]
    --r2_merge )  r2_merge="$2";   shift 2 ;; ## 合并临近窗口时的LD阈值
    --LD_maxD )   LD_maxD="$2";    shift 2 ;; ## plink参数，计算LD时，SNP对的最大距离kb
    --r2 )        r2="$2";         shift 2 ;; ## plink参数，计算LD时，输出结果文件中R2的阈值
    --inter )     inter="$2";      shift 2 ;; ## plink参数，计算LD时，间隔位点数多于inter的标记对不计算LD值
    --seed )      seed="$2";       shift 2 ;; ## MCMC抽样时的随机种子
    --iter )      iter="$2";       shift 2 ;; ## MCMC 总循环次数
    --burnin )    burnin="$2";     shift 2 ;; ## MCMC burn-in循环数
    --thin )      thin="$2";       shift 2 ;; ## MCMC抽样步长
    --report )    report_sep="$2"; shift 2 ;; ## MCMC报告间隔
    --code )      code="$2";       shift 2 ;; ## 脚本路径
    --thread )    thread="$2";     shift 2 ;; ## dmu任务并行数
    --alpha )     alpha="$2";      shift 2 ;; ## ssGBLUP relationship matrix alpha
    --vidf )      vidf="$2";       shift 2 ;; ## 验证群个体id文件名
    --out )       out="$2";        shift 2 ;; ## 输出准确性文件名前缀
    --software )  software="$2";   shift 2 ;; ## 可选值：JWAS/C
    --dirPre )    dirPre="$2";     shift 2 ;; ## JWAS输出文件夹增加的前缀
    --updatepri ) updatepri="$2";  shift 2 ;; ## 指定轮次更新prior方差尺度矩阵先验 [0]
    --dense )     dense=true;      shift   ;; ## 将DMU的ANALYSE中的method设为31，利用多线程计算方差组分
    --suffix )    suffix=true;     shift   ;; ## 在union/blend/multi等文件夹后添加品种名称后缀
    --debug )     debug=true;      shift   ;; ## 不跑DMU、gmatrix、bayes等时间长的步骤
    --rg_local )  rg_local=true;   shift   ;; ## 每个基因组分区的协方差不同，为局部ld、frq相关系数
    --res_const ) res_const=true;  shift   ;; ## JWAS性状间的残差效应约束为0
    --all_samps ) all_samps=true;  shift   ;; ## JWAS输出所有变量的样本
    --keep_all )  keep_all=true;   shift   ;; ## 保留所有Bayes结果文件
    --traceplot ) traceplot=true;  shift   ;; ## MCMC burn-in期后的trace plot
    --overwri )   overwrite=true;  shift   ;; ## 若存在已生成的准备文件(表型、id文件等)，则进行覆写
    -4 | --dmu4 ) dmu4=true;       shift   ;; ## 使用dmu4进行评估，不估计方差组分
    -h | --help ) echo "Open script $0 to view instructions" && exit 1 ;;
    -- ) shift; break ;;
    * ) break ;;
    esac
done

## 调用软件需要的模块
# module load GCC/11.3.0
# module load mkl

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## suffix
if [ "$suffix" = true ]; then
    suffix="_$(echo "${pops}" | tr ' ' '_')"
else
    suffix=""
fi

## 工作目录
workdir=$(pwd)

## 存放结果的目录
tpath=${workdir}/${type}${suffix}
mkdir -p ${tpath}

## 默认参数
nchr=${nchr:=30}
updatepri=${updatepri:=0}
LD_maxD=${LD_maxD:=10000}
r2=${r2:=0}
inter=${inter:=99999}
nsnp_win=${nsnp_win:=100}
win=${win:=50}
bin=${bin:=lava}
ref=${ref:=M}
bincol=${bincol:=1}
fold=${fold:=1}
rep=${rep:=1}
miss=${miss:=-99}
out=${out:=accuracy}
type=${type:=blend}
phereal=${phereal:=1}
h2=${h2:=0.5}
rg=${rg:=0.001}
re=${re:=0.001}
iter=${iter:=30000}
burnin=${burnin:=20000}
thin=${thin:=10}
DIR=${DIR:=${type}}
report_sep=${report_sep:=100}
GmatM=${GmatM:="single"}
invA=${invA:=1}
alpha=${alpha:=0.05}
method=${method:=GBLUP}
maf=${maf:=-0.01}
vidf=${vidf:=val.id}
software=${software:=C}
r2_merge=${r2_merge:=0.1}
add_rf=${add_rf:=1}
geno=${geno:=0.2}
mind=${mind:=0.2}
code=${code:=/BIGDATA2/cau_jfliu_2/liwn/code}
[[ ${tbvf} && ! ${tbv_col} ]] && tbv_col=2 ## 真实育种值文件存在，默认育种值在文件第二列，id在第一列

## 参数初始化
[[ ${res_const} ]] && res_const=" --constraint "
[[ ${all_samps} ]] && all_samps=" --all_samp_out "
[[ ${keep_all} ]] && keep_all=" --keep_all "
if [[ $(echo ${h2} | awk '{print NF}') -gt 1 ]]; then
  h2B=$(echo ${h2} | awk '{print $2}')
  h2=$(echo ${h2} | awk '{print $1}')
else
  h2B=${h2}
fi

## 随机数种子
if [[ ! ${seed} ]]; then
  seed=$RANDOM
  echo "$seed" >MCMC_seed.txt
fi

## 脚本
multiG=${code}/R/GS/multi_pop_Gmat.R
# keep_phe_gid=${code}/R/pheno_pre/pheno_all_genoid.R
pheMerge=${code}/R/pheno_pre/pheno_merge_pop.R
accur_cal=${code}/R/accuracy/accur_dmu_bayes_cal.R
bins_define=${code}/R/LD/bins_define.R
MCMC_polt=${code}/R/plot/MCMC_chain_plot.R
phe_order=${code}/R/pheno_pre/pheno_order.R
variance_prior=${code}/R/breeding/variance_prior.R
job_pool=${code}/shell/Func/job_pool.sh
bayesAS_JWAS=${code}/julia/MT-BayesAS/BayesAS.jl
ldbolck=${code}/shell/LD/bin_define_ldblock.sh
local_rg=${code}/R/LD/local_block_rg.R

## 并行作业数设置
if [[ ! ${thread} ]]; then
  if [[ ${SLURM_CPUS_ON_NODE} ]]; then thread=${SLURM_CPUS_ON_NODE}; else thread=$((rep * fold)); fi
fi

## 作业数池
# shellcheck source=/dev/null
source "${job_pool}"

## 作业池初始化
[[ ! ${debug} ]] && job_pool_init ${thread} 0

##################  自定义函数  ##############
#############################################
merge_plink() {
  ## 合并plink文件
  ## 参数1：plink文件列表，以空格分隔
  ## 参数2：输出文件名前缀
  ## 依赖：plink
  check_redirect() {
    ## 文件格式
    if [[ ! -f $1.fam ]]; then
      if [[ -f $1.map ]]; then
        plink \
          --file $1 \
          --chr-set ${nchr} \
          --make-bed \
          --out $1
      else
        echo "Error: $1.bed not found"
        return 1
      fi
    fi
    ## 合并文件列表
    echo "$1" >>$2
  }

  [[ -s $2.bed ]] && echo "WARN: plink binary file $2 already exists" && return 0

  local plink_list
  IF=" " read -r -a plink_list <<<"$1"

  ## 检查参数个数
  nfiles=${#plink_list[@]}
  if [[ ${nfiles} -lt 2 ]]; then
    echo "Error: at least two plink files are needed to merge"
    return 1
  fi

  ## 临时文件
  local seed_tmp=$RANDOM
  touch merge_${seed_tmp}.txt
  [[ -f merge_${seed_tmp}.txt ]] && echo "merge_${seed_tmp}.txt exist"

  ## 检查文件格式，并输出合并文件列表
  for i in "${!plink_list[@]}"; do
    if [[ ${i} -gt 0 ]]; then
      check_redirect ${plink_list[${i}]} merge_${seed_tmp}.txt
    else
      # touch merge_${seed_tmp}.txt
      # [[ -s merge_${seed_tmp}.txt ]] && cat merge_${seed_tmp}.txt
      [[ -s merge_${seed_tmp}.txt ]] && \
        plink \
          --bfile ${plink_list[0]} \
          --chr-set ${nchr} \
          --merge-list merge_${seed_tmp}.txt \
          --make-bed \
          --out $2
    fi
  done
  # rm merge_${seed_tmp}.txt
}

check_plink_bfile() {
  ## 检查提供的plink前缀对应的bfile是否存在
  ## 第1个参数为染色体数目
  ## 随后跟着任意个plink文件前缀，如/path/to/plinkA /path/to/plinkB ...
  ## 若只存在map/ped文件，将用plink转换为fam/bim/bed文件格式

  for i in $(seq 2 $#); do
    if [[ ! -f ${!i}.fam ]]; then
      if [[ -f ${!i}.map ]]; then
        plink \
          --bfile ${!i} \
          --chr-set $1 \
          --make-bed \
          --out ${!i}
      else
        echo "${!i}.map(.fam) not found! "
        exit 1
      fi
    fi
  done
}

##################  解析命令行参数  ##############
#################################################
np=$(echo ${pops} | awk '{print NF}')
IFS=" " read -r -a popN <<<"$pops"
IFS=" " read -r -a bfileWa <<<"$bfileW"

###############  检查品种内评估是否完成  #########
#################################################
for r in $(seq 1 ${rep}); do
  for f in $(seq 1 ${fold}); do
    for b in "${popN[@]}"; do
      for file in pheno.txt ${vidf}; do
        [[ ! -s ${workdir}/${b}/val${f}/rep${r}/${file} ]] &&
          echo "${workdir}/${b}/val${f}/rep${r}/${file} not found! " &&
          exit 1

        ## 获取表型文件中的整型和实型变量数
        if [[ ! ${num_int} ]]; then
          phef_within=${workdir}/${b}/val${f}/rep${r}/pheno.txt
          ncol=$(awk 'END{print NF}' ${phef_within})
          for i in $(seq 1 ${ncol}); do
            dot=$(awk -vl=${i} '{print $l}' ${phef_within} | grep -c "\.")
            [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
          done
        fi
        [[ ! ${num_real} ]] && num_real=$(($(awk 'END{print NF}' ${phef_within}) - num_int))
      done
    done
  done
done

##################  获取品种内模型效应设定  #############
########################################################
## 从其中一个品种子集中获取效应设定(不同品种效应相同，以第一个品种为准)
firstB=$(echo ${pops} | cut -d ' ' -f 1)
within_DIR=$(find ${workdir}/${firstB}/val1/rep1 -name "*.DIR" | head -n 1)
[[ ! -s ${within_DIR} ]] && echo "${within_DIR} not found! " && exit 1
model_line=$(grep -n "MODEL" ${within_DIR} | head -n 1 | cut -d ':' -f 1)
all_eff=$(sed -n "$((model_line + 3))p" ${within_DIR} | cut -d ' ' -f '4-')
ran_eff=$(sed -n "$((model_line + 4))p" ${within_DIR} | cut -d ' ' -f '2-')
## 效应个数
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')
## 增加品种固定效应
((nA++))
all_eff="$((num_int + 1)) ${all_eff}"

################  基因型文件  #############
##########################################
## 品种内的基因型文件路径
# IFS=" " read -r -a bfiles <<<"${bfileM} ${bfileW[*]}"
if [[ ! ${bfileW} ]]; then
  IF="\n" mapfile -t bfiles < <(printf "%s\n" "${popN[@]}" | xargs -I {} echo "${workdir}/{}/{}")
else
  IFS=" " read -r -a bfiles <<<"${bfileW}"
fi
## 检查fam/bim/bed文件是否存在
check_plink_bfile ${nchr} "${bfiles[@]}"

## 生成各品种合并后的plink文件
if [[ ${type} != "multi" && ! -s ${gmat} ]] || [[ ${type} == "multi" ]]; then
  if [[ ${bfileM} ]]; then
   check_plink_bfile ${nchr} ${bfileM}
  else
    ## 合并plink文件
    bfileM=${tpath}/merge
    merge_plink "${bfiles[*]}" ${bfileM}
  fi
fi

##################  参考群表型文件  ##############
#################################################
$pheMerge \
  --pops "${pops}" \
  --phef "${workdir}/#breed#/val#val#/rep#rep#/pheno.txt" \
  --rep ${rep} \
  --fold ${fold} \
  --nInt ${num_int} \
  --type ${type} \
  --overwri ${overwrite} \
  --pheCol ${phereal} \
  --fixPop \
  --out "${tpath}/val#val#/rep#rep#/pheno.txt"

#####################  系谱文件  #####################
#####################################################
if [[ ${pedf} ]]; then
  for pedi in ${pedf}; do
    if [[ ! -s ${pedi} ]]; then
      echo "${pedi} not found! " && exit 1
    else
      cat ${pedi} >pedi_merge.txt
    fi
  done
  echo "Number of individuals in pedigree: $(wc -l <pedi_merge.txt)"
else
  [[ -s pedi_merge.txt ]] && rm pedi_merge.txt
fi

####################  基因组区域确定  ##################
########################################################
cd ${tpath} || exit
## 基因组分区文件名
if [[ ! -s ${binf} ]]; then
  if [[ ${bin} == "fix" ]]; then
    bin_dir=${bin}_${nsnp_win}
  elif [[ ${bin} == "lava" ]]; then
    bin_dir=${bin}_${ref}_${nsnp_win}
  elif [[ ${bin} == "cubic" ]]; then
    bin_dir=${bin}_${ref}_${win}
  elif [[ ${bin} == "frq" || ${bin} == "ld" ]]; then
    bin_dir=${bin}_${r2_merge}_${nsnp_win}
  else
    echo "${bin} can only be fix, frq, ld or lava! "
    exit 1
  fi
  binf=${tpath}/${bin_dir}.txt
else
  # binf_base=$(basename "${binf}")
  # bin_dir="${binf_base%.*}"
  bin_dir=self_bin
  if [[ ${bincol} != 1 ]]; then
    awk -v col=${bincol} '{print $col}' ${binf} > bin_col1.txt
    binf=$(pwd)/bin_col1.txt
  fi
fi
## 生成文件
if [[ ${type} == 'multi' ]] && [[ ! -s ${binf} || ${overwrite} ]]; then
  echo 'Generating genome regions file...'

  ## 区间划分
  if [[ ${bin} == "fix" ]]; then
    if [[ ${nbin} ]]; then
      nsnp=$(wc -l <${bfileM}.bim)
      ## 约束最终划分的区间数与nbin相近
      nsnp_win=$((nsnp / nbin))
      echo "Number of SNPs per window set to: ${nsnp_win}"
    fi

    $bins_define \
      --win ${nsnp_win} \
      --map ${bfileM}.bim \
      --bin_merge ${bin} \
      --out ${binf}
  elif [[ ${bin} == "frq" || ${bin} == "ld" ]]; then
    ## 检查是否提供了每个品种的基因型文件 
    [[ ! "${bfileW[*]}" ]] && echo "Required parameter 'bfileW' is missing! " && exit 1

    ## 品种共有标记
    awk '{print $2}' ${bfileM}.bim >SNP_share_id.txt

    ## LD计算(只计算共有标记)
    for i in $(seq 0 $((np - 1))); do
      plink --bfile ${bfileWa[i]} \
        --chr-set ${nchr} \
        --extract SNP_share_id.txt \
        --r2 \
        --freq \
        --ld-window-kb ${LD_maxD} \
        --ld-window ${inter} \
        --ld-window-r2 ${r2} \
        --out ${popN[i]} >>${workdir}/plink.log
    done

    ## 根据要求定义区间
    $bins_define \
      --prefixs "${popN[*]}" \
      --bin_merge ${bin} \
      --win ${nsnp_win} \
      --seg ${r2_merge} \
      --map ${bfileM}.bim \
      --out ${binf}~

    ## 评估软件所需区间文件
    awk '{print $3}' ${binf}~ >${binf}
  elif [[ ${bin} == "lava" || ${bin} == "cubic" ]]; then
    bfile_block=bfile${ref}
    $ldbolck \
      --bfile ${!bfile_block} \
      --win ${win} \
      --maf ${maf} \
      --type ${bin} \
      --minSize ${nsnp_win} \
      --out ${binf}~ >${bin}_${ref}_block_${SLURM_JOB_ID}.log
    sed '1d' ${binf}~ | awk '{print $5}' >${binf}
  fi

  ## 检查是否出错
  if [[ ! -s ${binf} ]]; then
    echo 'error in creat bins file! '
    exit 1
  fi
fi

#####################  局部遗传相关  #####################
#########################################################
if [[ ${rg_local} && ${bin} == "lava" ]]; then
  summA=
  summB=
  if [[ -s ${summA}.assoc.txt && -s ${summB}.assoc.txt ]]; then
    echo "calculating local genetic correlations..."
    ## 软件待修改
    $local_rg \
      --pops ${pops} \
      --summ1 ${summA}.assoc.txt \
      --summ2 ${summB}.assoc.txt \
      --bfile ${!bfile_block} \
      --block ${binf}~ \
      --out ${binf}~ &>${bin}_local_rg_${SLURM_JOB_ID}.log
  else
    echo "Error: ${summA}.assoc.txt or ${summB}.assoc.txt ont found! "
    exit 1
  fi
fi

#####################  方差组分准备  #####################
#########################################################
## 方差组分文件名
if [[ ${priorVar} ]]; then
  parfA=$(basename "$(find ${workdir}/${A} -name "*PAROUT" | head -n 1)")
  parfB=$(basename "$(find ${workdir}/${B} -name "*PAROUT" | head -n 1)")

  ## 每个区间协方差先验不同
  [[ ${rg_local} ]] && rg_local=" --rg_local ${binf}~ "

  if [[ ${priorVar} == "pheno" ]]; then
    ## 根据表型方差推算遗传和残差方差
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/pheno.txt"
    varf_para2="${workdir}/${B}/val#val#/rep#rep#/pheno.txt"
  elif [[ ${priorVar} == "predict" ]]; then
    ## 用品种内评估的方差推算多品种遗传和残差方差
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/${parfA}"
    varf_para2="${workdir}/${B}/val#val#/rep#rep#/${parfB}"
  elif [[ ${priorVar} == "A" ]]; then
    ## 用A品种内估计的方差组分
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/${parfA}"
    varf_para2="null"
  elif [[ ${priorVar} == "B" ]]; then
    ## 用B品种内估计的方差组分
    varf_para1="${workdir}/${B}/val#val#/rep#rep#/${parfB}"
    varf_para2="null"
  fi

  ## 生成方差组分(初值)文件
  $variance_prior \
    --filea ${varf_para1} \
    --fileb ${varf_para2} \
    --rep ${rep} \
    --fold ${fold} \
    --h2 ${h2} \
    --h2B ${h2B} \
    --rg ${rg} \
    --re ${re} \
    --type ${type} \
    --var ${priorVar} \
    --norec \
    --add_rf ${add_rf} \
    --overwri ${overwrite} \
    --pcol $((num_int + phereal)) \
    ${rg_local} \
    --out ${tpath}/val#val#/rep#rep#/${type}_${dirPre}${bin}_prior.txt
fi

#####################  关系矩阵准备  #####################
##########################################################
if [[ ${type} == 'blend' || ${type} == 'union' ]]; then
  ## 生成基因型关系矩阵G阵
  if [[ -s ${gmat} ]]; then
    if [[ ${method} == "ssGBLUP" && ! -s ${gidf} ]]; then
      echo "${gidf} not found! " && exit 1
    else
      echo "Use user provided relationship matrix: ${gmat}"
    fi
  elif [[ ${GmatM} == 'multi' ]]; then
    ## 多品种关系矩阵
    [[ ! "${bfileW[*]}" ]] && echo "Required parameter 'bfileW' is missing! " && exit 1

    :> all_breed.id
    for i in $(seq 0 $((np - 1))); do
      awk '{print $2}' ${bfileWa[i]}.fam >${bfileW[i]}.ids
      cat ${bfileW[i]}.ids >>all_breed.id
      idf="${bfileW[i]}.ids ${idf}"
    done

    ## 格式转换，顺便剔除存在缺失的位点(只存在某个群体中的位点)
    plink \
      --bfile ${bfileM} \
      --chr-set ${nchr} \
      --geno 0 \
      --recode A \
      --out merge >>${workdir}/plink.log

    ## 多品种关系矩阵（待修改）
    $multiG \
      --rawf merge.raw \
      --idf ${idf} \
      --out merge
    mv merge.grm merge.agrm.id_fmt
  elif [[ ${GmatM} == 'single' ]]; then
    ## 生成基因型关系矩阵G阵
    if [[ ${method} == "GBLUP" ]]; then
      gmat_inv="--inv"
      [[ ! ${gmat} ]] && gmat=${tpath}/merge.agiv.id_fmt
    elif [[ ${method} == "ssGBLUP" ]]; then
      gmat_inv=""
      [[ ! ${gmat} ]] && gmat=${tpath}/merge.agrm.id_fmt
      [[ ! ${gidf} ]] && gidf=${tpath}/merge.id
    fi

    ## 构建基因组关系矩阵
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfileM} --grm agrm --out ${tpath}/merge ${gmat_inv}
    echo "G matrix created."

    ## 改名
    if [[ ${method} == "GBLUP" ]]; then
      ## 若指定了gmat文件名，则改名(同时可能移动文件)
      [[ ${gmat} && ! -s ${gmat} ]] && \
        mv merge.agiv.id_fmt ${gmat} && \
        echo "gamt created and has been rename to: ${gmat}"
    elif [[ ${method} == "ssGBLUP" ]]; then
      ## 若指定了gmat文件名，则改名(同时可能移动文件)
      [[ ${gmat} && ! -s ${gmat} ]] && \
        mv merge.agrm.id_fmt ${gmat} && \
        echo "gamt created and has been rename to: ${gmat}"
      [[ ${gidf} && ! -s ${gidf} ]] && mv ${tpath}/merge.id ${gidf}
    fi
  else
    echo "Gmat can only be multi or single! "
    exit 1
  fi
fi

###################  dmu参数卡  ####################
####################################################
cd ${tpath} || exit
if [[ ${type} == 'blend' || ${type} == 'union' ]]; then
  ## method
  if [[ ${method} == GBLUP ]]; then
    method_dmu="31"
  else
    method_dmu="1"
  fi

  ## ANALYSE
  [[ ${dmu4} ]] && ANALYSE="11 9 0 0" || ANALYSE="1 ${method_dmu} 0 0"

  ## MODEL
  if [[ ${type} == 'blend' ]]; then
    num_real=1
    MODEL="1\n0\n1 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
  else
    ABSORB="0"
    MODELS="1 0 ${nA} ${all_eff}"
    RANDOMS="${nR} ${ran_eff}"
    REGRES="0"
    NOCOV=""
    nNOCOV=0
    for i in $(seq 2 ${np}); do
      ABSORB="${ABSORB}\n0"
      MODELS="${MODELS}\n${i} 0 ${nA} ${all_eff}"
      RANDOMS="${RANDOMS}\n${nR} ${ran_eff}" # RANDOM为shell关键字，故改名
      REGRES="${REGRES}\n0"
      for j in $(seq 1 $((i -1))); do
        NOCOV="${NOCOV}\n${j} ${i}"
        ((nNOCOV++))
      done
    done

    num_real=${np}
    MODEL="${np}\\n${ABSORB}\\n${MODELS}\\n${RANDOMS}\\n${REGRES}\\n${nNOCOV}${NOCOV}"
  fi

  ## 写出参数卡
  [[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR has been overwrited! "
  {
    echo "\$COMMENT"
    echo "get EBV of individuals in validation population with reduced phenotypes"
    echo "\$ANALYSE ${ANALYSE}"
    echo "\$DATA  ASCII ($((num_int + 1)), ${num_real}, ${miss}) pheno.txt"
    echo -e "\$MODEL\n${MODEL}"
    echo "\$PRIOR %varf%"
  } >${DIR}.DIR

  ## 方差组分结构
  if [[ ${method} == "GBLUP" ]]; then
    echo "\$VAR_STR ${add_rf} GREL ASCII ${gmat}" >>${DIR}.DIR
    add_sol=3
  elif [[ ${method} == "ssGBLUP" ]]; then
    [[ ! -s pedi_merge.txt ]] && echo "pedi_merge.txt not found! " && exit 1
    echo "\$VAR_STR ${add_rf} PGMIX ${invA} ASCII pedi_merge.txt ${gidf} ${tpath}/merge.agrm.id_fmt ${alpha} G-ADJUST" >>${DIR}.DIR
    add_sol=4
  elif [[ ${method} == "PBLUP" ]]; then
    [[ ! -s pedi_merge.txt ]] && echo "pedi_merge.txt not found! " && exit 1
    echo "\$VAR_STR ${add_rf} PED ${invA} ASCII pedi_merge.txt" >>${DIR}.DIR
    add_sol=4
  else
    echo "method can only be PBLUP/GBLUP/ssGBLUP! "
    exit 1
  fi

  ## 输出效应值
  echo "\$SOLUTION" >>${DIR}.DIR
fi

#####################  子集处理和评估  ###################
##########################################################
for r in $(seq 1 ${rep}); do # r=1;f=1
  for f in $(seq 1 ${fold}); do
    ## 更改工作文件夹
    vali_path=${tpath}/val${f}/rep${r}
    cd ${vali_path} || exit

    ## 如果为调试模式则跳过估计育种值步骤
    [[ ${debug} ]] && continue

    #####  GBLUP模型  #####
    #######################
    if [[ ${type} == 'blend' || ${type} == 'union' ]]; then
      ## 拷贝DMU参数卡
      cp ${tpath}/${DIR}.DIR .

      ## 根据需要提供方差组分(或迭代初值)
      if [[ -s ${type}_${dirPre}${bin}_prior.txt ]]; then
        sed -i "s#%varf%#${type}_${dirPre}${bin}_prior.txt#" ${DIR}.DIR
      else
        sed -i '/$PRIOR/d' ${DIR}.DIR
      fi

      ## 育种值估计
      if [[ ${dmu4} ]]; then
        job_pool_run run_dmu4 ${DIR} ${vali_path}
      else
        job_pool_run run_dmuai ${DIR} ${vali_path}
      fi
    fi

    ##### MT-BayesAS模型 #####
    ############################
    if [[ ${type} == 'multi' ]]; then
      ## 固定效应
      fix_eff=${all_eff%" ${ran_eff}"}      

      ## JWAS软件
      if [[ ${software} == JWAS ]]; then
        echo "BayesAS is being implemented using JWAS..."
        ## JWAS软件
        job_pool_run $bayesAS_JWAS \
          --npop ${np} \
          --rawf ${bfileM}.raw \
          --phef ${vali_path}/pheno.txt \
          --output_dir ${vali_path}/${dirPre}${bin_dir} \
          --iter ${iter} \
          --burnin ${burnin} \
          --thin ${thin} \
          --binf ${binf} \
          --logf ${vali_path}/${dirPre}${bin_dir}_gibs_jwas${SLURM_JOB_ID}.log \
          --var_prior ${vali_path}/${type}_${dirPre}${bin}_prior.txt \
          --fix "${fix_eff}" \
          --rnd "${ran_eff}" \
          --update_priors "${updatepri}" \
          ${keep_all} \
          ${res_const} \
          ${all_samps}
      fi

      ## 自己开发软件
      if [[ ${software} == C ]]; then
        echo "BayesAS is being implemented using software written in C..."
        ## 自己编写的软件
        job_pool_run mbBayesAS \
          --bfile ${bfileM} \
          --phef ${vali_path}/pheno.txt \
          --fix "${fix_eff}" \
          --binf ${binf} \
          --iter ${iter} \
          --burnin ${burnin} \
          --thin ${thin} \
          --outdir ${vali_path} \
          --report ${report_sep} \
          --varOut var_${bin}_${r2_merge}_${nsnp_win}.txt \
          --effOut effect_${bin}_${r2_merge}_${nsnp_win}.txt \
          --gebvOut EBV_${bin} \
          --mcmcOut MCMC_process.txt_${bin} \
          --seed ${seed} \
          --logf ${bin}_${r2_merge}_${nsnp_win}_gibs_${SLURM_JOB_ID}.log &

        # ## 杂合度计算
        # plink --bfile ${bfileM} --chr-set ${nchr} --allow-no-sex --hardy >>${workdir}/plink.log
        # het=$(sed '1d' plink.hwe | awk '{ total += $7 } END { print total/NR }')

        # ## 参考群、验证群plink raw格式基因型文件
        # cat ${workdir}/${A}/val${f}/rep${r}/${vidf} >${vidf}
        # cat ${workdir}/${B}/val${f}/rep${r}/${vidf} >>${vidf}

        # [[ ! -s val.raw ]] || [[ -s val.raw && ${overwrite} ]] &&
        #   plink --bfile ${bfileM} --chr-set ${nchr} --keep ${vidf} --recode A --out val >>plink.log
        # [[ ! -s ref.raw ]] || [[ -s ref.raw && ${overwrite} ]] &&
        #   plink --bfile ${bfileM} --chr-set ${nchr} --remove ${vidf} --make-bed --recode A --out ref >>plink.log

        # ## 保证表型排序与raw文件中相同
        # $phe_order \
        #   --fam ref.fam \
        #   --keep geno \
        #   --pheno pheno.txt \
        #   --out pheno.txt
        # [[ $? -ne 0 ]] && exit 1

        # ## 只要表型列(倒数1、2列)
        # awk '{print $(NF-1),$NF}' pheno.txt >pheno_multi_C.txt
      fi
    fi
  done
done

## 等待后台程序运行完毕
[[ ! ${debug} ]] && job_pool_wait
## 释放线程
[[ ! ${debug} ]] && job_pool_shutdown

####################  准确性计算  ##################
###################################################
## 更换工作路径
cd ${tpath} || exit
## 准确性计算参数
[[ ${tbv_col} ]] && option=" --tbv_col ${tbv_col}"
[[ ${tbvf} ]] && option="${option} --tbvf ${tbvf} "
# [[ ${software} == 'C' ]] && option="${option} --famf ${tpath}/val#val#/rep#rep#/val.fam "
## 其他参数
if [[ ${type} == 'multi' ]]; then
  ## MT-Bayes模型
  if [[ ${software} == 'JWAS' ]]; then
    option="${option} --ebvf ${tpath}/val#val#/rep#rep#/${dirPre}${bin_dir}/EBV_y%i%.txt"
    ebv_col=2
  else
    option="${option} --ebvf ${tpath}/val#val#/rep#rep#/EBV_${bin}_y%i%.txt"
    ebv_col=2
  fi
else
  ## MT-GBLUP模型
  option="${option} --add_sol ${add_sol} --dir_val ${tpath}/val#val#/rep#rep#/${DIR}"
  if [[ ${type} == 'blend' ]]; then
    ebv_col=1
  else
    ebv_col=%i%
  fi
fi
## 计算准确性
for i in $(seq 0 $((np - 1))); do
  ip1=$((i+1))
  ## 输出文件名
  if [[ ${type} == 'multi' ]]; then
    outf=${out}_${dirPre}${software}_${bin}_${popN[${i}]}.txt
    [[ ${bin} == 'lava' ]] && outf=${out}_${ref}_${dirPre}${software}_${bin}_${popN[${i}]}.txt
  else
    outf=${out}_${popN[${i}]}.txt
  fi

  ## 计算准确性
  $accur_cal \
    --ebv_col ${ebv_col/\%i\%/${ip1}} \
    --val_idf ${workdir}/${popN[${i}]}/val#val#/rep#rep#/${vidf} \
    --fid ${popN[${i}]} \
    --rep ${rep} \
    --fold ${fold} \
    ${option/\%i\%/${ip1}} \
    --out ${tpath}/${outf}
done

#################  BayesAS 模型收敛作图  #################
#########################################################
## BayesAS 模型 收敛作图
if [[ ${type} == 'multi' ]]; then
  if [[ ${traceplot} ]]; then
    for r in $(seq 1 ${rep}); do
      ## Bayes MCMC链作图
      for f in $(seq 1 ${fold}); do
        $MCMC_polt \
          --files ${vali_path}/${bin}_MCMC_process.txt \
          --start ${burnin} \
          --end ${iter} \
          --thin ${report_sep} \
          --names "μ1 μ2 alpha11 alpha12" \
          --out "${vali_path}/${bin}_MCMC_process"
      done
    done
  fi
fi

###################  删除中间文件  #####################
########################################################
## 基因型文件
# [[ ${bfileM} && ${bfileM} =~ rmMiss ]] && rm ${bfileM}.*

# exit 0

## multi
# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Three/identical/AB2
pops="B C"
bfileM=/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Three/identical/AB2/merge
tbvf=/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Three/identical/AB2/pheno_sim.txt
tbv_col=5
fold=2
rep=5
phereal=4
type=multi
software=C
bin=fix
nsnp_win=100
nbin=279
rg=AB2
dirPre=equ_
thread=10
seed=12095
out=accur_bayes

# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Ratio/rep5/rand/identical/cor0.6
pops="breedA breedB"
bfileM=/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Ratio/rep5/rand/identical/cor0.6/pA_pB
tbvf=/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Ratio/rep5/rand/identical/cor0.6/breedA_breedB_pheno.txt
tbv_col=5
fold=5
rep=5
phereal=5
type=multi
# type=blend
ref=M
software=C
bin=lava
nsnp_win=100
thread=10
seed=12152
out=accur_bayes

# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Zhuang2019/LMA
pops="USA CAN"
bfileM=/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Zhuang2019/pA_pBb
tbvf=/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Zhuang2019/LMA/phe_adj_PBLUP.txt
fold=5
rep=5
type=multi
software=C
bin=fix
nsnp_win=100
thread=25
seed=4820
out=accur_bayes

# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/identical/cor0.2
pops="A B C D"
# bfileM=/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Lee2019/LIM_AAN
tbvf=/BIGDATA2/cau_jfliu_2/liwn/mbGS/QMSim/Four/rep1/identical/cor0.2/pheno_sim.txt
fold=5
rep=5
type=multi
software=C
bin=lava
nsnp_win=100
thread=25
seed=8031
phereal=5
out=accur_GBLUP
suffix=true

# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Xie2021/PFAI
pops="YY LL"
# bfileM=/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Xie2021/merge
tbvf=/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Xie2021/PFAI/phe_adj_PBLUP.txt
fold=5
rep=5
type=multi
phereal=1
software=C
bin=lava
nsnp_win=100
thread=8
seed=8031
out=accur_bayes
suffix=true
