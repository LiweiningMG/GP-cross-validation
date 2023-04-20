#!/usr/bin/bash
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

## 参数
while getopts "P:p:D:B:b:1:2:E:e:S:s:c:t:M:a:i:L:N:n:T:l:A:f:r:m:I:R:O:Vv4h" optname; do
  case "$optname" in
    "P") phefA="${OPTARG}";;      ## A群完整数据集表型文件
    "p") phefB="${OPTARG}";;      ## B群完整数据集表型文件
    "D") DIR="${OPTARG}";;        ## 参数卡模板
    "B") bfileA="${OPTARG}";;     ## plink二进制文件前缀
    "b") bfileB="${OPTARG}";;     ## plink二进制文件前缀
    "1") A="${OPTARG}";;          ## 群体A标识符，如YY
    "2") B="${OPTARG}";;          ## 群体A标识符，如DD
    "E") pedfA="${OPTARG}";;      ## 系谱文件(不提供则运行GBLUP)
    "e") pedfB="${OPTARG}";;      ## 系谱文件(不提供则运行GBLUP)
    "S") gpedA="${OPTARG}";;      ## 系谱文件(只用于划分参考-验证群)
    "s") gpedB="${OPTARG}";;      ## 系谱文件(只用于划分参考-验证群)
    "c") add_cor="${OPTARG}";;    ## 群体A和B的加性相关大小(两性状模型)
    "t") type="${OPTARG}";;       ## 群体合并方式 blend/union(两性状模型)
    "M") GmatM="${OPTARG}";;      ## 基因型矩阵构建方式[multi/single],详情看脚本~multi_pop_Gmat.R
    "a") invA="${OPTARG}";;       ## A逆构建方式(1/2/3/4/6)，1为考虑近交，2为不考虑近交，其他见DMU说明书
    "i") num_int="${OPTARG}";;    ## 整型列列数
    "L") all_eff="${OPTARG}";;    ## DIR中$MODEL第三行(定义所有效应)，前三位数不需要，只需要所有效应所在的列，如"2 3 1"
    "N") ran_eff="${OPTARG}";;    ## DIR中$MODEL第四行(随机效应分组)，第一位数不需要，只需要所有随机效应所在分组，如"1"
    "n") ran_effB="${OPTARG}";;   ## 多性状模型中第2性状的随机效应分组
    "T") tbv_col="${OPTARG}";;    ## 真实育种值在表型文件的第几列
    "l") preal="${OPTARG}";;      ## 表型在表型文件中实数列的位置
    "A") add_rf="${OPTARG}";;     ## 加性效应所在分组    
    "f") fold="${OPTARG}";;       ## 交叉验证倍数
    "r") rep="${OPTARG}";;        ## 重复计算次数
    "m") miss="${OPTARG}";;       ## 缺失表型表示符
    "I") name_int="${OPTARG}";;   ## 整型变量列列名
    "R") name_real="${OPTARG}";;  ## 实型变量列列名
    "O") out="${OPTARG}";;        ## 输出准确性文件名前缀
    "V") VarA=TRUE;;              ## 单性状联合评估中使用A群体的方差组分
    "v") VarB=TRUE;;              ## 单性状联合评估中使用B群体的方差组分
    "4") dmu4=TRUE;;              ## 使用dmu4进行评估，不估计方差组分
    "h") echo "Open script to view instructions" && exit 1;;
    ":") echo "No argument value for option $OPTARG";;
    "?") echo "Unknown option $OPTARG";;
    *) echo "Unknown error while processing options" ;;
  esac
done

####################
## 软件
module load R/4.1.0
module load PLINK/1.90
module load mkl  ## 宁超师兄的GMAT软件需要调用

## 检查必要参数是否提供
[[ ! -f ${phefA} || ! -f ${phefB} || ! -f ${DIR} || ! -f ${bfileA}.fam || ! -f ${bfileB}.fam || ! ${num_int} ]] && \
  echo "missing required parameters or files not found, please check! " && \
  exit 1

## 默认参数
[ ! ${fold} ] && fold=1
[ ! ${miss} ] && miss=-99
[ ! ${out} ] && out="accuracy"
[ ! ${type} ] && type="blend"
[ ! ${GmatM} ] && GmatM="single"
[ ! ${ran_effB} ] && ran_effB=${ran_eff}
[ ! ${invA} ] && invA=2

## 路径
code=/storage1/active/liwn/code
gmatrix=/storage1/active/liwn/Mysoftware/GMAT/gmatrix
multiG=${code}/R/GS/multi_pop_Gmat.R
pheMerge=${code}/R/pheno_pre/pheno_merge_pop2.R
workdir=$(pwd)

## 工作文件夹
mkdir -p ${workdir}/${type}
cd ${workdir}/${type} || exit

## 基因型合并
plink --bfile ${bfileA} --bmerge \
  ${bfileB}.bed \
  ${bfileB}.bim \
  ${bfileB}.fam \
  --out pA_pB

## 生成基因型关系矩阵G阵
if [[ ${GmatM} == 'multi' ]]; then
  ## 群体ID
  awk '{print $2}' ${bfileA}.fam > pA.ids
  awk '{print $2}' ${bfileB}.fam > pB.ids
  cat pA.ids > pA_pB.id
  cat pB.ids >> pA_pB.id
  
  ## 格式转换，顺便剔除存在缺失的位点
  plink --bfile pA_pB --geno 0 --recode A --out pA_pB
  
  ## 多品种关系矩阵
  $multiG \
    --rawf pA_pB.raw \
    --idAf pA.ids \
    --idBf pB.ids \
    --out pA_pB
  mv pA_pB.grm pA_pB.agrm.id_fmt
elif [[ ${GmatM} == 'single' ]]; then
  ## 生成基因型关系矩阵G阵
  $gmatrix --bfile pA_pB --grm agrm --out pA_pB --inv
else
  echo "Gmatm can only be multi or single! "
  exit 1
fi

for f in $(seq 1 ${fold}); do
  ## 检查表型文件是否存在
  [[ ! -d ${workdir}/${A}/val${f} || ! -d ${workdir}/${B}/val${f} ]] && \
    echo "validations folder does not exist! " && exit 1
  
  ## 合并参考群表型
  $pheMerge \
    --phefA ${workdir}/${A}/val${f}/val${f}_pheno.txt \
    --phefB ${workdir}/${B}/val${f}/val${f}_pheno.txt \
    --nInt ${num_int} \
    --method ${type} \
    --pheCol ${preal} \
    --out val${f}_pheno.txt
  
  ## 验证群id
  cat ${workdir}/${A}/val${f}/val${f}_id_2col.txt > val${f}_id_2col.txt
  cat ${workdir}/${B}/val${f}/val${f}_id_2col.txt >> val${f}_id_2col.txt
  
  ## 系谱
  if [[ ${pedfA} && ${pedfB} ]]; then
    cat ${pedfA} > pedi_${A}_${B}.txt
    cat ${pedfB} >> pedi_${A}_${B}.txt
  else
    [[ -f pedi_${A}_${B}.txt ]] && rm pedi_${A}_${B}.txt
  fi
done

## 表型文件列数
phecolA=$(awk 'END{print NF}' ${phefA})
num_real=$((phecolA - num_int))

## 效应个数
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')

## 方差组分
parfA=${workdir}/${A}/var.PAROUT
parfB=${workdir}/${B}/var.PAROUT
if [[ ${type} == 'union' ]]; then
  if [[ -s ${parfA} && -s ${parfB} ]]; then
    ## 从完整数据集中群体内估计结果提取方差组分
    add1=$(sed -n '1p' ${parfA} | awk '{print $4}')
    add2=$(sed -n '1p' ${parfB} | awk '{print $4}')
    res1=$(sed -n '2p' ${parfA} | awk '{print $4}')
    res2=$(sed -n '2p' ${parfB} | awk '{print $4}')
    
    ## 计算加性协方差，不提供则设为0
    if [ ${add_cor} ];then
      covar=$(echo | awk -v a=${add1} -v b=${add2} -v c=${add_cor} '{printf("%.10f",sqrt(a)*sqrt(b)*c)}')
      covar=${covar::-2}  ## 去掉末尾两位小数，防止相关系数大于1
    else
      covar=0
    fi
    
    ## 准备双性状模型方差组分文件
    {
      echo "1 1 1 ${add1}" 
      echo "1 2 1 ${covar}"
      echo "1 2 2 ${add2}" 
      echo "2 1 1 ${res1}" 
      echo "2 2 2 ${res2}" 
    } >bivar_var.PAROUT

    varf=${workdir}/union/bivar_var.PAROUT
  elif [ ${dmu4} ];then
    echo "${parfA} or ${parfB} not found."
    exit 1
  fi
elif [[ ${type} == 'blend' ]]; then
  ## 合并群单性状评估，用其中某个群体群体内估计的组分，不用合并数据进行组分估计
  if [[ ${dmu4} ]]; then
    if [[ ${VarB} ]]; then
      ## 使用B群方差组分
      if [[ -s ${parfB} ]];then
        varf=${parfB}/${B}/var.PAROUT
      else
        echo "${parfB} not found."
        exit 1
      fi
    else
      ## 使用A群方差组分
      if [[ -s ${parfA} ]];then
        varf=${parfA}/${A}/var.PAROUT
      else
        echo "${parfB} not found."
        exit 1
      fi
    fi
  fi
fi

## 合并验证群评估
for f in $(seq 1 ${fold}); do
  ## 多重交叉时创建单独文件夹，并行运算
  mkdir -p ${workdir}/${type}/val${f}
  cd ${workdir}/${type}/val${f} || exit
  
  ## 复制表型文件
  mv -f ${workdir}/${type}/val${f}_pheno.txt .
  mv -f ${workdir}/${type}/val${f}_id_2col.txt .
  
  ## 通用参数卡设置
  sed "s/num_int/$((num_int + 1))/" ${DIR} > val${f}.DIR
  sed -i "s/miss/${miss}/" val${f}.DIR
  sed -i "s/add_rf/${add_rf}/" val${f}.DIR
  sed -i "s/invA/${invA}/" val${f}.DIR  ## A逆构建方式
  sed -i "s/phef/val${f}_pheno.txt/" val${f}.DIR
  sed -i "s#G_matrix#${workdir}/${type}/pA_pB.agrm.id_fmt#" val${f}.DIR
  sed -i "s#gidf#${workdir}/${type}/pA_pB.id#" val${f}.DIR
  sed -i 's#$ANALYSE.*#$ANALYSE 1 1 0 0#g' val${f}.DIR
  sed -i '/$RESIDUALS.*/d' val${f}.DIR
  sed -i '/$DMU4.*/d' val${f}.DIR
  
  ## 效应及所在列设定
  if [[ ${type} == 'blend' ]]; then
    sed -i '/$PRIOR varf/d' val${f}.DIR
    sed -i "s/all_eff/1 0 $((nA + 1)) $((num_int + 1)) ${all_eff}/g" val${f}.DIR
    sed -i "s/ran_eff/${nR} ${ran_eff}/g" val${f}.DIR
    sed -i "s/num_real/1/" val${f}.DIR
  else
    sed -i "s/num_real/2/" val${f}.DIR
    sed -i "s/all_effA/1 0 ${nA} ${all_eff}/g" val${f}.DIR
    sed -i "s/all_effB/2 0 ${nA} ${all_eff}/g" val${f}.DIR
    sed -i "s/ran_effA/${nR} ${ran_eff}/" val${f}.DIR
    sed -i "s/ran_effB/${nR} ${ran_effB}/" val${f}.DIR
    
    ## 方差组分文件
    if [[ -s ${varf} && ${ran_eff} == "${ran_effB}" ]]; then
      sed -i "s#varf#${varf}#" val${f}.DIR
    else
      sed -i '/$PRIOR varf/d' val${f}.DIR  ## 不提供方差组分
    fi 
  fi
  
  ## 是否提供系谱(不提供则为GBLUP)
  if [[ -s ${workdir}/${type}/pedi_${A}_${B}.txt ]]; then
    sed -i "s#pedf#${workdir}/${type}/pedi_${A}_${B}.txt#" val${f}.DIR
  else
    sed -i "s#.*pedf.*#\$VAR_STR ${add_rf} GREL ASCII ${workdir}/${type}/pA_pB.agiv.id_fmt#" val${f}.DIR
    
    ## 只保留有基因型个体的表型
    mv val${f}_pheno.txt val${f}_pheno_ped.txt
    ${code}/R/pheno_pre/pheno_all_genoid.R \
      --phef val${f}_pheno_ped.txt \
      --gidf "${workdir}/${type}/pA_pB.id" \
      --num_int ${num_int} \
      --out val${f}_pheno.txt
  fi
  
  ## 整型和实型变量名，不提供则删除参数卡中$VARIABLE字段
  if [[ ${name_int} && ${name_real} ]]; then
    sed -i "s/name_int/${name_int}/g" val${f}.DIR
    sed -i "s/name_real/${name_real}/g" val${f}.DIR
  else
    sed -i '/$VARIABLE/d' val${f}.DIR
    sed -i '/name_real/d' val${f}.DIR
    sed -i '/name_int/d' val${f}.DIR
  fi
  
  ## 育种值估计
  if [[ ${dmu4} ]]; then
    sed -i 's#$ANALYSE.*#$ANALYSE 11 9 0 0#g' val${f}.DIR
    run_dmu4 val${f}
  else
    run_dmuai val${f}
  fi
  
  ## 检查ssGBLUP步骤是否出错
  if [[ $? -ne 0 && $? -ne 1000 ]]; then
    echo "Error in dmu4/ai, please check the log file: val${f}.lst"
    exit 1
  else
    echo "dmu of val${f} dataset completed."
  fi
  
  ## 复制保存文件
  [[ -f val${f}.DIR ]] && cp val${f}.DIR ${out}_val${f}.DIR
  [[ -f val${f}.SOL ]] && cp val${f}.SOL ${out}_val${f}.SOL
  [[ -f val${f}.lst ]] && cp val${f}.lst ${out}_val${f}.lst
  [[ -f val${f}.PAROUT ]] && cp val${f}.PAROUT ${out}_val${f}.PAROUT
done

cd ${workdir}/${type} || exit
## 可选参数
[ ! -s ${workdir}/${type}/pedi_${A}_${B}.txt ] && option='--add_rf 3' || option=' '
[ ${tbv_col} ] && option="${option} --tbv_col ${tbv_col}"
[[ ${type} == 'union' ]] && option="${option} --nTrait 2 "

## 准确性计算(A群)
${code}/R/accuracy/accur_dmu_cal.R \
  --phe_all ${phefA} \
  --dir_all ${workdir}/${A}/full \
  --traiti 1 \
  --fold ${fold} \
  --out ${workdir}/${type}/${out}_${A}.txt \
  ${option}
  
[[ ${type} == 'union' ]] && option="${option} --traiti 2"
## 准确性计算(B群)
${code}/R/accuracy/accur_dmu_cal.R \
  --phe_all ${phefB} \
  --dir_all ${workdir}/${B}/full \
  --fold ${fold} \
  --out ${workdir}/${type}/${out}_${B}.txt \
  ${option}

## debug
A=YY
B=DD
phefA=/storage1/active/liwn/tempWork/FCR_GWAS/dmu/${A}/pheno_dmu_${A}.txt
phefB=/storage1/active/liwn/tempWork/FCR_GWAS/dmu/${B}/pheno_dmu_${B}.txt
bfileA=/storage1/active/liwn/tempWork/FCR_GWAS/HIBLUP/${A}/${A}_numq
bfileB=/storage1/active/liwn/tempWork/FCR_GWAS/HIBLUP/${B}/${B}_numq
DIR=/storage1/active/liwn/SelectAndMate/DMU/ParaCard/single_trait_SS_GBLUP.DIR
pedfA=/storage1/active/liwn/tempWork/FCR_GWAS/dmu/${A}/pedi_dmu_${A}.txt
pedfB=/storage1/active/liwn/tempWork/FCR_GWAS/dmu/${B}/pedi_dmu_${B}.txt
num_int=4
type=blend
preal=1
miss=-99
fold=1
all_eff="2 3 1"
ran_eff="1"
add_rf=1
out=${A}_${B}_dmu_accuracy
