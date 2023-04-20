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
