#!/bin/sh
#$ -S /bin/sh

# suppress BLAS to use multiple CPUs
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

PATH=/home/imgishi/miniconda3/envs/ldsc/bin:$PATH

old_path=$1
new_path=$2
sig=$3
bim_path=$4


## 5. LD scores
bedlines=`ls ${old_path}/${sig} | cat | sort | uniq | cut -f 1 -d '.'`

for bedline in ${bedlines};do

 bedname=`echo ${bedline} | awk '{print $1}'`

 mkdir -p ${new_path}/${sig}/${bedname}

 for chr in $(seq 1 22);do
   
   #LD score
   /home/imgishi/tools/ldsc/ldsc.py \
      --l2  \
      --yes-really  \
      --bfile  ${bim_path}/1000G.EAS.QC.${chr}   \
      --ld-wind-cm  1  \
      --annot  ${old_path}/${sig}/${bedname}/${bedname}.${chr}.annot.gz  \
      --out  ${new_path}/${sig}/${bedname}/${bedname}.${chr}  \
      --print-snps  reference/LDSCORE/baseline_eas_v1.2/${chr}.snps 
 done
done



## 6. LDSC
# downloaded from LDSC page
# use baseline v1.2
weights2=reference/LDSCORE/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.
frqfile=reference/LDSCORE/1000G_Phase3_EAS_frq/1000G.EAS.QC.
baseline=reference/LDSCORE/baseline_eas_v1.2/baseline.

# SLE-EAS GWAS
trait=Lupus_ard2021
sumstats=reference/gwas_sumstats_ldsc/EAS/$trait/$trait.sumstats.gz

bedlines=`ls ${old_path}/${sig} | cat | sort | uniq | cut -f 1 -d '.'`

for bedline in ${bedlines};do

 bedname=`echo ${bedline} | awk '{print $1}'`

 mkdir -p ${new_path}/${sig}/${bedname}

 mainannot=${old_path}/${sig}/${bedname}/${bedname}.

 # background: union of all low-exp filter passed genes
 background=sclinker/LDscore/EAS/allgenes/allgenes/${bedname}/${bedname}.
   
 python  tools/ldsc/ldsc.py \
   --h2    $sumstats  \
   --ref-ld-chr  $mainannot,$background,$baseline \
   --w-ld-chr    $weights2  \
   --overlap-annot         \
   --frqfile-chr   $frqfile  \
   --out  ${new_path}/${sig}/${bedname}/${trait} \
   --print-coefficients
done






