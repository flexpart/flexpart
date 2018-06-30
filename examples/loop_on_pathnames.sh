#!/bin/bash

#OUTPUT="$(ls -1 pathnames*)"
OUTPUTdefault="$(ls -1 pathnames*)"

OUTPUT=${OUTPUTvar:-$OUTPUTdefault}

OUTPUTpwd="$(pwd)"

FLEXPARTdefault=FLEXPART_8d70e43
#FLEXPART=FLEXPART_8d70e43
FLEXPART=${FLEXPARTvar:-$FLEXPARTdefault}

echo $OUTPUT

for i in ${OUTPUT}
do

FP_slurm_batch_sl=slurm_batch_$i.sl
#FP_slurm_batch_sl=toto.sl

        echo  '#!/bin/bash' > $FP_slurm_batch_sl
        echo  '#SBATCH -J ' FP_$i  >> $FP_slurm_batch_sl
        echo  '#SBATCH -D ' ${OUTPUTpwd} '              # working directory of script'  >> $FP_slurm_batch_sl
        echo  '#SBATCH --mem=8092                       # total memory requirement for the node (in MB)'  >> $FP_slurm_batch_sl
        echo  '#SBATCH --mem-per-cpu=1024               # minimum amount of memory required pr. allocated CPU'  >> $FP_slurm_batch_sl
        echo  '#SBATCH -n 1                             # number of tasks, e.g. number of cores'  >> $FP_slurm_batch_sl
        echo  '#SBATCH -N 1                             # ensure all cores are on the same host/machine'  >> $FP_slurm_batch_sl
        echo  '#SBATCH --mail-type=ALL                  # when to send e-mail (valid options are: BEGIN,END,FAIL,REQUEUE,ALL)'  >> $FP_slurm_batch_sl
        echo  '#SBATCH --mail-user=ip@nilu.no           # who to send email to'  >> $FP_slurm_batch_sl
        echo  '#SBATCH -o output-%N-%j.out              # filename to send standard out to'  >> $FP_slurm_batch_sl
        echo  '#SBATCH -e error-%N-%j.err               # filename to send standard error to'  >> $FP_slurm_batch_sl
        echo  '                                          ' >> $FP_slurm_batch_sl
        echo  'srun -l ' $FLEXPART $i   >> $FP_slurm_batch_sl

	echo $FLEXPART  $i
done

