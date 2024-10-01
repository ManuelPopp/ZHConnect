#!/bin/bash
#
#SBATCH -A node
#SBATCH -p node
#SBATCH --job-name=Dami222
#SBATCH -e /home/poppman/shared/dami/log/%x%a.err
#SBATCH -o /home/poppman/shared/dami/log/%x%a.out
#
#SBATCH --qos normal
#SBATCH --time=0-00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3200
#SBATCH --array=1-2501%90
#
#SBATCH --mail-user=manuel.popp@wsl.ch
#SBATCH --mail-type=SUBMIT,END,FAIL

HABITAT=2_2_2.shp

RESIST="${HABITAT%.shp}"
THRESH=1.751681319215087

echo "Job ID: ${SLURM_JOB_ID}"
echo "Habitat type: ${HABITAT}"
echo "Resistance threshold: ${THRESH}"
echo "Resistance surface: ${RESIST}"

if [ -n "${SLURM_JOB_ID:-}" ] ; then
BSHPATH=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
else
BSHPATH=$(realpath "$0")
fi

MAINPATH=$(dirname "$(dirname "${BSHPATH}")")

if [ -n "${SLURM_JOB_ID:-}" ] ; then
  srun python3 $MAINPATH/py3/harary.py -hab ${HABITAT} -d ${THRESH} -r ${RESIST} -i ${SLURM_ARRAY_TASK_ID} -w /storage/poppman/shared/dami
else
	while getopts f:t: flag
	do
		case "${flag}" in
			f) FROM=${OPTARG};;
			t) TO=${OPTARG};;
		esac
	done
  for (( i=${FROM}; i<=${TO}; i++ ))
    do
      python3 $MAINPATH/py3/harary.py -hab ${HABITAT} -d ${THRESH} -r ${RESIST} -i ${i} -w /storage/poppman/shared/dami
    done
fi