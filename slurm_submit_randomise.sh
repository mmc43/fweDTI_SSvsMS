#!/bin/bash


#for f in {fw_SS_FA,fw_SS_MD,fw_SS_FW,fw_MS_FA,fw_MS_MD,fw_MS_FW}
for f in fw_MS_FW
do

sbatch --export=sub=$f randomise_script.sh

done

