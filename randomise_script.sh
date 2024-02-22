#!/bin/bash

module load fsl/6.0.5

dir="/imaging/correia/users/mc04/MRI_METHODS_TESTS/FreeWaterDiffusion/data/data_thr07/stats"

f=$sub

cd  $dir

randomise -i all_${f}_skeletonised.nii.gz -o results_${f} -m mean_FA_skeleton_mask -d design.mat -t design.con -n 5000 --T2 -V





