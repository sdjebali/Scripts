#!/bin/sh

     module load bioinfo/bcftools-1.9
     cd /work/project/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
     bcftools view 22.basic.vcf.gz | awk -v f=INFO -v m=0.5 -v M=1 -f /home/sdjebali/fragencode/tools/multi/Scripts/Awk/vcf_addbool_for_INFO_subfield.awk | bcftools view -O z -o 22.basic.INFOin0.5-1.vcf.gz

