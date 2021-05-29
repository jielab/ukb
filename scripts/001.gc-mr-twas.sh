#!/bin/bash

dir=/mnt/d/data/gwas/bbj 
ldsc_dir=/mnt/d/software_lin/ldsc
gcta=/mnt/d/software/gcta/gtca64
traits_x="bilirubin bmi.female bmi bmi.male calcium ck creatinine crp dbp hba1c hdl ht ldl menarche menopause monocyte na plt pulsep rbc sbp tg ua wbc"
traits_y="arrhythmia asthma breastcancer cad cataract chf copd lungcancer osteoporosis pad ra stroke t2d"


## 第一步：LDSC (https://github.com/bulik/ldsc) ##
conda activate ldsc
for trait in $traits_x $traits_y; do
	zcat $trait.gz | sed 's/\t\t/\tNA\t/g' | gzip -f > $trait.$dat.chk # 确保没有确实数据
	python2 $ldsc_dir/munge_sumstats.py --sumstats $dir/bbj.$trait.gz --chunksize 10000 --snp rsids --a1 alt --a2 ref --frq maf --p pval --N 200000 --signed-sumstats beta,0 --merge-alleles $ldsc_dir/hm3.snplist --out $trait
done
for trait in $traits_x $traits_y; do
	echo process $trait
    echo $trait $traits_x $traits_y | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % /mnt/d/software_lin/ldsc/ldsc.py --rg % --out $trait.rg --ref-ld-chr $ldsc_dir/eur_w_ld_chr/ --w-ld-chr $ldsc_dir/eur_w_ld_chr/
    awk '$1=="Summary" {printf NR}' $trait.rg.log | xargs -n1 -I % awk -v s=% 'FNR >=s' $trait.rg.log | sed 's/.sumstats.gz//g' > $trait.rg.txt
done


## 第二步：GSMR (https://cnsgenomics.com/software/gcta/#GSMR) ##
for trait in $traits_x $traits_y; do
	echo "SNP A1 A2 freq b se p N" > $trait.gcta.txt
	zcat $dir/bbj.$trait.gz | awk '{if(NR==1) print "SNP a"; else {if (arr[$5] !="Y") print $5, $4,$3,$10, $8,$9,$7, 200000; arr[$5]="Y"}' | fgrep -v NA  >> $trait.gcta.txt
done
for i in $traits_x; do
    echo "$i $i.gcta.txt" > $trait.exposure
	echo -e "
	source('/mnt/d/scripts/library/gsmr_plot.r')
	gsmr_data = read_gsmr_data('$i.eff_plot.gz')
	pdf('$i.pdf')
	" > $i.plot.R
    echo -n > $trait.outcome
    for j in $traits_y; do
        if [[ $i != $j ]]; then
            echo "$j $j.gcta.txt" >> $i.outcome
			echo "plot_gsmr_effect(gsmr_data, '$i', '$j', colors()[75])" >> $i.plot.R
        fi
    done
	echo "dev.off()" >> $i.plot.R
	echo run gsmr on $i vs. $traits_y
    gcta64 --bfile /mnt/d/data/hm3/hm3.b37 --gsmr-file $i.exposure $i.outcome --gsmr-direction 2 --diff-freq 0.3 --gwas-thresh 5e-8 --effect-plot --out $i
  	Rscript $i.plot.R
	awk -v t=$i 'NR==1 && $1==t' $i.gsmr > $i.way1.gsmr
	awk -v t=$i 'NR==1 && $2==t' $i.gsmr > $i.way2.gsmr
done
fgrep Error *log
cat *.gsmr | sort | uniq > 001.gsmr.txt
corrplot(M, method='color', col=col(200), type='upper', order='hclust', addCoef.col='black', tl.col='black', tl.srt=45, p.mat=p.mat, sig.level=0.01, diag=F)


## FUSION TWAS ##
## 安装 TWAS (http://gusevlab.org/projects/fusion/)
dir_tw=/mnt/d/data/twas_data
dir_gt=$dir_tw/GTEx_v7_multi_tissue
dir_ld=$dir_tw/LDREF
fusion=/mnt/d/software_lin/fusion_twas
for trait in RHR T2D; do
for tissue in `ls -d1 $dir_gt/*/ | sed 's/\/$//' | awk -F '/' '{print $NF}' | awk '{printf " "$1}'`; do
for chr in 7; do
    echo now process trait $trait, tissue $tissue, chr $chr
    Rscript $fusion/FUSION.assoc_test.R --sumstats $dir/summary/$trait.sumstats.gz --chr $chr --out $trait.$tissue.chr$chr.txt --weights $dir_gt/$tissue.P01.pos --weights_dir $dir_gt --ref_ld_chr $dir_ld/1000G.EUR.
done
done
done
