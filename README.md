# 表型数据提取以及GWAS分析流程. 

Author: Jie Huang, MD, PhD, Department of Global Health, Peking University School of Public Health



# #1.  提取一般表型数据（age, sex, race, bmi, etc.）

WINDOWS电脑建议安装系统自带的Ubuntu Linux系统，cd /mnt/d/。下载UKB小程序ukbunpack, unbconv, encoding.ukb 等。
苹果电脑，参考 https://github.com/spiros/docker-ukbiobank-utils。
打开ukbiobank.ac.uk, 点击中间的 Data Showcase 菜单。然后点击第一个“Essential Information”，阅读 Access and using your data。
写一个 VIP.fields.txt 文件，列出想提取的变量和对应的 data-field，比如 21022 age

```
awk '{print $1}' VIP.fields.txt > VIP.fields.ids
unkunpack ukb42156.enc 【数据密码】
ukbconv ukb42156.enc_ukb r -iVIP.fields.ids -oVIP
```

打开R ，用下面的几行代码，将上面生成的VIP.tab 数据读入，并且给每个变量赋予正确的名字。
下面的XXXX是文件路径，上述Linux 系统生成的 VIP.r文件，如果在Windows 系统里面运行R，需要修改文件的路径。

```
source("D:/XXXX/VIP.r")
pnames <- read.table("D:/XXXX/ukb.vip.fields", header=F)
pnames$V1 <- paste0("f.", pnames$V1, ".0.0")
phe <- subset(bd, select=grep("f.eid|\\.0\\.0", names(bd)))

```


# #2. 提取ICD数据（data field 42170）

http://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=131492。如果想自己动手来弄（不建议），可以用下面的通用代码来生成。

```
# ICD 这样的指标，包含了很多不同时间的时间点，量很大，建议分开来处理。
ukbconv ukb42156.enc_ukb r -s42170 -oicd
sed -i ‘s/”//g’ icd.tab

# 将 icd.tab 文件整合为两列，便于读入R。
cat icd.tab | sed -e 's/\tNA//g' -e 's/\t/,/2g' | awk '{ if(NR==1) print "IID\ticd10"; else if (NF==1) print $1 "\tNA"; else print $0 }' > icd.2cols

# 从 ICD.2cols 文件里面提取某一个变量，比如 bipolar（对应的ICD-10代码F31），用R读入数据后，生成一个 0/1/NA 变量。
phe$icd_bipolar = ifelse(“F31”, phe$icd10), 1, ifelse(“F”, phe$icd10), NA, 0))

# 如果需要批量处理很多ICD变量，先写一个 VIP.icd.txt 文件，第一列是ICD代码，第二列是相对应的变量的名字，比如I350|I35  stenosis，“|”表示“或者”。
# 这个文件第一行写上 codes names，然后用下面的R代码批量执行。
ICDnames <- read.table("ukb.vip.icd10", header=T)
for (i in 1:nrow(ICDnames)) { 
  phe[[paste0("icd_",ICDnames$names[i])]] = ifelse( grepl(ICDnames$codes[i], phe$icd10), 1, ifelse(grepl( substring(ICDnames$codes[i],1,1), phe$icd10), NA,0)) 
}

```


# #3. 提取ICD-Date 数据以及相对应的日期（data field 42180）


```
echo “41270\n41280” > vip.fields.txt
ukbconv ukb42156.enc_ukb r -ivip.fields.txt -oicd-date
sed -i ‘s/”//g’ icd-date.tab

# 提取单个ICD 的Date, 比如COPD  (代码J440，不是J44) 。
cnt=`head -1 icd-date.tab | awk '{printf NF}'` # 找出列数
awk -v cn=$cnt  -v co="J440" '{if (NR==1) print "IID", co; else {c=(cn-1)/2; printf $1;  for (i=2; i<=(c+1); i++) { if ($i==co) printf " "$(i+c) } printf "\n"  }}' icd-date.tab | awk ‘NF==2’ > icd-date.2cols

# 对于有一个不同的ICD-Date 的表型，比如 dementia 有5个ICD 代码“F00|F01|F02|F03|G30”，可以按照上述方法分别生成5个文件，比如 icdDate.F00.2cols, icdDate.F01.2cols，等。然后在R里面合并这些文件，并找出每人的最小的日期。  

```


# #4. 对表型数据进行 GWAS 运行之前的处理

提取需要研究的表型数据和相关的covariates，比如 age, sex, PCs。一般来说，quantitative的表型数据要 adjust for covariates 和转化成正态分布，这个可以在R里面用下面的命令来实现。
对于疾病的binary 表型，只需要把需要 adjust 的covarites 和表型数据放在同一个表型数据文件里面，然后在 GWAS里面的命令指明哪个是表型，哪些是 covariates。

```
trait_res = residuals(lm(trait ~ age+sex+PC1+PC2, na.action=na.exclude)
trait_inv = qnorm((rank(trait_res,na.last="keep")-0.5) / length(na.omit(trait_res)))
```


# #5. GWAS 运行
目前GWAS 由专人负责运行，以下链接可以随时下载公开的GWAS数据

```
Cardiovascular disease genomics http://www.broadcvdi.org/
fastgwa.info
```

# #6. GWAS 后续常规分析 


#6.1. 从千人基因组网站下载基因数据 https://www.internationalgenome.org/data 作为LD计算的参考。点击该页面 Phase 3 对应的VCF，下载所有以ALL开头的文件。可将下载后的VCF文件的名字改短为ALL.chr*.gz 这样的名字。然后用 PLINK 将 VCF格式转换为 PLINK格式 

```
for chr in {1..22}; do
  plink --vcf ALL.chr$chr.gz --make-bed --out chr$chr
done  
```


#6.2. 提取 significant 信号，添加简单的注释
	plink --annotate MY.gwas.txt NA attrib=snp129.attrib.txt ranges=glist-hg19 --border 10 --pfilter 5e-8 --out MY.gwas.top


#6.3. 提取 signifianct & independent 信号
用PLINK --clump 命令，如下。如果不考虑 LD, 只考虑距离，可以任意指定一个 LDfile同时设置  --clump-r2=0。

```
for chr in {1..22}; do
  plink --bfile chr$chr --clump Height.2018.txt --clump-p1 1e-08 --clump-p2 1e-8 --clump-kb 1000 --clump-r2 0 --out chr$chr
done  
```

#6.4. 用上述方法生成一个比较小的文件，用R里面的qqman package，或者我写的mhplot.R代码，绘制Manhattan plot，也可以绘制QQ plot（不太常用）。


#6.5. 从GWAS catalog (https://www.ebi.ac.uk/gwas) ) 寻找已知信号，通过R 的plot()来比较该 GWAS跟已经发表过的信号的EAF和BETA的一致性。


#6.6. 对于有统计显著性的重点locus，可以ZOOM 画图 (http://locuszoom.org)


 
# #7. GWAS的深度分析 


#7.1.	SNP频率和基本注解查询
```
GnomAD https://gnomad.broadinstitute.org
```

#7.2.	GWAS数据的功能性注释

```
post-GWAS analysis pipeline (github.com/Ensembl/postgap).
```

#7.3.	多个GWAS 之间的 genetic correlation 分析, LDSC (https://github.com/bulik/ldsc)

```
source activate ldsc
for trait in $traits; do
	zcat $dir/gwas/$trait.sumstats.gz | awk 'NF==5' | sed -e 's/\.000$//' -e 's/\t/ /g' | gzip -f > $trait.sumstats.gz 
done
for train in $traits; do
    echo $trait $traits | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % /mnt/d/software_lin/ldsc/ldsc.py --rg % --out $trait --ref-ld-chr $ld_dir --w-ld-chr $ld_dir
 #   awk '$1=="Summary" {printf NR}' $trait.log | xargs -n1 -I % awk -v s=% 'FNR >=s' *.log | sed 's/.sumstats.gz//g' > $trait.ldsc.txt
done
```


#7.4.	多基因风险评分PRS

```
PRSice: https://github.com/choishingwan/PRSice
LDpred2 https://privefl.github.io/bigsnpr/articles/LDpred2.html
```


#7.5.	因果分析 Mendelian Randomization，GSMR (https://cnsgenomics.com/software/gsmr)

```
for trait in RHR; do
    echo "$trait $dir/$trait.gcta.txt" > $trait.exposure
    echo -n > $trait.outcome
    for trait2 in $traits; do
        if [[ $trait2 != $trait ]]; then
            echo "$trait2 $dir/$trait2.gcta.txt" >> $trait.outcome
        fi
    done
   gcta64 --bfile g1k.EUR/chr$chr --gsmr-file $trait.exposure $trait.outcome --gsmr-direction 2 --gwas-thresh 5e-8 --effect-plot --out $trait
done
```


#7.6. TWAS (http://gusevlab.org/projects/fusion/)

```
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
```


# #8. 参考文献

```
2018. Adult height and risk of 50 diseases: a combined epidemiological and genetic analysis
 
2019. JACC. Genome-Wide Assessment for Resting Heart Rate and Shared Genetics With Cardiometabolic Traits and Type 2 Diabetes
Genome Wide Assessment of Shared Genetic Architecture Between Rheumatoid Arthritis and Cardiovascular Diseases Using the UK Biobank Data 

2019. JAMA. Association of Lifestyle and Genetic Risk With Incidence of Dementia
```

