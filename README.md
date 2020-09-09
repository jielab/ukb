# 表型数据提取以及GWAS分析流程. 

Author: Jie Huang, MD, PhD, Department of Global Health, Peking University School of Public Health



# #1.  提取一般表型数据（age, sex, race, bmi, etc.）

WINDOWS电脑建议安装系统自带的Ubuntu Linux系统，cd /mnt/d/。下载UKB小程序ukbunpack, unbconv, encoding.ukb 等。
苹果电脑，参考 https://github.com/spiros/docker-ukbiobank-utils。
打开ukbiobank.ac.uk, 点击中间的 Data Showcase 菜单。然后点击第一个“Essential Information”，阅读 Access and using your data。
写一个 VIP.fields.txt 文件，列出想提取的变量和对应的 data-field，比如 21022 age

```
unkunpack ukb42156.enc 【数据密码】
awk '{print $1}' ukb.vip.fields > ukb.vip.fields.id
sort ukb.vip.fields.id | uniq -d
ukbconv ukb42156.enc_ukb r -iukb.vip.fields.id -ovip
```

打开R ，用下面的几行代码，将上面生成的 vip.tab 数据读入，并且给每个变量赋予正确的名字。
下面的XXXX是文件路径，上述Linux 系统生成的 vip.r文件，如果在Windows 系统里面运行R，需要修改文件的路径。

```
source("D:/XXXX/vip.r")
pnames <- read.table("D:/XXXX/ukb.vip.fields", header=F)
pnames$V1 <- paste0("f.", pnames$V1, ".0.0")
phe <- subset(bd, select=grep("f.eid|\\.0\\.0", names(bd)))

```


# #2. 提取ICD数据（data field 42170）

http://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=131492。如果想自己动手来弄（不建议），可以用下面的通用代码来生成。

```
# ICD 这样的指标，包含了很多不同时间的时间点，量很大，建议分开来处理。
ukbconv ukb42156.enc_ukb r -s42170 -oicd
sed -i 's/"//g icd.tab

# 将 icd.tab 文件整合为两列，便于读入R。最后的 sed 命令将全部负数的人的记录替换成 "n"。
cat icd.tab | sed -e 's/\tNA//g' -e 's/\t/,/2g' | \
awk '{ if(NR==1) print "IID icd"; else if (NF==1) print $1 " NA"; else print $0"," }' | \
sed -e 's/-\w\+,/n/g' -e 's/n\+/n/g' > icd.2cols

# 从 ICD.2cols 文件里面提取某一个变量，比如 bipolar（对应的ICD-10代码F31），用R读入数据后，生成一个 0/1/NA 变量。
phe$icd_bipolar = ifelse("F31", phe$icd10), 1, ifelse(“F”, phe$icd10), NA, 0))

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
awk -v cn=$cnt -v co="J440" '{if (NR==1) print "IID", co; else {c=(cn-1)/2; printf $1;  
    for (i=2; i<=(c+1); i++) { if ($i==co) printf " "$(i+c) } printf "\n"  }}' icd-date.tab | awk ‘NF==2’ > icd-date.2cols

# 对于有一个不同的ICD-Date 的表型，比如 dementia 有5个ICD 代码“F00|F01|F02|F03|G30”，可以按照上述方法分别生成5个文件，比如 icdDate.F00.2cols, icdDate.F01.2cols，等。
然后在R里面合并这些文件，并找出每人的最小的日期。  

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


#6.1. 从千人基因组网站下载基因数据， 作为LD计算的参考。我们将千人基因组数据简称为 g1k.
之所有不简称为 1kg，是因为有的软件要求变量不能以数字开头。

```
打开 https://www.internationalgenome.org/data，在 Available data 下面，点击该页面 Phase 3 对应的 VCF 链接，可以看到以 “ALL.” 开头的文件，可以一个一个直接点击链接下载。
也可以用下面的命令下载, 并且随之将下载的VCF文件转换为PLINK格式

由于 chrX, chrY, chrMT 的文件名字跟其它染色体不同，用下面的命令下载的时候，里面的文件名字也需要相应调整。 

for chr in {1..22}; do
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
  plink --vcf ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out chr$chr
done  

除了下载上述页面上以 “ALL.” 开头的 VCF 文件，倒数第二个 integrated_call_samples_v3.20130502.ALL.panel 文件罗列了每一个样本的人群（pop）和人种 (super_pop)，以及性别。
根据这个文件，可以提取特定人种的样本，比如：
awk '$3=="EUR" {print $1,$1}' integrated_call_samples_v3.20130502.ALL.panel > g1k.EUR.keep
awk '$3=="EAS" {print $1,$1}' integrated_call_samples_v3.20130502.ALL.panel > g1k.EAS.keep

然后可以用 PLINK生成某一个特定人种的基因数据。下面的这个PLINK命令，可以写到上面的 PLINK 命令的下面，就不用写两次 "for chr in {1..22}" 了。
当然，也可以不用单独生成每个人种（比如EUR）的PLINK格式 基因数据，直接用一个包含所有样本的基因数据就可以了，只需要记住在相关的命令里面再写上 --keep g1k.EUR.keep（或对应的人种）。
for chr in {1..22}; do
  plink --bfile chr$chr --keep g1k.EUR.keep --make-bed --out EUR.chr$chr
done

不论是有所有2504个人基因数据的PLINK文件，还是只有某一个人种的PLINK文件，每个染色体都是单独的文件。
后续跑 GWAS 或提取 PRS 的时候，也是每条染色体的数据分开来跑，这样就可以进行并行计算（parallel computing）。
一般不建议把所有的染色体的数据合并成一个完整的单一的基因数据，毕竟8千多万个SNP，文件太大了，很多软件根本运行不了。

三个注意事项：
1. 其实，PLINK的网站 https://www.cog-genomics.org/plink/2.0/ 上也有千人基因组的数据，点击左下方菜单“1000 genomes phase3” 链接，按照操作下载和处理。
2. 不管是哪个方法得到的PLINK格式的数据，有的软件不允许 .bim 文件里面的 SNP 名字有重复，这个时候可以把原来的 .bim 文件备份一下，然后生成新的没有重复名字的 .bam 文件
	awk '{if(array[$2]=="Y") {i++; $2=$2".DUP"i}; print $0; array[$2]="Y"}' chr1.bim.COPY > chr1.bim 
3. g1k总共有8千多万个SNP，数据太大。如果 GSMR分析的GWAS文件总共才2百万个SNP，如果 GSMR --mbfile 命令太慢，可以考虑先用 plink --extract 命令生成一个只有那2百万个SNP的 g1k.new数据。 

```


#6.2. 提取 significant 信号，添加简单的注释

```
plink --annotate MY.gwas.txt NA attrib=snp129.attrib.txt ranges=glist-hg19 --border 10 --pfilter 5e-8 --out MY.gwas.top
```


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
    echo $trait $traits | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | \
    	xargs -n1 -I % /mnt/d/software_lin/ldsc/ldsc.py --rg % --out $trait --ref-ld-chr $ld_dir --w-ld-chr $ld_dir
 #   awk '$1=="Summary" {printf NR}' $trait.log | xargs -n1 -I % awk -v s=% 'FNR >=s' *.log | sed 's/.sumstats.gz//g' > $trait.ldsc.txt
done
```


#7.4.	多基因风险评分PRS

```
PRSice: https://github.com/choishingwan/PRSice
LDpred2 https://privefl.github.io/bigsnpr/articles/LDpred2.html
```


#7.5.	因果分析 Mendelian Randomization，GSMR (https://cnsgenomics.com/software/gcta/#GSMR)

GTCA 里面的 GSMR也需要用到上述提取的 g1k 基因数据作为 计算LD 的参考。
由于上述的 gak 数据是按照染色体分开的20多个数据，这个时候就需要用 --mbile （而不是 --bfile）来表明需要读取多个（multiple）bfile。
可以用这个命令生成一个 bfile.list 然后用到下面的命令里
```
seq 1 22 | xargs -n1 -I % echo chr% > bfile.list
```

```
for trait in RHR; do
    echo "$trait $dir/$trait.gcta.txt" > $trait.exposure
    echo -n > $trait.outcome
    for trait2 in $traits; do
        if [[ $trait2 != $trait ]]; then
            echo "$trait2 $dir/$trait2.gcta.txt" >> $trait.outcome
        fi
    done
   gcta64 --mbfile bfile.list --gsmr-file $trait.exposure $trait.outcome --gsmr-direction 2 --gwas-thresh 5e-8 --effect-plot --out $trait
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
    Rscript $fusion/FUSION.assoc_test.R --sumstats $dir/summary/$trait.sumstats.gz --chr $chr --out $trait.$tissue.chr$chr.txt \ 
    	--weights $dir_gt/$tissue.P01.pos --weights_dir $dir_gt --ref_ld_chr $dir_ld/1000G.EUR.
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

