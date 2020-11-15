setwd("C:/Users/黄捷/Desktop")
source("D:/scripts/compareP.f.R")

png("compareP.png", w=1200, h=400)
par(mfrow=c(1,3), mai=c(1,1,0.5,0.5))
compareP(
	f1="cad.txt",
	f1_name="Known", f1_snp="SNP", f1_ea="EA", f1_nea="NEA", f1_eaf="EAF", f1_beta="BETA", f1_p="P",
	f2="cvd.gwas.gz",
	f2_name="UKB", f2_snp="SNP", f2_ea="A1", f2_nea="A2", f2_eaf="A1_FREQ", f2_beta="BETA", f2_p="P"
	)
dev.off()







