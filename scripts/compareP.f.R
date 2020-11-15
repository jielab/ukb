compareP <-function(
	f1, f1_name, f1_snp, f1_ea, f1_nea, f1_eaf, f1_beta, f1_p,
	f2, f2_name, f2_snp, f2_ea, f2_nea, f2_eaf, f2_beta, f2_p
	){
	
	if(grepl(".gz$",f1)) dat1 <- read.table(gzfile(f1,'r'), header=T, as.is=T) else dat1 <- read.table(f1, header=T, as.is=T)
	if(grepl(".gz$",f2)) dat2 <- read.table(gzfile(f2,'r'), header=T, as.is=T) else dat2 <- read.table(f2, header=T, as.is=T)
	dat1$POS <- NULL; dat2$POS <- NULL; # avoid "POS" is used for "P"
	
	if (is.null(dat1$NEA)) dat1$NEA=NA;
	if (is.null(dat2$NEA)) dat2$NEA=NA;
	if (is.null(dat1$EAF)) dat1$EAF=0;
	if (is.null(dat2$EAF)) dat2$EAF=0;
	if (is.null(dat1$BETA) & !is.null(dat1$OR)) dat1$BETA=log(dat1$OR)
	if (is.null(dat2$BETA) & !is.null(dat2$OR)) dat2$BETA=log(dat2$OR)
	if (is.null(dat1$P)) dat1$P=NA; dat1$P <- -log10(dat1$P)
	if (is.null(dat2$P)) dat2$P=NA; dat2$P <- -log10(dat2$P)
	
	dat1 <- subset(dat1, select=c(f1_snp, f1_ea, f1_nea, f1_eaf, f1_beta, f1_p))
	dat2 <- subset(dat2, select=c(f2_snp, f2_ea, f2_nea, f2_eaf, f2_beta, f2_p))	
	names(dat1) <- paste0( c("SNP", "EA", "NEA", "EAF", "BETA", "P"), "_1")
	names(dat2) <- paste0( c("SNP", "EA", "NEA", "EAF", "BETA", "P"), "_2")
	dat <- merge(dat1, dat2, by.x="SNP_1", by.y="SNP_2")
	# dat <- na.omit(dat)
	for (var in c('EAF_1', 'EAF_2', 'BETA_1', 'BETA_2', 'P_1', 'P_2')) {
		dat[[var]] <- as.numeric(dat[[var]])
	}
	
	#dat <- subset(dat, (EA_1==EA_2 & NEA_1==NEA_2) | (EA_1==NEA_2 & NEA_1==EA_2))
	dat[which(dat$EA_1 != dat$EA_2), "EAF_1"]  <- 1- dat[which(dat$EA_1 != dat$EA_2), "EAF_1"]
	dat[which(dat$EA_1 != dat$EA_2), "BETA_1"] <- 0- dat[which(dat$EA_1 != dat$EA_2), "BETA_1"]
	# dat <- subset(dat, P_1 <0.05 & P_2<0.05); # subset(dat, abs(EAF_1 - EAF_2) <0.4)

	for (var in c("EAF", "BETA", "P")) { # unlist(strsplit(plots," "))
		dat1 <- subset(dat, select=c(paste0(var,'_1'), paste0(var,'_2')))
		dat1 <- na.omit(dat1)
		if (nrow(dat1)>0) {
			plot(dat1, main=paste0(var, " (N=", nrow(dat), ")"), xlab=f1_name, ylab=f2_name, cex=2, cex.lab=2, cex.axis=2, cex.main=2, font=2, font.lab=2)
		}
	}
}
