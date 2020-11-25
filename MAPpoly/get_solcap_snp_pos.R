#### Get potato infinium 8303 sequences ####
address <- "https://github.com/mmollina/SCRI_MAPpoly_training/raw/main/data/potato_8303SNPs_potato_dm_v4.03.gff3"
tempfl <- tempfile(pattern = "B2721_", fileext = ".gff3")
download.file(url = address, destfile = tempfl)
solcap.snp.pos <- ape::read.gff(file = tempfl)
head(solcap.snp.pos)
solcap.snp.pos$snp.names <- sapply(strsplit(solcap.snp.pos$attributes, split = ";|="), function(x) x[4])
solcap.snp.pos <- solcap.snp.pos[order(solcap.snp.pos$snp.names),]
solcap.snp.pos <- solcap.snp.pos[!duplicated(solcap.snp.pos$snp.names),]
solcap.snp.pos$ch <- sapply(strsplit(as.character(solcap.snp.pos$seqid), split = "chr"), function(x) as.numeric(x[2]))
solcap.snp.pos$ch[solcap.snp.pos$ch == 0] <- NA
solcap.snp.pos <- solcap.snp.pos[order(solcap.snp.pos[,"ch"],solcap.snp.pos[,"start"],decreasing=FALSE),]
solcap.snp.pos <- data.frame(chr = solcap.snp.pos$ch, pos = solcap.snp.pos$start, row.names = solcap.snp.pos$snp.names)
head(solcap.snp.pos)
write.csv(solcap.snp.pos, file = "~/repos/SCRI_MAPpoly_training/data/solcap_snp_pos_V4.03.csv")
