## Reading results from VCF2SM
require(mappoly)
x1<-read_vcf(file.in = "~/repos/SCRI/supermassa_vcf2sm/genotype_calling/sample_trifida_sm_chr1.vcf", 
             parent.1 = "PARENT1", 
             parent.2 = "PARENT2")
x2<-read_vcf(file.in = "~/repos/SCRI/supermassa_vcf2sm/genotype_calling/sample_trifida_sm_chr2.vcf", 
             parent.1 = "PARENT1", 
             parent.2 = "PARENT2")
x3<-read_vcf(file.in = "~/repos/SCRI/supermassa_vcf2sm/genotype_calling/sample_trifida_sm_chr3.vcf", 
             parent.1 = "PARENT1", 
             parent.2 = "PARENT2")

x<-merge_datasets(x1,x2)
x<-merge_datasets(x,x3)
plot(x)
print(x, detailed = TRUE)

# sx1<-read_vcf(file.in = "~/repos/SCRI/vcf2sm/sample_trifida_sm_chr1.vcf", parent.1 = "PARENT1", parent.2 = "PARENT2")
# sx2<-read_vcf(file.in = "~/repos/SCRI/vcf2sm/sample_trifida_sm_chr2.vcf", parent.1 = "PARENT1", parent.2 = "PARENT2")
# sx3<-read_vcf(file.in = "~/repos/SCRI/vcf2sm/sample_trifida_sm_chr3.vcf", parent.1 = "PARENT1", parent.2 = "PARENT2")
# 
# sx<-merge_datasets(sx1,sx2)
# sx<-merge_datasets(sx,sx3)
# plot(sx)
# print(sx, detailed = TRUE)
