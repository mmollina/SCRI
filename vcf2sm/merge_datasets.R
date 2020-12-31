x1<-read_vcf(file.in = "~/repos/SCRI/vcf2sm/trifida_sm_chr1.vcf", parent.1 = "PARENT1", parent.2 = "PARENT2")
x2<-read_vcf(file.in = "~/repos/SCRI/vcf2sm/trifida_sm_chr2.vcf", parent.1 = "PARENT1", parent.2 = "PARENT2")
x3<-read_vcf(file.in = "~/repos/SCRI/vcf2sm/trifida_sm_chr3.vcf", parent.1 = "PARENT1", parent.2 = "PARENT2")

x<-merge_datasets(x1,x2)
x<-merge_datasets(x,x3)
plot(x)
print(x, detailed = TRUE)
