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
polymapR.data<-export_data_to_polymapR(x)
polymapR.data[1:5, 1:10]
