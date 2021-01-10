grep "^#" chr1.vcf > header.vcf
grep "S3_" chr1.vcf > tmp.vcf
shuf -n 1000 tmp.vcf > output.vcf
cat header.vcf output.vcf > sample_ch1.vcf

grep "^#" chr2.vcf > header.vcf
grep "S9_" chr2.vcf > tmp.vcf
shuf -n 1000 tmp.vcf > output.vcf
cat header.vcf output.vcf > sample_ch2.vcf

grep "^#" chr3.vcf > header.vcf
grep "S12_" chr3.vcf > tmp.vcf
shuf -n 1000 tmp.vcf > output.vcf
cat header.vcf output.vcf > sample_ch3.vcf