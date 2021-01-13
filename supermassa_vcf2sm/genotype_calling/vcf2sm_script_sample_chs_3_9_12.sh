# take 2.5 minutes

# cd ~/repos/SCRI/supermassa_vcf2sm/genotype_calling/

# When using windows
# C:\Python27\python.exe VCF2SM.py -i sample_vcf/sample_ch+.vcf -o sample_trifida_sm_chr+.vcf --sF 1 --eF 3 -d 20 -D 2000 -g BT -1 Beauregard -2 Tanzania -S "../supermassa/src/SuperMASSA.py" -I f1 -M 2:6 -f 6 -p 0.80 -n 0.75 -c 0.75 -t 2

# python2.7 VCF2SM.py -i sample_vcf/sample_ch+.vcf -o sample_trifida_sm_chr+.vcf --sF 1 --eF 3 -d 20 -D 2000 -g BT -1 Beauregard -2 Tanzania -S "../supermassa/src/SuperMASSA.py" -I f1 -M 2:6 -f 6 -p 0.80 -n 0.75 -c 0.75 -t 2

python2.7 VCF2SM.py -i sample_vcf/sample_ch+.vcf \
                    -o sample_trifida_sm_chr+.vcf \
                    --sF 1 \
                    --eF 3 \
                    -d 20 \
                    -D 2000 \
                    -g BT \
                    -1 Beauregard \
                    -2 Tanzania \
                    -S ../supermassa/src/SuperMASSA.py \
                    -I hw \
                    -M 2:8:1 \
                    -f 6 \
                    -p 0.80 \
                    -n 0.75 \
                    -c 0.75 \
                    -t 15
