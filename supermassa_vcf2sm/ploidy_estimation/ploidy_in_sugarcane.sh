# cd ~/repos/SCRI/supermassa_vcf2sm/ploidy_estimation/
python2.7 ../supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 6:12 --naive_posterior_reporting_threshold 0.8 --file SugSNP_382.pro --f1_parent_data SugSNP_382.par > SugSNP_382.html > SugSNP_382.html

# --inference
# f1 or hw for respective biparental or Hardy-Weinberg model; 
# ploidy will use exclusively the ratio between abundances

# --naive_posterior_reporting_threshold
# In the example, we use -n 0.80 to only keep indiciduals with 
# associated probability of 0.80 or higher.

# Windows versison
# C:\Python27\python.exe ../supermassa/src/SuperMASSA.py  --print_genotypes --inference f1 --ploidy_range 6:12 --naive_posterior_reporting_threshold 0.8 --file SugSNP_382.pro --f1_parent_data SugSNP_382.par  > SugSNP_382.html

#python2.7 ../supermassa/src/SuperMASSA.py \
#       --print_genotypes \
#       --inference f1 \
#       --ploidy_range 6:12 \
#       --naive_posterior_reporting_threshold 0.8 \
#       --file SugSNP_382.pro \
#      --f1_parent_data SugSNP_382.par > SugSNP_382.html
       
#sh ploidy_in_sugarcane.sh 
