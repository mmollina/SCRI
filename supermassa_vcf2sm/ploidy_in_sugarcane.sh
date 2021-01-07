# mkdir ~/bin
# PATH=~/bin:$PATH
# ln -s /usr/bin/python2.7 ~/bin/python

# python supermassa/src/SuperMASSA.py --print_genotypes --inference f1 --ploidy_range 2:12 --naive_posterior_reporting_threshold 0.8 --file SugSNP_382.pro --f1_parent_data SugSNP_382.par > SugSNP_382.html

python supermassa/src/SuperMASSA.py \
       --print_genotypes \
       --inference f1 \
       --ploidy_range 2:12 \
       --naive_posterior_reporting_threshold 0.8 \
       --file SugSNP_382.pro \
       --f1_parent_data SugSNP_382.par 
       
#sh ploidy_in_sugarcane.sh > SugSNP_382.html

       