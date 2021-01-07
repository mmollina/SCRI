# mkdir ~/bin
# PATH=~/bin:$PATH
# ln -s /usr/bin/python2.7 ~/bin/python

python supermassa/src/SuperMASSA.py \
       --print_genotypes\
       --inference f1 \
       --ploidy_range 2:12 \
       --file sugarcane_data/SugSNP204_progeny \
       --f1_parent_data sugarcane_data/SugSNP204_parents 
       
       