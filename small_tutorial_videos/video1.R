require(mappoly)
setwd("~/repos/SCRI/small_tutorial_videos/")
dat <- read_fitpoly(file.in = "B2721_scores.dat", 
                    ploidy = 4,
                    parent1 = "Atlantic", 
                    parent2 = "B1829", 
                    verbose = TRUE)
plot(dat)
print(dat, detailed = TRUE)