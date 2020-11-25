require(fitpoly)
address <- "https://github.com/mmollina/B2721_map/raw/master/cluster_call/B2721_x_y.csv"
tempfl <- tempfile(pattern = "B2721_CC_", fileext = ".csv")
download.file(url = address, destfile = tempfl)
dat.raw<-read.csv(file = tempfl) ## Supplemenntary File S3
dat.raw[1:10,1:10]
dim(dat.raw)
# removing unused columns
datxy <- dat.raw[,-c(1,3,4)]
rownames(datxy) <- dat.raw$Name
# splitting data in two matrices, one containing x information and other containing y
X<-grep(pattern = "X", colnames(datxy))
Y<-grep(pattern = "Y", colnames(datxy))
X<-t(datxy[,X])
Y<-t(datxy[,Y])
rownames(X) <- stringr::str_remove(rownames(X), ".X")
rownames(Y) <- stringr::str_remove(rownames(Y), ".Y")
X[1:10, 1:10]
Y[1:10, 1:10]
n.ind <- nrow(Y)-4
# Formatting for fitPoly
Xm <- reshape2::melt(X, id.vars = c("MarkerName", "SampleName"))
Ym <- reshape2::melt(Y, id.vars = c("MarkerName", "SampleName"))
Xm <- Xm[,c(2,1,3)]
Ym <- Ym[,c(2,1,3)]
colnames(Xm) <- c("MarkerName", "SampleName", "X")
colnames(Ym) <- c("MarkerName", "SampleName", "Y")
mydata<-dplyr::left_join(Xm, Ym, by = c("MarkerName", "SampleName")) %>% mutate(R = X + Y, ratio = Y/(X+Y))
head(mydata)
mypop<-data.frame(SampleName = as.character(unique(mydata$SampleName)),
                  Population = c("pop1P1", "pop1P1", "pop1P2", "pop1P2", 
                                 rep("pop1", n.ind)))
mypar<-data.frame(popID = c("pop1", "pop1P1", "pop1P2"), 
                  parent1 = c("pop1P1",NA, NA),
                  parent2 = c("pop1P2", NA, NA))
# Running fitPoly
fitPoly::saveMarkerModels(ploidy=4,
                          data=mydata,
                          population=mypop,
                          pop.parents=mypar,
                          filePrefix = "B2721", 
                          ncores = 12)