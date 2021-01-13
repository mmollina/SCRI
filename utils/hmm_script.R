p <- 4
require(mappoly) 
n.mrk<-15
h.temp<-sim_homologous(m = 4, 
                       n.mrk = n.mrk, 
                       max.d = 2,
                       seed = 45637)

#h.temp$hom.allele.p <- lapply(h.temp$hom.allele.p, function(x){
#  if(any(x==0)) return(0)
#  return(1:length(x))
#})
#h.temp$hom.allele.q <- lapply(h.temp$hom.allele.q, function(x){
#  if(any(x==0)) return(0)
#  return(1:length(x))
#})

dat<-poly_cross_simulate(m = 4, 
                         rf.vec = mf_h(1), 
                         n.mrk=n.mrk,
                         n.ind=200, 
                         h.temp, 
                         seed = 653543, 
                         draw = TRUE,
                         file = "out_sim.pdf")
plot(dat)
s <- make_seq_mappoly(dat, "all")
tpt <- est_pairwise_rf(s)
map <- est_rf_hmm_sequential(input.seq = s,
                             thres.twopt = 5,
                             thres.hmm = 10,
                             extend.tail = 10,
                             tol = 0.1,
                             tol.final = 10e-5,
                             phase.number.limit = Inf,
                             twopt = tpt,
                             verbose = TRUE)
plot(map, mrk.names = TRUE)
map.err <- est_full_hmm_with_global_error(map, error = 0.01, verbose = T)

genoprob <- calc_genoprob_error(map.err, error = 0.01)

z<-round(genoprob$probs[,,1], 2)

#### Number of recombinant events ####
gm <- combn(1:p, p/2)
L <- matrix(NA, ncol(gm), ncol(gm))
for(i in 1:ncol(gm))
  for(j in 1:ncol(gm))
    L[i,j] <- sum(gm[,i]%in%gm[,j])
#image(L)  
#image(kronecker(L, L, FUN = "+"))

#### Transition function ####
transition <- function(r,p,L){
  T <- matrix(NA, ncol(gm), ncol(gm))
  for(i in 1:ncol(gm))
    for(j in 1:ncol(gm))
      T[i,j] <- ((1-r)^(p/2-L[i,j]) * r^L[i,j])/choose(p/2,L[i,j])
  kronecker(T,T)
}
#image(transition(0.5,p,L))

#### Emission ####
emission <- function(p1, p2, phaseP1, phaseP2){
  I <- NULL
  for(i in 1:length(p1)){
    x2<-x1<-numeric(p)
    if(p1[i]!=0)
    x1[phaseP1[[i]]] <- 1
    if(p2[i]!=0)
    x2[phaseP2[[i]]] <- 1
    I<-cbind(I, kronecker(apply(combn(x1, p/2), 2, sum),
                          apply(combn(x2, p/2), 2, sum), FUN = "+"))
  }
  I
}
p1 <- dat$dosage.p
p2 <- dat$dosage.q
I <- emission(p1, p2, map$maps[[1]]$seq.ph$P, map$maps[[1]]$seq.ph$Q)
g <- dat$geno.dose


ind<-1
#par(ask = T)
#for(ind in 1:50)
{
  cat(ind, "\n")
  par(bg = "black")
  plot(0, xlim = c(1,n.mrk), ylim = c(1,36), type = "n", 
       axes = F, xlab = "", ylab = "")
  for(k in 1:(n.mrk-1)){
    for(i in 1:36){
      for(j in 1:36){
        lines(x = c(k,k+1), y = c(i,j), 
              lwd = ifelse((I[i,k] == g[k, ind] &
                           I[j,k+1] == g[k+1, ind]), 0.5, 0.0), col = "white")
              #lwd = .2)
              #lwd = genoprob$probs[i,k,ind] * genoprob$probs[j,k+1,ind] * 3)
      }
    }
  }
  for(k in 1:n.mrk){
    points(rep(k,36), 1:36 ,
           pch = 19,
           #col = "red", cex = 2)
           col = ifelse(I[,k]==g[k, ind], "#e41a1c", "lightblue"), cex = 1)
  }
  
}
par(ask= T)
for(ind in 63:200){
  z<-round(genoprob$probs[,,ind], 1)
  image(z, main = ind)
}
dev.off()
# Make several .png images


ind<-64
g[,64]

z<-round(genoprob$probs[,,ind], 1)
system("rm *.png")
png(file="example%03d.png", width=1400, heigh=700)
cte <- 0
for(a in c(0.4972222, 0.15, rep(10,15)^-c(1:15))){
  a1<-0.5-a
  z1<-z
  #z1[]<-0.51
  z1[z1 == 0] <- a/34
  z1[z1 == 0.5] <- a1
  #par(bg = "black")
  plot(-0.5, xlim = c(1,n.mrk), ylim = c(1,36), type = "n", 
       axes = F, xlab = "", ylab = "", 
       main = paste0("EM iteration = ", cte),cex.main = 2.5)
  for(k in 1:(n.mrk-1)){
    for(i in 1:36){
      for(j in 1:36){
        lines(x = c(k,k+1), y = c(i,j), 
              lwd = ifelse((I[i,k] == g[k, ind] &
                              I[j,k+1] == g[k+1, ind]), 3, 0) * 1/abs(log(z1[i,k] * z1[j,k+1])), 
              col = "black") 
        #lwd = .2)
        #lwd = genoprob$probs[i,k,ind] * genoprob$probs[j,k+1,ind] * 3)
      }
    }
  }
  for(k in 1:n.mrk){
    points(rep(k,36), 1:36 ,
           pch = 19,
           #col = "red", cex = 2)
           col = ifelse(I[,k]==g[k, ind], "#e41a1c", "darkgray"), cex = 1)
    
  }
  text(x = rep(0.7,36), y = 1:36, labels = dimnames(z)[[1]], cex = 1.2)
  text(x = 1:15, y= rep(36.8, 15), labels = paste0("M",1:15))
  
  cte<-cte +1
}
dev.off()
system("convert -delay 30 *.png animated_hmm.gif")

png(file="fig0.png", width=1400, heigh=700)
plot(-0.5, xlim = c(1,n.mrk), ylim = c(1,36), type = "n", 
     axes = F, xlab = "", ylab = "", 
     #main = paste0("EM iteration = ", cte),cex.main = 2.5)
     main = "Initial states", cex.main = 2)
for(k in 1:(n.mrk-1)){
  for(i in 1:36){
    for(j in 1:36){
      lines(x = c(k,k+1), y = c(i,j), 
            lwd = ifelse((I[i,k] == g[k, ind] &
                            I[j,k+1] == g[k+1, ind]), 0.5, 0), #* 
            #1/abs(log(z1[i,k] * z1[j,k+1])), 
            col = "black") 
      #lwd = .2)
      #lwd = genoprob$probs[i,k,ind] * genoprob$probs[j,k+1,ind] * 3)
    }
  }
}
for(k in 1:n.mrk){
  points(rep(k,36), 1:36 ,
         pch = 19,
         #col = "red", cex = 2)
         col = ifelse(I[,k]==g[k, ind], "#e41a1c", "darkgray"), cex = 1)
  
}
text(x = rep(0.7,36), y = 1:36, labels = dimnames(z)[[1]], cex = 1.2)
text(x = 1:15, y= rep(36.8, 15), labels = paste0("M",1:15))
dev.off()

