######## simulation ###########

library(frbinom)
library(bbmle)

N=10  #N
num<-400 # n
m<-20 #repeat
x<-runif(num,-2,2)
a1_r<- -.5; a2_r<-2; #psi
b1_r<- 1 ; b2_r<--2; # eta
c1_r<-2; c2_r<-1;  # nu
p_r<-1/(1+exp(-(a1_r+x*a2_r )))
h_r<-1/(1+exp(-(b1_r+x*b2_r )))

c00_r<- apply(rbind(p_r,h_r) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
c_r<-c00_r/(1+exp(-(c1_r+x*c2_r )));

p0<-p1<-p2<-h0<-h1<-h2<-c0<-c1<-c2<-c00<-NULL
EST<-matrix(0, ncol=6, nrow=m)

DATA_simu0<-apply(rbind(p_r, h_r, c_r),2, function(x){
  rfrbinom( m, N, x[1],x[2],x[3] )})

for(i in 1:m){
  DATA_simu<-DATA_simu0[i,]
  likelihoodfunction_simu<-function(p0,p1,h0,h1,c0,c1){
    p_l<-1/(1+exp(-(p0+x*p1 )));
    h_l<-1/(1+exp(-(h0+ x*h1  )));
    c00<- apply(rbind(p_l,h_l) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
    c_l<-c00/(1+exp(-(c0+x*c1)));
    DAT<-rbind(DATA_simu,p_l,h_l,c_l )
    -sum(apply(DAT,  2 , function(x){log(dfrbinom(x[1],N,x[2],x[3],x[4] ))}))
  }

  est<-mle2(likelihoodfunction_simu, method="L-BFGS-B",
            start=list(  p0 = 0,p1 = 0,
                         h0 = 0, h1 = 0, c0 =0, c1 = 0) ,
            lower= c( -5,-5,-5,-5,-5,-5 ),
            upper=c(5,5,5,5,5,5)
  )

  EST[i,]<-est@coef
}



round(mean(EST[,1])-(-.5), digits=2 ) ; round(sd(EST[,1]), digits=2 )
round(mean(EST[,2])-(2), digits=2 ) ; round(sd(EST[,2]), digits=2 )
round(mean(EST[,3])-(1), digits=2 ) ; round(sd(EST[,3]), digits=2 )
round(mean(EST[,4])-(-2), digits=2 ) ; round(sd(EST[,4]), digits=2 )
round(mean(EST[,5])-(2), digits=2 ) ; round(sd(EST[,5]), digits=2 )
round(mean(EST[,6])-(1), digits=2 ) ; round(sd(EST[,6]), digits=2 )


###### Data analysis ########

##### apple data ####

install.packages("agridat")
library(agridat)
data(ridout.appleshoots)
dat <- ridout.appleshoots

pho<-rep(0,270)
pho[dat$photo==16]<-1
########################################################
###### photoperiod : categorical, BAP: numerical #########
########################################################

## fractional binomial regression ###
library(fbglm)
x0<-data.frame( pho=pho, bap=dat$bap)
fbglm(y=dat$roots, x=x0   )
## extended zero-inflated negative binomial regression ##
ZINB2( y=dat$roots, x=x0   )
## ZINB
library(pscl)
M3 <- zeroinfl(dat$roots~ pho+dat$bap  |
                 pho+dat$bap,
               dist = 'negbin')
summary(M3)
# ZIP
M4 <- zeroinfl(dat$roots~ pho+dat$bap  |
                 pho+dat$bap,
               dist = 'pois')
summary(M4)
## Vuong's test ##
# FB VS negative binomial #
test(  y=dat$roots, x=x0, model1="fbglm", model2="ZINB2"   )
# FB VS ZINB #
test(  y=dat$roots, x=x0, model1="fbglm", model2="ZINB"   )
#FB VS ZIP #
test(  y=dat$roots, x=x0, model1="fbglm", model2="ZIP"   )



########################################################
###### photoperiod : categorical, BAP: categorical #########
########################################################

bap1<-rep(0,270)
bap1[which(dat$bap==4.4) ]<-1
bap2<-rep(0,270)
bap2[which(dat$bap==8.8) ]<-1
bap3<-rep(0,270)
bap3[which(dat$bap==17.6) ]<-1
x<-data.frame(pho=pho, bap1=bap1, bap2=bap2, bap3=bap3)

## fractional binomial regression ###
library(fbglm)
fbglm(y=dat$roots, x=x   )
## extended zero-inflated negative binomial regression ##
ZINB2( y=dat$roots, x=x   )
## ZINB
library(pscl)
M3 <- zeroinfl(dat$roots~ pho+bap1+bap2+bap3  |
                 pho+bap1+bap2+bap3,
               dist = 'negbin')
summary(M3)
# ZIP
M4 <- zeroinfl(dat$roots~ pho+bap1+bap2+bap3 |
                 pho+bap1+bap2+bap3,
               dist = 'pois')
summary(M4)
## Vuong's test ##
# FB VS negative binomial #
test(  y=dat$roots, x=x, model1="fbglm", model2="ZINB2"   )
# FB VS ZINB #
test(  y=dat$roots, x=x, model1="fbglm", model2="ZINB"   )
#FB VS ZIP #
test(  y=dat$roots, x=x, model1="fbglm", model2="ZIP"   )
