# install.packages("devtools")
devtools::install_github("leejeo25/fbglm")
######## Section 4. Simulation ###########
###### subsection 4.1
library(frbinom)
library(bbmle)
# first simulation
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

# second simulation

N=20  #N
num<-400 # n
m<-20 #repeat
u1<-runif(num,-1,1)
u2<-runif(num,-1,1)
x1<-u1  # replace by x1<-runif(num,-1,1); x2<-runif(num,-1,1) for independent covariates
x2<-.5*u1+.5*u2
a1_r<- 1; a2_r<- -2; a3_r<-.5#psi
b1_r<- -.5 ; b2_r<-2; b3_r<-1 # eta
c1_r<-2; c2_r<-1 ; c3_r<-3 ;  # nu
p_r<-1/(1+exp(-(a1_r+x1*a2_r + x2*a3_r)))
h_r<-1/(1+exp(-(b1_r+x1*b2_r+ x2*b3_r )))

c00_r<- apply(rbind(p_r,h_r) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
c_r<-c00_r/(1+exp(-(c1_r+x1*c2_r+x2*c3_r )));

p0<-p1<-p2<-h0<-h1<-h2<-c0<-c1<-c2<-c00<-NULL
EST<-matrix(0, ncol=9, nrow=m)

DATA_simu0<-apply(rbind(p_r, h_r, c_r),2, function(x){
  rfrbinom( m, N, x[1],x[2],x[3] )})

for(i in 1:m){
  DATA_simu<-DATA_simu0[i,]
  likelihoodfunction_simu<-function(p0,p1,p2,h0,h1,h2,c0,c1,c2){
    p_l<-1/(1+exp(-(p0+x1*p1+x2*p2 )));
    h_l<-1/(1+exp(-(h0+ x1*h1+x2*h2  )));
    c00<- apply(rbind(p_l,h_l) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
    c_l<-c00/(1+exp(-(c0+x1*c1+x2*c2)));
    DAT<-rbind(DATA_simu,p_l,h_l,c_l )
    -sum(apply(DAT,  2 , function(x){log(dfrbinom(x[1],N,x[2],x[3],x[4] ))}))
  }
  
  est<-mle2(likelihoodfunction_simu, method="L-BFGS-B",
            start=list(  p0 = 0,p1 = 0,p2=0,
                         h0 = 0, h1 = 0, h2=0,c0 =0, c1 = 0,c2=0) ,
            lower= c( -5,-5,-5,-5,-5,-5,-5,-5,-5 ),
            upper=c(5,5,5,5,5,5,5,5,5)
  )
  
  EST[i,]<-est@coef
}

round(mean(EST[,1])-(a1_r), digits=2 ) ; round(sd(EST[,1]), digits=2 )
round(mean(EST[,2])-(a2_r), digits=2 ) ; round(sd(EST[,2]), digits=2 )
round(mean(EST[,3])-(a3_r), digits=2 ) ; round(sd(EST[,3]), digits=2 )
round(mean(EST[,4])-(b1_r), digits=2 ) ; round(sd(EST[,4]), digits=2 )
round(mean(EST[,5])-(b2_r), digits=2 ) ; round(sd(EST[,5]), digits=2 )
round(mean(EST[,6])-(b3_r), digits=2 ) ; round(sd(EST[,6]), digits=2 )
round(mean(EST[,7])-(c1_r), digits=2 ) ; round(sd(EST[,7]), digits=2 )
round(mean(EST[,8])-(c2_r), digits=2 ) ; round(sd(EST[,8]), digits=2 )
round(mean(EST[,8])-(c3_r), digits=2 ) ; round(sd(EST[,9]), digits=2 )

##### subsection 4.2
#if FB is the true model
N=10  #N
num<-200 # n
m<-20 #repeat
x1<-runif(num,0,1)
x2<-runif(num,0,1)
a1_r<- runif(m,-2,2); a2_r<-runif(m,0,1);  #psi
b1_r<- runif(m,0,2) ; b2_r<-runif(m,0,1);  # eta
c1_r<- runif(m,0,2); c2_r<-runif(m,0,1);   # nu
DATA_simu0<-matrix(0, nrow=m, ncol=num)
for(i in 1:m){
  p_r<-1/(1+exp(-(a1_r[i] +x1*a2_r[i] )))
  h_r<-1/(1+exp(-(b1_r[i] +x1*b2_r[i] )))
  c00_r<- apply(rbind(p_r,h_r) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
  c_r<-c00_r/(1+exp(-(c1_r[i] +x1*c2_r[i]  )));
  
DATA_simu0[i,]<-apply(rbind(p_r, h_r, c_r),2, function(x){
    rfrbinom( 1, N, x[1],x[2],x[3] )}) }

#if ZIP is the true model
num<-200 # n
m<-20 #repeat
x1<-runif(num,0,1)
x2<-runif(num,0,1)
b1_r<- runif(m,0.3, 1) ; b2_r<-runif(m,-.5,.5)  # eta
c1_r<- runif(m,-3,0); c2_r<-runif(m,-.5,.5)
DATA_simu0<-matrix(0, ncol=num , nrow=m)
for(i in 1:m) {th_r<-m_r<-pr<-c()
l_r<-exp((b1_r[i]+x1*b2_r[i] ))
p_0<-1/(1+exp(-(c1_r[i]+x1*c2_r[i] )))

DATA_simu0[i, ]<-apply(rbind(l_r, p_0),2, function(x){
  sample(seq(0,100,1),1,replace=TRUE,
         prob= c(x[2], rep(0,100))+(1-x[2])*dpois( seq(0,100,1), lambda=x[1],
                                                   log=FALSE) )})
}
#if ZINB is the true model
num<-200 # n
m<-20 #repeat
x1<-runif(num,-1,1)
x2<-runif(num,-1,1)
a1_r<- runif(m,-2,2); a2_r<-runif(m,-0.5,0.5) #psi
b1_r<- runif(m,0.3, 0.6) ; b2_r<-runif(m,-.5,.5)  # eta
c1_r<- runif(m,-3,0); c2_r<-runif(m,-.5,.5)
DATA_simu0<-matrix(0, ncol=num , nrow=m)
for(i in 1:m) {th_r<-m_r<-pr<-c()
th_r<-exp((a1_r[i] ))
m_r<-exp((b1_r[i]+x1*b2_r[i] ))
p_0<-1/(1+exp(-(c1_r[i]+x1*c2_r[i] )))
pr<-th_r/(m_r+th_r)

DATA_simu0[i, ]<-apply(rbind(th_r, pr, p_0),2, function(x){
  sample(seq(0,100,1),1,replace=TRUE,
         prob= c(x[3], rep(0,100))+(1-x[3])*dnbinom( seq(0,100,1), size=x[1],
                                                     prob=x[2], log=FALSE) )})
}
#if ZINB2 is the true model
num<-200 # n
m<-20 #repeat
x1<-runif(num,-1,1)
x2<-runif(num,-1,1)
a1_r<- runif(m,-2,2); a2_r<-runif(m,-0.5,0.5) #psi
b1_r<- runif(m,0.3, 0.6) ; b2_r<-runif(m,-.5,.5)  # eta
c1_r<- runif(m,-3,0); c2_r<-runif(m,-.5,.5)
DATA_simu0<-matrix(0, ncol=num , nrow=m)
for(i in 1:m) {th_r<-m_r<-pr<-c()
th_r<-exp((a1_r[i]+x1*a2_r[i] ))
m_r<-exp((b1_r[i]+x1*b2_r[i] ))
p_0<-1/(1+exp(-(c1_r[i]+x1*c2_r[i] )))
pr<-th_r/(m_r+th_r)

DATA_simu0[i, ]<-apply(rbind(th_r, pr, p_0),2, function(x){
  sample(seq(0,100,1),1,replace=TRUE,
         prob= c(x[3], rep(0,100))+(1-x[3])*dnbinom( seq(0,100,1), size=x[1],
                                                     prob=x[2], log=FALSE) )})
}
#if ZIB is the true model
N=10  #N
num<-200 # n
m<-20 #repeat
x1<-runif(num,0,1)
x2<-runif(num,0,1)
a1_r<- runif(m,-2,2); a2_r<-runif(m,0,1);  #psi
c1_r<- runif(m,-3,0); c2_r<-runif(m,-.5,.5)

DATA_simu0<-matrix(0, ncol=num , nrow=m)
for(i in 1:m) {p_0<-p_r<-c()
p_r<-1/(1+exp(-(a1_r[i] +x1*a2_r[i] )))

p_0<-1/(1+exp(-(c1_r[i]+x1*c2_r[i] )))

DATA_simu0[i, ]<-apply(rbind( p_r, p_0),2, function(x){
  sample(seq(0,100,1),1,replace=TRUE,
         prob= c(x[2], rep(0,100))+(1-x[2])*dbinom( seq(0,100,1), N,
                                                    x[1]) )})
}
## obtain mean AIC difference 

library(pscl)
AIC<-AIC_p<-AIC_nb<-AIC_nb2<-AIC_b<-c()
for(i in 1:20){
  y<-DATA_simu0[i,]
  
  M3 <- zeroinfl(y~ x1 | x1,
               dist = 'negbin')
  AIC_nb[i]<- AIC(M3)
  
  M4 <- zeroinfl(y~ x1  | x1,
                dist = 'pois')
  AIC_p[i] <-AIC(M4)
  AIC_nb2[i]<-ZINB2( y=y,x=data.frame(x1=x1) )$AIC
  AIC[i]<-fbglm( y=y,x=data.frame(x1=x1) )$AIC
  AIC_b[i]<-ZIB(y=y, x=data.frame(x1=x1))$AIC
}

mean(AIC_b-AIC)
length(which(AIC_b-AIC> 4))/20  
length(which(AIC_b-AIC< -4))/20  

mean(AIC_p-AIC)
length(which(AIC_p-AIC> 4))/20  
length(which(AIC_p-AIC< -4))/20

mean(AIC_nb-AIC)
length(which(AIC_nb-AIC> 4))/20  
length(which(AIC_nb-AIC< -4))/20

mean(AIC_nb2-AIC)
length(which(AIC_nb2-AIC> 4))/20  
length(which(AIC_nb2-AIC< -4))/20




###### Section 5. Data analysis ########

##### subsection 5.1. ####

install.packages("agridat")
library(agridat)
data(ridout.appleshoots)
dat <- ridout.appleshoots

pho<-rep(0,270)
pho[dat$photo==16]<-1
### fit the models #####
## FB ###
x0<-data.frame( pho=pho, bap=dat$bap)
fbglm(y=dat$roots, x=x0   )
## ZINB-2 ##
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
# ZIB
ZIB(y=dat$roots,x=x0  )

#### RQRs ####
#ZIB
ZIB2<-ZIB(y=dat$roots,x=data.frame(pho=pho)  )
coef_zib2<-c(ZIB2$coef  )
p_b<-1/(1+exp(-(coef_zib2[1] +pho*coef_zib2[2] )))
pi_b<-1/(1+exp(-(coef_zib2[3] +pho*coef_zib2[4] )))
r_b_2<-c()
for(i in 1:270){if(y[i]==0){r_b_2[i]<-runif(1,0,1)*(pi_b[i]+(1-pi_b[i])*dbinom( 0, N, prob=p_b[i] ))}
  else  { r_b_2[i]<-pi_b[i] +(1-pi_b[i])*pbinom( y[i]-1, N, prob=p_b[i])+runif(1,0,1)*(1-pi_b[i])*dbinom( y[i], N, prob=p_b[i] )   } } 

# FB
FB2<-fbglm(y=dat$roots,x=data.frame(pho=pho)  )
coef_fb2<-c(FB2$coef  )
p_r<-1/(1+exp(-(coef_fb2[1] +pho*coef_fb2[2] )))
h_r<-1/(1+exp(-(coef_fb2[3] +pho*coef_fb2[4] )))

c00_r<- apply(rbind(p_r,h_r) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
c_r<-c00_r/(1+exp(-( coef_fb2[5] +pho*coef_fb2[6]  )));
y<-dat$roots
N<-max(y)

r_fb2<-c()
for(i in 1:270){if(y[i]==0){r_fb2[i]<- runif(1,0,1)* dfrbinom( y[i], N, p_r[i],h_r[i],c_r[i] ) }   
  else{ r_fb2[i]<- sum(dfrbinom(seq(0,y[i]-1,1), N, p_r[i],h_r[i],c_r[i] )) +runif(1,0,1)* dfrbinom( y[i], N, p_r[i],h_r[i],c_r[i])}  } 

# ZIP
M6 <- zeroinfl(dat$roots~ pho  |
                 pho,
               dist = 'pois')
coef_p_2<-c(M6$coefficients$count, M6$coefficients$zero   )

l_r<-exp((coef_p_2[1]+pho*coef_p_2[2] ))
p_0_p<-1/(1+exp(-(coef_p_2[3]+pho*coef_p_2[4]   )))

r_p_2<-c()
for(i in 1:270){if(y[i]==0){r_p_2[i]<-runif(1,0,1)*(p_0_p[i]+(1-p_0_p[i])*dpois( 0, lambda=l_r[i], log=FALSE))}
  else  { r_p_2[i]<-p_0_p[i] +(1-p_0_p[i])*ppois( y[i]-1, lambda=l_r[i], log=FALSE )+runif(1,0,1)*(1-p_0_p[i])*dpois( y[i], lambda=l_r[i], log=FALSE )   } } 

#ZINB
M4 <- zeroinfl(dat$roots~ pho  |
                 pho  ,
               dist = 'negbin')

coef_nb_2<-c( M4$coefficients$count, M4$coefficients$zero ,2.51598)

th_r<-m_r<-pr<-c()
th_r<-exp(coef_nb_2[5])
m_r<-exp((coef_nb_2[1]+pho*coef_nb_2[2]))
p_0<-1/(1+exp(-(coef_nb_2[3]+pho*coef_nb_2[4] )))
pr<-th_r/(m_r+th_r)

r_nb_2<-c()
for(i in 1:270){if(y[i]==0){r_nb_2[i]<-runif(1,0,1)*(p_0[i]+(1-p_0[i])*dnbinom( 0, size=th_r, prob=pr[i] , log=FALSE))}
  else  { r_nb_2[i]<-p_0[i] +(1-p_0[i])*pnbinom( y[i]-1, size=th_r, prob=pr[i], log=FALSE )+runif(1,0,1)*(1-p_0[i])*dnbinom( y[i], size=th_r, prob=pr[i], log=FALSE )   } } 

#ZINB-2
Z2<-ZINB2(y=dat$roots, x=data.frame(pho=pho)   )
coef_nb2_2<- Z2$coef  

th_r<-m_r<-pr<-c()
th_r<-exp(coef_nb2_2[1]+pho*coef_nb2_2[2])
m_r2<-exp((coef_nb2_2[3]+pho*coef_nb2_2[4] ))
p_0<-1/(1+exp(-(coef_nb2_2[5]+pho*coef_nb2_2[6] )))
pr<-th_r/(m_r2+th_r)

r_nb2_2<-c()
for(i in 1:270){if(y[i]==0){r_nb2_2[i]<-runif(1,0,1)*(p_0[i]+(1-p_0[i])*dnbinom( 0, size=th_r[i], prob=pr[i] , log=FALSE))}
  else  { r_nb2_2[i]<-p_0[i] +(1-p_0[i])*pnbinom( y[i]-1, size=th_r[i], prob=pr[i], log=FALSE )+runif(1,0,1)*(1-p_0[i])*dnbinom( y[i], size=th_r[i], prob=pr[i], log=FALSE )   } } 

## QQplot
# FB
qqnorm(qnorm(r_fb2),  main="FB" ,cex.lab=2, cex.main=2.5, cex.axis=1.7 )
qqline(qnorm(r_fb2), col = "blue", lwd=2)
# ZIP
qqnorm(qnorm(r_p_2),  main="ZIP" ,cex.lab=2, cex.main=2.5, cex.axis=1.7)
qqline(qnorm(r_p_2), col = "blue", lwd=2)
# ZINB
 qqnorm(qnorm(r_nb_2),  main="ZINB" ,cex.lab=2, cex.main=2.5, cex.axis=1.7)
     qqline(qnorm(r_nb_2), col = "blue", lwd=2)
# ZINB-2
qqnorm(qnorm(r_nb2_2),  main="ZINB-2" ,cex.lab=2, cex.main=2.5, cex.axis=1.7)
     qqline(qnorm(r_nb2_2), col = "blue", lwd=2)
# ZIB
qqnorm(qnorm(r_b_2),  main="ZIB" ,cex.lab=2, cex.main=2.5, cex.axis=1.7)
     qqline(qnorm(r_b_2), col = "blue", lwd=2)

### SW test
#FB
coef_fb2<-c(FB2$coef  )
     p_r<-1/(1+exp(-(coef_fb2[1] +pho*coef_fb2[2] )))
     h_r<-1/(1+exp(-(coef_fb2[3] +pho*coef_fb2[4] )))
     
     c00_r<- apply(rbind(p_r,h_r) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
     c_r<-c00_r/(1+exp(-( coef_fb2[5] +pho*coef_fb2[6]  )));
     y<-dat$roots
     N<-max(y)
     p_fb2<-c()
    for(j in 1:100){ 
     r_fb2<-c()
     for(i in 1:270){if(y[i]==0){r_fb2[i]<- runif(1,0,1)* dfrbinom( y[i], N, p_r[i],h_r[i],c_r[i] ) }   
       else{ r_fb2[i]<- sum(dfrbinom(seq(0,y[i]-1,1), N, p_r[i],h_r[i],c_r[i] )) +runif(1,0,1)* dfrbinom( y[i], N, p_r[i],h_r[i],c_r[i])}  } 
     a<-shapiro.test(qnorm(r_fb2))
     p_fb2[j]<-a$p.value}
    mean(p_fb2)
     #  0.7527965
   
     #ZIP
     M6 <- zeroinfl(dat$roots~ pho  |
                      pho,
                    dist = 'pois')
   
     coef_p_2<-c(M6$coefficients$count, M6$coefficients$zero   )
     
     l_r<-exp((coef_p_2[1]+pho*coef_p_2[2] ))
     p_0_p<-1/(1+exp(-(coef_p_2[3]+pho*coef_p_2[4]   )))
     pval_p<-c()
     for(j in 1:100){
     r_p_2<-c()
     for(i in 1:270){if(y[i]==0){r_p_2[i]<-runif(1,0,1)*(p_0_p[i]+(1-p_0_p[i])*dpois( 0, lambda=l_r[i], log=FALSE))}
       else  { r_p_2[i]<-p_0_p[i] +(1-p_0_p[i])*ppois( y[i]-1, lambda=l_r[i], log=FALSE )+runif(1,0,1)*(1-p_0_p[i])*dpois( y[i], lambda=l_r[i], log=FALSE )   } } 
     
     a<-shapiro.test(qnorm(r_p_2))
     pval_p[j]<-a$p.value} 
     
     mean(pval_p)
     #  0.5190557
     #ZINB-2
     Z2<-ZINB2(y=dat$roots, x=data.frame(pho=pho)   )
     coef_nb2_2<- Z2$coef  
     
     th_r<-m_r<-pr<-c()
     th_r<-exp(coef_nb2_2[1]+pho*coef_nb2_2[2])
     m_r2<-exp((coef_nb2_2[3]+pho*coef_nb2_2[4] ))
     p_0<-1/(1+exp(-(coef_nb2_2[5]+pho*coef_nb2_2[6] )))
     pr<-th_r/(m_r2+th_r)
     pval_nb2<-c()
     for(j in 1:100){
     r_nb2_2<-c()
     for(i in 1:270){if(y[i]==0){r_nb2_2[i]<-runif(1,0,1)*(p_0[i]+(1-p_0[i])*dnbinom( 0, size=th_r[i], prob=pr[i] , log=FALSE))}
       else  { r_nb2_2[i]<-p_0[i] +(1-p_0[i])*pnbinom( y[i]-1, size=th_r[i], prob=pr[i], log=FALSE )+runif(1,0,1)*(1-p_0[i])*dnbinom( y[i], size=th_r[i], prob=pr[i], log=FALSE )   } } 
    a<-shapiro.test(qnorm(r_nb2_2))
     pval_nb2[j]<-a$p.value}
     
     mean(pval_nb2)
     # 0.4048557
     
     #ZINB
     M4 <- zeroinfl(dat$roots~ pho  |
                      pho  ,
                    dist = 'negbin')
 
     coef_nb_2<-c( M4$coefficients$count, M4$coefficients$zero ,2.51598)
 
     th_r<-m_r<-pr<-c()
     th_r<-exp(coef_nb_2[5])
     m_r<-exp((coef_nb_2[1]+pho*coef_nb_2[2]))
     p_0<-1/(1+exp(-(coef_nb_2[3]+pho*coef_nb_2[4] )))
     pr<-th_r/(m_r+th_r)
     pval_nb<-c()
    for(j in 1:100){ 
     r_nb_2<-c()
     for(i in 1:270){if(y[i]==0){r_nb_2[i]<-runif(1,0,1)*(p_0[i]+(1-p_0[i])*dnbinom( 0, size=th_r, prob=pr[i] , log=FALSE))}
       else  { r_nb_2[i]<-p_0[i] +(1-p_0[i])*pnbinom( y[i]-1, size=th_r, prob=pr[i], log=FALSE )+runif(1,0,1)*(1-p_0[i])*dnbinom( y[i], size=th_r, prob=pr[i], log=FALSE )   } } 
     a<-shapiro.test(qnorm(r_nb_2))
     pval_nb[j]<-a$p.value}
     
     mean(pval_nb)
     # 0.5157274
     
     #ZIB
   ZIB2<-ZIB(y=dat$roots,x=data.frame(pho=pho)  )
     coef_zib2<-c(ZIB2$coef  )
   p_b<-1/(1+exp(-(coef_zib2[1] +pho*coef_zib2[2] )))
     pi_b<-1/(1+exp(-(coef_zib2[3] +pho*coef_zib2[4] )))
     pval_b<-c()
   for(j in 1:100){  
     r_b_2<-c()
     for(i in 1:270){if(y[i]==0){r_b_2[i]<-runif(1,0,1)*(pi_b[i]+(1-pi_b[i])*dbinom( 0, N, prob=p_b[i] ))}
       else  { r_b_2[i]<-pi_b[i] +(1-pi_b[i])*pbinom( y[i]-1, N, prob=p_b[i])+runif(1,0,1)*(1-pi_b[i])*dbinom( y[i], N, prob=p_b[i] )   } } 
     a<-shapiro.test(qnorm(r_b_2))
     pval_b[j]<-a$p.value}
     mean(pval_b)
     #0.001589684

#####  subsection 5.2 #########

my_data <- read.delim("Deaths_5x10.txt", header = TRUE,sep = "")

Fe<-round(my_data$Female)
Ma<-round(my_data$Male )
Age<-rep( seq(0,23,1 )  , 10)
Year<-c()
for(i in 1:10){  Year<-c(Year, rep(i, 24))  }

D_data<-data.frame(Fe=Fe, MA=Ma, Age=Age, Year=Year )
for(i in 1:10){ a<-1+24*i; a4<-a-4; a3<-a-3; a2<-a-2;a1<-a-1 
Ma[a4]<-sum(Ma[ c(a4,a3, a2, a1) ]) 
Fe[a4]<-sum(Fe[ c(a4,a3, a2, a1) ]) 
}

D_data<-data.frame(Fe=Fe, MA=Ma, Age=Age, Year=Year )
A1<-which(Age==21)
A2<-which(Age==22)
A3<-which(Age==23)
D_data2<-D_data[-c(A1, A2, A3), ]

D_data3<-data.frame(Death=c(D_data2$MA, D_data2$Fe),
              Year=c(D_data2$Year, D_data2$Year) ,
              Age=c(D_data2$Age, D_data2$Age),
              Gender=c( rep("Male", 210), rep("Female",210)
                    ))
G<-rep(1,420)
G[which(D_data3$Gender=="Male")]<-0
XXX<-data.frame(Gender=G,Year=D_data3$Year )

data5<-data.frame(y=D_data3$Age, 
                  G=XXX$Gender, Year=D_data3$Year,
                  w=round(D_data3$Death/10000))
# ZINB
D_nb <- zeroinfl(y~ G+Year |
                   G+Year  ,data=data5, weights = w,
                 dist = 'negbin')
# ZIP
D_p <- zeroinfl(y~ G+Year  |
                  G+Year  ,data=data5, weights = w,
                dist  = 'pois')
# FB
FBB<-fbglm0(y=D_data3$Age, 
            x=XXX,
            w=data5$w)
# ZIB

ZIBB<-ZIB0(y=D_data3$Age, 
            x=XXX,
            w=data5$w)
## ZINB-2

ZINBB2<-ZINB20(y=D_data3$Age, 
           x=XXX,
           w=data5$w)

###### Graph ####

b<-which(XXX$Year==1 & XXX$Gender==0) #106 male&year 6

k<-min(b)

### ZIP
for(nnn in 1:1){
prob_pois<-c()

l_r<-exp((p_d[1]+G*p_d[2]+XXX$Year*p_d[3] ))
p_0_p<-1/(1+exp(-(p_d[4]+G*p_d[5]+XXX$Year*p_d[6]   )))
N<-20

for(i in 1:N){ 
  prob_pois[i]<-(1-p_0_p[k])*dpois( i, lambda=l_r[k], log=FALSE) 
}
prob_pois<-c(p_0_p[k] +(1-p_0_p[k])*dpois( 0, lambda=l_r[k] ,log=FALSE ),    prob_pois)


## ZINB

prob_nb<-c()
th_r<-m_r<-pr<-c()
th_r<-exp(nb_d[7])
m_r<-exp((nb_d[1]+G*nb_d[2]+ XXX$Year*nb_d[3]))
p_0<-1/(1+exp(-(nb_d[4]+G*nb_d[5]+ XXX$Year*nb_d[6] )))
pr<-th_r/(m_r+th_r)

for(i in 1:N){ 
  prob_nb[i]<-(1-p_0[k])*dnbinom( i, size=th_r[1], mu=m_r[k] , log=FALSE) 
}
prob_nb<-c(p_0[k] +(1-p_0[k])*pnbinom( 0, size=th_r[1], mu=m_r[k], log=FALSE ),    prob_nb)


### ZINB-2 ##
th_r<-m_r<-pr<-c()
th_r<-exp(nb2_d[1]+G*nb2_d[2]+XXX$Year*nb2_d[3])
m_r2<-exp((nb2_d[4]+G*nb2_d[5]+XXX$Year*nb2_d[6] ))
p_0<-1/(1+exp(-(nb2_d[7]+G*nb2_d[8]+XXX$Year*nb2_d[9] )))
pr<-th_r/(m_r2+th_r)
prob_nb2<-c()
for(i in 1:N){ 
  prob_nb2[i]<-(1-p_0[k])*dnbinom( i, size=th_r[k], prob=pr[k] , log=FALSE) 
}
prob_nb2<-c(p_0[k] +(1-p_0[k])*pnbinom( 0, size=th_r[k], prob=pr[k], log=FALSE ),    prob_nb2)

### ZIB  

p_b<-1/(1+exp(-(zib_d[1] +G*zib_d[2]+XXX$Year*zib_d[3] )))
pi_b<-1/(1+exp(-(zib_d[4] +G*zib_d[5]+XXX$Year*zib_d[6] )))
prob_b<-c()
for(i in 1:N){ 
  prob_b[i]<-(1-pi_b[k])*dbinom( i, N, prob=p_b[k] ) 
}
prob_b<-c(pi_b[k] +(1-pi_b[k])*pbinom( 0, N, prob=p_b[k]),    prob_b)


# FB ##

p_r<-1/(1+exp(-(fb_d[1] +G*fb_d[2]+XXX$Year*fb_d[3] )))
h_r<-1/(1+exp(-(fb_d[4] +G*fb_d[5]+XXX$Year*fb_d[6] )))
c00_r<- apply(rbind(p_r,h_r) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
c_r<-c00_r/(1+exp(-( fb_d[7] +G*fb_d[8]+XXX$Year*fb_d[9]  )));

prob_fb<-dfrbinom( seq(0,20,1),N, p_r[k], h_r[k], c_r[k] )

PR5<-data.frame(  x=seq(0,20,1),p00=prob_b,p1=prob_pois, p2=prob_nb, p3=prob_nb2, p4=prob_fb)


# start  ###
k20<-k+20
df.new3<-data.frame(x=seq(0,20,1), y=D_data3$Death[k:k20]/sum( D_data3$Death[k:k20] ) )

g5<-ggplot(df.new3, aes(x=x, y=y )) + 
  geom_col(   )

g5<-g5+geom_line( data=PR5, aes(x=x-.05, y=p00,color="ZIB"), size=1.5  )+
  geom_point(data=PR5, aes(x=x-.05, y=p00,color="ZIB"), size=1.5 )+
  geom_line( data=PR5, aes(x=x-.15, y=p3,color="ZINB-2"), size=1.5)+
  geom_point(data=PR5, aes(x=x-.15, y=p3,color="ZINB-2"), size=1.5)+
  geom_line( data=PR5, aes(x=x, y=p2,color="ZINB"), size=1.5)+
  geom_point( data=PR5, aes(x=x, y=p2,color="ZINB"), size=1.5)+
  geom_line( data=PR5, aes(x=x+.15, y=p1,color="ZIP"), size=1.5)+
  geom_point(data=PR5, aes(x=x+.15, y=p1,color="ZIP"), size=1.5 )+
  geom_line( data=PR5, aes(x=x+.05, y=p4,color="FB"), size=1.5  )+
  geom_point(data=PR5, aes(x=x+.05, y=p4,color="FB" ), size=1.5 )+
  scale_x_continuous( breaks=seq(0,20,1),labels=c(0, "1-4", "5-9","10-14",
                                                  "15-19", "20-24","25-29","30-34","35-39","40-44","45-49",
                                                  "50-54", "55-59", "60-64","65-69","70-74","75-79",
                                                  "80-84","85-89","90-94","95+"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20))+  
  scale_color_manual(values = colors)+
  labs(x="Age at Death", y="Probability", color="", alpha="",title="Male (1933-1939)")+
  theme(plot.title = element_text(size=27),
        axis.text.x = element_text(angle = 45, hjust = .5, vjust = 0.5)
        ,legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.text = element_text(size=20) ,legend.position=c(.25,.7) )
}
g5
