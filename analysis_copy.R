# Code to replicate the analysis done in 
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-018-2460-9
# The code is not very readable, but I am providing it in the interest of open research
# I promise I will come back and fix this when I have some time after my thesis is finished

# Load required packages

library(ggplot2)
library(reshape2)
library(dplyr)
library(magrittr)
library(dlnm)

# Read in data, available with the paper

# temp <- read.csv("~/Dropbox/PhD/Data/temp_rfe.csv") %>% 
#  mutate(pos=ceiling(ANCvisits*ANCprevalence),neg=ANCvisits-ceiling(ANCvisits*ANCprevalence))

temp <- read.csv("~/Downloads/12936_2018_2460_MOESM2_ESM.csv") %>% 
 mutate(pos=ceiling(ANCvisits*ANCprevalence),neg=ANCvisits-ceiling(ANCvisits*ANCprevalence))

 
# Unfortunately, the data I made available with the paper does not have dates, so we must add them manually

temp$date <- c(seq.Date(from=as.Date("2010-04-01"),to=as.Date("2015-12-01"),by = "month"), # Baraka
seq.Date(from=as.Date("2010-01-01"),to=as.Date("2014-12-01"),by = "month"), # Mweso
seq.Date(from=as.Date("2013-02-01"),to=as.Date("2014-12-01"),by = "month"), # Walikale
seq.Date(from=as.Date("2009-01-01"),to=as.Date("2014-12-01"),by = "month"), # Shamwana
seq.Date(from=as.Date("2012-01-01"),to=as.Date("2013-12-01"),by = "month"))

temp$year <- substr(temp$date,1,4)
 


# Hill function fitted later by DLNM model 

fhill <- function(x,theta,n){
  basis <- (x^n) / (theta^n + x^n)
  attributes(basis)$n <- n
  attributes(basis)$theta <- theta
  return(basis)
}

# Model LE
cb_le <- onebasis(temp$un5inc,fun="lin")
LE <- mgcv::gam(data=temp,cbind(pos,neg)~cb_le+factor(location)-1,
                  family = binomial,method="REML",group=location)
# Model LELL
cb_lell <- crossbasis(temp$un5inc,lag=3,arglag=list(fun="lin"),argvar=list(fun="lin"),group=temp$location)
LELL <- mgcv::gam(data=temp,cbind(pos,neg)~cb_lell+factor(location)-1,
                  family = binomial,method="REML",group=location)

# Model LENL
cb_lenl <-  crossbasis(temp$un5inc,lag=3,argvar=list(fun="lin"),arglag=list(fun="ps",df=8),group=temp$location)
cbpen_lenl <- cbPen(cb_lenl)
LENL <- mgcv::gam(cbind(pos,neg)~cb_lenl+factor(location)-1,family=binomial(link="logit"),temp,
                  paraPen=list(cb_lenl=cbpen_lenl), method="REML",group=location)

# Model NE
cb_ne <- onebasis(x=temp$un5inc,fun="fhill",theta=30.591975,n=0.587693)
NE <- mgcv::gam(cbind(pos,neg)~cb_ne+factor(location)-1,family=binomial(link="logit"),temp,
                method="REML",group=location)

# Model NELL
cb_nell <-  crossbasis(temp$un5inc,lag=3,argvar=list(fun="fhill",theta=30.591975,n=0.587693),arglag=list(fun="lin"),group=temp$location)
NELL <- mgcv::gam(cbind(pos,neg)~cb_nell+factor(location)-1,family=binomial(link="logit"),temp,
                   method="REML",group=location)

# Model NENL
cb_nenl <-  crossbasis(temp$un5inc,lag=3,arglag=list(fun="ps",df=8),argvar=list(fun="fhill",theta=30.591975,n=0.587693),group=temp$location)
cbpen_nenl <- cbPen(cb_nenl)
NENL <- mgcv::gam(cbind(pos,neg)~cb_nenl+factor(location)-1,family=binomial(link="logit"),temp,
                  paraPen=list(cb_nenl=cbpen_nenl), method="REML",group=location)


# AIC comparisons in Table 2

AIC(LE)
AIC(LELL)
AIC(LENL)
AIC(NE)
AIC(NELL)
AIC(NENL)

# Function to perform rolling origin cross-validation and work out the root mean squared error (rmse) for each model

rfcv <- function(loc,yr){
  train <- temp %>% dplyr::filter(location!=loc | (location==loc & year<yr))
  test <- temp %>% dplyr::filter(location!=loc | (location==loc & year<(yr+1)))

  ## LE
  cb_le <- onebasis(x=train$un5inc,fun="lin")
  LE_train <- mgcv::gam(data=train,cbind(pos,neg)~cb_le+factor(location),
                        family = binomial,method="REML",group=train$location)
  cb_le_t <- onebasis(x=test$un5inc,fun="lin")
  LE_pred <- predict(object=LE_train,newdata=list(cb_le=cb_le_t,location=test$location),type="response",se.fit=TRUE)

  test$y_LE <- LE_pred$fit
  test$y_LE_upr <- LE_pred$fit + (1.96*LE_pred$se.fit)
  test$y_LE_lwr <- LE_pred$fit - (1.96*LE_pred$se.fit)


  ## LELL
  cb_lell <- crossbasis(train$un5inc,lag=3,arglag=list(fun="lin"),argvar=list(fun="lin"),group=train$location)
  LELL_train <- mgcv::gam(data=train,cbind(pos,neg)~cb_lell+factor(location),
                          family = binomial,method="REML",group=location)
  cb_lell_t <- crossbasis(test$un5inc,lag=3,arglag=list(fun="lin"),argvar=list(fun="lin"),group=test$location)
  LELL_pred <- predict(object=LELL_train,newdata=list(cb_lell=cb_lell_t,location=test$location),type="response",se.fit=TRUE)

  test$y_LELL <- LELL_pred$fit
  test$y_LELL_upr <- LELL_pred$fit + (1.96*LELL_pred$se.fit)
  test$y_LELL_lwr <- LELL_pred$fit - (1.96*LELL_pred$se.fit)

  ## LENL
  cb_lenl <-  crossbasis(train$un5inc,lag=3,argvar=list(fun="lin"),arglag=list(fun="ps",df=8),group=train$location)
  cbpen_lenl <- cbPen(cb_lenl)
  LENL_train <- mgcv::gam(cbind(pos,neg)~cb_lenl+factor(location),family=binomial(link="logit"),train,
                          paraPen=list(cb_lenl=cbpen_lenl), method="REML",group=location)

  cb_lenl_t <- crossbasis(test$un5inc,lag=3,argvar=list(fun="lin"),arglag=list(fun="ps",df=8),group=test$location)
  LENL_pred <- predict(object=LENL_train,newdata=list(cb_lenl=cb_lenl_t,location=test$location),type="response",se.fit=TRUE)

  test$y_LENL <- LENL_pred$fit
  test$y_LENL_upr <- LENL_pred$fit + (1.96*LENL_pred$se.fit)
  test$y_LENL_lwr <- LENL_pred$fit - (1.96*LENL_pred$se.fit)

  ## NE
  cb_ne <- onebasis(x=train$un5inc,fun="fhill",theta=30.591975,n=0.587693)
  NE_train <- mgcv::gam(cbind(pos,neg)~cb_ne+factor(location),family=binomial(link="logit"),train,
                        method="REML",group=train$location)
  cb_ne_t <- crossbasis(x=test$un5inc,fun="fhill",theta=30.591975,n=0.587693)
  NE_pred <- predict(object=NE_train,newdata=list(cb_ne=cb_ne_t,location=test$location),type="response",se.fit=TRUE)


  test$y_NE <- NE_pred$fit
  test$y_NE_upr <- NE_pred$fit + (1.96*NE_pred$se.fit)
  test$y_NE_lwr <- NE_pred$fit - (1.96*NE_pred$se.fit)

  ## NELL
  cb_nell <- crossbasis(train$un5inc,lag=3,argvar=list(fun="fhill",theta=30.591975,n=0.587693),arglag=list(fun="lin"),group=train$location)
  NELL_train <- mgcv::gam(cbind(pos,neg)~cb_nell+factor(location),family=binomial(link="logit"),train,
                          method="REML",group=location)
  cb_nell_t <- crossbasis(test$un5inc,lag=3,argvar=list(fun="fhill",theta=30.591975,n=0.587693),arglag=list(fun="lin"),group=test$location)
  NELL_pred <- predict(object=NELL_train,newdata=list(cb_nell=cb_nell_t,location=test$location),se.fit=TRUE,type="response")

  test$y_NELL <- NELL_pred$fit
  test$y_NELL_upr <- NELL_pred$fit + (1.96*NELL_pred$se.fit)
  test$y_NELL_lwr <- NELL_pred$fit - (1.96*NELL_pred$se.fit)

  ## NENL
  cb_nenl <- crossbasis(train$un5inc,lag=3,arglag=list(fun="ps",df=8),argvar=list(fun="fhill",theta=30.591975,n=0.587693),group=train$location)
  cbpen_nenl <- cbPen(cb_nenl)
  NENL_train <- mgcv::gam(cbind(pos,neg)~cb_nenl+factor(location),family=binomial(link="logit"),train,
                          paraPen=list(cb_nenl=cbpen_nenl), method="REML",group=location)

  cb_nenl_t <- crossbasis(test$un5inc,lag=3,arglag=list(fun="ps",df=8),argvar=list(fun="fhill",theta=30.591975,n=0.587693),group=test$location)
  NENL_pred <- predict(object=NENL_train,newdata=list(cb_nenl=cb_nenl_t,location=test$location),type="response",se.fit=TRUE)

  test$y_NENL <- NENL_pred$fit
  test$y_NENL_upr <- NENL_pred$fit + (1.96*NENL_pred$se.fit)
  test$y_NENL_lwr <- NENL_pred$fit - (1.96*NENL_pred$se.fit)

  out <- test %>% dplyr::filter(location==loc & year==yr)
  return(out)
}

# Table of results for rolling origin cross validation

ndf <- rbind(rfcv("Baraka",2011),rfcv("Baraka",2012),rfcv("Baraka",2013),rfcv("Baraka",2014),rfcv("Baraka",2015),
             rfcv("Mweso",2011),rfcv("Mweso",2012),rfcv("Mweso",2013),rfcv("Mweso",2014),
             rfcv("Walikale",2014),
             rfcv("Lulimba",2013),
             rfcv("Shamwana",2010),rfcv("Shamwana",2011),rfcv("Shamwana",2012),rfcv("Shamwana",2013),rfcv("Shamwana",2014))

# RMSE value for each model, found in Table 2

ndf %>% summarise(rmse=sqrt(mean((ANCprevalence-y_LE)^2)))
ndf %>% summarise(rmse=sqrt(mean((ANCprevalence-y_LELL)^2)))
ndf %>% summarise(rmse=sqrt(mean((ANCprevalence-y_LENL)^2)))
ndf %>% summarise(rmse=sqrt(mean((ANCprevalence-y_NE)^2)))
ndf %>% summarise(rmse=sqrt(mean((ANCprevalence-y_NELL)^2)))
ndf %>% summarise(rmse=sqrt(mean((ANCprevalence-y_NENL)^2)))


# Code for Fig 3.

target <- c("Baraka","Mweso","Walikale","Shamwana","Lulimba")
temp$location <- gdata::reorder.factor(temp$location,new.order=target)


target <- c("Mweso","Walikale","Baraka","Lulimba","Shamwana")
prdf2 <- data.frame(date=temp$date,location=temp$location,type="ANC Prevalence",
                    CIlow=temp$ANCprevalence -(1.96*sqrt(1/temp$ANCvisits)*temp$ANCprevalence*(1-temp$ANCprevalence)),
                    CIhi=temp$ANCprevalence +(1.96*sqrt(1/temp$ANCvisits)*temp$ANCprevalence*(1-temp$ANCprevalence)),
                    obs=temp$ANCprevalence,val=NA)
incdf <- data.frame(date=temp$date,location=temp$location,type="Cases per child\n per year",CIlow=NA,CIhi=NA,obs=temp$un5inc,val=NA)




obsdf <- rbind(incdf,prdf2)
obsdf$location <- gdata::reorder.factor(obsdf$location,new.order=target)

# plot_pbs <- ggplot(data=obsdf) + facet_grid(type~location,scales="free",switch="y") + theme_bw() +
#   geom_vline(xintercept=as.numeric(seq.Date(from=as.Date("2009/07/01"),to=as.Date("2016/12/01"),by="year")),lty=2,alpha=0.25) +
#   geom_vline(xintercept=as.numeric(seq.Date(from=as.Date("2009/01/01"),to=as.Date("2016/12/01"),by="year")),alpha=0.35) +
#   geom_ribbon(aes(x=date,ymin=CIlow,ymax=CIhi,fill=location),alpha=0.5) +
#   geom_line(aes(x=date,y=obs,col=location),size=1.05) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.text.y=element_text(size=14,face="bold"),
#         strip.text.x=element_text(size=14,face="bold"),
#         strip.placement = "outside",
#         legend.position = "none",
#         axis.title.x = element_blank(),axis.title.y = element_blank(),
#         axis.text=element_text(size=16)) +
#   scale_x_date(date_breaks = "2 years", date_labels =  "%Y",expand=c(0.1,0))
# plot_pbs

max_ANC <- 0.6
min_ANC <- 0
max_inc <- 9.397008
# max_inc <- 10
min_inc <- 0

con_int <- obsdf %>% filter(type=="ANC Prevalence") %>%
 dplyr::select(CIhi,CIlow,location,date)

obsdf %>% dplyr::mutate(obs=ifelse(type!="ANC Prevalence",
                                   (max_ANC*obs)/(max_inc),obs)) %>%
 ggplot() + facet_wrap(~location,scales="free_x") + theme_bw() +
 geom_vline(xintercept=seq.Date(from=as.Date("2009/07/01"),to=as.Date("2016/12/01"),by="year"),lty=2,alpha=0.25) +
 geom_vline(xintercept=as.numeric(seq.Date(from=as.Date("2009/01/01"),to=as.Date("2016/12/01"),by="year")),alpha=0.35) +
 geom_ribbon(data=con_int,aes(x=date,ymin=CIlow,ymax=CIhi),fill= "red",alpha=0.5) +
 geom_line(aes(x=date,y=obs,lty=type)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       strip.text.y=element_text(size=14,face="bold"),
       strip.text.x=element_text(size=14,face="bold"),
       strip.placement = "outside",
       legend.position = "none",
       axis.title.x = element_blank(),axis.title.y = element_blank(),
       axis.text=element_text(size=16)) +
 scale_x_date(date_breaks = "2 years", date_labels =  "%Y",expand=c(0.1,0)) +
 scale_y_continuous(sec.axis=sec_axis(~.*max_inc/max_ANC),breaks=seq(0,0.6,0.2),labels=paste(seq(0,60,20),"%",sep=""))

# Code for Fig 4.

yearav <- temp %>% group_by(location,year) %>% summarise(x=mean(ANCprevalence),y=mean(un5inc))#,tem=mean(temp))
yearav$location <- gdata::reorder.factor(yearav$location,new.order=target)

temp2 <- temp
temp2$location <-  gdata::reorder.factor(temp2$location,new.order=target)


p_out2 <- function(loc,to1){
 
 x1 <- predict(NENL,newdata = list(cb_nenl=crossbasis(x=seq(0,to1,0.1),lag=3,arglag=list(fun="ps",df=8),argvar=list(fun="fhill",theta=30.591975,n=0.587693)),location=rep(loc,length(seq(0,to1,0.1)))),type="response",se=TRUE)
 return(data.frame(y=x1$fit,se=x1$se.fit,
                   x=seq(0,to1,0.1),
                   location=loc))
}

p_outs <- rbind(p_out2("Baraka",5),p_out2("Walikale",5),p_out2("Mweso",5),p_out2("Lulimba",5),p_out2("Shamwana",9))

plot_points <- temp2  %>% ggplot(aes(y=ANCprevalence,x=un5inc,col=location)) + geom_point(size=2) + theme_bw() + xlab("Cases per child per year") + ylab("Prevalence of infection in pregnant women") +
  scale_color_discrete(name="MSF hospital location:") +
  theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"),strip.text.x=element_text(size=14,face="bold"),
        strip.text.y=element_text(size=14,face="bold"),
        axis.text=element_text(size=12),axis.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14,face="bold"),legend.position = "bottom",legend.title=element_text(size=16,face="bold")) +
  geom_point(data=yearav,aes(y=x,x=y,fill=location),shape=22,col="black",size=4,show.legend = FALSE) +
 geom_ribbon(data=p_outs,inherit.aes=FALSE,aes(x=x,ymin=y-1.96*se,ymax=y+1.96*se,fill=location),alpha=0.2) +
 geom_line(data=p_outs,aes(x=x,y=y,col=location)) +
 scale_fill_discrete(guide="none") + coord_cartesian(ylim=c(0,0.7))


plot(plot_points)


# Code for Fig 6.


df7 <- data.frame("NENL" = NENL$fitted.values,
                  "date"=temp %>% arrange(location) %>% group_by(location) %>% slice(-c(1:3)) %>% ungroup %>% .$date,
                  "location"=temp %>% arrange(location) %>% group_by(location) %>% slice(-c(1:3)) %>% .$location,
                  "visits"=temp %>% arrange(location) %>% group_by(location) %>% slice(-c(1:3)) %>% .$ANCvisits,
                  "ANC"=temp %>%  arrange(location) %>% group_by(location) %>% slice(-c(1:3)) %>% .$ANCprevalence)
NENLdf <- full_join(data.frame(date=temp$date,location=temp$location,type="NENL",obs=temp$ANCprevalence),
                    ndf %>% dplyr::select(date,location,CIlow=y_NENL_lwr,CIhi=y_NENL_upr)) %>%
  full_join(temp %>% group_by(location) %>% slice(-c(1:3)) %>% dplyr::select(date,location) %>% ungroup) %>%
  full_join(df7 %>% dplyr::select(date,visits,location,"val"=NENL))
NENLdf$location <- gdata::reorder.factor(NENLdf$location,new.order=target)
rocv_plot <- NENLdf %>% mutate(ocil = obs - (1.96*sqrt(1/visits)*obs*(1-obs)),
                 ocih = obs + (1.96*sqrt(1/visits)*obs*(1-obs))) %>%
  ggplot() + facet_grid(location~type,scales="free") + theme_bw() +
  geom_vline(xintercept=as.numeric(seq.Date(from=as.Date("2009/07/01"),to=as.Date("2016/12/01"),by="year")),lty=2,alpha=0.25) +
  geom_vline(xintercept=as.numeric(seq.Date(from=as.Date("2009/01/01"),to=as.Date("2016/12/01"),by="year")),alpha=0.35) +
  geom_ribbon(aes(x=date,ymin=ocil,ymax=ocih,fill=location,alpha=0.5)) +
  geom_ribbon(aes(x=date,ymin=CIlow,ymax=CIhi),alpha=0.5) +
  geom_line(aes(x=date,y=obs,col=location),size=1.05) +
  geom_line(aes(x=date,y=val),size=1.05) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.x=element_blank(),
        strip.text.y=element_text(size=14,face="bold"),
        axis.title = element_text(size=14,face="bold"),
        legend.position = "none",
        axis.text=element_text(size=13)) +
  scale_x_date(date_breaks = "2 years", date_labels =  "%Y",expand=c(0.1,0)) +
  scale_y_continuous(expand=c(0,0.1)) +
  ylab("Prevalence of infection in pregnant women") + xlab("Date")

rocv_plot

ndf2 <- ndf
ndf2$location <- gdata::reorder.factor(ndf2$location,new.order=target)
ndf2 %>% ggplot() + geom_point(aes(x=ANCprevalence,y=y_NENL,col=location)) + theme_bw() +
  geom_abline(aes(intercept=0,slope=1),size=1.1) + 
  geom_errorbarh(aes(col=location,y=y_NENL,x=ANCprevalence,xmin=ANCprevalence - (1.96*sqrt(1/ANCvisits)*ANCprevalence*(1-ANCprevalence)),
                                            xmax=ANCprevalence + (1.96*sqrt(1/ANCvisits)*ANCprevalence*(1-ANCprevalence)))) +
  geom_errorbar(aes(col=location,ymin=y_NENL_lwr,ymax=y_NENL_upr,x=ANCprevalence)) +
  xlab("Observed ANC prevalence") + ylab("Cross-validation estimate of ANC prevalence by NENL model") +
  theme(axis.title = element_text(size=14,face="bold"),legend.text=element_text(size=14),legend.title=element_text(size=14,face="bold")) +
  scale_color_discrete(name="MSF location")

# Code for Fig 5.

cp_NELL <- crosspred(basis=cb_nell,NELL,from = 0.5,to=7.5,by = 0.1,cen = 1,cumul = TRUE)
cp_NENL <- crosspred(basis=cb_nenl,NENL,from = 0.5,to=7.5,by = 0.1,cen = 1,cumul = TRUE)

plot(cp_NENL,col="grey")
