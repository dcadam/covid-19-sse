###This script was modified from the below reference. 

############################################################################
#                                                                          #
#   THE ROLE OF SUPERSPREADING IN MIDDLE EAST RESPIRATORY                  #
#   SYNDROME CORONAVIRUS (MERS-COV) TRANSMISSION                           #
#                                                                          #
#   REFERENCE: KUCHARSKI AJ, ALTHAUS CL. EURO SURVEILL. 2015;20:pii=21167. #
#   WWW: http://www.eurosurveillance.org/ViewArticle.aspx?ArticleId=21167  #
#                                                                          #
#   Copyright (C) 2015 by Adam J. Kucharski (adam.kucharski@lshtm.ac.uk)   #
#   and Christian L. Althaus (christian.althaus@alumni.ethz.ch)            #
#                                                                          #
#   This program is free software; you can redistribute it and/or modify   #
#   it under the terms of the GNU General Public License as published by   #
#   the Free Software Foundation; either version 2 of the License, or      #
#   (at your option) any later version.                                    #
#                                                                          #
#   This program is distributed in the hope that it will be useful,        #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#   GNU General Public License for more details.                           #
#                                                                          #
#   You should have received a copy of the GNU General Public License      #
#   along with this program; if not, write to the                          #
#   Free Software Foundation, Inc.,                                        #
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              #
#                                                                          #
############################################################################

#Brancing Process Model for local clusters

#Data on Local Clusters sizes
chaindataSARS <- c(rep(1,46), rep(2,9), rep(3,3), rep(4,2), rep(5,1),rep(6,1),rep(9,3),19,22,106)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Functions to estimate R and k

# Outbreak size distribution function
hg.rj <- function(L,R,k) {
  p <- exp(lgamma(k*L+L-1) - lgamma(k*L) - lgamma(L+1) + log((R/k)^(L-1)/(1+R/k)^(k*L+L-1)))
  return(p)
}

# Likelihood function
likelihoodHG<- function(simdata,rstar,kstar){
  runs=length(simdata)
  liks=rep(1,runs)
  for(z in 1:runs){liks[z]=hg.rj(simdata[z],rstar,kstar)}
  sum(log(liks))
}

# Set up grid search
R0r1=0.10; R0r2=3.00; R0range=seq(R0r1,R0r2,by=.01)
kkr1=0.01; kkr2=55; kkrange=seq(kkr1,kkr2,by=0.01)

# Calculate likelihood surface
R0estimateHG<- function(simdata,R0range,kkrange){
  collect1 <- data.frame(matrix(NA, nrow=length(R0range),length(kkrange)))
  for(ii in 1:length(R0range)){
    for(jj in 1:length(kkrange)){
      rr=R0range[ii]
      kkr=kkrange[jj]
      collect1[ii,jj]=likelihoodHG(simdata,rr,kkr)
    }
  }
  collect1
}

liksurfSARS=R0estimateHG(chaindataSARS,R0range,kkrange)
likmSARS=(liksurfSARS==max(liksurfSARS))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate profile likelihoods

# For reference SARS 0.16 (90% confidence interval 0.11–0.64)
calc_profile<-function(liksurfSX,likmSX,chiV){
  
  prfk=apply(liksurfSX,2,function(x){max(x)})
  prfk2=kkrange[prfk-max(prfk)>-chiV]
  prfR=apply(liksurfSX,1,function(x){max(x)})
  prfR2=R0range[prfR-max(prfR)>-chiV]
  
  c(
    paste(R0range[sum(seq(1,length(R0range))%*%likmSX)]," (",min(prfR2),"-",max(prfR2),")",sep=""),
    paste(kkrange[sum(likmSX%*%seq(1,length(kkrange)))]," (",min(prfk2),"-",ifelse(max(prfk2)==max(kkrange),"Inf",max(prfk2)),")",sep="")
  )
  
}


prfk=apply(liksurfSARS,2,function(x){max(x)})
prfk2=kkrange[prfk-max(prfk)>-chiV]



#90% Confidence interval for R and K

chiV=qchisq(90/100, df=1)/2 
outputM1=calc_profile(liksurfSARS,likmSARS,chiV)

#95% cofidence interval for R and K
chiV=qchisq(95/100, df=1)/2 
outputM1=calc_profile(liksurfSARS,likmSARS,chiV)


output_parameters<-function(conf.interval){
  
  chiV=qchisq(conf.interval/100, df=1)/2 
  outputM1=calc_profile(liksurfSARS,likmSARS,chiV)
  
  outputdata=data.frame(outputM1)
  colnames(outputdata)=c(conf.interval)
  rownames(outputdata)=c('R', 'k')
  
  write.csv(outputdata,paste("param_estimates",conf.interval,".csv",sep=""))
  
}

output_parameters(90)
output_parameters(95)

#calculation of proportion of cases who do not spread to anyone 
dnbinom(0, size = 0.33, mu = 0.74)
dnbinom(0, size = 0.14, mu = 0.58)
dnbinom(0, size = 0.98, mu = 0.97)


####Supplementary Wave analysis WAVE ONE

#Cluster sizes by wave one (before March) and wave two (march onwards)
chaindataSARS_wave1 <- c(rep(1,11), rep(2,4), rep(3,2),5,9,19) #FIXED TO NEW PAIRS

liksurfSARS=R0estimateHG(chaindataSARS_wave1,R0range,kkrange)
likmSARS=(liksurfSARS==max(liksurfSARS))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate profile likelihoods

# For reference SARS 0.16 (90% confidence interval 0.11–0.64)
calc_profile<-function(liksurfSX,likmSX,chiV){
  
  prfk=apply(liksurfSX,2,function(x){max(x)})
  prfk2=kkrange[prfk-max(prfk)>-chiV]
  prfR=apply(liksurfSX,1,function(x){max(x)})
  prfR2=R0range[prfR-max(prfR)>-chiV]
  
  c(
    paste(R0range[sum(seq(1,length(R0range))%*%likmSX)]," (",min(prfR2),"-",max(prfR2),")",sep=""),
    paste(kkrange[sum(likmSX%*%seq(1,length(kkrange)))]," (",min(prfk2),"-",ifelse(max(prfk2)==max(kkrange),"Inf",max(prfk2)),")",sep="")
  )
  
}


prfk=apply(liksurfSARS,2,function(x){max(x)})
prfk2=kkrange[prfk-max(prfk)>-chiV]

#90% Confidence interval for R and K

chiV=qchisq(90/100, df=1)/2 
outputM1=calc_profile(liksurfSARS,likmSARS,chiV)

#95% cofidence interval for R and K
chiV=qchisq(95/100, df=1)/2 
outputM1=calc_profile(liksurfSARS,likmSARS,chiV)


output_parameters<-function(conf.interval){
  
  chiV=qchisq(conf.interval/100, df=1)/2 
  outputM1=calc_profile(liksurfSARS,likmSARS,chiV)
  
  outputdata=data.frame(outputM1)
  colnames(outputdata)=c(conf.interval)
  rownames(outputdata)=c('R', 'k')
  
  write.csv(outputdata,paste("param_estimates_wave_1",conf.interval,".csv",sep=""))
  
}

output_parameters(90)
output_parameters(95)

#calculation of proportion of cases who do not spread to anyone (wave 1) 
dnbinom(0, size = 2.31, mu = 0.66)
dnbinom(0, size = 0.22, mu = 0.44)
dnbinom(0, size = Inf, mu = 0.99)


####WAVE TWO

chaindataSARS_wave2 <- c(rep(1,35), rep(2,5), rep(3,1), rep(4,2),6,rep(9,2),22,106) #FIXED TO NEW PAIRS

liksurfSARS=R0estimateHG(chaindataSARS_wave2,R0range,kkrange)
likmSARS=(liksurfSARS==max(liksurfSARS))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate profile likelihoods

# For reference SARS 0.16 (90% confidence interval 0.11–0.64)
calc_profile<-function(liksurfSX,likmSX,chiV){
  
  prfk=apply(liksurfSX,2,function(x){max(x)})
  prfk2=kkrange[prfk-max(prfk)>-chiV]
  prfR=apply(liksurfSX,1,function(x){max(x)})
  prfR2=R0range[prfR-max(prfR)>-chiV]
  
  c(
    paste(R0range[sum(seq(1,length(R0range))%*%likmSX)]," (",min(prfR2),"-",max(prfR2),")",sep=""),
    paste(kkrange[sum(likmSX%*%seq(1,length(kkrange)))]," (",min(prfk2),"-",ifelse(max(prfk2)==max(kkrange),"Inf",max(prfk2)),")",sep="")
  )
  
}


prfk=apply(liksurfSARS,2,function(x){max(x)})
prfk2=kkrange[prfk-max(prfk)>-chiV]



#90% Confidence interval for R and K

chiV=qchisq(90/100, df=1)/2 
outputM1=calc_profile(liksurfSARS,likmSARS,chiV)

#95% cofidence interval for R and K
chiV=qchisq(95/100, df=1)/2 
outputM1=calc_profile(liksurfSARS,likmSARS,chiV)


output_parameters<-function(conf.interval){
  
  chiV=qchisq(conf.interval/100, df=1)/2 
  outputM1=calc_profile(liksurfSARS,likmSARS,chiV)
  
  outputdata=data.frame(outputM1)
  colnames(outputdata)=c(conf.interval)
  rownames(outputdata)=c('R', 'k')
  
  write.csv(outputdata,paste("param_estimates_wave2",conf.interval,".csv",sep=""))
  
}

output_parameters(90)
output_parameters(95)


#calculation of proportion of cases who do not spread to anyone (wave 2)
dnbinom(0, size = 0.2, mu = 0.77)
dnbinom(0, size = 0.08, mu = 0.54)
dnbinom(0, size = 0.63, mu = 1.13)
