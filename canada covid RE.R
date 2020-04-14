setwd('~/Dropbox/canada covid')
dat = read.csv('canada cases.csv',header=T)


ids = unique(dat$Case.identifier.number)

vals = as.character(unique(dat$Case.information))

month = day = sex = agegrp = trans = hosp = icu = c()
for (i in 1:length(ids)){
  month[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[3])]
  day[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[4])]
  sex[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[5])]
  agegrp[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[6])]
  trans[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[7])]
  hosp[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[8])]
  icu[i] = dat$VALUE[which(dat$Case.identifier.number==ids[i]&dat$Case.information==vals[9])]
  print(i)
}

month[month==99] = NA; day[day==99] = NA

caseDayName = paste(day,'-',month,'-','2020',sep='')
library(lubridate)
caseDay = dmy(caseDayName)
#which(is.na(caseDay))
#### impute for missing days the nearest available
caseDay[which(is.na(caseDay))] = caseDay[which(is.na(caseDay))-1]
caseDayNum = as.numeric(caseDay)
daynums = min(caseDayNum)+c(-30:100)
daynames = min(caseDay)+c(-30:100)

maxday = as.numeric(max(caseDay))
incl = which(daynums<=maxday)
daynums = daynums[incl]
daynames = daynames[incl]


sex[sex%in%c(7,9)] = NA; sex[sex==2] = 0
hosp[hosp==9] = NA; hosp[hosp==2] = 0 
icu[icu==9] = NA; icu[icu==2] = 0; icu[hosp==0] = 0
trans[trans==1] = 0; trans[trans==2] = 1; trans[trans==3] = NA

caseDat = cbind(caseDayNum,sex,agegrp,trans,hosp,icu)

newCas = c()
for (i in 1:length(daynums)){newCas[i] = sum(caseDayNum==daynums[i])}


t = 1:length(daynums)#1:(length(daynums)-1)
dens = dweibull(t,1.455,5.983) ### serial interval density
cumDens = pweibull(t,1.455,5.983)

daynames[107]
quantile(curShed[,107],c(0.5,0.025,0.975))
#############################################
#### for each of 1000 iterations ############
#############################################
### 1) sample hosp status (impute) ##########
### 2) sample infeciton time ################
### 3) sample total infections expected #####
### 4) compute Re ###########################
#############################################

agegrps = unique(agegrp)

hist(shedSx/(shedSx+rexp(1e5,1/9)))

set.seed(1)
infectShed = rweibull(1e5,1.455,5.983)
shedSx = rweibull(1e5,0.294,0.14)
sxHosp = rgamma(1e5,shape=5.0777220,rate=0.7645327)
totTime = infectShed + shedSx + sxHosp

#### generate probability of hosp|infect, take weighted average across ages for MISSING age group
load('pHospSx.Rdata')
ageprobs = table(agegrp[hosp==1])[1:8]/sum(table(agegrp[hosp==1])[1:8])
pHosp = array(NA,dim=c(1e3,9)); pHosp[,1:8] = pHospSx[,2:9]; pHosp[,9] = pHosp[,1:8]%*%ageprobs

set.seed(1)
totInfAge = array(NA,dim=c(1e3,9,length(daynums)))
totInf = array(NA,dim=c(1e3,length(daynums)))
for (k in 1:1e3){
  for (i in 1:length(agegrps)){
    sel = which(is.na(hosp)&agegrp==agegrps[i])
    hospImp[sel] = runif(length(sel),0,1)<mean(hosp[agegrp==agegrps[i]],na.rm=T)
    
    sel = which(hospImp==1&agegrp==agegrps[i]); tot = length(sel)
    dateInf = rep(NA,tot)
    dateInf = floor(caseDayNum[sel] - sample(timeTot,tot,replace=T))
    for (j in 1:length(daynums)){
      totInfAge[k,i,j] = (sum(dateInf==daynums[j],na.rm=T)/mean(totTime<=(maxday-daynums[j])))/pHosp[k,i]
    }
  }
  for (j in 1:length(daynums)){
    totInf[k,j] = sum(totInfAge[k,,j])
  }
  print(k)
}


Re = array(NA,dim=c(1e3,length(daynums)-7))
for (i in 1:1e3){
    cas = totInf[i,1:(length(daynums)-7)]
    p = matrix(0,length(cas),length(cas))
    for (j in 2:length(cas)){
      if (cas[j]>0){
        for (k in 1:(j-1)){
          if (cas[k]>0){
            p[j,k] = dens[j-k]/sum(dens[1:(j-1)]*cas[(j-1):1])
          }
        }
      }
    }
    rj = c()
    for (j in 1:length(cas)){
      rj[j] = (p[,j]%*%cas)/mean(timeIncub<=(maxday-daynums[j]))
    }
    Re[i,] = rj
  print(i)
}

q95fn = function(x){return(quantile(x,c(0.5,0.025,0.975),na.rm=T))}
reQ = apply(Re,2,q95fn)



set.seed(1)
timeIncub = rweibull(1e5,1.455,5.983)
timeShed = rexp(1e5,1/9)
timeTot = timeIncub+timeShed

curShedAge = curRecovAge = cumInfectAge = array(NA,dim=dim(totInfAge))
curShed = curRecov = cumInfect = array(NA,dim=dim(totInf))
pRecov = c(); for (i in 0:500){pRecov[i] = mean(timeTot<i)}
pIncub = c(); for (i in 0:500){pIncub[i] = mean(timeIncub<i)}
pShed = pIncub - pRecov
for (i in 1:1e3){
  for (j in 1:9){
    for (k in 1:length(daynums)){
      curRecovAge[i,j,k] = totInfAge[i,j,1:k]%*%pRecov[k:1]
      curShedAge[i,j,k] = totInfAge[i,j,1:k]%*%pShed[k:1]
      cumInfectAge[i,j,k] = sum(totInfAge[i,j,1:k])
    }
  }
  for (k in 1:length(daynums)){
    curRecov[i,k] = sum(curRecovAge[i,,k])
    curShed[i,k] = sum(curShedAge[i,,k])
    cumInfect[i,k] = sum(cumInfectAge[i,,k])
  }
  print(i)
}
curIncubAge = cumInfectAge - (curRecovAge+curShedAge)
curIncub = cumInfect - (curRecov+curShed)

save(curShedAge,file='curShedAge.Rdata'); save(curShed,file='curShed.Rdata')
save(curIncubAge,file='curIncubAge.Rdata'); save(curIncub,file='curIncub.Rdata')
save(cumInfectAge,file='cumInfectAge.Rdata'); save(cumInfect,file='cumInfect.Rdata')
save(curRecovAge,file='curRecovAge.Rdata'); save(curRecov,file='curRecov.Rdata')
save(totInfAge,file='totInfAge.Rdata'); save(totInf,file='totInf.Rdata')

infQ = apply(cumInfect,2,q95fn)/(37.59e6)


reSel = 48:100
infSel = 48:107

months = month(daynames)
days =day(daynames)

monthlabs = c('Jan','Feb','Mar','Apr','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
daychar = paste(monthlabs[months],' ',days,sep='')

setwd('~/Dropbox/canada covid')
save(daynames,file='daynames.Rdata')


pdf(file='rePlot.pdf',width=6,height=2.5)
par(mfrow=c(1,2))
par(lwd=0.5)
yub = 4
par(mar=c(3.25,2,1,1)); par(mgp=c(3,0.35,0)); par(tck=-0.02)
plot(reQ[1,reSel],x=daynums[reSel],type='n',ylim=c(0,yub),axes=F,ann=F)
mtext(side=3,'A',cex=0.75,font=2,adj=0)
abline(h=1,lwd=0.5,col='darkgoldenrod3')
polygon(y=c(reQ[2,reSel],rev(reQ[3,reSel])),x=daynums[c(reSel,rev(reSel))],col=rgb(1,0,0,0.25),lty=0)
lines(reQ[2,reSel],x=daynums[reSel],lwd=0.25,col='darkred')
lines(reQ[3,reSel],x=daynums[reSel],lwd=0.25,col='darkred')
lines(reQ[1,reSel],x=daynums[reSel],lwd=0.75,col='darkred')
text(x=seq(min(daynums[reSel]),max(daynums[reSel]),7),
     y=-0.1*yub,
     daychar[seq(min(reSel),max(reSel),7)],xpd=T,srt=45,adj=1,cex=0.65)
box(bty='l')
axis(side=2,at=0:yub,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
axis(side=1,at=seq(min(daynums[reSel]),max(daynums[reSel]),7),labels=NA,lwd=0,lwd.ticks=0.5)
mtext(side=1,'Date infection acquired',cex=0.75,line=2.25)
mtext(side=2,expression(italic(R)[italic(E)]),cex=0.75,line=1)


yub = 0.01
plot(infQ[1,infSel],x=daynums[infSel],type='n',ylim=c(0,yub),axes=F,ann=F)
mtext(side=3,'B',cex=0.75,font=2,adj=0)
polygon(y=c(infQ[2,infSel],rev(infQ[3,infSel])),x=daynums[c(infSel,rev(infSel))],col=rgb(0,0,1,0.25),lty=0)
lines(infQ[2,infSel],x=daynums[infSel],lwd=0.25,col='darkblue')
lines(infQ[3,infSel],x=daynums[infSel],lwd=0.25,col='darkblue')
lines(infQ[1,infSel],x=daynums[infSel],lwd=0.75,col='darkblue')
text(x=seq(min(daynums[infSel]),max(daynums[infSel]),7),
     y=-0.1*yub,
     daychar[seq(min(infSel),max(infSel),7)],xpd=T,srt=45,adj=1,cex=0.65)
box(bty='l')
axis(side=2,at=seq(0,yub,0.002),labels=seq(0,yub,0.002)*100,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
axis(side=1,at=seq(min(daynums[infSel]),max(daynums[infSel]),7),labels=NA,lwd=0,lwd.ticks=0.5)
mtext(side=1,'Date',cex=0.75,line=2.25)
mtext(side=2,'Cumulative infection prevalence (%)',cex=0.75,line=1.25)
dev.off()


