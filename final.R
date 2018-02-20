library(MotifDb)
library(MotIV)
library(seqLogo)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)

#Gather Data
genome = BSgenome.Hsapiens.UCSC.hg19
sox = query(query(query(MotifDb,'sapien'),'Sox'),'jolma')
soxx=c(sox[1],sox[2],sox[7],sox[10],sox[13],sox[16],sox[20],sox[27],sox[30])
unlist(soxx)
seqLogo(as.list(soxx)[[1]])

#Transcription Factors PWM
sox9 = (matchPWM(unlist(soxx[1]),genome$chr1)) #9
sox10 = (matchPWM(unlist(soxx[2]),genome$chr1)) #15
sox14 = (matchPWM(unlist(soxx[3]),genome$chr1)) #12
sox15 = (matchPWM(unlist(soxx[4]),genome$chr1)) #15
sox18 = (matchPWM(unlist(soxx[5]),genome$chr1)) #15
sox21 = (matchPWM(unlist(soxx[6]),genome$chr1)) #15
sox2 = (matchPWM(unlist(soxx[7]),genome$chr1)) #17
sox7 = (matchPWM(unlist(soxx[8]),genome$chr1)) #16
sox8 = (matchPWM(unlist(soxx[9]),genome$chr1)) #16

#Short Binding Sites
smallbs = c(start(sox9),start(sox14))
hist(smallbs)
plot(density(smallbs))
head(smallbs)
tail(smallbs)
smallbins = seq(0,2.5*10^8,by=1*10^7)
smallbins = smallbins[1:25]
smallbs.cut = cut(smallbs,smallbins,right=FALSE)
a = table(smallbs.cut)
freqsmallbs = c(17895,19682,20572,22840,25517,28849,30167,35574,34555,32751,35236,29104,3732,0,17685,25993,28926,30252,31651,37543,23364,31494,27571,29111,27572)
plot(smallbins,freqsmallbs)
smallmodel <- lm(freqsmallbs ~ poly(smallbsmallins,2))
#predsmall = smallmodel$coefficients[1]+smallmodel$coefficients[2]*smallbins+smallmodel$coefficients[3]*smallbins^2
#points(smallbins, predict(smallmodel), type="l", col="red", lwd=2)
dfsmall = data.frame(Position=smallbins,Frequency=freqsmallbs)
smallplot = ggplot(dfsmall,aes(x=Position,y=Frequency))+geom_point()+geom_smooth()+labs(title="Short Binding Sites")
dfsmallfront = data.frame(Position=smallbins[1:8],Frequency=freqsmallbs[1:8])
smallfrontmodel = lm(dfsmallfront$Frequency~dfsmallfront$Position)
dfsmallback = data.frame(Position=smallbins[20:25],Frequency=freqsmallbs[20:25])
smallbackmodel = lm(dfsmallback$Frequency~dfsmallback$Position)

#Long Binding Sites
largebs = c(start(sox2),start(sox7),start(sox8))
hist(largebs)
head(largebs)
tail(largebs)
plot(density(largebs))
largebins = seq(0,2.5*10^8,by=2*10^7)
largebins = largebins[1:12]
largebs.cut = cut(largebs,largebins,right=FALSE)
b = table(largebs.cut)                                                                                       
freqlargebs = c(916,1087,1473,2196,2198,2137,106,1119,1880,2504,1737,1676)
dflarge = data.frame(Position=largebins,Frequency=freqlargebs)
largeplot = ggplot(dflarge,aes(x=Position,y=Frequency))+geom_point()+geom_smooth()+labs(title="Large Binding Sites")
dflargefront = data.frame(Position=largebins[1:5],Frequency=freqlargebs[1:5])
largefrontmodel = lm(dflargefront$Frequency~dflargefront$Position)
dflargeback = data.frame(Position=largebins[10:12],Frequency=freqlargebs[10:12])
largebackmodel = lm(dflargeback$Frequency~dflargeback$Position)

#Analyze
anova(smallbackmodel,largebackmodel)
anova(smallfrontmodel,largefrontmodel) #significant
anova(smallfrontmodel,largebackmodel) #significant
anova(smallbackmodel,largefrontmodel)
