#A = NNO
#B = N15NO
#C = 15NNO
#D = NN18O

#read in necessary packages for spline fitting 
library("splines")
library("MASS")
library("ConSpline")
library("cgam")

#first correcting the A data
#read in data
library("readxl")
caldata<- read_excel("Desktop/Dissertation Chapters/Calibration model/N2O analyzer calibration dataset_github.xlsx", sheet = "standards w indicators_B&C")
View(caldata)

#extracting observed N2O concentrations and true N2O concentrations
obs<-caldata$`[14N14N16O] = A` 
tru<-caldata$"TRUEA"
tru<-tru+0.0426 #the 0.0426 accounts for the trace N2O contamination in the zero-grade air

#plotting the true vs. observed [N2O]
plot(tru,obs)
abline(0,1)

#we want to do the polynomial fit on the log scale to better observe the differences in [N2O] and improve the fit 
logobs<-log(obs) 
logtru<-log(tru)
plot(logtru,obs)
plot(logtru,logobs)

#first we need to correct the A data in small, medium, and large chunks so that we can estimate
#accurrate B/A ratios later on to then determine the secondary standard alpha value. 

###########################
#####SMALL A CHUNK#########
#correcting A < or = to 1.49ppm N2O 
logobs.4<-logobs[which(logobs<=log(1.5))] #separating out the standards that are less than or equal to 1.49ppm 
logtru.4<-logtru[which(tru <=1.565)] #had to set it at 1.565 of non-logged true so the number of samples
#would match up for logobs vs. logtru (added in 0.0426 for contaminant) 
m1s<-lm(logobs.4~logtru.4+I(logtru.4^2)) #m1s = model 1 small NNO range
summary(m1s) ##summary tells us that X and X^2 are both significantly different from zero 
plot(logtru.4,logobs.4)
lines(sort(logtru.4), m1s$fitted.values[order(logtru.4)], lwd = 2, col = "red")
residm1s<-m1s$residuals #note that plots of residuals for each NNO range are now more homoscedastic 
plot(logtru.4,residm1s)

#solving quadratic equation to calculate A for 0-1.49ppm N2O   
summary(m1s) #quadratic equation formula: ax^2 + bx +c 
a1<-m1s$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1s$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1s$coefficients[1]-logobs.4

#solving the quadratic equation - when A < or = to 1.49ppm
top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the corrected small A values
Acorrsm<-exp(xhat1)
AcorrsmUSGS<-Acorrsm[c(17:20)] #separates USGS A vals
AcorrsmTANK<-Acorrsm[c(1:16)] #separates secondary standard A vals

############################
#####MEDIUM A CHUNK#########
#correcting A >1.49ppm and < or = to 33.12ppm N2O 
logobs33<-logobs[which(logobs>log(1.5)&logobs<=log(33.18))] #separating out the standards >1.49 and < or =33.12ppm N2O
logtru33<-log(tru[which(tru>1.565&tru<=33.185)]) #again, needed to tweak values a bit so numbers would match up
m1m<-lm(logobs33~logtru33+I(logtru33^2)) #m1m = model 1 medium NNO range
summary(m1m) ##summary tells us that X and X^2 are both significantly different from zero 
plot(logtru33,logobs33)
lines(sort(logtru33), m1m$fitted.values[order(logtru33)], lwd = 2, col = "red")
residm1m<-m1m$residuals
plot(logtru33,residm1m)

#solving quadratic equation to calculate A for >1.49ppm and < or = to 33.12ppm N2O  
summary(m1m) #quadratic equation formula: ax^2 + bx +c 
a1<-m1m$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1m$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1m$coefficients[1]-logobs33

#solving the quadratic equation when A >1.49ppm and < or = to 33.12ppm N2O 
top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the corrected medium A values
Acorrmed<-exp(xhat1)
AcorrmedUSGS<-Acorrmed[41:50] 
AcorrmedTANK<-Acorrmed[1:40]

##########################
#####LARGE A CHUNK#########
#correcting A >33.12 and < or = to 300ppm N2O 
logobs300<-logobs[which(logobs>log(33.18)&logobs<=max(logobs))] #separating out the standards >33.12ppm N2O
logtru300<-log(tru[which(tru>33.185&tru<=300.065)]) #again, kind of tweaked numbers to account for zero-air contaminant
m1l<-lm(logobs300~logtru300+I(logtru300^2)) #m1l = model 1 large NNO range
summary(m1l) ##summary tells us that X and X^2 are both significantly different from zero 
plot(logtru300,logobs300)
lines(sort(logtru300), m1l$fitted.values[order(logtru300)], lwd = 2, col = "red")
residm1l<-m1l$residuals
plot(logtru300,residm1l)

#solving quadratic equation to calculate A for >33.12 up to 300ppm N2O   
summary(m1l) #quadratic equation formula: ax^2 + bx +c 
a1<-m1l$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1l$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1l$coefficients[1]-logobs300

#solving the quadratic equation when A >33.12 up to 300ppm N2O 
top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the corrected large A values
Acorrlg<-exp(xhat1)
AcorrlgUSGS<-Acorrlg[c(17:21)] 
AcorrlgTANK<-Acorrlg[c(1:16)] 

#########################################################################################
#########################################################################################
#########################################################################################
#Now correcting the USGS52 B data (not in chunks) across the entire [N2O] range because this 
#correction function will be used subsequently to correct the unknown secondary standard alpha 
#values across the entire [N2O] range. 
Bobs<-caldata$`[14N15N16O] = B` #extracting B values from dataset 
AcorrUSGS<-c(AcorrsmUSGS,AcorrmedUSGS,AcorrlgUSGS) #corrected A values come from the chunked
#m1's (m1s, m1m, m1l)

#now separating out USGS only data b/c we'll use the known USGS alpha value as our correction 
usgs <- which(caldata$Source == "USGS52 stock")
Bobs1=Bobs[usgs] #B values, USGS only 

#now calculating true B values from USGS data 
#we did this by rearranging eqn delta=((B/A-std)/std)*1000 b/c if we know the delta val =13.52, and we
#also, since we corrected the A values, we can plug in the known delta (13.52) and 
#the known As and solve for the concentration of each B 
std<-0.0036765
#mixing model to account for contaminant in zero-air in standards 
Btru1<-(((13.52/1000)*std)+std)*(AcorrUSGS-0.0426)+(((0/1000)*std)+std)*(0.0426)#this is the rearranged eqn: B = (((13.52/1000)*std)+std)*Acorr
Btru1 #true B concentrations
plot(Btru1,Bobs1)
#we are log transforming the B values for the same reasons we log transfored A - it improves 
#the polynomial fit when we use a log scale!
logBobs1<-log(Bobs1) 
logBtru1<-log(Btru1)
plot(logBtru1,logBobs1)

#m2 = model to correct all the B concentrations for the USGS52 N2O standard
m2<-lm(logBobs1~ logBtru1 + I(logBtru1^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logBtru1), m2$fitted.values[order(logBtru1)], lwd = 2, col = "red")
resid<-m2$residuals
plot(logBtru1,resid)

#this is solving for the USGS-only B coefficients 
summary(m2) #tells us if x and x^2 are significantly different from zero 
a2<-m2$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2$coefficients[1]-logBobs1 # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale USGS B concentrations

exp(xhat2) #these are the corrected non-log scale USGS B concentrations 
BcorrUSGS<-exp(xhat2)

#now calculating d15Nalpha values for the USGS52 data across the entire [N2O] range 
ba1<-BcorrUSGS/AcorrUSGS #the corrected B []s divided by the corrected USGS A []s
ba1 #these are the B/A  ratios 
deltaB1<-(((ba1-std))/std)*1000 #we use the B/A ratios to calculate the B delta values
deltaB1
mean(deltaB1)
plot(AcorrUSGS,deltaB1) #plot of corrected deltas vs. corrected [N2O]s
plot(log(AcorrUSGS),deltaB1) #log-scale helps us better visualize the correction

######################################################################################################
######################################################################################################
######################################################################################################
#Now fitting m2 to the 500ppm tank (secondary standard) data to correct the B values from the 
#secondary standard, since the d15Nalpha for the secondary standard was not known until we 
#estimated it. We did this by using the calibrated d15Nalpha values from the USGS52 standard 

#separating out secondary standard data  to then calibrate using the model fitted 
#from the USGS B data (m2) 
stock <- which(caldata$Source == "500ppm N2O stock")
Bobs2=Bobs[stock] #B []s from secondary standard only - these have not been corrected yet 
logBobs2<-log(Bobs2)   

#this is solving for the secondary standard B coefficients
summary(m2) #tells us if x and x^2 are significantly different from zero 
a2<-m2$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2$coefficients[1]-logBobs2 # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale secondary standard B concentrations

exp(xhat2) #these are the corrected non-log scale secondary standard B concentrations 
BcorrTANK<-exp(xhat2)

AcorrTANK<-c(AcorrsmTANK,AcorrmedTANK,AcorrlgTANK) #corrected secondary standard A values 

#now we'll estimate the delta values for the secondary standard. We will use these values to 
#generate a spline fit and THEN we'll add the coefficients from that spline to the corrected 
#mean USGS d15Nalpha value. This subsequent vector of deltas will then be averaged, and that 
#final value will be the secondary standard d15Nalpha. 
ba2<-exp(xhat2)/AcorrTANK 
ba2 
deltaB2<-(((ba2-std))/std)*1000
deltaB2
mean(deltaB2)
plot(AcorrTANK,deltaB2)
plot(log(AcorrTANK),deltaB2)

######################################################################################################
######################################################################################################
######################################################################################################
#Now using m3 to estimate the delta B value of the secondary standard

#CONCATONATING THE USGS AND SECONDARY STANDARD DATA 
x1<-c(AcorrTANK,AcorrUSGS) #these are all of the corrected A concentrations for USGS and the secondary standard
y<-c(deltaB2,deltaB1) #these are all of the corrected B delta values for USGS and the secondary standard
x2<-c(caldata$Indicator)  #this is an indicator term, 1=if tank, 0=if USGS

logx1<- log(x1) #need to log scale the corrected A values 

#setting up model 
#model m3 = spline fit for estimating coefficient to add to primary standard to then 
#estimate secondary standard 
myknots = quantile(logx1, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) #these knots were determined by a statistician on our team 
m3 <- lm( y ~ bs(logx1, knots = myknots)+x2)
summary(m3)
plot(log(AcorrTANK),deltaB2)
points(log(AcorrUSGS),deltaB1, col="red")
lines(sort(logx1), m3$fitted.values[order(logx1)], lwd = 2, col = "red")
#check the residuals of the fit 
residm3<-m3$residuals
plot(logx1,residm3)

#now we can estimate the secondary standard true d15Nalpha value
tankB<-mean(deltaB1) + m3$coefficients[10] 
tankB#this is the estimated mean secondary standard d15Nalpha

#mean(deltaB1) = the averaged corrected USGS52 delta B values
#m3$coefficients[10] = the coefficient from m4 that we added to the USGS mean d15Nalpha because 
#that's the offset between the two models.

#note that the tankB value calculated here differs slightly from the tankB value written in
#Stuchiner et al. (2020). This is because the final value used in our manuscript included ~20
#additional standard runs that we added into the dataset to compute the final, mean, d15Nalpha 
#value. The methodology for estimating the value still stands exactly the same. 
#########################################################
#########################################################
#########################################################








