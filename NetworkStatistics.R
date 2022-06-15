##1 Getting Started
##Statnet has recently been updated and this code is based on the current latest version of R (4.1.0) and Statnet (4.0.0)
library(statnet)

##The sna and network packages are part of statnet, but might need to be loaded separately to attach data
#library(network)
#library(sna)

#make sure you're in the correct file directory or change the following command
load("SESWorkshopData.Rdata")



###########################################################
###########################################################
#2.1 Graph correlation and bivariate QAP
# Remember the Florentine families data ## good to use also with attribute data
data(florentine)
gplot(flobusiness) # Examine business ties
gplot(flomarriage) # Examine marriage ties

# Could the two be related?
par(mfrow=c(1,2))
coords<-plot(flobusiness,label=flobusiness%v%"vertex.names",label.cex=.5,pad=1)
title("Business Ties")
plot(flomarriage,label=flomarriage%v%"vertex.names",label.cex=.5,pad=1,coord=coords)
title("Marriage Ties")
par(mfrow=c(1,1))

# Let's try a graph correlation
gcor(flobusiness,flomarriage)

par(mfrow=c(1,2))
gplot(flobusiness[,],label=flobusiness%v%"vertex.names",label.cex=.5,pad=1,coord=coords)
title("Business Ties")
gplot(rmperm(flobusiness),label=flobusiness%v%"vertex.names",label.cex=.5,pad=1,coord=coords)
title("Permuted Business Ties")
par(mfrow=c(1,1))


#and use it
qt<-qaptest(list(flobusiness,flomarriage),gcor,g1=1,g2=2)
summary(qt) # Examine the results
plot(qt) # Plot the QAP distribution

# Testing against covariate effects
wealth<-sapply(flomarriage%v%"wealth",rep,network.size(flomarriage))
wealth
wealthdiff<-abs(outer(flomarriage%v%"wealth",flomarriage%v%"wealth","-"))
wealthdiff
qt1<-qaptest(list(flomarriage,wealth),gcor,g1=1,g2=2)
qt2<-qaptest(list(flomarriage,wealthdiff),gcor,g1=1,g2=2)
summary(qt1) # Do wealthy families have more ties?
summary(qt2) # Is there a wealth difference effect?

# For more information....
?qaptest
?gcor
?outer
?sapply
?rep

###########################################################
#2.2 Network Regression

# Combine the previous tests (takes a while to perform QAP test)
marriagefit<-netlm(flomarriage,list(flobusiness,wealth,wealthdiff))
summary(marriagefit) # Examine the results

# Another example: we begin by preparing the response variable. We will use the Cheyenne
# EMON in valued form, but need to recode the frequency data
data(emon)
Y<-as.sociomatrix(emon[[1]], "Frequency") # Extract frequencies
Y[Y>0]<-5-Y[Y>0] # Now, higher -> more frequent

# Extract some vertex attributes
crk<-emon[[1]]%v% "Command.Rank.Score" # Command rank
spon<-emon[[1]]%v%"Sponsorship" # Type of organization

# Form some predictor matrices (X variables)
Xcr<-sapply(crk,rep,length(crk)) # Column effects for command rank
Xcr
Xsp<-outer(spon,spon,"!=") # Dyadic effect for type difference
Xsp

# Fit the model (takes a while to perform QAP test)
cmfit<-netlm(Y,list(Xcr,Xsp))
summary(cmfit) # Examine the results


cmfitB<-netlm(Y,list(Xcr,Xsp),nullhyp="classical")
summary(cmfitB) # In this example, pretty similar

# For more information....
?outer
?netlm



###########################################################
###########################################################
##3 Simple univariate conditional uniform graph tests
# The cug.test function provides a simple interface for univariate CUG tests.
# Let's try testing some data on trade in complex manufactured goods to see if overall
# activity (density) is greater then would be expected from size alone.
cug.test(ctrade,gden) # Call cug.test with the gden (density) function

# Is there more reciprocity than density would suggest? Let's see.
cug.test(ctrade,grecip,cmode="edges") # Conditioning on edges, calling grecip

# Given biases in density and reciprocity, we might look to see if there is a
# bias in transitivity, too. This time, let's condition on all of the above.
cug.test(ctrade,gtrans,cmode="dyad") # Conditioning on dyad census

# We are not limited to simple commands. Let's try indegree centralization:
ct<-cug.test(ctrade,centralization,cmode="dyad",FUN.arg=list(FUN=degree, cmode="indegree")) 
# Note that we here must pass not only arguments to
# cug.test, but also to centralization and degree!

ct # Print the result
plot(ct) # Can also plot it!

#Using CUG to produce reference distributions

mod1<-netlm(ctrade,list(trade[1,,],trade[3,,],trade[4,,],trade[5,,]),nullhyp = "classical")

mod2<-netlm(ctrade,list(trade[1,,],trade[3,,],trade[4,,],trade[5,,]),nullhyp = "qap")

mod3<-netlm(ctrade,list(trade[1,,],trade[3,,],trade[4,,],trade[5,,]),nullhyp = "cugden")

