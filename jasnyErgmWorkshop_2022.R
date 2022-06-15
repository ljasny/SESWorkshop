##################
## 0. Load Packages
library(statnet)

##################
## 1. Network data structure
## first, set your working directory to where the edgelist file is stored
## setwd("C:/Users/lorien/etc")

#load and plot data
edgelistdata<-read.csv("edgeList.csv",header=T,stringsAsFactors=F)
edgelistnet<-network(edgelistdata,matrix.type="edgelist")
edgelistnet
plot(edgelistnet,displaylabels=T) #what's missing?

#weights
plot(edgelistnet,displaylabels=T,edge.lwd=10*edgelistnet%e%"Weight")
edgelistnet[,]
as.sociomatrix.sna(edgelistnet,"Weight")

#attribute data
edgelistnet%v%"gender"<-c("M","F","F","M","M","M")
edgelistnet%v%"age"<-1:4
edgelistnet%v%"age" #what happened?

##################
##2.1 Fitting A Basic ERG Model

data(package="ergm") # tells us the datasets in our packages
data(florentine) # loads flomarriage & flobusiness data
flomarriage # Let's look at the flomarriage data
plot(flomarriage, displaylabels = TRUE) # Let's view the flomarriage network

flomodel.01 <- ergm(flomarriage ~ edges) # fit model
plogis(-1.6094) #compare to gden(flomarriage)
gden(flomarriage) #the edges term controls for density

flomodel.01 # look at the model
summary(flomodel.01) # look in more depth

#Let's add a term often thought to be a measure of "clustering" -- the number of completed triangles
set.seed(1664) #Seeding the RNG for demo purposes; otherwise, your output will vary
flomodel.02 <- ergm(flomarriage ~ edges + triangle)
summary(flomodel.02)

#Let's take a closer look at the ergm object itself:
class(flomodel.02) # this has the class ergm
names(flomodel.02) # let's look straight at the ERGM obj.
flomodel.02$formula # the $ allows you to pull an element out from
coef(flomodel.02) 

# Must check MCMC chains
mcmc.diagnostics(flomodel.02)

# Must also check model fit
gof.02 = gof(flomodel.02)
plot(gof.02)

# additional goodness of fit
gof.02.2 = gof(flomodel.02~triadcensus)
plot(gof.02.2)

# Including node attributes in the model
wealth <- flomarriage %v% 'wealth'
wealth 
plot(flomarriage, vertex.cex = wealth/20, displaylabels = TRUE) 
# network plot with vertex size proportional to wealth

#We can test whether edge probabilities are a function of wealth:
flomodel.03 <- ergm(flomarriage ~ edges + nodecov('wealth'))
summary(flomodel.03)

# Can test edge probabilities as a function of edge attributes too
wealth.diff <- abs(outer(wealth, wealth, FUN="-"))
wealth.diff

flomodel.04 <- ergm(flomarriage ~ edges + edgecov(wealth.diff))
summary(flomodel.04)

flomodel.04.2 <- ergm(flomarriage ~ edges + absdiff('wealth'))
summary(flomodel.04.2)

# How do you find the available model terms?
?ergm.terms

#####################################################################
#####################################################################
##2.2 Examples with Directed Networks and More Attributes

data(samplk) # Let's try a model or two on
samplk3  # directed data: Sampson's Monks (liking)
plot(samplk3, displaylabels = TRUE)
sampmodel.01 <- ergm(samplk3 ~ edges + mutual)# Is there a statistically significant
summary(sampmodel.01) # tendency for ties to be reciprocated ("mutuality")?

data(faux.mesa.high) # Let's try a larger network
mesa <- faux.mesa.high
mesa
plot(mesa, vertex.col = 'Grade')
legend('bottomleft', fill=7:12, legend=paste('Grade', 7:12), cex=0.75)

fauxmodel.01 <- ergm(mesa ~ 
                       edges + 
                       nodefactor ('Grade') +
                       nodefactor ('Race') +
                       nodematch('Grade',diff=T) + 
                       nodematch('Race',diff=T)
)
summary(fauxmodel.01) #Note that two of the coefficients are estimated as -Inf (the nodematch coefficients for race Black and Other). Why is this?
table(mesa %v% "Race")
mixingmatrix(mesa, "Race")

fauxmodel.02 <- ergm(mesa ~ 
                       edges +
                       nodefactor('Grade') +
                       nodematch('Grade',diff=T) + 
                       nodefactor('Race') +
                       nodematch('Race',diff=T,levels=c(2,3,5))
)

summary(fauxmodel.02)

fauxmodel.03 <- ergm(mesa ~ 
                       edges +
                       triangle +
                       nodefactor('Grade') +
                       nodematch('Grade',diff=T) + 
                       nodefactor('Race') +
                       nodematch('Race',diff=T,levels=c(2,3,5))
)

##not run for time
# fauxmodel.04 <- ergm(mesa ~ 
#                        gwesp(0.5, fixed=T) +
#                        triangle +
#                        nodefactor('Grade') +
#                        nodematch('Grade',diff=T) + 
#                        nodefactor('Race') +
#                        nodematch('Race',diff=T,levels=c(2,3,5))
# )
# 
# summary(fauxmodel.04)

summary(mesa~gwesp(0.5))

#########################################################
#########################################################
##3 Diagnostics: Troubleshooting and Checking for Model Degeneracy
fit <- ergm(flobusiness ~ edges + degree(1), 
            control=control.ergm(MCMC.interval=1, MCMC.burnin=1, seed=1))
mcmc.diagnostics(fit)

fit <- ergm(flobusiness ~ edges + degree(1),
            control=control.ergm(MCMC.interval=100000, MCMC.burnin=100000, seed=1))
mcmc.diagnostics(fit)

data('faux.magnolia.high')
?faux.magnolia.high
magnolia <- faux.magnolia.high
plot(magnolia)

fit <- ergm(magnolia ~ edges + triangle, 
            control = control.ergm(seed=1))

fit <- ergm(magnolia ~ edges + triangle,
            control = control.ergm(seed=1,
                                   MCMC.burnin=100000,
                                   MCMC.interval=500000))

## not run for time
#fit <- ergm(magnolia ~ edges + gwesp(0.5, fixed = TRUE),
#            control = control.ergm(seed = 1))
#mcmc.diagnostics(fit)
#########################################################
#########################################################
##4 Degree Constraints
##fixing maximum outdegree
fit.2 <- ergm(magnolia ~ edges + gwesp(0.5, fixed = TRUE),
            control = control.ergm(seed = 1),constraints=~bd(maxout=3))

max(degree(magnolia))
fit.2 <- ergm(magnolia ~ edges + gwesp(0.5, fixed = TRUE),
              control = control.ergm(seed = 1),constraints=~bd(maxout=16))
mcmc.diagnostics(fit.2)

##setting who can and cannot send a tie (for example, if some individuals were interviewed and some were not)
isRespondent<-sample(c(TRUE,FALSE),network.size(magnolia),replace=T)

resp_maxOut<-vector(length=network.size(magnolia))
resp_maxOut[]<-NA
resp_maxOut[!isRespondent]<-0
table(resp_maxOut)

fit.3 <- ergm(magnolia ~ edges + gwesp(0.5, fixed = TRUE),
              control = control.ergm(seed = 1),constraints=~bd(maxout=resp_maxOut))

#########################################################
#########################################################
##5 Bipartite Data

#look at the network (it's a matrix)
gplot(clanByPatch,displaylabels = T,usearrows = F,vertex.col=c(rep("blue",8),rep("green",14)))

#make it a network object
CPNet<-as.network(clanByPatch,bipartite=T)

#let's make some random attributes
#we'll start with how many people are in each clan
people<-sample(10:100,CPNet%n%"bipartite",replace=T)
people

CPNet%v%"people"<-people
CPNet%v%"people" ##what happened???

CPNet%v%"people"<-c(people,rep(NA,CPNet%n%"n"-CPNet%n%"bipartite"))
CPNet%v%"people"

##do the same for attributes on the 2nd mode
CPNet%v%"public"<-c(rep(NA,CPNet%n%"bipartite"),sample(0:1,(CPNet%n%"n"-CPNet%n%"bipartite"),replace=T))

CPMod1<-ergm(CPNet~edges+b1cov("people")+b2factor("public"))
CPMod2<-ergm(CPNet~edges+b1cov("people")
             +b2factor("public")
             +gwb2dsp(.5,fixed=T)
             +gwb1degree(.5,fixed=T))
#########################################################
#########################################################
##6 Multilevel models

plot(madagascar,vertex.col=ifelse(madagascar%v%"social","blue","green"))

##extract social-only network
madagascar_social<-network(madagascar[,][1:8,1:8])

#look at the density
mad_soc_mod1<-ergm(madagascar_social~edges)
summary(mad_soc_mod1)

##density for social nodes should should be identical to density of the social only network
mad_mod1<-ergm(madagascar~nodemix("social"))
summary(mad_mod1)

##what about switching out the different terms so it reads cleaner
mad_mod1<-ergm(madagascar~nodemix("social",levels2=c(1,2,4)))
summary(mad_mod1)

##next,we'll think about writing a term to study socio-ecological triads with either 1 or 2 social nodes (and 2 or 1 ecological)
#we'll use the recently released 'Filter' and 'Sum' terms for ERGM under Statnet 4
#these make writing ERGM terms more like it's own language so anything you can think of (pretty much) can be written
#to do this, it can be easier (at least it is for me) to work with a tiny demo network to make sure my code for a new term is right

plot(demo,vertex.col=ifelse(demo%v%"Social","blue","green"),vertex.cex=2,displaylabels=T)

##How many triangles are there?
summary(demo~triangles)

##How many ecological triangles are there?
summary(demo~F(~triangles, 
               ~nodematch(~Social!=1)))

##Check to make sure
demo[8,7]
demo[8,7]<-0
summary(demo~F(~triangles, 
               ~nodematch(~Social!=1)))
demo[8,7]<-1 ##reset the network




summary(demo~Sum(list(+1~nodemix(~Social,levels2=2),
                +1~nodemix(~Social,levels2=3),
                +1~nodemix(~Social,levels2=4)),label="SSE_triangle"))

summary(demo~F(~triangles,~Sum(list(+1~nodemix(~Social,levels2=2),
                      +1~nodemix(~Social,levels2=3),
                      +1~nodemix(~Social,levels2=4)),label="SSE_triangle")))

demo[5,6]
demo[5,6]<-0
summary(demo~F(~triangles,~Sum(list(+1~nodemix(~Social,levels2=2),
                                    +1~nodemix(~Social,levels2=3),
                                    +1~nodemix(~Social,levels2=4)),label="SSE_triangle")))
demo[5,6]<-1


mod_soceco<-ergm(madagascar~nodemix("social",levels2=c(1,2,4))+F(~triangles,~Sum(list(+1~nodemix(~social,levels2=2),
                                                                                      +1~nodemix(~social,levels2=3),
                                                                                      +1~nodemix(~social,levels2=4)),label="SSE_triangle")))

summary(mod_soceco)

#########################################################
#########################################################
##7 Network Simulation
flomodel.01
names(flomodel.01)
flomodel.01$network
plot(flomodel.01$network, displaylabels = TRUE)

flomodel.01.sim <- simulate(flomodel.01, nsim=10)
class(flomodel.01.sim)
length(flomodel.01.sim)
# Take a look at the first simulated network
flomodel.01.sim[[1]]
plot(flomodel.01.sim[[1]], displaylabels = TRUE)

# Extract arbitrary statistics from the simulations
flo.sim.triad.census <- lapply(flomodel.01.sim, triad.census, mode = "graph")
flo.sim.triad.census
# Take those list elements and put them into a matrix
tris = do.call(rbind, flo.sim.triad.census)
tris

#compare simulations and empirical
par(mfrow=c(2,2))
for(i in 1:4){
  hist(tris[,i], main = paste(colnames(tris)[i], "ties"))
  # Or density curves:
  #  plot(density(tris[,i]), main = paste(colnames(tris)[i], "ties"))
  abline(v=triad.census(nflo,mode="graph")[i],col="red")
}

# What do the same statistics look like for simulations from our other models?
# Write a function to do the same analyses using any model
simTriadsCompare = function(model, 
                            nsims = 10, 
                            directed = model$network %n% "directed")
  # This is a function to compare tridic censuses (censi?) on networks simulated
  # from an ERGM versus the empirical network.
  # Arguments:
  # model is an exponential random graph model (class = "ergm")
  # nsims is the number of networks to simulate from the ERGM
  # Returns:
  # Density plots for each of the possible triadic arrangements with a 
  # vertical line drawn at the observed count
  
{
  # Run the simulations
  sims = simulate(model, nsim = nsims)
  
  # Extract graph type to pass to triad.census
  graphType = ifelse(directed, "digraph", "graph")
  
  # Perform triad census on each simulated graph and rbind into a matrix
  tris = lapply(sims, triad.census, mode = graphType)
  tris = do.call(rbind, tris)
  
  # Extract the empirical triadic census values
  empTris = triad.census(model$network, mode = graphType)
  
  # Output density plots
  par(mfrow=c(2,2))
  for(i in 1:ncol(tris)){
    plot(density(tris[,i]), main = paste(colnames(tris)[i], "ties"))
    abline(v=empTris[i],col="red")
  }
}

# Now call that function with each of our models
simTriadsCompare(model = flomodel.01)
simTriadsCompare(model = flomodel.02)
# What was that second model again?
flomodel.02$formula
simTriadsCompare(model = flomodel.03)
simTriadsCompare(model = flomodel.04)

