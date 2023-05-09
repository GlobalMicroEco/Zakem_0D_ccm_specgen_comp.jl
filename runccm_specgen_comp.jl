#5/9/23: Send to Lee. Tinkered with some parameter values (randomized vmax, for example). Replaced \rho with v. Changed name of .nc output to "out..."

#8/10/22: Clean up a bit for Liang

#11/19/20: Upload to Github

#Apr 19, 2020 by ejz
#For input into function runccm() within ccm.jl
#Parameters are assembled into type Prms (using "struct"), which is defined in ccm.jl 
#runccm then requires prms as input, i.e. a specific instance of the type Prms

#output file name (to which date and .nc will be added): 

fsave = "out_specgen_comp"

##############################################
#Parameters

using Statistics
using StatsBase

T = 365*1 #days -- run time
nd = 10 #n OM pools
nb = nd*10 #n microbial pops
Pt = 1 #mmol C/m3 total OM supply

nrec = 500 # of timepoints to record (will be nrec + 1 with t=0)

#Consumption matrix

#####################################################################################
# #Multiple specialists for each D
# 
# #spec:
# #CM = [i == j for i = 1:nd, j=1:nb] #diagonal only (specialists)
# CMp = [i == j for i = 1:nd, j=1:nd] #diagonal only (specialists)
# #specgen:
# #CM = [i >= j for i = 1:nd, j=1:nb] #fills the lower diagonal (spec to gen gradient)
# #random:
# #CM = rand(nd,nb) .> 0.8 #randomly fills above specified value
# 
# #now put multiple specialist CMs together:
# CM = hcat([CMp for i = 1:(nb ÷ nd)]...) #the ... turns the elements of the list into arguments for the function

#####################################################################################
#Make two matrices and combine into one: the first with just specialists (to make sure they are there), and the second with a range of generalist ability 

#FIRST HALF: Specialists only
CMone = zeros(nd,nb÷2)
#CM = [i == j for i = 1:nd, j=1:nb÷2] #diagonal only (specialists)
#CM[i == j for i = 1:nd, j=1:nb÷2]
for j = 1:nd
    CMone[j,j] = 1
end

##############################################################################
#SECOND HALF: random
##############################################################################
#METHOD 1: assign nup first

#1.1. nup = # of substrates consumed by each pop
#lin space:
#nupr = rand(1:nd,nb÷2) # n substrates consumed by each, nb long

#log space: 
xhi = log10(nd)
xlo = log10(1)
nupr = 10 .^ ((rand(nb÷2).-0.5).*(xhi-xlo).+mean([xhi xlo]))
nupr = round.(Int,nupr) #round to integers:

#nupr[:] .= nd/10 #to test only impact of # consumers 
#now for $ieach B, assign its substrates randomly: rand(1:nd,nup[j]) but UNIQUELY
#uniquely using sample: sample(1:nd,nb,replace=false) from StatsBase

#food=[sample(1:nd,nupr[j], replace = false) for j=1:nb] #an array of indices

#1.2. assign a likelihood of consumption for each substrate
w = [1/nd:1/nd:1; ] #ordered
wv = StatsBase.ProbabilityWeights(w)
#sample(1:nd,wv,10,replace=false) #over range 1:nd, with weight wv, give 10 samples, all unique

food=[sample(1:nd,wv,nupr[j], replace = false) for j=1:nb÷2] #an array of indices

CMtwo = zeros(nd,nb÷2)
for j = 1:nb÷2
    CMtwo[food[j],j] .= 1
end

##############################################################################
#METHOD 2: assign cms (# of consumers) first

# #1.1. nup = # of substrates consumed by each pop
# #lin space:
# #cmsr = rand(1:nb÷2,nd)
# 
# #log sapce:
# xhi = log10(nb÷2)
# xlo = log10(1)
# cmsr = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
# cmsr = round.(Int,cmsr) #round to integers
# 
# #1.2. assign a likelihood of generality for each population (how many substrates does it eat?)
# #w = [1/nd:1/nd:1; ] #ordered
# w = [1:nb÷2; ] #ordered -- same as above
# #w[:] .= 1
# #w = w .^ 4 #if this is increased, then things are a bit more evenly distributed (makes some MUCH less likely)
# wv = StatsBase.ProbabilityWeights(w)
# 
# #1.3. each is a list of the consumers for each substrate
# consumers=[sample(1:nb÷2,wv,cmsr[j], replace = false) for j=1:nd] #an array of indices
# 
# CMtwo = zeros(nd,nb÷2)
# for j = 1:nd
#     CMtwo[j,consumers[j]] .= 1
# end

#####################################
#concatenate

CM = hcat(CMone,CMtwo)

#####################################

#convert to boolean:
CM = convert(Array{Bool}, CM .== 1)

nup = sum(CM,dims=1)[1,:]

#penalty for generalists:
pen = 1 ./ nup

#presence absence for each population (sets overall probability):
PA = ones(nb) #length nb
#PA = rand(nb)

#supply weight for each pool (sets overall probability)
#Cw = rand(nd) #length nd 
Cw = ones(nd) #equal supply

#Dilution rate for additional sink (for D and B)
D = 0 #1/d 

###################################################
#SUBSTRATE TRAITS ("averages")

#1. Max uptake rate (1/d)
#Random assignment:
xhi = log10(1e1)
xlo = log10(1e-1) 
vmax_i = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
#Ordered assignment:
#vmax_i = collect(10 .^ range(-2,stop=1,length=nd)) #like logspace

#2. affinity (m3/mmol/d) = pmax/km
# xhi = log10(100)
# xlo = log10(1)
# #aff = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
# #km_i = vmax_i./aff #half-sat constant (mmol/m3)
# #constant affinity:
km_i = vmax_i./10 #if km_ij is set by tradeoff below, this represents a "default" affinity in the case of no tradeoff

#yield
y_i = rand(nd)*0.5  
#y_i = ones(nd)*0.3 #mol B/mol OM

####################################################
#MICROBIAL TRAITS

#Now assign each population some variation of the average traits (pmax, affinity, and y)

#Tradeoff between growth rate and affinity:
#Consider a fraction of the proteome devoted to growth vs. affinity

#Constant: F_g + F_a 
F_g = rand(nb)
F_a = 1. .- F_g

#now assign traits based on F. Easiest way: pmax scales with F 
vmax_ij = zeros(nd,nb)
for i = 1:nd
    vmax_ij[i,:] = vmax_i[i]*F_g./0.5.*CM[i,:]
end

#How does F_a relate to km? through affinity 
aff = zeros(nd,nb)
km_ij = zeros(nd,nb)
for i = 1:nd
    #old aff function used for previous runs (old Vmax vs affinity)
    #aff[i,:] = exp.(2*F_b/F_bA).*CM[i,:]
    #new: scales same as y and Vmax
    aff[i,:] = 10*F_a./0.5.*CM[i,:]
    km_ij[i,:] = vmax_ij[i,:]./aff[i,:]
end

# #Constant yield 
y_ij = broadcast(*,y_i,CM)

#Mortality 
#mq = rand(nb)*1. #quadratic mort: m3/mmol/d
#mlin = rand(nb)*1e-2 #linear mort: 1/d
mq = ones(nb) #quadratic mort: m3/mmol/d
mlin = ones(nb)*1e-2 #linear mort: 1/d

##############################################
#Initial conditions:

#random:
#xhi = log10(10)
#xlo = log10(0.1)
#dIC = 10 .^ ((rand(nd).-0.5).*(xhi-xlo).+mean([xhi xlo]))
dIC = ones(nd)*0.1

bIC = ones(nb)*0.1

##############################################
#load model including Prms struct:
include("ccm_comp.jl")

##############################################
#assemble into type Prms (defined in ccm.jl)
params = Prms(fsave,T,Pt,nd,nb,CM,pen,PA,Cw,D,vmax_i,km_i,y_i,vmax_ij,km_ij,y_ij,mq,mlin,dIC,bIC,nrec);

##############################################
#run model and save into fsave_*date*.nc:
d, b, timet, Ctot, v, prms = runccm(params);




