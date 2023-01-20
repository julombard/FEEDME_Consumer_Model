import Parameters_common
global beta, d, gamma, v ,m, e, p, theta, Taillepop, nbpatches

beta = Parameters_common.beta  #Infectious contact rate
d = Parameters_common.d #Per capita natural death rate
gamma = Parameters_common.gamma  #Proportion of vertical transmission
v = Parameters_common.v # Virulence
m = Parameters_common.m # Dispersal propensity
e = Parameters_common.e # Ressource Encounter rate
p = Parameters_common.p # Profitability (conversion ressource -> reproduction)
theta = Parameters_common.theta # Medium supply

class Site(): #Site object containing (non explicit) individuals
    def __init__(self,effectifS,effectifI): #First try with arg way to implement feature unsure
        self.effectifS = effectifS #S density
        self.effectifI = effectifI #I density


class Event():
    def __init__(self,name, propensity, Schange, Ichange, Rchange):
        self.name = name # Event name in letter and not in memory address, handful to identify what's happening
        self.S = 0 #Has to take density values to understand the maths
        self.I = 0
        self.R = 0
        self.formula = propensity # The unique formule (str) given by model construction
        self.propensity = eval(self.formula)#Convert string in maths instruction, very useful to externalise model building
        self.Ichange = eval(Ichange) # State Change due to event, Typically -1, 0, 1
        self.Schange = eval(Schange)
        self.Rchange = eval(Rchange)


    def UpdatePropensity(self, S, I, R): # Class method to compute propensities without creating new objects
        #print(self.name)
        self.S = S
        self.I = I
        self.R = R
        self.propensity = eval(self.formula) # Changes propensity values while keeping formula untouched
        return self.propensity