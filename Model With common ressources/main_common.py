import Functions_common
import numpy as np
import pandas as pd
import multiprocessing
import Parameters_common
import Classes_common
import itertools
from copy import deepcopy
#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
#Including New features for metapopulations modelling, designed by Massol F., Lion S. and bibi
#Damn efficient compared to previous try (on lab laptop 100 sites simulated over 16K iterations took 2'30)
#Parameters extracted from param files


#beta = Parameters.beta  #Infectious contact rate
#d = Parameters.d #Per capita natural death rate
#gamma = Parameters.gamma  #Proportion of vertical transmission
#v = Parameters.v # Virulence
#m = Parameters.m # Dispersal propensity
#e = Parameters.e # Ressource Encounter rate
#p = Parameters.p # Profitability (conversion ressource -> reproduction)
#theta = Parameters.theta # Medium supply

def RunModel(seed,IDsim, vecparam) :
    global beta, d, gamma, v, m, e, p, theta, Taillepop, nbpatches
    #Simulation parameters
    print(IDsim, "ID")
    print(vecparam,'Params')

    beta = vecparam[0]  # Infectious contact rate, because parasite property
    gamma = vecparam[1]  # Proportion of vertical transmission, because parasite property
    v = vecparam[2]  # Virulence, because parasite property
    m = vecparam[3]  # Dispersal propensity, because parasite property
    m_write = deepcopy(m) # Because the value of m can change during the run (due to closing, opening)
    theta = vecparam[4]  # Medium supply, Because why not

    print('Current m parameter running', beta, gamma, v, m, theta)
    #Update parameters needed in other files
    Classes_common.beta = beta

    Functions_common.beta = beta
    Classes_common.m =m
    Functions_common.m =m
    Classes_common.gamma = vecparam[1]  # Proportion of vertical transmission, because parasite property
    Classes_common.v = vecparam[2]  # Virulence, because parasite property
    Classes_common.m = vecparam[3]  # Dispersal propensity, because parasite property
    m_write = deepcopy(m)  # Because the value of m can change during the run (due to closing, opening)
    Classes_common.theta = vecparam[4]  # Medium supply, Because why not

    #Initialize Medium
    Medium_ressource = Parameters_common.Initial_Medium_Quality


    # Stocker les sorties dans un dictionnaire
    dico_densities_df = {}

    # Simulation parameters
    np.random.seed(seed) #Set seed for reproducibility
    nb_iterations = 0 #Store the number of interations to define times that are saved (later)
    sim_time = 0 # Simulation time (model time, not an iteration number)
    vectime = [0] # to keep t variable
    tmax = 300 # Ending time
    nbpatches = Parameters_common.nbpatches # Number of patches
    Taillepop = Parameters_common.Taillepop # Initial local population sizes

    #Set the landscape
    ListSites = Functions_common.SetMetapop(nbpatches, Taillepop)

    # Stocker les sorties dans un dictionnaire
    dico_densities_df = {}
    dico_propensities_df ={}

    #Event definition
    MediumRefreshing = Classes_common.Event(name='MediumRefreshing',propensity='theta', Schange='0', Ichange='0', Rchange ='1')
    ReproductionS = Classes_common.Event(name='Reproduction S',propensity='e * p * self.R * self.S + e * (1 - v)* (1-gamma) * p *self.R * self.I ', Schange='1', Ichange='0',Rchange ='-1')
    ReproductionI = Classes_common.Event(name='Reproduction I',propensity='e * (1 - v)* gamma * p *self.R * self.I', Schange='0', Ichange='1',Rchange ='-1') # Integrate the probability of vertical transmission separately
    DeathS = Classes_common.Event(name='Death S',propensity='d*self.S', Schange='-1', Ichange='0',Rchange ='0')
    DispersalS = Classes_common.Event(name='Dispersal S',propensity='m*self.S', Schange='-1', Ichange='0',Rchange ='0')
    DispersalI = Classes_common.Event(name='Dispersal I',propensity='m*self.I', Schange='0', Ichange='-1',Rchange ='0')
    Infection = Classes_common.Event(name='Infection',propensity='beta *self.S*self.I', Schange='-1', Ichange='1',Rchange ='0')
    DeathI = Classes_common.Event(name='Death I',propensity='d *self.I', Schange='0', Ichange='-1',Rchange ='0')

    #Event vector, cause tidying up things is always nice
    Events = [MediumRefreshing,ReproductionS,ReproductionI,DeathS,DispersalS,DispersalI,Infection,DeathI]

    #Initializing outputs storage
    Densities_out = [] # Collect densities outputs
    Propensities_out =[] # Collect propensities outputs
    IsTrackPropensites = False #Set false/ true to not track/track propensities
    for i in range(len(Events)):
        Propensities_out.append([0]) #Store times series of propensities

    #Get initial lists and values in outputs
    #We want to get one list per Sx(t) and Ix(t) to store them easily in a dataframe at the end
    for index, i in enumerate(ListSites):
        dico_densities_df[f"S{index}"] = [i.effectifS]
        dico_densities_df[f"I{index}"] = [i.effectifI]

    #Model main Loop
    while sim_time < tmax :
        print('We have currently passed', sim_time,'time in the simulation') # Kind of a loading bar
        vectime.append(sim_time) #Update time vector

        #Defnie the landscape configuration at given time -> Are the gates opened ?
        #Let's say that each X time, gates are openened during one time step
        if round(sim_time,1)%20 == 0.0 : # X=5
            if round(sim_time,0)==0 :
                m = 0
                Classes_common.m = 0
            else :
                m = vecparam[3]
                Classes_common.m = vecparam[3]
        else :
            m = 0
            Classes_common.m = 0
        #print('VOICI M', m, Classes.m)
        #Compute the propensities
        Propensities, Sum_propensities = Functions_common.GetPropensites(ListSites, Events, Medium_ressource) # Get a vector of propensities ordered by event and by sites
        #print("Propensities", Propensities)
        #print("Sum Propensities", Sum_propensities)
        #print("Events", Events)
        SumS, SumI = Functions_common.SumDensities(ListSites) # Get total densities of species

        Tau_candidates = {}
        events_indexes = {}
        # Defining the time increment according to Gillespie algorithm and First Reaction Method
        for i in range(len(Events)) :
            name_event = Events[i].name
            r1 = np.random.uniform(0, 1, 1)  # Draw random numbers
            aj = Sum_propensities[i]  # Propensity of the event
            if aj==0:
                Candidate =[0]
            else :
                Candidate = (1 / aj) * np.log(1 / r1)  # Time increment candidate
                Tau_candidates[name_event] = Candidate[0] # Store the event and its tentative reaction
                events_indexes[name_event] = i

        Next_event = min(Tau_candidates, key=Tau_candidates.get) # Gives the name of the next event
        Tau = Tau_candidates[Next_event] # Gives the time increment
        #print(Tau)
        #Probas = [i / a0 for i in Sum_propensities]  # Individual propensity
        #NextReaction = list(np.random.multinomial(1, Probas))  # List to get index easily
        NextReactionIndex = events_indexes[Next_event]
        #print(NextReactionIndex)
        #print('The next reaction ', Next_event)

        # Determine where the reaction is going to happen
        #print('Next Reaction', Next_event)
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        #print(props, "les props")
        sumprop = sum(props)
        proba_site = [i / sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1) # Get the index of the patch in which the event takes place
        # print('Next places', NextPlaceIndex)

        # This part apply the effect of events in site populations
        event = Events[NextReactionIndex]
        site = ListSites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS += event.Schange
            site.effectifI += event.Ichange
            #nbmigrants = 1
            # Here we distribute successful migrants among neighboring sites

            # Determine which site will receive the dispersing individual
            # Get Index of the current site, and choose with bernouilli trial if it goes forward or backward
            receiving_sites = [NextPlaceIndex -1,NextPlaceIndex + 1]
            for i in receiving_sites :
                if i < 0 or i > nbpatches-1 : # If we should allow dispersal to a site that does not exist
                    receiving_sites.remove(i) # Delete the coordinate
            site_destination = np.random.choice(receiving_sites)  # destination is a site object
            # print('The destination is', site_destination)

            # add individual to destination
            if abs(event.Schange) > 0:  # if S are dispersers
                ListSites[site_destination].effectifS += 1
            elif abs(event.Ichange) > 0:
                ListSites[site_destination].effectifI += 1
            else:
                pass
                # print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            site.effectifS += event.Schange
            site.effectifI += event.Ichange
            Medium_ressource += event.Rchange


        # Update time
        sim_time += Tau
        print('Effectifs', SumS, SumI)
        print('Medium', Medium_ressource)
        # print('time increment', Tau)
        # print('Le temps passe si vite',sim_time)

        # 1. Densities
        indexlist = 0
        for index, i in enumerate(ListSites):
            if i.effectifS < 0:  # Avoid negative population in the "big fat brute" way
                i.effectifS = 0
            if nb_iterations % 10 == 0:
                #print('Saving...')
                dico_densities_df[f"S{index}"].append(i.effectifS)
            indexlist += 1
            if i.effectifI < 0:
                i.effectifI = 0
            if nb_iterations % 10 == 0:
                dico_densities_df[f"I{index}"].append(i.effectifI)
            indexlist += 1
        # 2. Propensities
        if IsTrackPropensites == True:
            # Propensities of each event in a list sorted by event
            for index, propensitiy in enumerate(Sum_propensities):
                Propensities_out[index].append(propensitiy)

        # Structuring outputs to get a .csv file event if the loop has broken
        if IsTrackPropensites == True:  # if we track propensities
            # Creating propensities dataframe
            dataprop = pd.DataFrame(columns=['t'])
            for event in Events:
                EventName = event.name
                dataprop[EventName] = []
            # Filling the dataframe
            for index, colname in enumerate(dataprop):
                if index == 0:
                    dataprop[colname] = vectime
                else:
                    dataprop[colname] = Propensities_out[index - 1]
            # Saving into .csv file
            dataprop.to_csv('Propensities_outputs_' + str(seed) + '.csv')

    #Creating the time series dataframe
    datadensity = pd.DataFrame.from_dict(data=dico_densities_df)
    #print("YOOOOOOOOOOOOOOUHOUUUUUUUUUUUUUUUUU", len(datadensity))
    VectimeDf = pd.DataFrame(data=vectime)
    datadensity.insert(0, "Time", VectimeDf, allow_duplicates=False)
    #Complete the DF with parameters values
    datadensity.insert(0, "beta", beta)
    datadensity.insert(0, "gamma", gamma)
    datadensity.insert(0, "v", v)
    datadensity.insert(0, "m", m_write)
    datadensity.insert(0, "theta", theta)


    datadensity.to_csv('Sim_outputs_'+IDsim+ "_" + str(seed) + '.csv')

################### MULTIPROCESSING PART ###########


# Paramètres de multiprocessing
list_seeds = [1,2,3,4,5]
nbsims = len(list_seeds)

#Create a list of parameters
#Sorting matters ! Follow this one : beta, gamma, v , m , theta
dict_param = {}

# Create parameters combinations
beta_levels = [0.005,0.007,0.009]
gamma_levels = [0.5]
v_levels =[0.3]
m_levels = [0.8]
theta_levels = [20]

Superlist = [beta_levels,gamma_levels,v_levels,m_levels,theta_levels] # List of list of levels, listception is my leitmotiv
my_combinations = list(itertools.product(*Superlist)) # Compute the set of combination between lists and return it as a set of tuples


dico_param = {} # Enable an empty dictionary
for i in range(len(my_combinations)) :
    dict_param[str(i)] = list(my_combinations[i]) # Assign each tuple (converted into list of param) to a specific ID
print(dico_param)


# In the end, list_params must contain each parameters combinations
#Launch a batch of nbsim simulations
#Do not f*ck with that guy otherwise the whole computer can crash
if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count()-2 #Carefully count the number of available threads, and leave two of them alone (in order to do something else when simulating)
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb) # I don't really know what this is doing, probably creating some kind of logical space for thread assignation
    for key in dict_param.keys():
        for i in range(nbsims):
            pool.apply_async(RunModel, args=(list_seeds[i],key,dict_param[key])) # Launch Nbsim simulation, beware because that makes you loose error messages
            #RunModel(list_seeds[i],key,dict_param[key]) #Launch sims one by one, by makes the error messages reappear (useful for debugging)
    pool.close() # Ends something
    pool.join() # Do something