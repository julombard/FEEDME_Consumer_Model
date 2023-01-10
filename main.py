import Functions
import numpy as np
import pandas as pd
import multiprocessing
import Classes
import Parameters
#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
#Including New features for metapopulations modelling, designed by Massol F., Lion S. and bibi
#Damn efficient compared to previous try (on lab laptop 100 sites simulated over 16K iterations took 2'30)
#Parameters extracted from param files


beta = Parameters.beta  #Infectious contact rate
d = Parameters.d #Per capita natural death rate
gamma = Parameters.gamma  #Proportion of vertical transmission
v = Parameters.v # Virulence
m = Parameters.m # Dispersal propensity
e = Parameters.e # Ressource Encounter rate
p = Parameters.p # Profitability (conversion ressource -> reproduction)
theta = Parameters.theta # Medium supply

def RunModel(seed, param) :
    #Simulation parameters
    # This part changes parameters value for long autonomous runs
    m = param
    print('Current m parameter running', m)
    Classes.m = param
    Functions.m = param

    # Stocker les sorties dans un dictionnaire
    dico_densities_df = {}

    # Simulation parameters
    np.random.seed(seed) #Set seed for reproducibility
    nb_iterations = 0 #Store the number of interations to define times that are saved (later)
    sim_time = 0.0 # Simulation time (model time, not an iteration number)
    vectime = [0] # to keep t variable
    tmax = 10 # Ending time
    nbpatches = Parameters.nbpatches # Number of patches
    Taillepop = Parameters.Taillepop # Initial local population sizes

    #Set the landscape
    ListSites = Functions.SetMetapop(nbpatches, Taillepop)

    # Stocker les sorties dans un dictionnaire
    dico_densities_df = {}
    dico_propensities_df ={}

    #Event definition
    MediumRefreshing = Classes.Event(name='MediumRefreshing',propensity='theta', Schange='0', Ichange='0', Rchange ='1')
    ReproductionS = Classes.Event(name='Reproduction S',propensity='e * p * self.R * self.S ', Schange='1', Ichange='0',Rchange ='-1')
    ReproductionI = Classes.Event(name='Reproduction I',propensity='e * (1 - v)  * p *self.R * self.I', Schange='1', Ichange='0',Rchange ='-1') # Integrate the probability of vertical transmission separately
    DeathS = Classes.Event(name='Death S',propensity='d*self.S', Schange='-1', Ichange='0',Rchange ='0')
    DispersalS = Classes.Event(name='Dispersal S',propensity='m*self.S', Schange='-1', Ichange='0',Rchange ='0')
    DispersalI = Classes.Event(name='Dispersal I',propensity='m*self.I', Schange='0', Ichange='-1',Rchange ='0')
    Infection = Classes.Event(name='Infection',propensity='beta *self.S*self.I', Schange='-1', Ichange='1',Rchange ='0')
    DeathI = Classes.Event(name='Death I',propensity='d *self.I', Schange='0', Ichange='-1',Rchange ='0')

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
        if sim_time%5 == 0 : # X=5
            m = Parameters.m
        else :
            m = 0

        #Compute the propensities
        Propensities, Sum_propensities = Functions.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites
        #print("Propensities", Propensities)
        SumS, SumI, SumR = Functions.SumDensities(ListSites) # Get total densities of species

        # Defining the time increment according to Gillespie algorithm and First Reaction Method
        r1 = np.random.uniform(0, 1, 1)  # Draw random numbers
        a0 = sum(Sum_propensities)  # Overall propensity

        Tau = (1 / a0) * np.log(1 / r1)  # Time increment
        Probas = [i / a0 for i in Sum_propensities]  # Individual propensity
        NextReaction = list(np.random.multinomial(1, Probas))  # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        # print('The next reaction Index', NextReactionIndex)

        # Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
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
            site.effectifR += event.Rchange
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
            site.effectifR += event.Rchange


        # Update time
        sim_time += Tau[0]
        #print('Effectifs', SumS, SumI, SumR)
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
            dataprop.to_csv('Propensities_outputs_m'+str(param)+ "_" + str(seed) + '.csv')

    #Creating the time series dataframe
    datadensity = pd.DataFrame.from_dict(data=dico_densities_df)
    print("YOOOOOOOOOOOOOOUHOUUUUUUUUUUUUUUUUU", len(datadensity))
    VectimeDf = pd.DataFrame(data=vectime)
    datadensity.insert(0, "Time", VectimeDf, allow_duplicates=False)
    datadensity.to_csv('Sim_outputs_m'+str(param)+ "_" + str(seed) + '.csv')

################### MULTIPROCESSING PART ###########


# Paramètres de multiprocessing
list_seeds = [1,2,3,4,5,6]
list_params =[0.1,0.2,0.3,0.4,0.5]
nbsims = len(list_seeds)


# In the end, list_params must contain each parameters combinations
#Launch a batch of nbsim simulations
#Do not f*ck with that guy otherwise the whole computer can crash
if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count()-2 #Carefully count the number of available threads, and leave two of them alone (in order to do something else when simulating)
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb) # I don't really know what this is doing, probably creating some kind of logical space for thread assignation
    for j in range(len(list_params)) :
        for i in range(nbsims):
            pool.apply_async(RunModel, args=(list_seeds[i],list_params[j])) # Launch Nbsim simulation, beware because that makes you loose error messages
            #RunModel(list_seeds[i],list_params[j]) #Launch sims one by one, by makes the error messages reappear (useful for debugging)
    pool.close() # Ends something
    pool.join() # Do something