import numpy as np
from copy import  deepcopy
import Classes


def DoDirectMethod(Propensities, Sum_propensities, exactsteps, events, sites):
    rho = 0  # Not very convenient to put it here but it has to....
    print(Propensities, 'HELLOOOOOO')
    for i in range(exactsteps):
        r1 = np.random.uniform(0,1,1) #Draw random numbers
        a0 = sum(Sum_propensities) # Overall propensity

        Tau = (1/a0) * np.log(1/r1) #Time increment
        Probas = [i / a0 for i in Sum_propensities] #Individual propensity
        NextReaction = list(np.random.multinomial(1, Probas)) # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        #print('The next reaction Index', NextReactionIndex)

        #Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        sumprop = sum(props)
        proba_site = [i/sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1)
        #print('Next places', NextPlaceIndex)


        # This part apply the effect of events in site populations
        event = events[NextReactionIndex]
        site = sites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
            nbmigrants = 1
            # Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
            SuccessfulMigrants = 0

            roll4urlife = np.random.uniform(0, 1, 1)
            if roll4urlife > rho: SuccessfulMigrants += 1
            # Here we distribute successful migrants among neighboring sites
            # This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
            if SuccessfulMigrants == 1:
                # Determine which site will receive the dispersing individual
                receiving_sites = deepcopy(sites)  # Create a working copy of sites
                # print('receivers', receiving_sites)
                del receiving_sites[NextPlaceIndex]  # removing departure site from the copy
                # print('receivers post suppression', receiving_sites)
                site_destination = np.random.choice(receiving_sites)  # destination is a site object
                # print('The destination is', site_destination)

                # add individual to destination
                if abs(event.Schange) > 0:  # if S are dispersers
                    site_destination.effectifS += 1
                elif abs(event.Ichange) > 0:
                    site_destination.effectifI += 1
                else:
                    pass
                    #print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
        return Tau

def SetMetapop(nbpatches, taillepop): #Creates sites objects containing populations
    ListSites=[] # List that will contain all sites
    for i in range(nbpatches): # Creates sites, the 1st will always contain one infected and the other 0
        if i == 0:
            newsite = Classes.Site(effectifS=taillepop-10, effectifI=10, effectifR=100)
            ListSites.append(newsite)
        else:
            newsite = Classes.Site(effectifS=0, effectifI=0, effectifR=100)
            ListSites.append(newsite)
    print(ListSites)
    return ListSites

def GetPropensites (Sites, Events): # Compute the propensities
    Propensities = []
    for i in Events: # For each event
        PropEvent =[]
        for j in Sites : # And each site
            S, I, R = j.effectifS, j.effectifI, j.effectifR # Get the xi
            Prop = i.UpdatePropensity(S,I,R) #Compute propensity
            PropEvent.append(Prop)
        Propensities.append(PropEvent)
    sumpropensities = []
    for i in Propensities :
        sumpropensities.append(sum(i))
    return Propensities, sumpropensities

def GetPropensites_Per_sites (Sites, Events): # Compute the propensities
    Propensities_per_sites = []
    for i in Sites: # For each event
        PropSite =[]
        S, I = i.effectifS, i.effectifI # get the xi
        for j in Events : # And each site
            Prop = j.UpdatePropensity(S,I) #Compute propensity
            PropSite.append(Prop)
        Propensities_per_sites.append(PropSite)
    sumpropensities = []
    for i in Propensities_per_sites :
        sumpropensities.append(sum(i))
    return Propensities_per_sites, sumpropensities

def SumDensities(Sites) : # Count all S and I individual
    SumS = 0
    SumI = 0
    SumR = 0
    for i in Sites:
        SumS += i.effectifS
        SumI += i.effectifI
        SumR += i.effectifR
    return SumS, SumI, SumR