#Michael Craig
#August 30, 2021
#Import transmission and regional data

import pandas as pd, numpy as np, geopandas as gpd, os, sys

################################################################################
#Define regions with inter-regional transmission constraints
#Output: dict of zone:[p-regions from REEDS] (zone is same as p-region in WECC & ERCOT since no aggregation)
def defineTransmissionRegions(interconn,balAuths):
    transRegions = dict()
    if balAuths == 'full': #running full interconnection
        if interconn == 'ERCOT' or interconn == 'WECC':
            regions = pd.read_csv(os.path.join('Data','REEDS','regions_default.csv'))
            csvRegions = {'ERCOT':'texas','WECC':'western','EI':'eastern'}
            regionRows = regions.loc[regions['interconnect']==csvRegions[interconn]]
            pRegions = regionRows['p'].unique()
            for p in pRegions: transRegions[p] = [p] #want value to be list so consistent format across interconns 
        elif interconn == 'EI':
            transRegions = createEIGroupings(transRegions)
    else:
        sys.exit('defineTransmissionRegions not set up for non-full run!')
    return transRegions

def createEIGroupings(transRegions):
    transRegions['SERC'] = list(range(87,99)) + list(range(101,103))
    transRegions['NY'] = [127,128]
    transRegions['NE'] = list(range(129,135))
    transRegions['MISO'] = [37] + list(range(42,47)) + list(range(68,87)) + [56,58,66] + list(range(103,109))
    transRegions['PJM'] = list(range(109,127)) + [99,100]
    transRegions['SPP'] = [35,36] + list(range(38,42)) + list(range(47,56)) + [57]
    for r,p in transRegions.items():
        transRegions[r] = ['p' + str(i) for i in p]
    return transRegions
################################################################################

################################################################################
def setupTransmissionAndZones(genFleetFull,transRegions,interconn,balAuths):
    #Read regions and map to zones
    pRegionShapes = loadRegions(transRegions)
    #Assign generators to zones (through p-regions)
    genFleet = assignGensToPRegions(genFleetFull.copy(),pRegionShapes)
    if interconn == 'ERCOT': genFleet = pd.concat([genERCOTRegionReassign(genFleetFull.copy()),genFleet])
    checkZones(genFleet,transRegions) 
    #Get transmission constraints between zones
    limits,costs,dists = getTransmissionData(interconn,transRegions,pRegionShapes)
    return genFleet,transRegions,limits,dists,costs,pRegionShapes

############################
def loadRegions(transRegions):
    pRegionShapes = importPRegions()
    pRegions = getAllPRegions(transRegions)
    pRegionShapes = pRegionShapes.loc[pRegionShapes['PCA_Code'].isin(pRegions)]
    transRegionsReversed = reverseTransRegions(transRegions)
    pRegionShapes['region'] = pRegionShapes['PCA_Code'].map(transRegionsReversed)
    return pRegionShapes

def importPRegions():
    return gpd.read_file(os.path.join('Data','REEDS','Shapefiles','PCAs.shp'))

def getAllPRegions(transRegions):
    allPRegions = list()
    for r,pRegions in transRegions.items(): allPRegions += pRegions
    return allPRegions

def reverseTransRegions(transRegions):
    transRegionsReversed = dict()
    for zone,pRegions in transRegions.items(): 
        for p in pRegions: 
            transRegionsReversed[p] = zone
    return transRegionsReversed

############################
# Assign generators to appropriate regions
def assignGensToPRegions(genFleet,pRegionShapes):
    genFleetGdf = gpd.GeoDataFrame(genFleet, geometry=gpd.points_from_xy(genFleet.Longitude, genFleet.Latitude))
    genFleetPCA = gpd.sjoin(genFleetGdf, pRegionShapes, how="inner", op='intersects')
    genFleet = pd.DataFrame(genFleetPCA)
    return genFleet

############################
#REEDS regions don't cover many ERCOT generators, so push ERCOT gens outside of 
#REEDS ERCOT bounds into ERCOT per regionSubs.
def genERCOTRegionReassign(genFleet):
    regionSubs = {'p57':'p63','p48':'p60','p66':'p64'}
    pRegionShapes = importPRegions()
    pRegionShapes = pRegionShapes.loc[pRegionShapes['PCA_Code'].isin([k for k in regionSubs])]
    pRegionShapes['region'] = pRegionShapes['PCA_Code'].map(regionSubs)
    genFleet = assignGensToPRegions(genFleet,pRegionShapes)
    return genFleet

###########################
def checkZones(genFleet,transRegions):
    genZones = genFleet['region'].unique()
    transZones = [z for z in transRegions]
    for z in genZones:
        if z not in transZones:
            sys.exit(('Generator zone ' + z + ' not in transmission zones ',transZones))

##########################
def getTransmissionData(interconn,transRegions,pRegionShapes):
    limits,costs,dists = importTransmissionData(interconn,transRegions,pRegionShapes)
    dists,limits,costs = expandTransmissionData(limits,costs,dists)
    #Should have distances, costs, and limits for all possible lines; if not, need to align (likely add zero existing capacities). 
    assert((dists.shape[0] == limits.shape[0]) and (dists.shape[0] == costs.shape[0]))
    return limits,costs,dists

def importTransmissionData(interconn,transRegions,pRegionShapes):
    #Read in data
    limits = pd.read_csv(os.path.join('Data','REEDS','transmission_capacity_initial.csv'),header=0)
    costs = pd.read_csv(os.path.join('Data','REEDS','transmission_line_cost.csv'),names=['r','cost($/mw-mile)'])
    #Filter dfs to only include lines that start/end within p-regions of interest
    pRegions = getAllPRegions(transRegions)
    limits = limits.loc[limits['r'].isin(pRegions)]
    limits = limits.loc[limits['rr'].isin(pRegions)]
    costs = costs.loc[costs['r'].isin(pRegions)]
    #If have multi-p-region zones, aggregate transmission data (median of costs & distances; sum of initial capacity)
    if interconn == 'EI':
        #Get regional (instead of p-region) centroids for distance calculations
        regionShapes = pRegionShapes.dissolve(by='region')
        centroids = regionShapes.centroid
        #Map p-regions to zones for aggregation
        transRegionsReversed = reverseTransRegions(transRegions)
        limits['rZone'] = limits['r'].map(transRegionsReversed)
        limits['rrZone'] = limits['rr'].map(transRegionsReversed)
        limits = limits.loc[limits['rZone']!=limits['rrZone']]
        totalLimits,medianCosts,dists = list(),list(),list()
        #For each pair of regions, get inter-regional costs, distances, capacity
        for r in transRegions:
            for rr in transRegions:
                if r != rr:
                    #Get inter-regional lines (if any)
                    interZoneLimits = limits.loc[(limits['rZone']==r) & (limits['rrZone']==rr)]
                    if interZoneLimits.shape[0]>0:
                        #Get limits as sum of interregional capacity
                        totalLimitAC,totalLimitDC = interZoneLimits['AC'].sum(),interZoneLimits['DC'].sum()
                        totalLimits.append(pd.Series({'r':r,'rr':rr,'AC':totalLimitAC,'DC':totalLimitDC}))
                        #Get costs for p-regions in both BAs, take median by BA, then average across BAs
                        interZonePRegions = interZoneLimits['r'].unique()
                        interZoneCostsR = costs.loc[costs['r'].isin(interZonePRegions)]
                        interZonePRegions = interZoneLimits['rr'].unique()
                        interZoneCostsRR = costs.loc[costs['r'].isin(interZonePRegions)]
                        medianR = np.median(interZoneCostsR['cost($/mw-mile)'])
                        medianRR = np.median(interZoneCostsRR['cost($/mw-mile)'])
                        medianCosts.append(pd.Series({'r':r,'rr':rr,'cost($/mw-mile)':((medianR + medianRR)/2)}))
                        #Get distance b/wn region centroids
                        centR,centRR = centroids[r],centroids[rr]
                        dist = haversine(centR,centRR)
                        dists.append(pd.Series({'r':r,'rr':rr,'dist(mile)':haversine(centR,centRR)}))
        limits = pd.concat(totalLimits,axis=1).T
        costs = pd.concat(medianCosts,axis=1).T
        dists = pd.concat(dists,axis=1).T
    else: #load distances if not aggregating data
        dists = pd.read_csv(os.path.join('Data','REEDS','transmission_distance.csv'),header=0)
        dists = dists.loc[dists['r'].isin(pRegions)]
        dists = dists.loc[dists['rr'].isin(pRegions)]
    return limits,costs,dists

def haversine(centR,centRR,kmToMi=0.621371):
    latR,lonR = np.radians(centR.y),np.radians(centR.x)
    latRR,lonRR = np.radians(centRR.y),np.radians(centRR.x)
    R = 6378.137 #km; mean Earth radius
    h = (np.sin((latR-latRR)/2))**2 + np.cos(latR) * np.cos(latRR) * (np.sin((lonR-lonRR)/2))**2
    c = 2 * np.arcsin(np.sqrt(h))
    return R*c*kmToMi #dist in mi

#Expand transmission data by (1) reversing source-sink line limits & (2) mapping costs to lines (not regions)
#& (3) add GAMS symbols as r + rr. 
def expandTransmissionData(limits,costs,dists):
    #Reverse source-sink limits
    limitsReversed = limits.rename(columns={'r':'rr','rr':'r'})
    limits = pd.concat([limits,limitsReversed])
    limits.reset_index(inplace=True,drop=True)
    #Map costs to lines
    costs = costs.set_index(costs['r'].values).to_dict()['cost($/mw-mile)']
    costsAllLines = dists[['r','rr']].copy()
    costsAllLines['rCost'],costsAllLines['rrCost'] = costsAllLines['r'].map(costs),costsAllLines['rr'].map(costs)
    costsAllLines['Line Cost ($/mw-mile)'] = (costsAllLines['rCost'] + costsAllLines['rrCost'])/2
    #Add GAMS symbols
    for df in [dists,limits,costsAllLines]: df['GAMS Symbol'] = df['r'] + df['rr']
    return dists,limits,costsAllLines
################################################################################

#Unused distance code in getTransmissionData
# dists = pd.read_csv(os.path.join('Data','REEDS','transmission_distance.csv'),header=0)
# dists = dists.loc[dists['r'].isin(pRegions)]
# dists = dists.loc[dists['rr'].isin(pRegions)]
# dists['rZone'] = dists['r'].map(transRegionsReversed)
# dists['rrZone'] = dists['rr'].map(transRegionsReversed)
# dists = dists.loc[dists['rZone']!=dists['rrZone']]
# #Get distances
# interZoneDists = dists.loc[(dists['rZone']==r) & (dists['rrZone']==rr)]
# medDistAC,medDistDC = np.median(interZoneDists['AC'].values),np.median(interZoneDists['DC'].values)
# medianDists.append(pd.Series({'r':r,'rr':rr,'AC':medDistAC,'DC':medDistDC}))