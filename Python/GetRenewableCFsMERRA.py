import os, copy, datetime
from os import path
from statistics import mode
from AuxFuncs import *
import pandas as pd
import datetime as dt
import numpy as np
from netCDF4 import Dataset

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getREGen(genFleet,tgtTz,reYear,currYear):
    #Isolate wind & solar units
    windUnits,solarUnits = getREInFleet('Wind',genFleet),getREInFleet('Solar',genFleet)
    #Get list of wind / solar sites in region
    lats,lons,cf = loadMerraData(reYear)
    #Match to CFs
    get_cf_index(windUnits,lats,lons),get_cf_index(solarUnits,lats,lons)
    #Get hourly generation (8760 x n df, n = num generators). Use given met year data but set dt index to currYear.
    windGen = get_hourly_RE_impl(windUnits,cf['wind'],currYear)
    solarGen = get_hourly_RE_impl(solarUnits,cf['solar'],currYear)
    #Shift into right tz
    windGen,solarGen = shiftTz(windGen,tgtTz,currYear,'wind'),shiftTz(solarGen,tgtTz,currYear,'solar')
    #Combine by region and fill in missing regions with zeros
    windGenByRegion = windGen.groupby(level='region',axis=1).sum()
    solarGenByRegion = solarGen.groupby(level='region',axis=1).sum()
    for genByRegion in [windGenByRegion,solarGenByRegion]:
        regionsNoGen = [r for r in genFleet['region'].unique() if r not in genByRegion.columns]
        for r in regionsNoGen: genByRegion[r] = 0
    return windGen,solarGen,windGenByRegion,solarGenByRegion

def getREInFleet(reType,genFleet):
    reUnits = genFleet.loc[genFleet['FuelType']==reType]
    return reUnits

# Get all necessary information from powGen netCDF files, VRE capacity factors and lat/lons
# Outputs: numpy arrays of lats and lons, and then a dictionary w/ wind and solar cfs
# as an np array of axbxc, where a/b/c = # lats/# longs/# hours in year
def loadMerraData(reYear):
    #File and dir
    dataDir = 'Data\\MERRA'
    solarFile,windFile = path.join(dataDir,str(reYear) + '_solar_generation_cf.nc'),path.join(dataDir,str(reYear) + '_wind_generation_cf.nc')
    # Error Handling
    if not (path.exists(solarFile) and path.exists(windFile)):
        error_message = 'Renewable Generation files not available:\n\t'+solarFile+'\n\t'+windFile
        raise RuntimeError(error_message)
    #Load data
    solarPowGen = Dataset(solarFile)
    windPowGen = Dataset(windFile) #assume solar and wind cover same geographic region
    #Get lat and lons for both datasets
    lats,lons = np.array(solarPowGen.variables['lat'][:]), np.array(solarPowGen.variables['lon'][:])
    #Store data
    cf = dict()
    cf["solar"] = np.array(solarPowGen.variables['cf'][:]) 
    cf["wind"] = np.array(windPowGen.variables['cf'][:])
    solarPowGen.close(),windPowGen.close()
    # Error Handling
    if cf['solar'].shape != (lats.size, lons.size, 8760):
        print("powGen Error. Expected array of shape",lats.size,lons.size,8760,"Found:",cf['solar'].shape)
        return -1
    return lats,lons,cf

# Convert the latitude and longitude of the vg into indices for capacity factor matrix,
#then save that index into df
# More detail: The simulated capacity factor maps are of limited resolution. This function
#               identifies the nearest simulated location for renewable energy generators
#               and replaces those generators' latitudes and longitudes with indices for 
#               for the nearest simulated location in the capacity factor maps
def get_cf_index(RE_generators, powGen_lats, powGen_lons):
    RE_generators.loc[:,"lat idx"] = find_nearest_impl(RE_generators["Latitude"].astype(float), powGen_lats).astype(int)
    RE_generators.loc[:,"lon idx"] = find_nearest_impl(RE_generators["Longitude"].astype(float), powGen_lons).astype(int)

# Find index of nearest coordinate. Implementation of get_RE_index
def find_nearest_impl(actual_coordinates, discrete_coordinates):
    indices = []
    for coord in actual_coordinates:
        indices.append((np.abs(coord-discrete_coordinates)).argmin())
    return np.array(indices)

# Find expected hourly capacity for RE generators. Of shape (8760 hrs, num generators)
def get_hourly_RE_impl(RE_generators,cf,yr):
    # combine summer and winter capacities
    RE_nameplate = np.tile(RE_generators["Capacity (MW)"].astype(float),(8760,1))
    # multiply by variable hourly capacity factor
    hours = np.tile(np.arange(8760),(RE_generators['Capacity (MW)'].size,1)).T # shape(8760 hrs, num generators)
    RE_capacity = np.multiply(RE_nameplate, cf[RE_generators["lat idx"], RE_generators["lon idx"], hours])
    #hours (no leap day!)
    yr = str(yr)
    idx = pd.date_range('1/1/'+yr + ' 0:00','12/31/' + yr + ' 23:00',freq='H')
    idx = idx.drop(idx[(idx.month==2) & (idx.day ==29)])
    reGen = pd.DataFrame(RE_capacity,index=idx,columns=[RE_generators['GAMS Symbol'].values,RE_generators['region'].values])
    reGen.columns.names = ['GAMS Symbol','region']
    return reGen

#shift tz (MERRA in UTC)
def shiftTz(reGen,tz,yr,reType):
    origIdx = reGen.index
    tzOffsetDict = {'CST': -6,'EST': -5}
    reGen.index = reGen.index.shift(tzOffsetDict[tz],freq='H')
    reGen = reGen[reGen.index.year==yr]
    reGen = reGen.append([reGen.iloc[-1]]*abs(tzOffsetDict[tz]),ignore_index=True)
    if reType=='solar': reGen.iloc[-5:] = 0 #set nighttime hours to 0
    reGen.index=origIdx
    return reGen