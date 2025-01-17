#Michael Craig
#October 4, 2016
#Functions update input fleets' fuel prices for given year

from SetupGeneratorFleet import *

def updateFuelPricesAndCosts(fleet,currYear,fuelPrices,regCostFrac):
    fleet = addFuelPrices(fleet,currYear,fuelPrices)
    fleet = calcOpCost(fleet)
    fleet['RegOfferCost($/MW)'] = regCostFrac*fleet['OpCost($/MWh)']*fleet['RegOfferElig']
    return fleet