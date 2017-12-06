import netCDF4 as nc
import os
import numpy as np
import pandas as pd
import datetime
## Require for multi-core simulations
import multiprocessing as mp
import shapefile
##

### Missing value#####
MV = -9999

### Generic functions#######
def generateDateRange(startTime, endTime):
  times = []
  dayCount = 0
  while startTime + datetime.timedelta(days=dayCount) <= endTime:
    times.append(startTime + datetime.timedelta(days=dayCount))
    dayCount += 1
  return np.array(times)

def addMonthsToDate(startTime, months):
  times = []
  year = startTime.year
  month = startTime.month + months
  if month > 12: year += 1; month -= 12
  return datetime.datetime(year, month, 1)

#### Functions to read forecast data (adjust for each modelling group!!!) #######

def makeRefFileName(hyd_mod):
    if hyd_mod == 'mHM':
        referenceFile = "/data/edge/data/processed/mhm_output/edge_domain/seasonal_forecast_2017/mhm_ref_run_eobs/reference_EOBS_forced_%s_monthly_streamflow_categorical_values_scii_08_011993_122011.nc" %(hyd_mod)
    elif hyd_mod == 'noah-mp':
        referenceFile = "/data/edge/data/processed/noahmp_output/edge_domain/eobs/reference_EOBS_forced_%s_monthly_streamflow_categorical_values_scii_08_011993_122011.nc" %(hyd_mod)
    elif hyd_mod == 'PCRGLOBWB':
        referenceFile = "/home/land1/niko/EDgE/reference/routing/mRM_Fluxes_States_1993_2011.nc"
    elif hyd_mod == 'VIC':
        referenceFile = "/home/land1/mpan/for_stephan/E-OBS/mRM/mRM_Fluxes_States_199301-201412.nc"
    print "Reading %s" %(referenceFile)
    return referenceFile

def makeOutputNameReference(uncertaintyOutput, GHM="PCRGLOBWB"):
    outputFileName = "reference-%s-discharge.nc" %(GHM)
    outputFile = "%s/%s" %(uncertaintyOutput, outputFileName)
    return outputFile


def readReferenceData(obsLon, obsLat, startTime, endTime, varName = "scii_category", GHM="PCRGLOBWB"):
  ## Read the data from the reference NetCDF file and report back as Numpy Array
  if GHM == "PCRGLOBWB":
    dates = generateDateRange(datetime.datetime(1993,1,1), datetime.datetime(2011,12,31))
  if GHM == "VIC":
    dates = generateDateRange(datetime.datetime(1993,1,1), datetime.datetime(2014,12,31))
  ## Select start and end of simulation time period
  start = np.argmin(dates < startTime)
  end = np.argmax(dates > endTime)
  ## Create reference file name
  referenceFile = makeRefFileName(GHM)
  ## Read reference data
  refFile = nc.Dataset(referenceFile)
  dataPoints = np.zeros((len(range(start, end)), len(obsLat), len(obsLon)))
  dataCount = 0
  for time in range(start, end):
    print time
    refData = refFile.variables[varName][time,:,:]
    dataPoints[dataCount,:,:] = refData[obsLat, obsLon]
    dataCount += 1
  lat = refFile.variables["northing"][:]
  lon = refFile.variables["easting"][:]
  refFile.close()
  return dataPoints

def readCoordinatesData(GHM="PCRGLOBWB"):
  ## Read coordinates and time properties from the reference file
  referenceFile = makeRefFileName(GHM)
  refFile = nc.Dataset(referenceFile)
  time = refFile.variables["time"][:]
  lat = refFile.variables["northing"][:]
  lon = refFile.variables["easting"][:]
  refFile.close()
  return time, lat, lon

#### Functions to read observation data#######

def readObservationsStations(fileName):
  f = open(fileName)
  lines=f.readlines()
  ## Read first line
  tempLine = np.array(lines[0].split("\t"))[1:]
  ## Extract station Numbers
  stationNums = []
  for line in tempLine:
    stationNums.append(int(line))
  f.close()
  return np.array(stationNums)


def readObservationsData(fileName, matchOrder, stationNums, startTime, endTime):
  ## Function to read observations from txt file and conver to Numpy array
  f = open(fileName)
  lines=f.readlines()
  ## Extract observations and times
  numLines = len(lines) - 1
  obsData = np.zeros((numLines, len(stationNums))) - 999.
  obsTimes = []
  lineCount = 0
  ## Loop over line numbers
  for line in range(numLines):
    tempLine = lines[line+1].split("\t")
    tempTime = str(tempLine[0]).split('-')
    currDate = datetime.datetime(int(tempTime[0]),int(tempTime[1]),int(tempTime[2]))
    ## Save values only if date is within the range of the modelled time period
    if currDate >= startTime and currDate <= endTime:
      obsTimes.append(currDate)
      tempData = np.zeros((len(stationNums)))
      tempData[:] = tempLine[1:]
      tempData[tempData < 0.0] = np.nan
      for s in range(len(matchOrder)):
        if matchOrder[s] != -999:
          obsData[lineCount,s] = tempData[matchOrder[s]]
        else:
		  obsData[lineCount,s] = np.nan
      lineCount += 1
  f.close()
  return obsData[:lineCount,:], np.array(obsTimes)

######## NetCDF functions ##################3

def createNetCDFnoTime(ncFileName, varName, varUnits, stations, lat, lon,\
                                      longName = None):
	
    rootgrp= nc.Dataset(ncFileName,'w', format='NETCDF4')
    #-create dimensions - time is unlimited, others are fixed
    rootgrp.createDimension('station',len(stations))
    rootgrp.createDimension('time',None)
    
    date_time= rootgrp.createVariable('time','f4',('time',))
    date_time.standard_name= 'time'
    date_time.long_name= 'Days since 1901-01-01'
    date_time.units= 'Days since 1901-01-01' 
    date_time.calendar= 'standard'
		
    station= rootgrp.createVariable('station','f4',('station',))
    station.standard_name= 'FID Discharge'
    station.long_name= 'FID Discharge station, derived from SMHI database'
    station.units= '-'
    
    lats= rootgrp.createVariable("lat",'f4',('station',) ,fill_value=MV,zlib=True)
    lats.standard_name= 'Latitude'
    lats.long_name= 'Latitude of stations'
    lats[:] = lat

    lons= rootgrp.createVariable("lon",'f4',('station',) ,fill_value=MV,zlib=True)
    lons.standard_name= 'Longitude'
    lons.long_name= 'Longitude of stations'
    lons[:] = lon
	
    station[:]= stations
    
    rootgrp.description = 'Memory analysis for EDgE-SWICCA intercomparison'
    rootgrp.history = 'Created by Niko Wanders (n.wanders@uu.nl, nwanders@princeton.edu) and Stephan Thober (stephan.thober@ufz.de) on Dec 2017'
    rootgrp.source = 'Comparison of observed discharge and forecasted discharge over the EDgE-SWICCA EU domain'
    
    shortVarName = varName
    longVarName  = varName
    if longName != None: longVarName = longName
    var= rootgrp.createVariable(shortVarName,'f4',('time','station',) ,fill_value=MV,zlib=True)
    var.standard_name = varName
    var.long_name = longVarName
    var.units = varUnits
        
    rootgrp.sync()
    rootgrp.close()

def data2NetCDF(ncFile,varName,varField,timeStamp,posCnt = None):
  #-write data to netCDF
  rootgrp= nc.Dataset(ncFile,'a')    
  
  shortVarName= varName        
  
  date_time= rootgrp.variables['time']
  if posCnt == None: posCnt = len(date_time)
  
  date_time[posCnt]= nc.date2num(timeStamp,date_time.units,date_time.calendar)
  rootgrp.variables[shortVarName][posCnt,:]= (varField)
  
  rootgrp.sync()
  rootgrp.close()

def writeMetricsToFile(outputFile, data, varName, startTime, endTime):
  ## write metrics to NetCDF output file
  timeRange = generateDateRange(startTime, endTime)
  pos =0
  for time in timeRange:
    data2NetCDF(outputFile,varName, data[pos,:], timeStamp = time)
    pos += 1

########## Matching functions #####################
def matchObsPointsModel(lat, lon, shapeFileName = "MIP_Stations_LAEA/MIP_Stations_LAEA.shp"):
  ## Match coordinate of observations to model's mask. Used to extract the correct modelled timeseries
  sf = shapefile.Reader(shapeFileName)
  points = sf.shapes()
  fields = sf.records()
  obsLon = []
  obsLat = []
  Lon = []
  Lat = []
  stationS = []
  ## Loop over points
  for point in points:
    ## Find point with minimal distance to observation location
    selY = np.argmin(np.abs(lat - point.points[0][1]))
    selX = np.argmin(np.abs(lon - point.points[0][0]))
    Lat.append(point.points[0][1])
    Lon.append(point.points[0][0])
    obsLon.append(selX)
    obsLat.append(selY)
  for field in fields:
    stationS.append(field[9])
  return obsLat, obsLon, Lat, Lon, stationS

def matchObsPointsObservation(stationNums, shapeFileName = "MIP_Stations_LAEA/MIP_Stations_LAEA.shp"):
  ## Matches the Observation points to a shape file ID
  sf = shapefile.Reader(shapeFileName)
  points = sf.shapes()
  fields = sf.records()
  dataMatch = []
  for field in fields:
    match = np.argmax(np.array(stationNums) == int(field[9]))
    if int(field[9]) != stationNums[match]:
      dataMatch.append(-999)
    else:
      dataMatch.append(match)
  return dataMatch

#### Main functions to handle data processing ###########

def extractTimeseriesReference(GHM, startTime, endTime):
    ## Reading and Matching observation data
    stationNums = readObservationsStations("MIP_Q_series.txt")
    matchOrder = matchObsPointsObservation(stationNums)
    data, time = readObservationsData("MIP_Q_series.txt", matchOrder, stationNums, startTime = startTime, endTime=endTime)
    
    ## Point to output directories of seasonal forecasts
    if GHM == 'PCRGLOBWB':
        uncertaintyOutput = "/home/land1/niko/EDgE/EDgE_SWICCA_intercomparison/output"
    elif GHM == 'VIC':
        uncertaintyOutput = "/home/land1/niko/EDgE/EDgE_SWICCA_intercomparison/output"
    else:
        uncertaintyOutput = "/home/land1/niko/EDgE/EDgE_SWICCA_intercomparison/output"
    try:
      os.mkdir(uncertaintyOutput)
    except:
      pass
    ## Point to location  of reference (baseline) simulation
    referenceFile = makeRefFileName(GHM)
    outputFile = makeOutputNameReference(uncertaintyOutput, GHM)
    print outputFile
    ## Read coordinates of stations
    time, lat, lon = readCoordinatesData(GHM)
    obsLat, obsLon, stationLat, stationLon, stations = matchObsPointsModel(lat=lat, lon=lon, shapeFileName = "MIP_Stations_LAEA/MIP_Stations_LAEA.shp")
    print data.shape
    ## Create output file
    createNetCDFnoTime(outputFile, "discharge", "-", stations, stationLat, stationLon, longName= "discharge")
    if GHM == "PCRGLOBWB":
      createNetCDFnoTime("/home/land1/niko/EDgE/EDgE_SWICCA_intercomparison/output/observationalData.nc", "discharge", "-", stations, stationLat, stationLon, longName= "discharge")

    refData = readReferenceData(obsLon, obsLat, startTime, endTime, varName = "Qrouted", GHM=GHM)
    print refData[0,0]
    writeMetricsToFile(outputFile, refData, "discharge", startTime, endTime)
    if GHM == "PCRGLOBWB":
      writeMetricsToFile("/home/land1/niko/EDgE/EDgE_SWICCA_intercomparison/output/observationalData.nc", data, "discharge", startTime, endTime)

def getData(varName, GHM, startTime = datetime.datetime(1993,1,1), endTime = datetime.datetime(2009,12,31)):
  extractTimeseriesReference(GHM, startTime, endTime)
  return "Computation discharge finished for GHM %s" %(GHM)


####### Run single instance seasonal ##########

getData("discharge",'PCRGLOBWB', datetime.datetime(1994,1,1), datetime.datetime(2009,12,31))
getData("discharge",'VIC', datetime.datetime(1994,1,1), datetime.datetime(2009,12,31))
# computeUncertainty("streamflow", 1, 'noah-mp', datetime.datetime(1995,1,1), datetime.datetime(1996,12,31))
# computeUncertainty("streamflow", 1, 'mHM', datetime.datetime(1995,1,1), datetime.datetime(1996,12,31))

####### Run parralel #################

varS = []
modS = []
ghmS = []

for var in ["streamflow"]:
  for mod in range(5):
    varS.append(var)
    modS.append(mod)
    ghmS.append("PCRGLOBWB")

nTot = len(varS)
pool = mp.Pool(processes=5)
#results = [pool.apply_async(computeUncertainty,args=(varS[num], modS[num], ghmS[num])) for num in range(nTot)]
#output = [p.get() for p in results]

#print output
