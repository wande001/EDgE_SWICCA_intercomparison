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
  year = startTime.year
  month = startTime.month
  while datetime.datetime(year, month, 1) < endTime:
    times.append(datetime.datetime(year, month, 1))
    month += 1
    if month > 12: year += 1; month -= 12
  return np.array(times)

def addMonthsToDate(startTime, months):
  times = []
  year = startTime.year
  month = startTime.month + months
  if month > 12: year += 1; month -= 12
  return datetime.datetime(year, month, 1)

def aggregateToMonth(data, times):
  dt, dx = data.shape
  nextMonth = addMonthsToDate(datetime.datetime(times[0].year, times[0].month, 1),1) - datetime.timedelta(days=1)
  numMonth = 1
  for t in range(dt):
    if times[t] == nextMonth and t != dt-1:
      numMonth += 1
      nextMonth = addMonthsToDate(times[t+1], 1) - datetime.timedelta(days=1)
  monthCount = 0
  monthData = np.zeros((numMonth,dx))
  monthTimes = []
  curMonth = 0
  nextMonth = addMonthsToDate(datetime.datetime(times[0].year, times[0].month, 1),1) - datetime.timedelta(days=1)
  for t in range(dt):
    monthData[curMonth,:] += data[t,:]
    monthCount += 1
    if times[t] == nextMonth and t != dt-1:
      monthData[curMonth,:] = monthData[curMonth,:]/monthCount
      curMonth += 1
      monthCount = 0
      monthTimes = datetime.datetime(times[t].year, times[t].month, 1)
      nextMonth = addMonthsToDate(times[t+1], 1) - datetime.timedelta(days=1)
    elif t == dt-1:
	  monthData[curMonth,:] = monthData[curMonth,:]/monthCount
  return monthData, monthTimes
  
def computeQuintiles(data, month=True):
  if month:
    dt, dx = data.shape
    quintData = np.ones((dt,dx))
    for x in range(dx):
      for m in range(12):
        tSel = np.arange(m, dt, 12)
        dataMonthTemp = data[tSel,x]
        limits = np.nanpercentile(dataMonthTemp, [20,40,60,80,100], axis=0)
        quintTemp = np.ones((len(dataMonthTemp)))
        for l in range(4):
          quintTemp[dataMonthTemp > limits[l]] = l+2
        quintData[tSel,x] = quintTemp
  else:
    dt, dx = data.shape
    quintData = np.ones((dt,dx))
    for x in range(dx):
      dataMonthTemp = data[:,x]
      limits = np.nanpercentile(dataMonthTemp, [20,40,60,80,100], axis=0)
      quintTemp = np.ones((len(dataMonthTemp)))
      for l in range(4):
        quintTemp[dataMonthTemp > limits[l]] = l+2
      quintData[:,x] = quintTemp
  quintData[np.isnan(data)] = np.nan
  return quintData
    
#### Functions to read forecast data (adjust for each modelling group!!!) #######

modelS = ["ECMWF-S4","LFPW","CanCM4","FLOR", "ESP"]

def makeInputNameSeasonal(var, mod, hyd_mod):
    model = modelS[mod]
    if hyd_mod == 'mHM':
        sciiInput = "/data/edge/data/processed/SCIIs/2017/pan_eu/SCII-8/categorical_values"
    if hyd_mod == 'noah-mp':
        sciiInput = "/data/edge/data/processed/noahmp_output/edge_domain/seasonal_forecast/SCII-8/categorical_values"
    if hyd_mod == 'PCRGLOBWB':
        sciiInput = "/home/land1/niko/EDgE/outputSCII/seasonal_new/categorical_values_scii_8"
    inputFileName = "%s_%s_%s-probabilistic-quintile-distribution_monthly_1993_01_2012_05.nc" %(model, hyd_mod, var)
    inputFile = "%s/%s" %(sciiInput, inputFileName)
    print "Reading %s" %(inputFile)
    return inputFile

def makeRefFileName(hyd_mod):
    if hyd_mod == 'mHM':
        referenceFile = "/data/edge/data/processed/mhm_output/edge_domain/seasonal_forecast_2017/mhm_ref_run_eobs/reference_EOBS_forced_%s_monthly_streamflow_categorical_values_scii_08_011993_122011.nc" %(hyd_mod)
    if hyd_mod == 'noah-mp':
        referenceFile = "/data/edge/data/processed/noahmp_output/edge_domain/eobs/reference_EOBS_forced_%s_monthly_streamflow_categorical_values_scii_08_011993_122011.nc" %(hyd_mod)
    if hyd_mod == 'PCRGLOBWB':
        referenceFile = "/home/land1/niko/EDgE/reference/routing/reference_EOBS_forced_PCRGLOBWB_monthly_streamflow_categorical_values_scii_08_011993_122011.nc"
    print "Reading %s" %(referenceFile)
    return referenceFile

def makeOutputNameSeasonal(var, mod, uncertaintyOutput, GHM="PCRGLOBWB"):
    model = modelS[mod]
    outputFileName = "%s_%s_seasonal-forecast-%s-uncertainty.nc" %(var, model, GHM)
    outputFile = "%s/%s" %(uncertaintyOutput, outputFileName)
    return outputFile


def readForecastData(var, mod, obsLon, obsLat, startTime, endTime, varName = "scii_category", GHM="PCRGLOBWB", lead=0, quant=0):
  dates = generateDateRange(datetime.datetime(1993,1,1), addMonthsToDate(datetime.datetime(2011,12,31),lead))
  ## Select start and end of forecast simulation time period
  start = np.argmin(dates < startTime)
  end = np.argmax(dates > endTime)
  ## Create seasonal input file name
  inputFile = makeInputNameSeasonal(var, mod, GHM)
  ## Read seasonal forecast data for start to end, for quant (between 0 and 5) and for lead month (between 0-6)
  ncFile = nc.Dataset(inputFile)
  data = ncFile.variables[varName][range(start, end),quant,lead,:,:]
  ncFile.close()
  dataPoints = data[:,obsLat,obsLon]
  return dataPoints

def readReferenceData(obsLon, obsLat, startTime, endTime, varName = "scii_category", GHM="PCRGLOBWB", lead=0):
  ## Read the data from the reference NetCDF file and report back as Numpy Array
  dates = generateDateRange(addMonthsToDate(datetime.datetime(1993,1,1),lead), datetime.datetime(2011,12,31))
  ## Select start and end of simulation time period
  start = np.argmin(dates < startTime)
  end = np.argmax(dates > endTime)
  ## Create reference file name
  referenceFile = makeRefFileName(GHM)
  ## Read reference data
  refFile = nc.Dataset(referenceFile)
  refData = refFile.variables[varName][range(start, end),:,:]
  lat = refFile.variables["y"][:]
  lon = refFile.variables["x"][:]
  refFile.close()
  dataPoints = refData[:,obsLat,obsLon]
  return dataPoints

def readCoordinatesData(GHM="PCRGLOBWB"):
  ## Read coordinates and time properties from the reference file
  referenceFile = makeRefFileName(GHM)
  refFile = nc.Dataset(referenceFile)
  time = refFile.variables["time"][:]
  lat = refFile.variables["y"][:]
  lon = refFile.variables["x"][:]
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
                                      var4D, longName = None, loop=False, quantDim = 5):
	
    rootgrp= nc.Dataset(ncFileName,'w', format='NETCDF4')
    #-create dimensions - time is unlimited, others are fixed
    rootgrp.createDimension('station',len(stations))
    rootgrp.createDimension('lead_time',6)
    rootgrp.createDimension('quantile',quantDim)
		
    station= rootgrp.createVariable('station','f4',('station',))
    station.standard_name= 'FID Discharge'
    station.long_name= 'FID Discharge station, derived from SMHI database'
    station.units= '-'
    
    lead= rootgrp.createVariable('lead_time','f4',('lead_time',))
    lead.standard_name= 'lead time'
    lead.long_name= 'lead time from the forecasting date'
    lead.units= 'month-id'
    
    quant= rootgrp.createVariable('quant','f4',('quantile',))
    quant.standard_name= 'quantile'
    quant.long_name= '1:Q_forecast <= Qref_quantile_level_20; 2: Qref_quantile_level_20 < Q_forecast <= Qref_quantile_level_40; 3: Qref_quantile_level_40 < Q_forecast <= Qref_quantile_level_60; 4: Qref_quantile_level_60 < Q_forecast <= Qref_quantile_level_80; 5: Q_forecast > Qref_quantile_level_80'
    quant.units= '-'
    quant[:] = range(1,6) #[0.1,0.3,0.5,0.7,0.9]

    lats= rootgrp.createVariable("lat",'f4',('station',) ,fill_value=MV,zlib=True)
    lats.standard_name= 'Latitude'
    lats.long_name= 'Latitude of stations'
    lats[:] = lat

    lons= rootgrp.createVariable("lon",'f4',('station',) ,fill_value=MV,zlib=True)
    lons.standard_name= 'Longitude'
    lons.long_name= 'Longitude of stations'
    lons[:] = lon
	
    station[:]= stations
    
    lead[:] = range(1,7)
        
    rootgrp.description = 'Uncertainty and skill analysis for EDgE-SWICCA intercomparison'
    rootgrp.history = 'Created by Niko Wanders (n.wanders@uu.nl, nwanders@princeton.edu) and Stephan Thober (stephan.thober@ufz.de) on May 2017'
    rootgrp.source = 'Comparison of observed discharge and forecasted discharge for five quintiles over the EDgE-SWICCA EU domain'
    
    if loop:
        for i in range(len(varName)):
            shortVarName = varName[i]
            longVarName  = longName[i]
            if longName != None: longVarName = longName
            if var4D[i]:
              var= rootgrp.createVariable(shortVarName,'f4',('quantile','lead_time','station',) ,fill_value=MV,zlib=True)
            else:
              var= rootgrp.createVariable(shortVarName,'f4',('lead_time','station',) ,fill_value=MV,zlib=True)
            var.standard_name = varName[i]
            var.long_name = longVarName[i]
            var.units = varUnits[i]
    else:
        shortVarName = varName
        longVarName  = varName
        if longName != None: longVarName = longName
        if var4D[i]:
          var= rootgrp.createVariable(shortVarName,'f4',('quantile','lead_time','station',) ,fill_value=MV,zlib=True)
        else:
          var= rootgrp.createVariable(shortVarName,'f4',('lead_time','station',) ,fill_value=MV,zlib=True)
        var.standard_name = varName
        var.long_name = longVarName
        var.units = varUnits
        
    rootgrp.sync()
    rootgrp.close()

    
def data2NetCDFnoTime(ncFile,varName,var4D,varField, quant=0, lead=0):
  #-write data to netCDF
  rootgrp= nc.Dataset(ncFile,'a')
  
  shortVarName= varName
  if var4D:
    rootgrp.variables[shortVarName][quant,lead,:]= (varField)
  else:
    rootgrp.variables[shortVarName][lead,:]= (varField)
  rootgrp.sync()
  rootgrp.close()
    
def writeMetricsToFile(outputFile, metrics, varNames, quant, lead, var4D):
  ## write metrics to NetCDF output file
  for var in range(len(varNames)):
    data2NetCDFnoTime(outputFile,varNames[var],var4D[var], metrics[var], quant=quant, lead=lead)	

def writeAcToFile(outputFile, metrics, varNames, quant, lead, var4D):
  ## Write AC to NetCDF output file
  data2NetCDFnoTime(outputFile,varNames,var4D, metrics, quant=quant, lead=lead)

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

########### Metric functions ########################3

def brier(obs, pred, quant=1, limit=0.2):
  out = np.mean(((obs == quant) - (pred > limit))**2, axis=0)
  return out

def correlation2D(obs, pred):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = obs - obs.mean(axis=0)
    B_mB = pred - pred.mean(axis=0)
    # Sum of squares across rows
    ssA = (A_mA**2).sum(axis=0)
    ssB = (B_mB**2).sum(axis=0)
    # Finally get corr coeff
    return np.sum(A_mA * B_mB, axis=0)/np.sqrt(ssA * ssB)

def PODandFAR(obs, pred, quant = 1, limit=0.2):
  dt, dx = obs.shape
  events = obs.flatten() == float(quant)
  hits = ((events) & (pred.flatten() > limit)).reshape((dt,dx))
  falseAlarm = ((events == False) & (pred.flatten() > limit)).reshape((dt,dx))
  misses = ((events) & (pred.flatten() <= limit)).reshape((dt,dx))
  correctNegative = ((events == False) & (pred.flatten() <= limit)).reshape((dt,dx))
  numEvents = np.sum(events.reshape((dt,dx)), axis=0, dtype=np.float)
  hits = np.sum(hits, axis=0, dtype=np.float)
  misses = np.sum(misses, axis=0, dtype=np.float)
  falseAlarm = np.sum(falseAlarm, axis=0, dtype=np.float)
  correctNegative = np.sum(correctNegative, axis=0, dtype=np.float)
  PODout = hits/numEvents
  FARout = falseAlarm/(hits+falseAlarm)
  POFDout = falseAlarm/(correctNegative+falseAlarm)
  TSout = hits/(hits+misses+falseAlarm)
  return PODout, FARout, POFDout, TSout

def forecastProbability(pred):
  forecastProb = np.mean(pred, axis=0)
  return forecastProb

def computeMetrics(refData, forecastData, quant, limit=0.2):
  ## Main function to compute metrics, will call other functions for the individual metrics, returns an object with all metrics
  brierout = brier(refData, forecastData, quant=quant+1, limit=limit)
  PODout, FARout, POFDout, TSout = PODandFAR(refData, forecastData, quant=quant+1, limit=limit)
  forecastProbout = forecastProbability(forecastData)  
  metrics = [brierout, PODout, FARout, POFDout, TSout, forecastProbout]
  return metrics

def computePercentiles(forecastData, probs = [10,25,50,75,90]):
  metrics = np.percentile(forecastData, probs, axis=0)
  return metrics

def catchmentStatistics(inputData, catchmentSWICCA):
  result = np.ones((inputData.shape))
  if len(inputData.shape) == 4:
    quant, lead, lat, lon = inputData.shape
    for q in range(quant):
      for l in range(lead):
        df = pd.DataFrame({ 'catchment' : catchmentSWICCA.reshape((lat*lon)), 'data' : inputData[q,l,:,:].reshape(lat*lon)})
        catchAvg = df.groupby(['catchment']).mean()
        locs = catchAvg.index.values
        catchData = np.array(catchAvg["data"])[:]
        for c in range(len(locs)):
          resultsMask = catchmentSWICCA==locs[c]
          result[q,l,:,:][resultsMask] = catchData[c]
        resultsMask = catchmentSWICCA == -999
        result[q,l,:,:][resultsMask] = MV
  return result

#### Main functions to handle data processing ###########

def seasonalUncertainty(var, mod, GHM, startTime, endTime):
    ## Reading and Matching observation data
    stationNums = readObservationsStations("MIP_Q_series.txt")
    matchOrder = matchObsPointsObservation(stationNums)
    ## Vars for output
    varNames = ["brier", "POD", "FAR", "POFD", "TS","forecastProb","AC"]
    longNameS = ["Brier Score", "Probability of detection", "False Alarm Rate", "Probability of False Detection", "Threat Score", "Forecast probability for quantile", "Anomaly correlation"]
    varUnits = ["-", "-", "-", "-", "-", "-", "-"]
    ## Define if output variable is 4D, dimensions, lon, lat, initialization time and quantile
    var4D = [True, True, True, True, True, True, False]
    ## Point to output directories of seasonal forecasts
    if GHM == 'PCRGLOBWB':
        uncertaintyOutput = "/home/land1/niko/EDGE/"
    else:
        uncertaintyOutput = "/work/thober/tmp/"
    quantDim = 5
    try:
      os.mkdir(uncertaintyOutput)
    except:
      pass
    ## Point to location  of reference (baseline) simulation
    referenceFile = makeRefFileName(GHM)
    outputFile = makeOutputNameSeasonal(var, mod, uncertaintyOutput, GHM)
    print outputFile
    ## Read coordinates of stations
    time, lat, lon = readCoordinatesData(GHM)
    obsLat, obsLon, stationLat, stationLon, stations = matchObsPointsModel(lat=lat, lon=lon, shapeFileName = "MIP_Stations_LAEA/MIP_Stations_LAEA.shp")
    ## Create output file
    createNetCDFnoTime(outputFile, varNames, varUnits, stations, stationLat, stationLon, var4D, longName= longNameS, loop=True, quantDim = quantDim)

    for lead in range(6):
      ## Read Observational time series
      data, time = readObservationsData("MIP_Q_series.txt", matchOrder, stationNums, startTime = addMonthsToDate(startTime,lead), endTime=addMonthsToDate(endTime, lead))
      ## Aggregate data to monthly values
      monthData, monthTimes = aggregateToMonth(data, time)
      ## Convert to quintile data
      quintileData = computeQuintiles(monthData)
      ## Loop over the quantiles, input NetCDF files have a 4 dimensional structure or lon, lat, time, quantile for Reference and 5D for forecast, lon, lat time, quantile, Lead
      for quant in range(quantDim):
        print var, modelS[mod], lead, quant
        ## Read forecast data for certaint GHM, quantile and lead
        forecastData = readForecastData(var, mod, obsLon, obsLat, startTime, endTime, varName = "prob_quintile_dist", GHM=GHM, quant=quant-1, lead=lead)
        ## For first quantile output matrix should be created
        if quant == 0:
          ensForecast = forecastData * (quant-1)
        ## New data added to the existing output matrix
        else:
          ensForecast += forecastData * (quant-1)
        refData = readReferenceData(obsLon, obsLat, startTime, endTime, varName = "scii_category", GHM=GHM, lead=lead)
        ## Forecast validation against reference data
        metrics = computeMetrics(refData, forecastData, quant)
        ## Forecast validation against observations, disabled for now
        # metrics = computeMetrics(quintileData, forecastData, quant)
        writeMetricsToFile(outputFile,metrics, varNames[:-1], quant=quant, lead=lead, var4D=var4D[:-1])
        ## If last quantile write AC correlation for all quantiles
        if quant == 4:
          acMetrics = correlation2D(quintileData, ensForecast)
          writeAcToFile(outputFile, acMetrics, varNames[-1], quant=quant, lead=lead, var4D=var4D[-1])       


def computeUncertainty(varName, mod, GHM, startTime = datetime.datetime(1993,1,1), endTime = datetime.datetime(2009,12,31)):
  seasonalUncertainty(varName, mod, GHM, startTime, endTime)
  return "Computation SCII %s model %s finished for GHM %s" %(varName, modelS[mod], GHM)


####### Run single instance seasonal ##########

computeUncertainty("streamflow", 1, 'PCRGLOBWB', datetime.datetime(1995,1,1), datetime.datetime(1996,12,31))
# computeUncertainty("streamflow", 1, 'noah-mp', datetime.datetime(1995,1,1), datetime.datetime(1996,12,31))
# computeUncertainty("streamflow", 1, 'mHM', datetime.datetime(1995,1,1), datetime.datetime(1996,12,31))

####### Run parralel #################

#varS = []
#modS = []

#for var in ["discharge"]:
#  for mod in range(5):
#    varS.append(scii)
#    modS.append(mod)

#nTot = len(sciiS)
#pool = mp.Pool(processes=5)
#results = [pool.apply_async(computeUncertainty,args=(sciiS[num], modS[num])) for num in range(nTot)]
#output = [p.get() for p in results]

#print output
