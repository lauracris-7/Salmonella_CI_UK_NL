#---- CONDITIONAL INCIDENCE ANALYSIS in England and Wales ----#
    # The code looks at how the risk of salmonellosis in humans depends on weather factors


# This script assumes the data has gone through some manipulation and preparation to fit this code
# Contents:
  # detrended conditional incidence analysis for:
    # - individual weather factors analysis, 
    # - 2-way factor analysis, 
    # - 3-way factor analysis. 
  # plots with narrow ribbon displaying only the yearly variability, rahter than the weather variability for all the catchment areas across the country.


  
# Load required packages

library(ISOweek)
library(lubridate)
library(ggplot2)
require(MASS)
library(scales)
require(pheno)
library(timeDate)
library(pastecs)
library(stringi)
library(timeSeries)
library(readr)
library(data.table)
library(xts)
library(MALDIquant)
library(plyr)
library(dplyr)
library(tools)
library(wesanderson)
library(zoo)
library(slider)
library(stats)
library(utils)
library(dplR)
library(WaveletComp)
library(gdata)
library (scales)


#### 1. Set global parameters and read input data #### 

################### Global parameters ###################
# Set the length of the desired time-lag
width <-7
width_char <-paste(width)


# Set years of analysis
first_year <- 2000
last_year <- 2016


################### Incidence data ###################
  # Read Simulation of Salmonella cases averaged for an estimated day of infection. 
  # Reminder that these cases are only for the PC with residents information provided by Public Health England.

Env_Pathogen_data_all2 <-read_csv(paste("/Code/Simulated_Salmonella_environment_",width_char,".csv",sep=""))
Env_Pathogen_data_all2 <- as.data.frame(Env_Pathogen_data_all2)

Env_Pathogen_data_all2 <-Env_Pathogen_data_all2[,-1]
colnames(Env_Pathogen_data_all2) <-c("PostCode","Date","Cases",
                                     "Global_radiation",
                                     "Maximum_air_temperature",
                                     "Mean_air_temperature",
                                     "Minimum_air_temperature",
                                     "Mean_cloud_cover",
                                     "Mean_dewpoint_temperature",
                                     "Mean_soil_temperature_1m",
                                     "Mean_soil_temperature_10_cm",
                                     "Mean_soil_temperature_20_cm",
                                     "Mean_soil_temperature_30_cm",
                                     "Mean_soil_temperature_50_cm",
                                     "Mean_surface_air_pressure",
                                     "Mean_visibility",
                                     "Mean_wind_speed",
                                     "Cumul_Precipitation",
                                     "Mean_Precipitation",
                                     "Relative_humidity",
                                     "Snow_depth_9_h",
                                     "Sunshine_duration",
                                     "daylength",
                                     "residents")

Env_Pathogen_data_all2 <- distinct(Env_Pathogen_data_all2)

Env_Pathogen_data_all2$Temperature_amplitude <- Env_Pathogen_data_all2$Maximum_air_temperature - Env_Pathogen_data_all2$Minimum_air_temperature


  # Indoor weather conditions:   
  # Calculation of the outdoor absolute humidity in kPa based on the mean oudoor temperature using the formula from [70] Mander P. 2012 (used in the paper of Verheyen 2022)
  Env_Pathogen_data_all2$Absolute_humidity <- 6.112*((exp((17.625*Env_Pathogen_data_all2$Mean_air_temperature)/(243.04+Env_Pathogen_data_all2$Mean_air_temperature))*
                                                        Env_Pathogen_data_all2$Relative_humidity*2.1674)/(273.16+Env_Pathogen_data_all2$Mean_air_temperature))
  
  Env_Pathogen_data_all2$Temperature_indoor <- ifelse (Env_Pathogen_data_all2$Mean_air_temperature < 21, 21, Env_Pathogen_data_all2$Mean_air_temperature)
  
    # Calculation of saturation vapour density (VDsat) in (gm^-3) from Nave R. 2016 (used in the paper of Verheyen 2022)
  Env_Pathogen_data_all2$Vapour_density_sat <- 5.02+0.32*Env_Pathogen_data_all2$Temperature_indoor +
    8.18*1E-3*Env_Pathogen_data_all2$Temperature_indoor^2 +
    3.12*1E-4*Env_Pathogen_data_all2$Temperature_indoor^3
  
  Env_Pathogen_data_all2$Relative_humidity_indoor <- Env_Pathogen_data_all2$Absolute_humidity/Env_Pathogen_data_all2$Vapour_density_sat


  
################### Weather data ###################
  # Read the simulation of the weather factors for every Post Code averaged for the estimated delay in the past

Env_laboratory_data2 <-read_csv(paste("/Code/Simulated_Laboratory_",width_char,".csv",sep=""))
Env_laboratory_data2 <- as.data.frame (Env_laboratory_data2)
Env_laboratory_data2 <-Env_laboratory_data2[,-1]
colnames(Env_laboratory_data2) <-c("PostCode","Date",
                                   "Global_radiation",
                                   "Maximum_air_temperature",
                                   "Mean_air_temperature",
                                   "Minimum_air_temperature",
                                   "Mean_cloud_cover",
                                   "Mean_dewpoint_temperature",
                                   "Mean_soil_temperature_1m",
                                   "Mean_soil_temperature_10_cm",
                                   "Mean_soil_temperature_20_cm",
                                   "Mean_soil_temperature_30_cm",
                                   "Mean_soil_temperature_50_cm",
                                   "Mean_surface_air_pressure",
                                   "Mean_visibility",
                                   "Mean_wind_speed",
                                   "Cumul_Precipitation",
                                   "Mean_Precipitation",
                                   "Relative_humidity",
                                   "Snow_depth_9_h",
                                   "Sunshine_duration",
                                   "daylength",
                                   "residents")

Env_laboratory_data2 <-distinct(Env_laboratory_data2) 

Env_laboratory_data2$Temperature_amplitude <- Env_laboratory_data2$Maximum_air_temperature - Env_laboratory_data2$Minimum_air_temperature


  # Indoor weather conditions:  
  # Calculation of the outdoor absolute humidity in kPa based on the mean oudoor temperature using the formula from Mander P. 2012 (used in the paper of Verheyen 2022)
  Env_laboratory_data2$Absolute_humidity <- 6.112*((exp((17.625*Env_laboratory_data2$Mean_air_temperature)/(243.04+Env_laboratory_data2$Mean_air_temperature))*
                                                      Env_laboratory_data2$Relative_humidity*2.1674)/(273.16+Env_laboratory_data2$Mean_air_temperature))
  
  Env_laboratory_data2$Temperature_indoor <- ifelse (Env_laboratory_data2$Mean_air_temperature < 21, 21, Env_laboratory_data2$Mean_air_temperature)
  
  # Calculation of saturation vapour density (VDsat) in (gm^-3) from Nave R. 2016 (used in the paper of Verheyen 2022)
  Env_laboratory_data2$Vapour_density_sat <- 5.02+0.32*Env_laboratory_data2$Temperature_indoor +
    8.18*1E-3*Env_laboratory_data2$Temperature_indoor^2 +
    3.12*1E-4*Env_laboratory_data2$Temperature_indoor^3
  
  Env_laboratory_data2$Relative_humidity_indoor <- Env_laboratory_data2$Absolute_humidity/Env_laboratory_data2$Vapour_density_sat



  # Select the records for the years of interest. We remove the rows for which there is no information for any of the weather variables. Dates match available data for residents catchment areas
  Env_Pathogen_data_all2 <-(subset(Env_Pathogen_data_all2, year(as.Date(Env_Pathogen_data_all2$Date))>=first_year & year(as.Date(Env_Pathogen_data_all2$Date))<=last_year))
  Env_laboratory_data_all2 <-(subset(Env_laboratory_data2, year(as.Date(Env_laboratory_data2$Date))>=first_year & year(as.Date(Env_laboratory_data2$Date))<=last_year))




################### Catchment areas ###################

catchment_population_df <-read.csv(paste("../Data/Catchment areas/Sum_ByLab_1987_2016.csv",sep=""))

colnames(catchment_population_df) <-c("PostCode","col2","col3",
                                      "residents_1987",
                                      "residents_1988",
                                      "residents_1989",
                                      "residents_1990",
                                      "residents_1991",
                                      "residents_1992",
                                      "residents_1993",
                                      "residents_1994",
                                      "residents_1995",
                                      "residents_1996",
                                      "residents_1997",
                                      "residents_1998",
                                      "residents_1999",
                                      "residents_2000",
                                      "residents_2001",
                                      "residents_2002",
                                      "residents_2003",
                                      "residents_2004",
                                      "residents_2005",
                                      "residents_2006",
                                      "residents_2007",
                                      "residents_2008",
                                      "residents_2009",
                                      "residents_2010",
                                      "residents_2011",
                                      "residents_2012",
                                      "residents_2013",
                                      "residents_2014",
                                      "residents_2015",
                                      "residents_2016")


  # Reshape the residents information to make a data frame with the postcode, year and number of residents
  years_col <-c()
  residents_col <-c()
  
  for (h in c(1987:2016)){
    years_col <- rbind(years_col, cbind(rep(as.character(h), times=length(catchment_population_df$PostCode))))
    h_internal <-h-1987+4 # this takes the residents info column of the year of interest
    residents_col <-rbind(residents_col, cbind(catchment_population_df[, h_internal]))
  }
  
  population_df <-data.frame(as.character(rep(catchment_population_df$PostCode, times=30)),c(years_col),c(residents_col))
  colnames(population_df) <-c("PostCode","Year","residents_Lab")
  
  # Result: population_df data frame with the number of residents per year per PC for comparison with the model



#### 2. Calculate the detrend value of cases to minimize the inter-annual variability ####

Env_Pathogen_data0 <- ddply (Env_Pathogen_data_all2, ~Date, summarise, cases=sum(Cases)) # cases at national level per day

salmonella_data_national<-Env_Pathogen_data0[order(as.Date(Env_Pathogen_data0$Date)),]       # sort from least recent to most recent 
dt0<-1

salmonella_data_national$Date <- as.Date (salmonella_data_national$Date, "%Y-%m-%d")
colnames(salmonella_data_national) [1] <- "date"




################### Detrend incidence data ###################

salmonella_data_national$cases<-as.numeric(salmonella_data_national$cases)
sal <- salmonella_data_national

    # get the average number of cases for the last 5 years of study (2012-2016, inclusive), to be used in the reconstruction
    sal2 <- sal[sal[,1]>="2012-01-01",]
    mean_cases_2011_16 <- mean(sal2$cases) # this has the risk of not accounting for days without cases = 15.39245
    
    national_detrended <-detrend.series(y=sal$cases, y.name = "Raw reported cases", make.plot = TRUE,
                                        method = c("Spline"),return.info=TRUE)
    

sal_df <-data.frame(sal$date, sal$cases, national_detrended$series, national_detrended$curves)
colnames(sal_df)<-c("date","cases","detrend","curve")
    
salmonella_cases<- Env_Pathogen_data_all2 # renaming to match previous code
colnames(salmonella_cases)[2]<-"date"


salmonella_cases_pc<-merge(salmonella_cases, sal_df, by="date")  # !careful with this dataframe. Cases corresponds to the cases per PC, and cases to the sum of Cases per day.

# Assumption: applying the national curve data to a smaller PC area is similar to calculating the curve of the specific PC and takes less computing power.
salmonella_cases_pc$Cases_det<-salmonella_cases_pc$Cases/salmonella_cases_pc$curve # correct to use Cases here
salmonella_cases_pc$Cases_det_adjusted <- salmonella_cases_pc$Cases_det * mean_cases_2011_16 # mean daily cases equals #15.39245

salmonella_cases_detrend<-salmonella_cases_pc




#### 3.Conditional probability Uniform#####
# Analysis was done following an uniform division of the range of the environmental variables, independently of the number of observations

##### Select relevant weather combinations: ######

# --Options:
# "Maximum_air_temperature"
# "Mean_air_temperature"
# "Mean_dewpoint_temperature"
# "Temperature_amplitude"
# "Temperature_indoor"

# "Relative_humidity"
# "Mean_Precipitation"
# "Cumul_Precipitation"
# "Relative_humidity_indoor"
 
# "daylength"
# "Global_radiation"
# "Sunshine_duration"
# "Mean_surface_air_pressure"

variable_1 <- "Mean_air_temperature"
variable_2 <- "Relative_humidity"
variable_3 <- "daylength"


      # Select the variables of interest:
      
      Env_Pathogen_data <- salmonella_cases_detrend [,c("PostCode","date","Cases_det_adjusted", variable_1, variable_3, variable_2,"residents")] # UK adjusted!
      Env_Pathogen_data <- Env_Pathogen_data[complete.cases(Env_Pathogen_data),] # Remove NA. Removes instances with no weather records                                                    
      
      Env_laboratory_data <- Env_laboratory_data_all2[,c("PostCode","Date", variable_1, variable_3, variable_2, "residents")]
      Env_laboratory_data <- Env_laboratory_data[complete.cases(Env_laboratory_data),] # Remove NA. Removes instances with no weather records
      


      ################### Divide the domain of weather variable in bins of delta size for classifying the weather in uniform intervals
      
      if (variable_3=="daylength" | variable_3== "Sunshine_duration"){
        
        delta_var3 <- 1
        breaks_var3 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var3), 0)
        
      } else if (variable_3=="Mean_surface_air_pressure") {
        
        delta_var3 <- 5
        breaks_var3 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var3), 0)
        
      } else if (variable_3=="Mean_Precipitation") {
        
        delta_var3 <- 0.5
        breaks_var3 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var3), 2)
        
      } else if (variable_3=="Mean_air_temperature") {
        
        delta_var3 <-1 # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature
        breaks_var3 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                                   max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 0)
        
      } else if (variable_3=="Relative_humidity") {
        
        delta_var3 <- 1 # permutation checks. idem for all the deltas below for relative humidity
        breaks_var3 <-round(seq(min((Env_laboratory_data[[variable_2]])), 
                               max(na.omit(Env_laboratory_data[[variable_2]])), by=delta_var2), 0) 
        
      } else {
        delta_var3 <- 1000 # global radiation. The values of radiation go between and seem to change more or less in 1000 units.
        breaks_var3 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])), by=delta_var3), 0)
        
      }
      
      
      
      if (variable_1=="Temperature_indoor"){
        
        delta_var1 <-0.5 
        breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                                   max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 2)
        
      } else if (variable_1=="Mean_surface_air_pressure") {
        
        delta_var1 <- 5
        breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var1), 0)
        
      } else if (variable_1=="Mean_air_temperature") {
        
        delta_var1 <-1 # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature
        breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                                   max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 0)
        
      }else if (variable_1=="daylength") {
        
        delta_var1 <- 1
        breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var3), 0)
        
      }else { # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature
      
      delta_var1 <-1 # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature
      breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                                 max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 0)
      }
      
      
      
      if (variable_2=="Relative_humidity"){
        
        delta_var2 <-1
        breaks_var2 <-round(seq(min((Env_laboratory_data[[variable_2]])), 
                               max(na.omit(Env_laboratory_data[[variable_2]])), by=delta_var2), 0) 
        
      } else if (variable_2=="Mean_Precipitation" | variable_2=="Cumul_Precipitation") {
        
        delta_var2 <- 0.5
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_2]])), 
                               max(na.omit(Env_laboratory_data[[variable_2]]) ), by=delta_var2), 2) 
        
      } else if (variable_2== "Relative_humidity_indoor") {
        
        delta_var2 <- 0.1
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_2]])), 
                               max(na.omit(Env_laboratory_data[[variable_2]]) ), by=delta_var2), 2) 
        
      } else if (variable_2=="Mean_dewpoint_temperature") {
        
        delta_var2 <- 1
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_2]])), 
                               max(na.omit(Env_laboratory_data[[variable_2]]) ), by=delta_var2), 0) 
        
      } else if (variable_2=="Sunshine_duration") {
        
        delta_var2 <- 1
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var2), 0)
      } else if (variable_2=="Mean_surface_air_pressure") {
        
        delta_var2 <- 5
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var2), 0)
      } else if (variable_2=="Mean_air_temperature") {
        
        delta_var2 <-1 # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                                   max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 0)
        
      }else if (variable_2=="daylength") {
        
        delta_var2 <- 1
        breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                                     max(na.omit(Env_laboratory_data[[variable_3]])) , by=delta_var3), 0)
      }
      
      
      # Create a data frame to fill with
      estratification_weather_cases <-data.frame(character(), character(), character(), numeric(), numeric(), numeric())
      colnames(estratification_weather_cases) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")
      
      
      
      ################### Conditional incidence stratification:
      inner_select3_min <-1 
      inner_select3_max <-length(breaks_var3)
      
      
      #gives the index of the categorical "bin" that the value of the variable belongs to for every nested weather variable:
      
      for (inner_select3 in c(inner_select3_min:inner_select3_max)){
        
        # below I round the values of the variables of interest in Env_Pathogen_data to the closest break value and get the index value of the nearest break
        
        wt <- match.closest (Env_Pathogen_data[[variable_3]], breaks_var3) # new way to round to the nearest and get the index
        ww <-which(wt==inner_select3)                                                                      # number of events when the variable coincides with the category
        Env_Pathogen_data_z <-Env_Pathogen_data[ww, ]                                            # subset of all cases where daylength is between the category above
        
        # same for laboratory data set
        wt_t <- match.closest (Env_laboratory_data[[variable_3]], breaks_var3)
        ww_t <-which(wt_t==inner_select3)
        Env_laboratory_data_z <-Env_laboratory_data[ww_t, ]
                                                                                                 # subsetting the same as above including humidity delimitation included
        {
          {
            inner_select2_min <- 1
            inner_select2_max <-length(breaks_var2)
            
            for (inner_select2 in c(inner_select2_min:inner_select2_max)){ # from the bins resulting above, find another variable stratification
              
              wt_y <- match.closest (Env_Pathogen_data_z[[variable_2]], breaks_var2)
              ww_y <-which(wt_y==inner_select2)
              Env_Pathogen_data_y <-Env_Pathogen_data_z[ww_y, ]
              
              wt_y_t <- match.closest (Env_laboratory_data_z[[variable_2]], breaks_var2)
              ww_y_t <-which(wt_y_t==inner_select2)
              Env_laboratory_data_y <-Env_laboratory_data_z[ww_y_t, ] 
              
             
              {
                inner_select1_min <-1
                inner_select1_max <- length(breaks_var1)
                
                for (inner_select1 in c(inner_select1_min:inner_select1_max)){
                  
                  wt_x <- match.closest (Env_Pathogen_data_y[[variable_1]], breaks_var1)
                  ww_x <-which(wt_x==inner_select1)
                  Pathogen_stratified <-Env_Pathogen_data_y[ww_x,] 
                  
                  
                  wt_x_t <- match.closest (Env_laboratory_data_y[[variable_1]], breaks_var1)
                  ww_x_t <-which(wt_x_t==inner_select1)
                  Laboratory_stratified <-Env_laboratory_data_y[ww_x_t,] 
                  
                  Total_cases <-sum((as.numeric(na.omit(Pathogen_stratified$Cases_det_adjusted))))
                  residents <-sum((as.numeric(na.omit(Pathogen_stratified$residents)))) # number of people for areas with a recorded case for specific values of the 3 weather variables.
                  residents_tot <-sum((as.numeric(na.omit(Laboratory_stratified$residents)))) #number of people exposed to these weather conditions, irrespectively whether there was a case found or not. Will be used to calculate the incidence
                  
                  
                  data_df <-data.frame (
                    breaks_var3[inner_select3],
                    breaks_var2[inner_select2],
                    breaks_var1[inner_select1],
                    Total_cases,
                    residents,
                    residents_tot
                  )
                  
                  colnames(data_df) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")
                  estratification_weather_cases <-rbind(estratification_weather_cases, data_df) 
                  print(c(inner_select1, inner_select2,inner_select3, Total_cases)) 
                }
              }
            }
          }
        }
      } # output: estratification_weather_cases
      
      write.table(estratification_weather_cases, paste(variable_3,"_",variable_2,"_",variable_1,"_",width_char,
                                      "_Simulated_for_rec_UNIFORM_narrowedvalues_DETREND.csv",sep=""), col.names = NA, row.names = TRUE , sep = ",",eol = "\n") # will be used in code No.2
      
      ############ Reconstruction from Conditional probability
      # The code does look at how the risk of salmonellosis in humans depends on environmental variables
      Env_Pathogen_data$PostCode <-as.factor(Env_Pathogen_data$PostCode)
      Env_laboratory_data$year <-year(as.Date(Env_laboratory_data$Date))
      
      
      modelled_cases<-c()
      
      
      #----insert of breaks ---

      dt_var3 = data.table(breaks_var3, val = breaks_var3)
      dt_var3 <- na.omit (dt_var3)
      setattr(dt_var3, "sorted", "breaks_var3")
      dt_var3.2 <- dt_var3[J(Env_laboratory_data[[variable_3]]), roll = "nearest"] 
      
      weather_bins <- as.data.frame(Env_laboratory_data[,c("PostCode","Date", "residents")])
      weather_bins$break_var3<- dt_var3.2$val
      
      
      # Relative humidity
      dt_var2 = data.table(breaks_var2, val = breaks_var2)
      dt_var2 <- na.omit (dt_var2)
      setattr(dt_var2, "sorted", "breaks_var2")
      dt_var2.2 <- dt_var2[J(Env_laboratory_data[[variable_2]]), roll = "nearest"]
      
      
      weather_bins$break_var2<- dt_var2.2$val
      
      
      
      # Maximum air temperature
      dt_var1 = data.table(breaks_var1, val = breaks_var3)
      dt_var1 <- na.omit (dt_var1)
      setattr(dt_var1, "sorted", "breaks_var1")
      dt_var1.2  <- dt_var1[J(Env_laboratory_data[[variable_3]]), roll = "nearest"]
      
      weather_bins$break_var1<- dt_var1.2$val
      colnames(weather_bins) <- c("PostCode","Date", "residents", variable_3, variable_2, variable_1)
      
      
      CI_per_weather_bin <- left_join (weather_bins,estratification_weather_cases, by= c(variable_2, variable_1, variable_3))
      CI_per_weather_bin <- na.omit(CI_per_weather_bin)
      CI_per_weather_bin <- distinct(CI_per_weather_bin)
      CI_per_weather_bin <-CI_per_weather_bin[order(as.Date(CI_per_weather_bin$Date)),]
      
      #colnames(variable_df_dis)[6] ="residents"
      colnames(CI_per_weather_bin)[8] ="residents" # not used. No. of people exposed to the same weather conditions in an area with at least one case.
      colnames(CI_per_weather_bin)[9] ="residents_total" 
      
      # calculation of all the potential incidence based on the weather conditions of the area. The counts used as numerator correspond to the simulation of cases per weather condition.
      CI_per_weather_bin$incidence <-CI_per_weather_bin$counts/CI_per_weather_bin$residents_total 
      lambda <-CI_per_weather_bin$incidence
      comp_cases <-lambda*CI_per_weather_bin$residents.x # estimation of the "real" incidence based on the possible cases due to weather and the residents number in the catchment area
      
      modelled_cases <-data.frame(CI_per_weather_bin$Date, comp_cases, lambda, CI_per_weather_bin$PostCode)
      colnames(modelled_cases) <-c("Date","Cases","Lambda","Lab")

      
      
      ####### Plot simulated reconstruction
      real_cases_all <-salmonella_data_national_det 
      colnames(real_cases_all) <- c("date", "Cases_det_adjusted")
      
      real_cases_all$date <- as.Date(real_cases_all$date, format="%Y-%m-%d")
      

      
      modelled_cases$Date <- as.Date(modelled_cases$Date)
      colnames(modelled_cases) <-c("Date","Cases","Lambda","PostCode")
      modelled_cases$Year <-year(modelled_cases$Date)
      modelled_cases <-modelled_cases[order(modelled_cases$Date, na.last = NA),]
      
      modelled_cases <-merge(modelled_cases, population_df, by=c('Year','PostCode')) # Add residents information to the time series.
      modelled_cases <-na.omit(modelled_cases)
      
      
      ### sum cases per PC 
      modelled_cases <- aggregate (Cases ~ Date, modelled_cases, sum) 
      
      ts_roll_mean2 <- modelled_cases %>% mutate (rolling_mean2= slider::slide_mean(Cases, before=3 , after=3 )) %>% ungroup()
      modelled_cases$rolling_mean <-ts_roll_mean2$rolling_mean2  
      
      modelled_cases_summary <-data.frame(modelled_cases$Date, modelled_cases$rolling_mean, rep("Model",times=length(modelled_cases[,1]))) 
      colnames(modelled_cases_summary)<-c("Date","Mean","source")
      
      ### Real cases
      ## add all cases per date
      real_cases_national <- aggregate (Cases_det_adjusted ~ date, real_cases_all, sum) 
      
      real_cases_mean <-real_cases_national %>% mutate (rolling_mean= slider::slide_mean(Cases_det_adjusted, before=3 , after=3 ))
      real_cases_national$rolling_mean <-real_cases_mean$rolling_mean 
      
      real_cases_summary<-cbind(real_cases_national, rep("Cases",times=length(real_cases_national[,1])))
      real_cases_summary <- real_cases_summary[,-2] # keep rolling mean of the detrended cases
      colnames(real_cases_summary) <-c("Date","Mean","source")              
      
      real_cases_summary$Date<-as.factor(real_cases_summary$Date)
      modelled_cases_summary$Date <- as.factor(modelled_cases_summary$Date)
      
      modelled_cases_all <-rbind(modelled_cases_summary, real_cases_summary)
      
      
      
      ######## Plots #####
      

      modelled_cases_all$Date <- as.Date(modelled_cases_all$Date)
      
      modelled_cases_plot <-ggplot(modelled_cases_all, aes(x=Date, y=Mean, colour=source, group=2))+
        
                                geom_line(linewidth=0.75) + 
                                xlab("Date") + ylab("Salmonellosis Cases")
                                  
      
      # Save the plot

      
      
      ############## Average per day of the year 
      modelled_cases$yday <-as.factor(yday(modelled_cases$Date))
      
      modelled_cases_average1 <-ddply(modelled_cases,~yday, summarise, mean=mean(rolling_mean)) # mean of all cases for the same day of the year for the 17 years of study
      modelled_cases_quantile1 <-ddply(modelled_cases,~yday, function (x) quantile(x$rolling_mean , c(.25,.5,.75))) # quantiles of cases for all the years for the same day of the year
      modelled_cases_average2 <-cbind(modelled_cases_average1,modelled_cases_quantile1[,-1])
      
      
      real_cases_national$yday <-as.factor(yday(as.Date(as.character(real_cases_national$date))))
      real_cases_average1 <-ddply(real_cases_national, ~yday, summarise, mean=mean(rolling_mean))
      real_cases_quantile1 <-ddply(real_cases_national, ~yday,function (x) quantile(x$rolling_mean, c(.25,.5,.75)))
      real_cases_average2 <-cbind(real_cases_average1, real_cases_quantile1[,-1])

      real_cases_average <- data.frame(real_cases_average2[,1], real_cases_average2[,c(2:5)])
      
      
      df1 <-data.frame(real_cases_average,rep("Cases",times=length(real_cases_average[,1])))
      colnames(df1) <-c("Day","Mean","f_quant","median","s_quant","source")
      df2 <-data.frame(modelled_cases_average2,rep("Model",times=length(modelled_cases_average2[,1])))
      colnames(df2) <-c("Day","Mean","f_quant","median","s_quant","source")
      
      average_data <-rbind(df1,df2)
      average_data$Day <-as.numeric(average_data$Day)
      
 
      
      variable_1.char <- stri_replace_all_fixed(variable_1, "_", " ")
      variable_2.char <- stri_replace_all_fixed(variable_2, "_", " ")
      variable_3.char <- stri_replace_all_fixed(variable_3, "_", " ")
      if(variable_3 == "daylength"){
        variable_3.char <- "Day length"
      }
      
      yearly_average <-ggplot(average_data, aes(x=Day, y=Mean, colour=source))+
        geom_line(linewidth=2)+ theme_bw(20)+
        geom_ribbon(aes(ymin=f_quant, ymax=s_quant, fill=source), alpha=0.15) +
        xlab("Date") + ylab("Salmonellosis Cases") +
        

        ggtitle(paste0(variable_1.char,", ", variable_2.char,", ", variable_3.char))
      
      

      
      # Save the plot



#### 4. Conditional probability Quantile ####
# The code looks at how the risk of salmonellosis in humans on a quantile division of the range of the environmental variables


############# Manually do the weather combinations of interest identified in the reconstruction of the uniform quantiles #############
# Change here for doing other combinations

# 1. Mean_dewpoint_temperature, Mean_Precipitation, daylength
# 2. Mean_air_temperature, Mean_Precipitation, daylength
# 3. Mean_dewpoint_temperature, Relative_humidity, daylength
# 4. Mean_air_temperature, Mean_wind_speed, daylength
# 5. Mean_dewpoint_temperature, Mean_surface_air_pressure, Mean_Precipitation
# 6. Maximum_air_temperature, Mean_Precipitation, Global_radiation
# 7. "Mean_air_temperature", "Mean_Precipitation", "Mean_wind_speed"  Additional combination that did NOT work well to see the CI plot.

variable_1 <- "Mean_air_temperature"
variable_2 <- "Relative_humidity"
variable_3 <- "daylength"
      
# Select the variables of interest:
      
Env_Pathogen_data <- salmonella_cases_detrend [,c("PostCode","date","Cases_det_adjusted", variable_1, variable_3, variable_2,"residents")]
Env_Pathogen_data <- Env_Pathogen_data[complete.cases(Env_Pathogen_data),] # Remove NA. Removes instances with no weather records                                                    
      
Env_laboratory_data <- Env_laboratory_data_all2[,c("PostCode","Date",variable_1, variable_3, variable_2, "residents")]
Env_laboratory_data <- Env_laboratory_data[complete.cases(Env_laboratory_data),] # Remove NA. Removes instances with no weather records
      


############# Create the breaks for a quantile division of the bins: detrended cases #############
  # The bins collect the value of bins_size_var3 (e.g. 25% of the data). The size of the bins is different but the number of observations is the same. The bins depend one from another weather variable to keep the number of observations per bin
  # I created 4 quantiles where there is not much diversity of values (i.e. daylength, RH...), and smaller for higher diversity (Radiation, Temperature..)

bins_size_var1 <-0.1 # temp, dewpoint
bins_size_var2 <-0.25 # precip, widspeed, air pressure
bins_size_var3 <-0.25 # DL, global radiation, windspeed

#bins_size_var2 <- 1 #for removing precipitation from the analysis

# Divisions for weather factor 3

breaks_var3 <-function (variable_3, bins_size_var3)
{
  breaks_var3 <-as.numeric(quantile(na.omit(Env_laboratory_data[[variable_3]]), probs=seq(0, 1, by=bins_size_var3), na.rm=TRUE))
  breaks_var3[length(breaks_var3)] <-ceiling(as.numeric(quantile(na.omit(Env_laboratory_data[[variable_3]]), probs=seq(0,1, by=bins_size_var3), na.rm=TRUE)))[length(breaks_var3)]
  breaks_var3[1] <-floor(as.numeric(quantile(na.omit(Env_laboratory_data[[variable_3]]), probs=seq(0,1, by=bins_size_var3), na.rm=TRUE)))[1]
  
  return(breaks_var3)
}


#	Divisions for weather factor 2

breaks_var2 <-function (variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)
{
  wt <-(findInterval(Env_laboratory_data[[variable_3]], breaks_var3(variable_3, bins_size_var3)))
  ww <-which(wt==inner_select3)
  Env_laboratory_data_some <-Env_laboratory_data[ww, ]
  
  
  if (length(Env_laboratory_data_some[ ,1])!=0) {
    
    breaks_var2 <-as.numeric(quantile(na.omit(Env_laboratory_data_some[[variable_2]]), probs=seq(0,1, by=bins_size_var2), na.rm=TRUE))
    breaks_var2[length(breaks_var2)] <-ceiling(as.numeric(quantile(na.omit(Env_laboratory_data_some[[variable_2]]), probs=seq(0, 1, by=bins_size_var2), na.rm=TRUE)))[length(breaks_var2)]
    breaks_var2[1] <-floor(as.numeric(quantile(na.omit(Env_laboratory_data_some[[variable_2]]), probs=seq(0, 1, by=bins_size_var2), na.rm=TRUE)))[1]
    
  }else{
    
    breaks_var2 <-c()
  }
  
  return(breaks_var2)
}

#	Divisions for weather factor 1

breaks_var1 <-function(variable_3, variable_2, variable_1, bins_size_var3, bins_size_var2, bins_size_var1, inner_select3, inner_select2)
{
  wt <-(findInterval(Env_laboratory_data[[variable_3]],breaks_var3(variable_3, bins_size_var3)))
  ww <-which(wt==inner_select3)
  Env_laboratory_data_some <-Env_laboratory_data[ww, ]
  
  if(is.na(breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)[inner_select2])=='FALSE')
  {
    wt2 <-(findInterval(Env_laboratory_data_some[[variable_2]], breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)))
    ww2 <-which(wt2==inner_select2) 
    Env_laboratory_data_some2 <-Env_laboratory_data_some[ww2, ]
    
    if (length(Env_laboratory_data_some2[ ,1])!=0) 
    {
      
      breaks_var1 <-as.numeric(quantile(na.omit(Env_laboratory_data_some2[[variable_1]]), probs=seq(0,1, by=bins_size_var1), na.rm=TRUE))
      breaks_var1[length(breaks_var1)] <-ceiling(as.numeric(quantile(na.omit(Env_laboratory_data_some2[[variable_1]]), probs=seq(0, 1, by=bins_size_var1), na.rm=TRUE)))[length(breaks_var1)]
      breaks_var1[1] <-floor(as.numeric(quantile(na.omit(Env_laboratory_data_some2[[variable_1]]), probs=seq(0,1, by=bins_size_var1), na.rm=TRUE)))[1]
      
    } 
    else 
    { 
      breaks_var1 <-c() 
    }
  }  
  else 
  {
    breaks_var1<-c()
  }
  return(breaks_var1)
}




################# Create a data frame to fill:

Conditional_incidence_quantiles <-data.frame(character(), character(),character(),numeric(),numeric(),numeric())
colnames(Conditional_incidence_quantiles) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")

residents_i_var <-0
residents_universal <-0

################### Divide the domains of the variables in bins according to quantiles
  #Grouping the values according to the breaks of each weather variable (1,2,3)

inner_select3_min <-1 
inner_select3_max <-length(breaks_var3(variable_3, bins_size_var3))-1

# Grouping the values according to the breaks of each weather variable (1,2,3)

for (inner_select3 in c(inner_select3_min:inner_select3_max)){ 
 
  wt <- findInterval (Env_Pathogen_data[[variable_3]], breaks_var3(variable_3, bins_size_var3))
  ww <-which(wt==inner_select3)                                                                      # number of events when the variable coincides with the category
  Env_Pathogen_data_z <-Env_Pathogen_data[ww, ]                                            # subset of all cases where daylength is between the category above
  
  
  wt_t <- findInterval (Env_laboratory_data[[variable_3]], breaks_var3(variable_3, bins_size_var3))
  ww_t <-which(wt_t==inner_select3)
  Env_laboratory_data_z <-Env_laboratory_data[ww_t, ]
  

    if (length(breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))!=0)
    {
      inner_select2_min <- 1
      inner_select2_max <-length(breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))-1
      
      for (inner_select2 in c(inner_select2_min:inner_select2_max)){ # from the bins resulting above, find another variable stratification
        
        
        wt_y <- findInterval (Env_Pathogen_data_z[[variable_2]], breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))
        ww_y <-which(wt_y==inner_select2)
        Env_Pathogen_data_y <-Env_Pathogen_data_z[ww_y, ]
        
        
        wt_y_t <- findInterval (Env_laboratory_data_z[[variable_2]], breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))  
        ww_y_t <-which(wt_y_t==inner_select2)
        Env_laboratory_data_y <-Env_laboratory_data_z[ww_y_t, ] 
        
        
        if (length(breaks_var1(variable_3, variable_2, variable_1, bins_size_var3, bins_size_var2, bins_size_var1, inner_select3, inner_select2))!=0)
        {
          inner_select1_min <-1
          inner_select1_max <- length(breaks_var1(variable_3,variable_2,variable_1,bins_size_var3,bins_size_var2,bins_size_var1,inner_select3,inner_select2))-1
          
          for (inner_select1 in c(inner_select1_min:inner_select1_max)){

            wt_x <- findInterval (Env_Pathogen_data_y[[variable_1]], breaks_var1(variable_3,variable_2,variable_1,bins_size_var3,bins_size_var2,bins_size_var1,inner_select3,inner_select2))
            ww_x <-which(wt_x==inner_select1)
            Yt1 <-Env_Pathogen_data_y[ww_x, ] #where there are cases for specific values of the 3 weather variables. residents is the number of people exposed to these situations during the entire duration of the dataset. It is the number of times the people is exposed to this caracteristics, only when a case is found. It can count a same catchment area twice. Residents total is the same irrespectively of the presence of a case

            
            wt_x_t <- findInterval (Env_laboratory_data_y[[variable_1]], breaks_var1(variable_3,variable_2,variable_1,bins_size_var3,bins_size_var2,bins_size_var1,inner_select3,inner_select2))
            ww_x_t <-which(wt_x_t==inner_select1)
            Y_tot <-Env_laboratory_data_y[ww_x_t, ] #NA in residents at Env_laboratory_data. we do not care if there is a case reported
            
            #Yt1 <- merge(Y_tot, salmonella_cases_detrend)
            
            Total_cases <-sum((as.numeric(na.omit(Yt1$Cases_det_adjusted))))
            residents <-sum((as.numeric(na.omit(Yt1$residents)))) # we may not use.
            residents_tot <-sum((as.numeric(na.omit(Y_tot$residents)))) #number of people exposed to these weather conditions. Will be used to calculate the incidence
            
            
            data_df <-data.frame (
              breaks_var3(variable_3, bins_size_var3)[inner_select3],
              breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)[inner_select2],
              breaks_var1(variable_3, variable_2, variable_1, bins_size_var3, bins_size_var2, bins_size_var1, inner_select3, inner_select2)[inner_select1],
              Total_cases,
              residents,
              residents_tot
            )
            
            colnames(data_df) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")
            Conditional_incidence_quantiles <-rbind(Conditional_incidence_quantiles, data_df) 
            print(c(inner_select1, inner_select2,inner_select3, Total_cases)) # check in the output log file if x y and z cover all the conditions established in breaks()
          }
        }
      }
    }
} # output: Conditional_incidence_quantiles



############# Plots CI Quantile #############

# In x, y, z order:
# (1.) Mean_dewpoint_temperature, Mean_Precipitation, daylength
# 1. Mean_air_temperature, Mean_Precipitation, daylength
# 2. Mean_dewpoint_temperature, Relative_humidity, daylength
# 3. Mean_air_temperature, Mean_wind_speed, daylength
# 4. Mean_dewpoint_temperature, Mean_surface_air_pressure, Mean_Precipitation
# 5. Maximum_air_temperature, Mean_Precipitation, Global_radiation
# 6. "Mean_air_temperature", "Mean_Precipitation", "Mean_wind_speed"  Additional combination that did NOT work well to see the CI plot.


variable_1 <- "Maximum_air_temperature"
variable_2 <- "Mean_Precipitation"
variable_3 <- "Global_radiation"




Conditional_incidence_quantiles$incidence <-Conditional_incidence_quantiles$counts/Conditional_incidence_quantiles$residents_tot

# Normalized incidence per 10,000,000 inhabitants to normalize our plots later
norm <- 1/1e7 # time steps included in the data itself, so no need to divide the number of days
labelling_char<-c()



#### Plots 3 variables quantiles ####


  for (i_light in c(1:length(unique(Conditional_incidence_quantiles[[variable_3]])))){

    light_label <-round(as.numeric(unique(Conditional_incidence_quantiles[[variable_3]])))
    
    if (i_light < length(light_label)){ 
      
      light_label_char <-paste(variable_3.char," [",light_label[i_light],"-",light_label[i_light+1],") ",units.z, sep="")
    } else {
      light_label_char<-paste(variable_3.char," >",light_label[i_light]," ",units.z, sep="")
    }
      
    
    wt <-which(Conditional_incidence_quantiles[[variable_3]]==unique(Conditional_incidence_quantiles[[variable_3]])[i_light])
    Conditional_incidence_df <-Conditional_incidence_quantiles[wt,] # Here we have selected the rows of the df corresponding with each of the 4 daylength bins (n=40 rows)
    Conditional_incidence_df <-Conditional_incidence_df[order(Conditional_incidence_df[[variable_2]]),]     
    

      Conditional_incidence_df$incidence <-Conditional_incidence_df$incidence/norm
      
      y=Conditional_incidence_df$counts 
      n=Conditional_incidence_df$residents_tot # input data as frequencies
      
      pw=(y + 2)/(n + 4) # Wilson's point estimate
      Wald_adj_lower <- (pw-2*sqrt(pw*(1-pw)/n))/norm # Lower conf. limit. by using 2 we are "approximating". Opted to use this method by suggestion of the textbook
      Wald_adj_upper <- (pw+2*sqrt(pw*(1-pw)/n))/norm  # Upper conf. limit 

      
      Conditional_incidence_df$conf_minus <-Wald_adj_lower
      Conditional_incidence_df$conf_plus <-Wald_adj_upper
      Conditional_incidence_df <-na.omit(Conditional_incidence_df)
      
      
      Conditional_incidence_df[[variable_2]] <-round(Conditional_incidence_df[[variable_2]])
      

      #Plot
      
      Conditional_incidence_plot  <- ggplot(Conditional_incidence_df, 
                                aes(x=.data[[variable_1]], y=incidence, color=factor(.data[[variable_2]]), fill=factor(.data[[variable_2]])))+
        
        geom_ribbon(aes(ymin=conf_minus, ymax=conf_plus), alpha=0.35)+
        
        ylab("Conditional Incidence of salmonellosis Cases (per 10 million)")+
        xlab(variable_1.char)+ 
        ggtitle(light_label_char)
          
      # Save plot

    }


####2.4 Plots 2 variables quantiles####
variable_1 <- "Mean_air_temperature"
#variable_2 <- "Mean_Precipitation" #not shown
variable_3 <- "daylength"

# Select the variables of interest:

Env_Pathogen_data <- salmonella_cases_detrend [,c("PostCode","date","Cases_det_adjusted", variable_1, variable_3, variable_2,"residents")]
Env_Pathogen_data <- Env_Pathogen_data[complete.cases(Env_Pathogen_data),] # Remove NA. Removes instances with no weather records                                                    

Env_laboratory_data <- Env_laboratory_data_all2[,c("PostCode","Date",variable_1, variable_3, variable_2, "residents")]
Env_laboratory_data <- Env_laboratory_data[complete.cases(Env_laboratory_data),] # Remove NA. Removes instances with no weather records



######################### Create the breaks for a quantile division of the bins: detrended cases
# The bins collect the value of bins_size_var3 (e.g. 25% of the data). The size of the bins is different but the number of observations is the same. The bins depend one from another weather variable to keep the number of observations per bin
# I created 4 quantiles where there is not much diversity of values (i.e. daylength, RH...), and smaller for higher diversity (Radiation, Temperature..)

bins_size_var1 <-0.1 # temp, dewpoint
#bins_size_var2 <- 1 #for removing precipitation from the analysis
bins_size_var3 <-0.25 # DL, global radiation



# Divisions for weather factor 3

breaks_var3 <-function (variable_3, bins_size_var3)
{
  breaks_var3 <-as.numeric(quantile(na.omit(Env_laboratory_data[[variable_3]]), probs=seq(0, 1, by=bins_size_var3), na.rm=TRUE))
  breaks_var3[length(breaks_var3)] <-ceiling(as.numeric(quantile(na.omit(Env_laboratory_data[[variable_3]]), probs=seq(0,1, by=bins_size_var3), na.rm=TRUE)))[length(breaks_var3)]
  breaks_var3[1] <-floor(as.numeric(quantile(na.omit(Env_laboratory_data[[variable_3]]), probs=seq(0,1, by=bins_size_var3), na.rm=TRUE)))[1]
  
  return(breaks_var3)
}



#	Divisions for weather factor 2

breaks_var1 <-function(variable_3, variable_1, bins_size_var3, bins_size_var1, inner_select3)
{
  wt <- match.closest(Env_laboratory_data[[variable_3]],breaks_var3(variable_3, bins_size_var3))
  ww <-which(wt==inner_select3)
  Env_laboratory_data_some <-Env_laboratory_data[ww, ]
 
  breaks_var1 <-as.numeric(quantile(na.omit(Env_laboratory_data_some[[variable_1]]), probs=seq(0,1, by=bins_size_var1), na.rm=TRUE))
  breaks_var1[length(breaks_var1)] <-ceiling(as.numeric(quantile(na.omit(Env_laboratory_data_some[[variable_1]]), probs=seq(0, 1, by=bins_size_var1), na.rm=TRUE)))[length(breaks_var1)]
  breaks_var1[1] <-floor(as.numeric(quantile(na.omit(Env_laboratory_data_some[[variable_1]]), probs=seq(0,1, by=bins_size_var1), na.rm=TRUE)))[1]
    
  return(breaks_var1)
}




################# Create a data frame to fill:

Conditional_incidence_quantiles <-data.frame(character(), character(),numeric(),numeric(),numeric())
colnames(Conditional_incidence_quantiles) <-c(variable_3, variable_1,"counts","residents","residents_tot")

residents_i_var <-0
residents_universal <-0

################### Divide the domains of the variables in bins according to quantiles
#Grouping the values according to the breaks of each weather variable (1,2,3)

inner_select3_min <-1 
inner_select3_max <-length(breaks_var3(variable_3, bins_size_var3))#-1


# Grouping the values according to the breaks of each weather variable (1,2,3)

for (inner_select3 in c(inner_select3_min:inner_select3_max)){ 
  
  wt <- findInterval (Env_Pathogen_data[[variable_3]], breaks_var3(variable_3, bins_size_var3))
  ww <-which(wt==inner_select3)                                                                      # number of events when the variable coincides with the category
  Env_Pathogen_data_z <-Env_Pathogen_data[ww, ]                                            # subset of all cases where daylength is between the category above
  
  
  wt_t <- findInterval (Env_laboratory_data[[variable_3]], breaks_var3(variable_3, bins_size_var3))
  ww_t <-which(wt_t==inner_select3)
  Env_laboratory_data_z <-Env_laboratory_data[ww_t, ]
  

      if (length(breaks_var1(variable_3, variable_1, bins_size_var3, bins_size_var1, inner_select3))!=0)
      {
        inner_select1_min <-1
        inner_select1_max <- length(breaks_var1(variable_3,variable_1,bins_size_var3,bins_size_var1,inner_select3))#-1
        
        for (inner_select1 in c(inner_select1_min:inner_select1_max)){
          
          wt_x <- findInterval (Env_Pathogen_data_z[[variable_1]], breaks_var1(variable_3,variable_1,bins_size_var3,bins_size_var1,inner_select3))
          ww_x <-which(wt_x==inner_select1)
          Yt1 <-Env_Pathogen_data_z[ww_x, ] #where there are cases for specific values of the 3 weather variables. residents is the number of people exposed to these situations during the entire duration of the dataset. It is the number of times the people is exposed to this caracteristics, only when a case is found. It can count a same catchment area twice. Residents total is the same irrespectively of the presence of a case
          
          
          wt_x_t <- findInterval (Env_laboratory_data_z[[variable_1]], breaks_var1(variable_3,variable_1,bins_size_var3,bins_size_var1,inner_select3))
          ww_x_t <-which(wt_x_t==inner_select1)
          Y_tot <-Env_laboratory_data_z[ww_x_t, ] #NA in residents at Env_laboratory_data. we do not care if there is a case reported
          
          #Yt1 <- merge(Y_tot, salmonella_cases_detrend)
          
          Total_cases <-sum((as.numeric(na.omit(Yt1$Cases_det_adjusted))))
          residents <-sum((as.numeric(na.omit(Yt1$residents)))) #we may not use.
          residents_tot <-sum((as.numeric(na.omit(Y_tot$residents)))) #number of people exposed to these weather conditions. Will be used to calculate the incidence
          
          
          data_df <-data.frame (
            breaks_var3(variable_3, bins_size_var3)[inner_select3],
            breaks_var1(variable_3, variable_1, bins_size_var3, bins_size_var1, inner_select3)[inner_select1],
            Total_cases,
            residents,
            residents_tot
          )
          
          colnames(data_df) <-c(variable_3, variable_1,"counts","residents","residents_tot")
          Conditional_incidence_quantiles <-rbind(Conditional_incidence_quantiles, data_df) 
          print(c(inner_select1, inner_select3, Total_cases)) # check in the output log file if x y and z cover all the conditions established in breaks()
        }
      }
   }

Conditional_incidence_quantiles$incidence <-Conditional_incidence_quantiles$counts/Conditional_incidence_quantiles$residents_tot




# Labels
  for (i_light in c(1:length(unique(Conditional_incidence_quantiles[[variable_3]])))){
    
    light_label <-round(as.numeric(unique(Conditional_incidence_quantiles[[variable_3]])))
    
    if (i_light < length(light_label)){ 
      
      light_label_char <-paste(variable_3.char," [",light_label[i_light],"-",light_label[i_light+1],") ",units.z, sep="")
    } else {
      light_label_char<-paste(variable_3.char," >",light_label[i_light]," ",units.z, sep="")
    }
    
    
    wt <-which(Conditional_incidence_quantiles[[variable_3]]==unique(Conditional_incidence_quantiles[[variable_3]])[i_light])
    Conditional_incidence_df <-Conditional_incidence_quantiles[wt,] # Here we have selected the rows of the df corresponding with each of the 4 daylength bins (n=40 rows)
 
      
      Conditional_incidence_df$incidence <-Conditional_incidence_df$incidence/norm
      
      
      y=Conditional_incidence_df$counts 
      n=Conditional_incidence_df$residents_tot # input data as frequencies
      
      
      pw=(y + 2)/(n + 4) # Wilson's point estimate
      Wald_adj_lower <- (pw-2*sqrt(pw*(1-pw)/n))/norm # Lower conf. limit. by using 2 we are "approximating". Opted to use this method by suggestion of the textbook
      Wald_adj_upper <- (pw+2*sqrt(pw*(1-pw)/n))/norm  # Upper conf. limit 
      
      
      Conditional_incidence_df$conf_minus <-Wald_adj_lower
      Conditional_incidence_df$conf_plus <-Wald_adj_upper
      Conditional_incidence_df <-na.omit(Conditional_incidence_df)
      

      
      #Plot
      
      Conditional_incidence_plot_2var <- ggplot(Conditional_incidence_df, 
                                            aes(x=.data[[variable_1]], y=incidence))+

                                            geom_ribbon(aes(ymin=conf_minus, ymax=conf_plus), alpha=0.35)+
        
                                            ylab("Conditional Incidence of salmonellosis Cases (per 10 million)")+
                                            xlab(variable_1.char)+
                                            ggtitle(light_label_char) 
        
      
      # Save plot  

  }

   
####2.5 Plots 1 variable quantiles#### 
variable_1 <- "Mean_air_temperature"


# Select the variables of interest:

Env_Pathogen_data <- salmonella_cases_detrend [,c("PostCode","date","Cases_det_adjusted", variable_1, variable_3, variable_2,"residents")]
Env_Pathogen_data <- Env_Pathogen_data[complete.cases(Env_Pathogen_data),] # Remove NA. Removes instances with no weather records                                                    

Env_laboratory_data <- Env_laboratory_data_all2[,c("PostCode","Date",variable_1, variable_3, variable_2, "residents")]
Env_laboratory_data <- Env_laboratory_data[complete.cases(Env_laboratory_data),] # Remove NA. Removes instances with no weather records



######################### stratify just based on one breaks of values

breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                     max(na.omit(Env_laboratory_data[[variable_1]])), by=1), 0)



################# Create a data frame to fill:

Conditional_incidence_quantiles <-data.frame(character(), numeric(),numeric(),numeric())
colnames(Conditional_incidence_quantiles) <-c(variable_1,"counts","residents","residents_tot")

residents_i_var <-0
residents_universal <-0

################### Divide the domains of the variables in bins according to quantiles
#Grouping the values according to the breaks of each weather variable (1,2,3)

inner_select3_min <-1 
inner_select3_max <-length(breaks_var1)

# Grouping the values according to the breaks of each weather variable (1,2,3)

for (inner_select3 in c(inner_select3_min:inner_select3_max)){ 
 
  wt <- match.closest (Env_Pathogen_data[[variable_1]], breaks_var1)
  ww <-which(wt==inner_select3)                                            # number of events when the variable coincides with the category
  Yt1 <-Env_Pathogen_data[ww, ]                                            # subset of all cases where day length is between the category above
  
  
  wt_t <- match.closest (Env_laboratory_data[[variable_1]], breaks_var1)
  ww_t <-which(wt_t==inner_select3)
  Y_tot <-Env_laboratory_data[ww_t, ]
  
         Total_cases <-sum((as.numeric(na.omit(Yt1$Cases_det_adjusted))))
          residents <-sum((as.numeric(na.omit(Yt1$residents)))) #we may not use.
          residents_tot <-sum((as.numeric(na.omit(Y_tot$residents)))) #number of people exposed to these weather conditions. Will be used to calculate the incidence
          
          
          data_df <-data.frame (
           breaks_var1 [inner_select3],
            Total_cases,
            residents,
            residents_tot
          )
          
          colnames(data_df) <-c(variable_1,"counts","residents","residents_tot")
          Conditional_incidence_quantiles <-rbind(Conditional_incidence_quantiles, data_df) 
          print(c(inner_select3, Total_cases)) # check in the output log file if x y and z cover all the conditions established in breaks()
        }
      
    

norm <- 1e-07 # set value for visualization purposes

Conditional_incidence_quantiles$incidence <-Conditional_incidence_quantiles$counts/Conditional_incidence_quantiles$residents_tot
Conditional_incidence_quantiles$incidence <-Conditional_incidence_quantiles$incidence/norm
  
  
  y=Conditional_incidence_quantiles$counts 
  n=Conditional_incidence_quantiles$residents_tot # input data as frequencies
  
  
  pw=(y + 2)/(n + 4) # Wilson's point estimate
  Wald_adj_lower <- (pw-2*sqrt(pw*(1-pw)/n))/norm # Lower conf. limit. by using 2 we are "approximating". Opted to use this method by suggestion of the textbook
  Wald_adj_upper <- (pw+2*sqrt(pw*(1-pw)/n))/norm  # Upper conf. limit 
  
  
  Conditional_incidence_quantiles$conf_minus <-Wald_adj_lower
  Conditional_incidence_quantiles$conf_plus <-Wald_adj_upper
  Conditional_incidence_quantiles <-na.omit(Conditional_incidence_quantiles)
  
  
 
  #Plot
  
  Conditional_incidence_plot_1var <- ggplot(Conditional_incidence_quantiles, 
                                            aes(x=.data[[variable_1]], y=incidence))+
                                            #aes(x=.data[[variable_3]], y=incidence))+

                                            geom_ribbon(aes(ymin=conf_minus, ymax=conf_plus), alpha=0.35)+
                                            ylab("Conditional Incidence of salmonellosis Cases (per 10 million)")+
                                            xlab(variable_1.char) 
    
  #Save plot
