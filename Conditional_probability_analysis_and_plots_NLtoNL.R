
### ---- Conditional incidence derived from Dutch data and applied to the NL itself---- ###


#### 1.Conditional probability Uniform#####

# The code looks at how the risk of salmonellosis in humans depends on environmental variables
# Analysis was done following an uniform division of the range of the environmental variables, independently of the number of observations

rm(list=ls(all=TRUE)) 
# 
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




# Set the length of the desired time-lag

width <-7
width_char <-paste(width)


# Set years of analysis
first_year <- 2015
last_year <- 2019

# Selection of the variables of interest. 
#Change HERE depending on the variable object of study, by copy-pasting from the list below

variable_1 <-"Mean_air_temperature" 
variable_2 <-"Relative_humidity"
variable_3 <-"daylength"

################### Read dataframes ###################

# Read the Salmonella cases for RIVM
salmonella_data_all_RIVM <- read_csv(paste0("/Simulated_Salmonella_environment_",width_char,".csv"))
salmonella_data_all_RIVM <- as.data.frame(salmonella_data_all_RIVM [,-1])
salmonella_data_all_RIVM$Date <- as.Date(salmonella_data_all_RIVM$Date)

salmonella_data_RIVM <- salmonella_data_all_RIVM 
salmonella_data_RIVM$Cases <- as.numeric(salmonella_data_RIVM$Cases)
salmonella_data_RIVM$Cases[is.na(salmonella_data_RIVM$Cases)] <- 0
salmonella_data_RIVM <- salmonella_data_RIVM[complete.cases(salmonella_data_RIVM),] # remove instances where climate may not be available

colnames(salmonella_data_RIVM) <-c("PostCode","Date", "Cases",
                                   variable_1,variable_2,
                                   variable_3, "residents" )


# Read the weather averaged for 7 days for the Netherlands
Env_laboratory_data_all_RIVM <- read_csv(paste0("/Simulated_Laboratory_",width_char,"_tmean,rh,day.csv"))
Env_laboratory_data_all_RIVM <- as.data.frame(Env_laboratory_data_all_RIVM[,-1])
Env_laboratory_data_all_RIVM$Date <- as.Date(Env_laboratory_data_all_RIVM$Date)

colnames(Env_laboratory_data_all_RIVM) <-c("PostCode","Date",
                                           variable_1,variable_2,
                                           variable_3, "residents" )



# Select the records for the years of interest. We remove the rows for which there is no information for any of the weather variables. Dates match available data for residents catchment areas
Env_Pathogen_data_all2 <-(subset(salmonella_data_RIVM, year(as.Date(salmonella_data_RIVM$Date))>=first_year & year(as.Date(salmonella_data_RIVM$Date))<=last_year))
Env_laboratory_data_all2 <-(subset(Env_laboratory_data_all_RIVM, year(as.Date(Env_laboratory_data_all_RIVM$Date))>=first_year & year(as.Date(Env_laboratory_data_all_RIVM$Date))<=last_year))



################### Calculate the detrend value of cases to lower the behaviour component of the data ###################

salmonella_data_RIVM0 <- ddply (Env_Pathogen_data_all2, ~Date, summarise, cases=sum(Cases)) # cases at national level per day

salmonella_data_national<-salmonella_data_RIVM0[order(as.Date(salmonella_data_RIVM0$Date)),]       # sort from least recent to most recent 
dt0<-1

salmonella_data_national$Date <- as.Date (salmonella_data_national$Date, "%Y-%m-%d")
colnames(salmonella_data_national) [1] <- "date"



################### Detrend incidence data ###################

salmonella_data_national$cases<-as.numeric(salmonella_data_national$cases)
sal <- salmonella_data_national

# get the average number of cases for the last 5 years of study (2012-2016, inclusive)
mean_cases_2011_16 <- mean(sal$cases) # this has the risk of not accounting for days without cases! #15.39245

national_detrended <-detrend.series(y=sal$cases, y.name = "Raw reported cases", make.plot = TRUE,
                                    method = c("Spline"),return.info=TRUE)

sal_df <-data.frame(sal$date, sal$cases, national_detrended$series, national_detrended$curves)
colnames(sal_df)<-c("date","cases","detrend","curve")

salmonella_cases<- Env_Pathogen_data_all2 # renaming to match GLI code
colnames(salmonella_cases)[2]<-"date"


salmonella_cases_pc<-merge(salmonella_cases, sal_df, by="date")  # !careful with this dataframe, Cases corresponds to the cases per PC, and cases to the sum of Cases per day.
# we assumed that applying the national curve data to a smaller PC area is similar to calculating the curve of the specific PC.
salmonella_cases_pc$Cases_det<-salmonella_cases_pc$Cases/salmonella_cases_pc$curve # correct to use Cases here
salmonella_cases_pc$Cases_det_adjusted <- salmonella_cases_pc$Cases_det * mean_cases_2011_16 # mean daily cases equals #15.39245


salmonella_cases_detrend<-salmonella_cases_pc
salmonella_data_national_det<-ddply(salmonella_cases_detrend, ~date, summarise, cases_detrended=sum(Cases_det_adjusted))

#sort from least recent to most recent 
salmonella_data_national_det<-salmonella_data_national_det[order(as.Date(salmonella_data_national_det$date)),]


Env_laboratory_data_all2$year <- year(Env_laboratory_data_all2$Date)
population_df <- Env_laboratory_data_all2 [,c("PostCode","year", "residents")]
population_df <- unique(population_df)
colnames(population_df) [2] <- "Year"


 # Select the variables of interest:
      
Env_Pathogen_data <- salmonella_cases_detrend [,c("PostCode","date","Cases_det_adjusted", variable_1, variable_3, variable_2,"residents")] 

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
  
}else { # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature

delta_var1 <-1 # Same breaks for Tmax, Tvar, Tmean and dewpoint temperature
breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                           max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 0)
}



if (variable_2=="Relative_humidity"){
  
  delta_var2 <-5
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
}


# Create a data frame to fill with
estratification_weather_cases <-data.frame(character(), character(), character(), numeric(), numeric(), numeric())
colnames(estratification_weather_cases) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")



################### Conditional incidence stratification:
inner_select3_min <-1 
inner_select3_max <-length(breaks_var3)#-1


#gives the index of the categorical "bin" that the value of the variable belongs to for every nested weather variable:

for (inner_select3 in c(inner_select3_min:inner_select3_max)){
  
  # below I round the values of the variables of interest in Env_Pathogen_data to the closest break value and get the index value of the nearest break
  
  wt <- match.closest (Env_Pathogen_data[[variable_3]], breaks_var3) # new way to round to the nearest and get the index
  ww <-which(wt==inner_select3)                                                                      # number of events when the variable coincides with the category
  Env_Pathogen_data_z <-Env_Pathogen_data[ww, ]                                            # subset of all cases where daylength is between the category above
  
  #same for lab dataset
  wt_t <- match.closest (Env_laboratory_data[[variable_3]], breaks_var3)
  ww_t <-which(wt_t==inner_select3)
  Env_laboratory_data_z <-Env_laboratory_data[ww_t, ]
                                                                                           # subsetting the same as above including humidity delimitation included
  {
    {
      inner_select3_min <- 1
      inner_select3_max <-length(breaks_var2)#-1
      
      for (inner_select3 in c(inner_select3_min:inner_select3_max)){ # from the bins resulting above, find another variable stratification
        
        wt_y <- match.closest (Env_Pathogen_data_z[[variable_2]], breaks_var2)
        ww_y <-which(wt_y==inner_select3)
        Env_Pathogen_data_y <-Env_Pathogen_data_z[ww_y, ]
        
        wt_y_t <- match.closest (Env_laboratory_data_z[[variable_2]], breaks_var2)
        ww_y_t <-which(wt_y_t==inner_select3)
        Env_laboratory_data_y <-Env_laboratory_data_z[ww_y_t, ] 
        
       
        {
          inner_select1_min <-1
          inner_select1_max <- length(breaks_var1)#-1
          
          for (inner_select1 in c(inner_select1_min:inner_select1_max)){
            
            wt_x <- match.closest (Env_Pathogen_data_y[[variable_1]], breaks_var1)
            ww_x <-which(wt_x==inner_select1)
            Pathogen_stratified <-Env_Pathogen_data_y[ww_x,] #
            
            
            wt_x_t <- match.closest (Env_laboratory_data_y[[variable_1]], breaks_var1)
            ww_x_t <-which(wt_x_t==inner_select1)
            Laboratory_stratified <-Env_laboratory_data_y[ww_x_t,] 
            
            Total_cases <-sum((as.numeric(na.omit(Pathogen_stratified$Cases_det_adjusted))))
            residents <-sum((as.numeric(na.omit(Pathogen_stratified$residents)))) # number of people for areas with a recorded case for specific values of the 3 weather variables.
            residents_tot <-sum((as.numeric(na.omit(Laboratory_stratified$residents)))) #number of people exposed to these weather conditions, irrespectively whether there was a case found or not. Will be used to calculate the incidence
            
            
            data_df <-data.frame (
              breaks_var3[inner_select3],
              breaks_var2[inner_select3],
              breaks_var1[inner_select1],
              Total_cases,
              residents,
              residents_tot
            )
            
            colnames(data_df) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")
            estratification_weather_cases <-rbind(estratification_weather_cases, data_df) 
            print(c(inner_select1, inner_select3,inner_select3, Total_cases)) 
          }
        }
      }
    }
  }
} # output: estratification_weather_cases

write.table(estratification_weather_cases, paste(variable_3,"_",variable_2,"_",variable_1,"_",width_char,
                                "_Simulated_for_rec_UNIFORM_narrowedvalues_DETREND_NL.csv",sep=""), 
            col.names = NA, row.names = TRUE , sep = ",",eol = "\n")

############ Reconstruction from Conditional probability
# The code does look at how the risk of salmonellosis in humans depends on environmental variables
Env_Pathogen_data$PostCode <-as.factor(Env_Pathogen_data$PostCode)
Env_laboratory_data$year <-year(as.Date(Env_laboratory_data$Date))


modelled_cases<-c()


#----insert of breaks ---
# Weather Variable 3

dt_var3 = data.table(breaks_var3, val = breaks_var3)
dt_var3 <- na.omit (dt_var3)
setattr(dt_var3, "sorted", "breaks_var3")
dt_var3.2 <- dt_var3[J(Env_laboratory_data[[variable_3]]), roll = "nearest"] 

weather_bins <- as.data.frame(Env_laboratory_data[,c("PostCode","Date", "residents")])
weather_bins$break_var3<- dt_var3.2$val


# Weather Variable 2
dt_var2 = data.table(breaks_var2, val = breaks_var2)
dt_var2 <- na.omit (dt_var2)
setattr(dt_var2, "sorted", "breaks_var2")
dt_var2.2 <- dt_var2[J(Env_laboratory_data[[variable_2]]), roll = "nearest"]


weather_bins$break_var2<- dt_var2.2$val



# Weather variable 1
dt_var1 = data.table(breaks_var1, val = breaks_var1)
dt_var1 <- na.omit (dt_var1)
setattr(dt_var1, "sorted", "breaks_var1")
dt_var1.2  <- dt_var1[J(Env_laboratory_data[[variable_1]]), roll = "nearest"]

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

data_table_1 = data.table(modelled_cases, key= c('Year','PostCode')) 
data_table_2 = data.table(population_df, key=c('Year','PostCode')) 
modelled_cases <- data_table_1 [data_table_2, residents:=residents]

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

real_cases_mean <-real_cases_national %>% mutate (rolling_mean= slider::slide_mean(Cases_det_adjusted, before=3 , after=3 )) # %>% ungroup()
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


yearly_average <-ggplot(average_data, aes(x=Day, y=Mean, colour=source))+

  geom_ribbon(aes(ymin=f_quant, ymax=s_quant, fill=source), alpha=0.15) +
  xlab("Date") + ylab("Salmonellosis Cases") +
  
  ggtitle(paste0(variable_1.char,", ", variable_2.char,", ", variable_3.char))



# Save the plot




#### 2. Conditional probability Quantile ####

# The code does look at how the risk of salmonellosis in humans on a quantile division of the range of the environmental variables
# load("Datasets_detrend.RData")


# Select the variables of interest:
      
Env_Pathogen_data <- salmonella_cases_detrend [,c("PostCode","date","Cases_det_adjusted", variable_1, variable_3, variable_2,"residents")]
Env_Pathogen_data <- Env_Pathogen_data[complete.cases(Env_Pathogen_data),] # Remove NA. Removes instances with no weather records                                                    
      
Env_laboratory_data <- Env_laboratory_data_all2[,c("PostCode","Date",variable_1, variable_3, variable_2, "residents")]
Env_laboratory_data <- Env_laboratory_data[complete.cases(Env_laboratory_data),] # Remove NA. Removes instances with no weather records
      


######################### Create the breaks for a quantile division of the bins: detrended cases
  # The bins collect the value of bins_size_var3 (e.g. 25% of the data). The size of the bins is different but the number of observations is the same. The bins depend one from another weather variable to keep the number of observations per bin
  # I created 4 quantiles where there is not much diversity of values (i.e. daylength, RH...), and smaller for higher diversity (Radiation, Temperature..)

bins_size_var1 <-0.1 # temp, dewpoint
bins_size_var2 <-0.25 # precip, widspeed, air pressure
bins_size_var3 <-0.25 # DL, global radiation, windspeed



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

breaks_var1 <-function(variable_3, variable_2, variable_1, bins_size_var3, bins_size_var2, bins_size_var1, inner_select3, inner_select3)
{
  wt <-(findInterval(Env_laboratory_data[[variable_3]],breaks_var3(variable_3, bins_size_var3)))
  ww <-which(wt==inner_select3)
  Env_laboratory_data_some <-Env_laboratory_data[ww, ]
  
  if(is.na(breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)[inner_select3])=='FALSE')
  {
    wt2 <-(findInterval(Env_laboratory_data_some[[variable_2]], breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)))
    ww2 <-which(wt2==inner_select3) 
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
  #Grouping the values according to the breaks of each weather variable

inner_select3_min <-1 
inner_select3_max <-length(breaks_var3(variable_3, bins_size_var3))-1

# Grouping the values according to the breaks of each weather variable

for (inner_select3 in c(inner_select3_min:inner_select3_max)){ 
 
  wt <- findInterval (Env_Pathogen_data[[variable_3]], breaks_var3(variable_3, bins_size_var3))
  ww <-which(wt==inner_select3)                                                                      # number of events when the variable coincides with the category
  Env_Pathogen_data_z <-Env_Pathogen_data[ww, ]                                            # subset of all cases where daylength is between the category above
  
  
  wt_t <- findInterval (Env_laboratory_data[[variable_3]], breaks_var3(variable_3, bins_size_var3))
  ww_t <-which(wt_t==inner_select3)
  Env_laboratory_data_z <-Env_laboratory_data[ww_t, ]
  

    if (length(breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))!=0)
    {
      inner_select3_min <- 1
      inner_select3_max <-length(breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))-1
      
      for (inner_select3 in c(inner_select3_min:inner_select3_max)){ # from the bins resulting above, find another variable stratification
        
        
        wt_y <- findInterval (Env_Pathogen_data_z[[variable_2]], breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))
        ww_y <-which(wt_y==inner_select3)
        Env_Pathogen_data_y <-Env_Pathogen_data_z[ww_y, ]
        
        
        wt_y_t <- findInterval (Env_laboratory_data_z[[variable_2]], breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3))  
        ww_y_t <-which(wt_y_t==inner_select3)
        Env_laboratory_data_y <-Env_laboratory_data_z[ww_y_t, ] 
        
        
        if (length(breaks_var1(variable_3, variable_2, variable_1, bins_size_var3, bins_size_var2, bins_size_var1, inner_select3, inner_select3))!=0)
        {
          inner_select1_min <-1
          inner_select1_max <- length(breaks_var1(variable_3,variable_2,variable_1,bins_size_var3,bins_size_var2,bins_size_var1,inner_select3,inner_select3))-1
          
          for (inner_select1 in c(inner_select1_min:inner_select1_max)){

            wt_x <- findInterval (Env_Pathogen_data_y[[variable_1]], breaks_var1(variable_3,variable_2,variable_1,bins_size_var3,bins_size_var2,bins_size_var1,inner_select3,inner_select3))
            ww_x <-which(wt_x==inner_select1)
            Yt1 <-Env_Pathogen_data_y[ww_x, ] #where there are cases for specific values of the 3 weather variables. residents is the number of people exposed to these situations during the entire duration of the dataset. It is the number of times the people is exposed to this caracteristics, only when a case is found. It can count a same catchment area twice. Residents total is the same irrespectively of the presence of a case

            
            wt_x_t <- findInterval (Env_laboratory_data_y[[variable_1]], breaks_var1(variable_3,variable_2,variable_1,bins_size_var3,bins_size_var2,bins_size_var1,inner_select3,inner_select3))
            ww_x_t <-which(wt_x_t==inner_select1)
            Y_tot <-Env_laboratory_data_y[ww_x_t, ] #NA in residents at Env_laboratory_data. we do not care if there is a case reported
            
            #Yt1 <- merge(Y_tot, salmonella_cases_detrend)
            
            Total_cases <-sum((as.numeric(na.omit(Yt1$Cases_det_adjusted))))
            residents <-sum((as.numeric(na.omit(Yt1$residents)))) #we may not use.
            residents_tot <-sum((as.numeric(na.omit(Y_tot$residents)))) #number of people exposed to these weather conditions. Will be used to calculate the incidence
            
            
            data_df <-data.frame (
              breaks_var3(variable_3, bins_size_var3)[inner_select3],
              breaks_var2(variable_3, variable_2, bins_size_var3, bins_size_var2, inner_select3)[inner_select3],
              breaks_var1(variable_3, variable_2, variable_1, bins_size_var3, bins_size_var2, bins_size_var1, inner_select3, inner_select3)[inner_select1],
              Total_cases,
              residents,
              residents_tot
            )
            
            colnames(data_df) <-c(variable_3, variable_2, variable_1,"counts","residents","residents_tot")
            Conditional_incidence_quantiles <-rbind(Conditional_incidence_quantiles, data_df) 
            print(c(inner_select1, inner_select3,inner_select3, Total_cases)) # check in the output log file if x y and z cover all the conditions established in breaks()
          }
        }
      }
    }
} # output: Conditional_incidence_quantiles



#### Plots CI Quantile ####

# Read data frame if needed
Conditional_incidence_quantiles <-read.csv(paste0("Conditional_probability_",variable_3,"_",variable_2,"_",variable_1,"_",width_char,"_quantile_detrend_2000-2016.csv"))
Conditional_incidence_quantiles <-Conditional_incidence_quantiles[,-1]
Conditional_incidence_quantiles$incidence <-Conditional_incidence_quantiles$counts/Conditional_incidence_quantiles$residents_tot

# Normalized incidence per 1,000,000 inhabitants to normalize our plots later
norm <- 1/1e7 # changed to 10 million on 23/05/23 after adding the subset of cases



#### Plots 3 variables quantiles####


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

