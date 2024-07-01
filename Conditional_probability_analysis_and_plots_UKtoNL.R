
###----Conditional incidence analysis applying the conditional probability calculated for England and Wales to the Netherlands ----###

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



# Set the length of the desired time-lag

width <-7
width_char <-paste(width)


# Set years of analysis
first_year <- 2000
last_year <- 2016

# Selection of the variables of interest. 

variable_1 <-"Mean_air_temperature" 
variable_2 <-"Relative_humidity"
variable_3 <-"daylength"

################### Read data frames ###################
  # Read the Salmonella cases for RIVM
salmonella_data_all_RIVM <- read.csv(paste0("/Simulated_Salmonella_environment_",width_char,".csv"), na.strings = NA)
salmonella_data_all_RIVM <- salmonella_data_all_RIVM [,-1]

salmonella_data_RIVM <- salmonella_data_all_RIVM [year(salmonella_data_all_RIVM$Date) %in% c(2015:2019),]  
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


# Select years of interest: 2015 to 2019

Env_laboratory_data_RIVM <- Env_laboratory_data_all_RIVM [year(Env_laboratory_data_all_RIVM$Date) %in% c(2015:2019),]  
Env_laboratory_data_RIVM <- Env_laboratory_data_RIVM[complete.cases(Env_laboratory_data_RIVM),]


# Load dataframes with detrend applied in England and Wales:
 load("/Datasets_detrend.RData")

Env_laboratory_data <- Env_laboratory_data_all2 [,c("PostCode","Date", variable_1, variable_3, variable_2, "residents")]
Env_laboratory_data <- Env_laboratory_data[complete.cases(Env_laboratory_data),]

# Read CI stratification from England and Wales:
estratification_weather_cases <- read.csv(paste0(variable_3,"_",variable_2,"_",variable_1,"_",width_char,"_Simulated_for_rec_UNIFORM_narrowedvalues_DETREND.csv"))
estratification_weather_cases <- estratification_weather_cases [,-1]

# Remove the correcting factor from the England and Wales
mean_cases_2011_16 <- 15.39245 # England and Wales
mean_cases_RIVM <- 2.405965*59642000/17243913 # NL country adjusted for the latest 
estratification_weather_cases$counts <- estratification_weather_cases$counts/mean_cases_2011_16*mean_cases_RIVM


###### Reconstruction of NL #####
      
# Repeated values as per England and Wales:
delta_var2 <-5 
delta_var1 <-1
delta_var3 <-1


#Quantile analysis is good for plotting, but for the reconstruction we do for now from uniform. Here we are defining the size of the uniform bins for the data to be clasified afterwards:
### We use the lab info because we are applying the conditional incidence to the whole data set for all the PC even if there was no case
# to match the England and Wales breaks
breaks_var2 <-round(seq(min(na.omit(Env_laboratory_data[[variable_2]])), 
                       max(na.omit(Env_laboratory_data[[variable_2]])), by=delta_var2), 0) 

breaks_var1 <-round(seq(min(na.omit(Env_laboratory_data[[variable_1]])),   
                                           max(na.omit(Env_laboratory_data[[variable_1]])), by=delta_var1), 0)

breaks_var3 <-round(seq(min(na.omit(Env_laboratory_data[[variable_3]])), 
                            max(na.omit(Env_laboratory_data[[variable_3]])), by=delta_var3),0)


modelled_cases<-c()



#----Round weather values to the closest break ---
# remove values that are +-delta_var2 from the England and Wales values. 6434 instances removed

Env_laboratory_data_RIVM <- Env_laboratory_data_RIVM [Env_laboratory_data_RIVM[[variable_2]]>= min(Env_laboratory_data_RIVM[[variable_2]]-delta_var2)&
                                                        Env_laboratory_data_RIVM[[variable_2]]<= max(Env_laboratory_data_RIVM[[variable_2]]+delta_var2) ,] 
Env_laboratory_data_RIVM <- Env_laboratory_data_RIVM [Env_laboratory_data_RIVM[[variable_1]]>= min(Env_laboratory_data_RIVM[[variable_1]]-delta_var1)&
                                                        Env_laboratory_data_RIVM[[variable_1]]<= max(Env_laboratory_data_RIVM[[variable_1]]+delta_var1) ,] 



# -- Daylength 
dt_var3 = data.table(breaks_var3, val = breaks_var3)
dt_var3 <- na.omit (dt_var3)
setattr(dt_var3, "sorted", "breaks_var3")
dt_var3.2 <- dt_var3[.(val = Env_laboratory_data_RIVM[[variable_3]]), on = "val", roll = "nearest"]

weather_bins <- as.data.frame(Env_laboratory_data_RIVM[,c("PostCode","Date", "residents")])
weather_bins$break_var3<- dt_var3.2$breaks_var3


# -- Relative humidity
dt_var2 = data.table(breaks_var2, val = breaks_var2)
dt_var2 <- na.omit (dt_var2)
setattr(dt_var2, "sorted", "delta_var2")

dt_var2.2 <- dt_var2[.(val = Env_laboratory_data_RIVM[[variable_2]]), on = "val", roll = "nearest"]

weather_bins$break_var2<- dt_var2.2$breaks_var2



# -- Temperature
dt_var1 = data.table(breaks_var1, val = breaks_var1)
dt_var1 <- na.omit (dt_var1)
setattr(dt_var1, "sorted", "breaks_var1")
dt_var1.2 <- dt_var1[.(val = Env_laboratory_data_RIVM[[variable_1]]), on = "val", roll = "nearest"]


weather_bins$break_var1<- dt_var1.2$breaks_var1
colnames(weather_bins) <- c("PostCode","Date", "residents", variable_3, variable_2, variable_1)


CI_per_weather_bin <- left_join (weather_bins, estratification_weather_cases, by = c(variable_3, variable_2, variable_1)) # inner join, keep only matching values



# calculation of all the potential incidence based on the weather conditions of the area. The counts used as numerator correspond to the simulation of cases per weather condition.
CI_per_weather_bin$incidence <-CI_per_weather_bin$counts/CI_per_weather_bin$residents_tot 
CI_per_weather_bin$incidence[is.na(CI_per_weather_bin$incidence)] <- 0 

comp_cases <-CI_per_weather_bin$incidence*CI_per_weather_bin$residents.x # estimation of the "real" incidence based on the possible cases due to weather and the residents number in the grid area


modelled_cases <-data.frame(CI_per_weather_bin$Date, comp_cases, CI_per_weather_bin$incidence, CI_per_weather_bin$PostCode) # modelled time series
colnames(modelled_cases) <-c("Date","Cases","Lambda","Lab")

modelled_cases <- distinct(modelled_cases)



      
####### Plot simulated reconstruction ####
real_cases_all <-salmonella_data_RIVM 
real_cases_all$Date <- as.Date(real_cases_all$Date, format="%Y-%m-%d")


modelled_cases$Date <- as.Date(modelled_cases$Date)
colnames(modelled_cases) <-c("Date","Cases","Lambda","PostCode")
modelled_cases$Year <-year(modelled_cases$Date)
modelled_cases <-modelled_cases[order(modelled_cases$Date, na.last = NA),]




### sum cases per PC 
modelled_cases_national <- aggregate (Cases ~ Date, modelled_cases, sum) 

modelled_cases_summary <-cbind(modelled_cases_national, rep("Model",times=length(modelled_cases_national[,1]))) 
colnames(modelled_cases_summary)<-c("Date","Mean","source")

### Real cases
## add all cases per date
real_cases_national <- aggregate (Cases ~ Date, real_cases_all, sum) 

real_cases_summary <-cbind(real_cases_national, rep("Cases",times=length(real_cases_national[,1])))
colnames(real_cases_summary) <-c("Date","Mean","source")              

real_cases_summary$Date<-as.factor(real_cases_summary$Date)
modelled_cases_summary$Date <- as.factor(modelled_cases_summary$Date)

modelled_cases_all <-rbind(real_cases_summary, modelled_cases_summary)



######## Plots #####

modelled_cases_all$Date <- as.Date(modelled_cases_all$Date)
modelled_cases_all$Date <- as.character(modelled_cases_all$Date)


modelled_cases_plot <-ggplot(modelled_cases_all, aes(x=Date, y=Mean, colour=source, group=year(Date)))+
  
                          geom_line(linewidth=0.75) + 
                          xlab("Date") + ylab("Salmonellosis Cases")


# Save the plot




############## Average per day of the year 

modelled_cases_national$yday <-as.factor(yday(modelled_cases_national$Date))

modelled_cases_average1 <-ddply(modelled_cases_national,~yday, summarise, mean=mean(Cases)) # mean of all cases for the same day of the year for the 17 years of study
modelled_cases_quantile1 <-ddply(modelled_cases_national,~yday, function (x) quantile(x$Cases , c(.25,.5,.75))) # quantiles of cases for all the years for the same day of the year
modelled_cases_average2 <-cbind(modelled_cases_average1,modelled_cases_quantile1[,-1])


real_cases_national$yday <-as.factor(yday(as.Date(as.character(real_cases_national$Date))))
real_cases_average1 <-ddply(real_cases_national, ~yday, summarise, mean=mean(Cases))
real_cases_quantile1 <-ddply(real_cases_national, ~yday,function (x) quantile(x$Cases, c(.25,.5,.75)))
real_cases_average2 <-cbind(real_cases_average1, real_cases_quantile1[,-1])



df1 <-data.frame(real_cases_average2,rep("Cases",times=length(real_cases_average2[,1])))
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


