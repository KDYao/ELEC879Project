##############################
# Ref:https://github.com/CatalystCode/sportssensor with modifications, lots of errors in original script:( 
# New feature introduced from our appraoch: pitch, yaw, roll

################################################################
# This R script transforms the "long" sensory data to be "wide". 
################################################################

# Usually, the sensory data is saved in a "long" way, meaning that at each sampling point, 
# the reading of a single sensor is saved as a row in the dataset. If each sensor is measuring 10 variables, each row will have 10 columns for these 
# 10 variables. If we have 5 sensors which are sampling at the same time stamp, we will have 5 rows like this at each time stamp. 
# Such kind of data storage format makes data recording easy, but presents some challenge for data analysis and modeling. For data analysis and modeling,
# we prefer that at each time stamp, for each experimental subject, there is only one record, and each column represents a variable measured by a sensor. 
# In the previous example, if each sensor is measuring 10 variables, and we have 5 sensors, data analysis and modeling prefer to have 1 row with 50 columns
# for each sampling time stamp. We call this format of data as "Wide". 

Sys.setlocale("LC_TIME", "C")
Sys.setlocale("LC_COLLATE", "C") 
Sys.setlocale("LC_TIME", "English")

##install libraries if not installed yet
list.of.packages <- c('reshape2', 'dplyr', 'zoo', 'ggplot2', 'base64enc')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]

if(length(new.packages))
  install.packages(new.packages)

library(base64enc)
library(dplyr)
library(reshape2) 
library(ggplot2)
library(zoo)

# Read the 1.35GB ski sensory data from Azure blob storage SAS url. Depending on the speed of internet and the configuration of your computer,
# this may take a while (even longer than a couple of hours). Alternatively, you can first download a copy of the data to your local machine,
# and load the data from this local copy to R workspace, might speed up the data loading process.
#raw_sensor_file <- "http://publicdatarepo.blob.core.windows.net/sportssensor/ski_sensor_long_anonymized.csv?sv=2014-02-14&sr=c&sig=E0%2BRNOt6%2FqiKEnVGOmV5Uu7rFYJCih9NDJYmx7wW4TU%3D&st=2017-06-22T07%3A00%3A00Z&se=2020-01-01T08%3A00%3A00Z&sp=r"
setwd("/Users/yaokundi/Documents/Courses/Queens/ELEC879/Project/Ski")
exp_dataset <- read.csv("ski_sensor_long_anonymized.csv", header=T, sep=",")
exp_dataset$index<-rownames(exp_dataset)

features<-c(
  "ExperimentDate",
  "ExperimentTime",
  "ExperimentTimeFraction",
  "experimentId",
  "activityTypeId",
  "subjectId",
  "tag",
  "x",
  "y",
  "z",
  "aX",
  "aY",
  "aZ",
  "qW",
  "qX",
  "qY",
  "qZ",
  "speed"
)
exp_dataset<-exp_dataset[features]

melted <- melt(exp_dataset, id.vars = c("ExperimentDate", "ExperimentTime","ExperimentTimeFraction","experimentId","activityTypeId","subjectId","tag"))
melted$variable<-paste(melted$tag,melted$variable)
melted$tag<-NULL 

widedata<- dcast(melted, ExperimentDate + ExperimentTime + ExperimentTimeFraction + experimentId + activityTypeId + subjectId  ~ variable)

# Replace spaces in column names to provide convenience to further analysis
col_names <- colnames(widedata)
# Replace "Left aX" to "Left_aX"
col_names <- gsub(' ', '_', col_names)
colnames(widedata) <- col_names

output_file <- "wide_data.csv" # provide a local file name to save the "wide" data
write.csv(widedata, file=output_file, row.names=FALSE, sep=",", quote=FALSE)



##########################
#Step2:
# This script enhances the sensor data by deriving more columns to better describe the 
# activities of subjects.
# Inputs:
# It takes two inputs:
# 1. A data file with all sensor readings and 
# some columns for subject_id, experiment_id, timestamps, activity id, etc. 
# 2. A csv file describing what extra columns you want to add to better describe 
#    the subject activities. Each row should have the following columns:
#    2.1 Feature_Type: what type of feature it is. Currently support four types: Distance,
#        CenterPoint, NormalLinetoPlane, and AngleCosine
#    2.2 Output_ColName: name of the derived column added to the data
#    2.3 Input_ColName1: names of the column representing point 1
#    2.4 Input_ColName2: names of the column representing point 2
#    2.5 Input_ColName3: names of the column representing point 3
# Note: in 2.3-2.5, if multiple columns are needed to define each point, the names of these 
#    columns should be separated by :, where each column is a coordinate. 
#    For instance, for feature type NormalLinetoPlane, three points are needed to define a plane in 3-D space.
#    Each point should be defined by 3 columns: x, y, and z. Then, the names of these 3 coordinates of a point
#    such as left shoulder should be concatenated by : like leftshoulder_x:leftshoulder_y:leftshoulder_z. 
# Outputs: 
# It outputs an enhanced dataset with extra columns which are derived 
# from the original sensor data


#datafile <- "C:\\Projects\\SensorReport\\final_wide_data_0422.csv"
datafile <- "wide_data.csv"

dataset <- read.csv(datafile, sep=",", header=TRUE)

#derivedSensorFile <-  "C:\\Projects\\SensorReport\\Step2_Config_DerivedSensors.csv"
derivedSensorFile <- "Rscripts/Step2_Config_DerivedSensors.csv"
dataset2 <- read.csv(derivedSensorFile, sep=",", header=TRUE, stringsAsFactors = F)
# Remove the leading and tailing spaces around Output_ColName
# Remove extra space at the beginning of each string
dataset2[,'Output_ColName'] <- gsub("^\\s+|\\s+$", "", dataset2[,'Output_ColName'])

#TODO: REMOVE THIS 
#dataset1 <- dataset[1:1000, ]

# Function to extract column names from a string, where multiple column names are 
# separated by '<sep>'
# remove extra space
extractColNames <- function(input_cols, sep=':'){
  return(gsub("^\\s+|\\s+$", "", strsplit(input_cols, sep)[[1]]))
}

# This is the function which calculates the normal line to a plane, defined by 3 points in 3-D space. 
# To define the plane, we need to first find two directions on this plane. 
# In this function, direction 1 (d1) is from p3->p1, direction 2 (d2) is from p3->p2
# Based on right-hand rule (see https://en.wikipedia.org/wiki/Right-hand_rule for details). 
# The orthogonal direction is also a 3-D vector, and normalized such that the norm (sum of square) of this
# vector is 1.
orthogonal_direction <- function(data, input_point1, input_point2, input_point3, output_name){
  nrows <- nrow(data) 
  input_point1_array <- extractColNames(input_point1)
  input_point2_array <- extractColNames(input_point2)
  input_point3_array <- extractColNames(input_point3)
  input_point1_wrt_p3 <- data[,input_point1_array] - data[,input_point3_array]
  input_point2_wrt_p3 <- data[,input_point2_array] - data[,input_point3_array]
  
  orth_dir <- as.data.frame(array(0, dim=c(nrows,3)))
  a <- as.data.frame(array(0, dim=c(2,3)))
  for (i in 1:nrows){
    a[1,] <- input_point1_wrt_p3[i,]
    a[2,] <- input_point2_wrt_p3[i,]
    a1 <- as.matrix(a)
    orth_dir[i,1] <- det(a1[,2:3])
    orth_dir[i,2] <- det(a1[,c(1,3)])
    orth_dir[i,3] <- det(a1[,1:2])
    norm <- sqrt(sum(orth_dir[i,]^2))
    orth_dir[i,] <- orth_dir[i,]/norm
  }
  orth_dir <- data.frame(orth_dir) 
  colnames(orth_dir) <- paste(output_name, c('x','y','z'), sep='_')
  data <- data.frame(data, orth_dir)
  return(data)
}

# Function to calculate the Euclidean distance between two points. 
# These two points can be 1-D, 2-D, or 3-D. 
# These two points should be using the same dimensions. 
distance <- function(data, input_coord1, input_coord2, output_name){
  input_coord1_array <- extractColNames(input_coord1)
  input_coord2_array <- extractColNames(input_coord2)
  num_dim <- length(input_coord1_array)
  dist <- NULL
  if (num_dim == 1){
    dist <- abs(data[[input_coord1_array]] - data[[input_coord2_array]])
  } else{
    for (i in 1:num_dim){
      if (i == 1) {
        dist <- (data[[input_coord1_array[i]]] - data[[input_coord2_array[i]]])^2
      } else{
        dist <- dist + (data[[input_coord1_array[i]]] - data[[input_coord2_array[i]]])^2
      }
      
    }
    dist <- sqrt(dist)
  }
  data[[output_name]] <- dist
  return(data)
}

# This function calculates the center point between two points. 
# These two points can be 1-D, 2-D, or 3-D. 
# These two points should be using the same dimensions.
center_point <- function(data, input_coord1, input_coord2, output_name){
  input_coord1 <- extractColNames(input_coord1)
  input_coord2 <- extractColNames(input_coord2)
  data[[output_name]] <- (data[[input_coord1]] + data[[input_coord2]])/2
  return (data)
}

# This function calculates the cosine of two vectors. 
# For instance, to calculate the angle between two planes, you can calculate
# the angle between the normal directions of these two planes. 
# Function orthogonal_direction can be used to calculate the normal direction to a plane
angleCosine <- function(data, input_coord1, input_coord2, output_name){
  input_coord1_array <- extractColNames(input_coord1)
  input_coord2_array <- extractColNames(input_coord2)
  data[[output_name]] <- rowSums(data[,input_coord1_array]*data[,input_coord2_array])
  return(data)
}


################# Yaw,pitch and roll 
#This function calculate the rotation matrix, it takes five inputs: qw, qx, qy, qz outputname
gyro <- function(data, input_coord1, input_coord2,input_coord3, input_coord4, output_name){
  input_coord1_array <- extractColNames(input_coord1)
  input_coord2_array <- extractColNames(input_coord2)
  input_coord3_array <- extractColNames(input_coord3)
  input_coord4_array <- extractColNames(input_coord4)
  
  out = as.data.frame(array(0, dim=c(nrow(data),3)))
  
  for (i in 1:nrow(data)){
    qw = data[,input_coord1_array[1]][i]
    qx = data[,input_coord2_array[1]][i]
    qy = data[,input_coord3_array[1]][i]
    qz = data[,input_coord4_array[1]][i]
    
    
    # create a rotation matrix values
    m00 = 1- 2*qy^2 - 2*qz^2
    m01 = 2*qx*qy - 2*qz*qw
    m02 = 2*qx*qz + 2*qy*qw
    m10 = 2*qx*qy + 2*qz*qw
    m11 = 1- 2*qx^2 - 2*qz^2
    m12 = 2*qy*qz - 2*qx*qw
    m20 = 2*qx*qz - 2*qy*qw
    m21 = 2*qy*qz + 2*qx*qw
    m22 = 1- 2*qx^2 - 2*qy^2
    
    yaw = atan2(m20, m21)
    patch = acos(m22)
    roll = -atan2(m02,m12)
    
    
    out[i,1] = yaw
    out[i,2] = patch
    out[i,3] = roll
  }
  
  colnames(out) <- paste(output_name, c('yaw','pitch','roll'), sep='_')
  
  
  print(head(out))
  
  data <- data.frame(data, out)
  return(data)
}



dataset1<-dataset
# This is the main execution part of this script. 
# It calls the functions defined above to derive the extra columns.
num_extra_features <- nrow(dataset2)
for (i in 1:num_extra_features){
  feature_type_i <- dataset2[i,'Feature_Type']
  if (feature_type_i == 'Distance'){
    dataset1 <- distance(dataset1, dataset2[i, 'Input_ColName1'], dataset2[i, 'Input_ColName2'], dataset2[i, 'Output_ColName'])
  } else if (feature_type_i == 'CenterPoint'){
    dataset1 <- center_point(dataset1, dataset2[i, 'Input_ColName1'], dataset2[i, 'Input_ColName2'], dataset2[i, 'Output_ColName'])
  } else if (feature_type_i == 'NormalLinetoPlane') {
    dataset1 <- orthogonal_direction(dataset1, dataset2[i, 'Input_ColName1'], dataset2[i, 'Input_ColName2'], dataset2[i, 'Input_ColName3'], dataset2[i, 'Output_ColName'])
  } else if (feature_type_i == 'AngleCosine') {
    dataset1 <- angleCosine(dataset1, dataset2[i, 'Input_ColName1'], dataset2[i, 'Input_ColName2'], dataset2[i, 'Output_ColName'])
  } else if (feature_type_i == 'Gyro'){
    dataset1 <- gyro(dataset1, dataset2[i, 'Input_ColName1'], dataset2[i, 'Input_ColName2'], dataset2[i, 'Input_ColName3'], dataset2[i, 'Input_ColName4'], dataset2[i, 'Output_ColName'])
  }
}

# You might want to write your enhanced data to a local csv file for further analysis, 
# feature engineering and machine learning modeling. 
# You can directly use the following lines to write to a local csv file. 
# Depending on the size of the data to be written to the destination file, it may take a while.

# outputfile <- <the path and file name to a destination csv file>
# write.csv(dataset1, file=outputfile, sep=",", row.names=FALSE, quote=FALSE)

outputfile <- "step2_output.csv"
write.csv(dataset1, file=outputfile, sep=",", row.names=FALSE, quote=FALSE)





###########################Step3
# This script does the feature engineering on sensor data. 
# For each sensor column in the input dataset, it calculates the descriptive statistics
# in each segment in an experiment of a subject. The width of the segment is defined 
# by the Window_Size in the parameter json file. 
# It also calculates the correlation between two sensors to understand how different body parts
# coordinates. 
# Inputs:
# It takes two inputs:
# 1. A data file with all sensor readings and 
# some columns for subject_id, experiment_id, timestamps, activity id, etc. 
# 2. A csv file with content in json format describing the features to generate, how to 
# slice and dice the data, how to segment the data, etc. 
# Outputs: 
# It outputs a feature set that can be used to train the machine learning models, and to analyze 
# what are the major differentiators between groups of subjects

datafile <- "step2_output.csv"
#datafile<- "wide_data.csv"
data <- read.csv(datafile, sep=",", header=TRUE, stringsAsFactors=FALSE)

# Read the text file which defines the parameters for feature engineering
# in json format
parameter_file <- "Rscripts/Step3_Config_Feature_Dict.csv"
dataset2 <- read.csv(parameter_file, sep="\t", header=F, quote="'")
library(rjson)
library(stringr)
# dataset2 will be read in as rows, if the original json text is in multiple rows. 
# Therefore, we need to concatenate these multiple rows into 
# a single string with the right json structure
parameters_json <- paste(dataset2[['V1']], collapse='')
parameters <- fromJSON(parameters_json)


# Take arguments WindowSize from command line
# For fast processing 
args = commandArgs(trailingOnly=TRUE)

# Revise windows size
if (length(args)==1) {
  # Window size
  parameters$Window_Size = as.integer(args[1])
} else if (length(args)>1) {
  stop("Only one argument (window size integer) would be accepted", call.=FALSE)
}

# Get the numeric indices of sensor columns from the parameters object
sensor_columns <- NULL
sensorColStr <- parameters$Sensor_Columns
sensorColStr_array <- strsplit(sensorColStr, ',')[[1]]
num_col_segs <- length(sensorColStr_array)
for (i in 1:num_col_segs){
  seg_i <- sensorColStr_array[i]
  seg_i <- strsplit(seg_i, '-')[[1]]
  if (length(seg_i) == 2){
    sensor_columns <- c(sensor_columns, c(as.numeric(seg_i[1]):as.numeric(seg_i[2])))
  } else{
    sensor_columns <- c(sensor_columns, as.numeric(seg_i[1]))
  }
}

# Get the list of column names from a string,
# where column names are separated by '<sep>'.
extract_columns <- function(inputstring, sep=','){
  outputArray <- strsplit(inputstring, sep)[[1]]
  outputArray <- gsub("^\\s+|\\s+$", "", outputArray)
  return(outputArray)
}

# Get the ID columns from the parameters object. These ID columns should be used to differentiate the rows in the feature set that is going to 
# be generated by this R script. These ID columns should be helpful to identify the the subject and the time when each row of feature is about
Id_Columns <- extract_columns(parameters$Id_Columns)

# Get the Segment Column names from the parameters object. These columns will be used to slice and dice the raw data when we 
# aggregate the orignal sensory data into segments for each unique combination of Segment Columns. For instance, if the 
# Segment Columns are SubjectId and ExperimentId, then the data of each unique (SubjectId, ExperimentId) is split incto small windows
# where window size is defined by parameters$Window_Size, and descriptive statistics and other features of interests in each window 
# are extracted to represent the activity in that window
Seg_Columns <- extract_columns(parameters$Segmented_By)

# Function to calculate the descriptive statistics from a segment of a single sensor signal.
# The size of the segment is defined by parameters$Window_Size. The sensor singal data passed
# to this function is just sensor readings in a single segment. 
# There are two sets of statistics: staistics in the time domain, including median, 
# standardivation, max, min, 1st, 3rd quantiles;
# and statistics in the frequency domain after Fourier transformation, 
# including the constant power, the average powers in the lower, median, and 
# higher frequency bands. 
# The lower frequency band is the first 1/3 of half of the sampling frequency, 
# the mid band is the second 1/3, and the high band is the last 1/3.
segment_signal <- function(sensor_signal_seg, sampling_freq){
  signal_statistics <- rep(0,10) #median, standardivation, max, min, 1st quantile, 3rd quantile, constant power, low-band avg power, mid-band avg power, and high-band avg power
  signal_quantile <- as.numeric(quantile(sensor_signal_seg, na.rm=TRUE))
  signal_statistics[1:6] <- c(signal_quantile[3], sd(sensor_signal_seg), signal_quantile[5], signal_quantile[1], signal_quantile[2], signal_quantile[4])
  fourierComponents <- fft(sensor_signal_seg);
  #get the absolute value of the coefficients  
  fourierCoefficients <- abs(fourierComponents);
  normalizedFourierComponents = fourierCoefficients / (sampling_freq/2);#normalize
  lowerband = round(sampling_freq/2/3)
  midband = round(sampling_freq/3)
  highband = round(sampling_freq/2)
  signal_statistics[7:10] <- c(normalizedFourierComponents[1], mean(normalizedFourierComponents[2:lowerband]), 
                               mean(normalizedFourierComponents[(lowerband+1):midband]), mean(normalizedFourierComponents[(midband+1):highband]))
  return(signal_statistics)
}

# This function splits the readings of a sensor in an experiment into segments, 
# whose width is defined by parameters$Window_Size.
# Then it calls segment_signal() to calculate the statistics of this sensor in this segment.
sensor_aggregation <- function(timestamp, sensor_signal, window_size, sampling_freq, signal_name, var_count){
  nrows <- length(sensor_signal)
  num_windows <- floor(nrows/window_size)
  seg_index <- seq(from=window_size, to=nrows, by=window_size)
  seg_index <- c(0, seg_index)
  num_index <- length(seg_index)
  feature_set <- NULL
  for (i in c(1:(num_index-1))){
    start_index <- seg_index[i] + 1
    end_index <- seg_index[i+1]
    feature_set_i <- segment_signal(sensor_signal[start_index:end_index], sampling_freq)
    if (var_count == 1){
      feature_set_i <- c(i, feature_set_i)
    }
    feature_set <- rbind(feature_set, feature_set_i)
  }
  feature_set_entire <- segment_signal(sensor_signal, sampling_freq)
  if (var_count == 1){
    feature_set_entire <- c(0, feature_set_entire)
  }
  feature_set <- rbind(feature_set, feature_set_entire)
  timestamps <- timestamp[c(seg_index[2:num_index], nrows),]
  id_names <- colnames(timestamp)
  if (var_count == 1){
    feature_set <- cbind(timestamps, feature_set)
    col_names <- c(id_names, 'seg', 'median','std','max','min','first_q','third_q','constant','lowband','midband','highband')
  } else{
    col_names <- c('median','std','max','min','first_q','third_q','constant','lowband','midband','highband')
  }
  if (var_count == 1){
    col_names[(length(id_names)+2):length(col_names)] <- paste(signal_name, col_names[(length(id_names)+2):length(col_names)], sep='_')
  } else{
    col_names <- paste(signal_name, col_names, sep='_')
  }
  feature_set <- data.frame(feature_set)
  colnames(feature_set) <- col_names
  return(feature_set)
}

# This function calculates the correlation between readings of two sensors in each segment. 
# This helps understand how two different body parts are coordinating along with each other 
# during the activity
cross_sensor_aggregation <- function(sensor_signal1, sensor_signal2, window_size, signal_name){
  nrows <- length(sensor_signal1)
  num_windows <- floor(nrows/window_size)
  seg_index <- seq(from=window_size, to=nrows, by=window_size)
  seg_index <- c(0, seg_index)
  num_index <- length(seg_index)
  feature_set <- NULL
  for (i in c(1:(num_index-1))){
    start_index <- seg_index[i] + 1
    end_index <- seg_index[i+1]
    feature_set_i <- cor(sensor_signal1[start_index:end_index], sensor_signal2[start_index:end_index])
    feature_set <- rbind(feature_set, feature_set_i)
  }
  feature_set_entire <- cor(sensor_signal1, sensor_signal2)
  feature_set <- rbind(feature_set, feature_set_entire)
  return(feature_set)
}

# Seg_Columns defines how you would like to slice and dice the data. For instance, 
# we can slice and dice the data by subject_id and experiment_id, and then for each 
# chunk of data defined by a unique (subject_id, experiment_id), we segment the data 
# into multiple segments, and calculate statistics in each segment
unique_seg_ids <- unique(data[,Seg_Columns])
num_unique_seg_ids <- nrow(unique_seg_ids)
num_seg_ids <- length(Seg_Columns)
num_rows <- nrow(data)
num_correlations <- length(parameters$Correlations)
correlation_names <- names(parameters$Correlations)
col_names <- colnames(data)
all_sub_signals <- NULL
for (i in 1:num_unique_seg_ids){
  row_index <- rep(TRUE, num_rows)
  for (j in 1:num_seg_ids) {
    row_index_i <- data[,Seg_Columns[j]] == unique_seg_ids[i, j]
    row_index <- row_index & row_index_i
  }
  exp_data <- data[row_index, ]
  exp_signals <- NULL
  var_count <- 1
  for (col_index in sensor_columns){
    signal <- sensor_aggregation(exp_data[, Id_Columns], exp_data[,col_index], as.numeric(parameters$Window_Size), 
                                 as.numeric(parameters$Sampling_Frequency), col_names[col_index], var_count)
    if (var_count == 1){
      exp_signals <- signal
    } else{
      exp_signals <- cbind(exp_signals, signal)
    }
    var_count <- var_count + 1
  }
  nrows <- nrow(exp_signals)
  exp_id <- NULL
  for (j in 1:num_seg_ids){
    exp_id <- cbind(exp_id, rep(unique_seg_ids[i,j], nrows))
  }
  exp_id <- as.data.frame(exp_id)
  colnames(exp_id) <- Seg_Columns
  col_names1 <- colnames(exp_signals)
  exp_signals <- data.frame(experimentId=exp_id, exp_signals)
  colnames(exp_signals) <- c(Seg_Columns, col_names1)
  
  # Calculate the correlations defined in parameters object
  for (j in 1:num_correlations){
    name_i <- correlation_names[j]
    var1_name = parameters$Correlations[[name_i]]$var1
    var2_name = parameters$Correlations[[name_i]]$var2
    correlation_col <- cross_sensor_aggregation(exp_data[,var1_name], exp_data[,var2_name], as.numeric(parameters$Window_Size), name_i)
    colnames_now <- colnames(exp_signals)
    exp_signals <- cbind(exp_signals, correlation_col)
    colnames(exp_signals) <- c(colnames_now, name_i)
  }
  all_sub_signals <- rbind(all_sub_signals, exp_signals)
}


# You might want to write your feature set to a local csv file for further analysis 
# and machine learning modeling. 
# You can directly use the following lines to write to a local csv file. 

outputfile_nosuffix <- 'dataset/step3_output'
outputfile = paste(outputfile_nosuffix,"_",parameters$Window_Size, ".csv", sep = '')
write.csv(all_sub_signals, file=outputfile, sep=",", row.names=FALSE, quote=FALSE)
# outputfile <- '<the path and file name to a destination csv file>'
# write.csv(all_sub_signals, file=outputfile, sep=",", row.names=FALSE, quote=FALSE)
