## Final Project - Hubway trip analysis

# Michael Patterson - Data Science 350
# Summer Quarter

# To run this script, you need the following 4 files, are on Google Drive:
# hubway_stations.csv, boston_weather.csv, hubway_trips.csv, and IMG_ELEV2005.tif
# https://drive.google.com/folderview?id=0By0Yeqq6Gk4rMVI5a0tVRXhUWGs&usp=sharing

# import packages
library(psych)
library(data.table)
library(raster)
library(ggplot2)
library(logging)
library(glmnet)
library(boot)
library(igraph)


# load trip info from csv
load_trips <- function(trip_filename = 'hubway_trips.csv') {
  if( !(trip_filename %in% dir()) ){
    setwd('C:\\Users\\Me\\Google Drive\\Science\\Data Science\\2nd quarter\\Final project')
  }
  trips <- read.csv(trip_filename, stringsAsFactors = FALSE)
  trips$subsc_type <- as.factor(trips$subsc_type)
  
  # get rid of weird trip lengths
  short_trips <- trips$duration < 60 # remove trips shorter than 1 minute, as not likely a real trip
  trips <- trips[!short_trips,]
  long_trips <- trips$duration > 12 * 60 * 60
  trips <- trips[!long_trips,] # get rid of trips longer than 12 hours
  
  return(data.table(trips))
}

# turn trips datetime strings into separate date and hour columns
clean_datetimes_trips <- function(trips) {
  
  # split the datetime into date and time, then assign each part to correct column
  split_datetimes <- matrix(unlist(strsplit(as.character(trips$start_date), ' ') ), ncol = 2, byrow = TRUE)
  trips$start_date <- as.Date(split_datetimes[,1], format="%m/%d/%Y")
  trips$start_hour <- as.factor( substr(split_datetimes[,2], 1, 2  ) )
  
  # repeat for end_date
  split_datetimes <- matrix(unlist(strsplit(as.character(trips$end_date), ' ') ), ncol = 2, byrow = TRUE)
  trips$end_date <- as.Date(split_datetimes[,1], format="%m/%d/%Y")
  trips$end_hour <- as.factor( substr(split_datetimes[,2], 1, 2  ) )
  
  return(trips)
}

# function used by load_stations to get elevation for each station
get_elevation <- function(long, lat, elevation_raster) {
  # convert longitude and latitude to raster point format
  point <- SpatialPoints(cbind( long, lat), proj4string = CRS('+proj=longlat' ) )
  # return the elevation of the point as looked up on the raster
  return(extract(elevation_raster, point))
}

## load bike station info from csv
load_stations <- function(stations_filename = 'hubway_stations.csv') {  
  if( !(stations_filename %in% dir()) ){
    setwd('C:\\Users\\Me\\Google Drive\\Science\\Data Science\\2nd quarter\\Final project')
  }
  stations <- read.csv(stations_filename)
  
  # load elevation data for each station
  elevation <- raster('IMG_ELEV2005.tif')
  plot(elevation)
  
  # add elevation column to stations info
  stations <- transform(stations, elev = get_elevation(lng, lat, elevation))
  stations$id <- as.factor(stations$id)
  return(data.table(stations))
}

## load weather information scraped from Weather Underground
load_weather <- function(weather_filename = 'boston_weather.csv') {
  # load the csv
  if( !(weather_filename %in% dir()) ){
    setwd('C:\\Users\\Me\\Google Drive\\Science\\Data Science\\2nd quarter\\Final project')
  }
  weather_data <- read.csv(weather_filename)
  weather_data$X = NULL
  
  # extract the date and Hour from the datetime
  weather_data$datetime = paste(weather_data$date,weather_data$time)
  weather_data$datetime = strptime(weather_data$datetime, format="%Y-%m-%d %I:%M %p")
  weather_data$date <- as.Date(weather_data$date)
  weather_data$start_hour = as.factor(format(round(weather_data$datetime, units="hours"), format="%H"))
  # POSIX not accepted in data.tables
  weather_data$datetime = NULL
  weather_data$time = NULL
  
  # rename 'date' for compatability with other data.tables
  setnames(weather_data, 'date', 'start_date')
  
  # remove nonsense data
  weather_data <- weather_data[weather_data$temp > -100,]
  
  # remove spaces from column names (this was for use in formulas later, but never panned out)
  old_levels <- levels(weather_data$conditions)
  levels(weather_data$conditions) <- gsub( ' ', '_', old_levels)
  
  # merge rare conditions like Thunderstorms using regex!
  # e.g. Thunderstorms will match "Heavy Thunderstorm" and "Light Thunderstorms"
  weather_data$conditions[grepl('Thunderstorms', weather_data$conditions)] <- as.factor('Thunderstorm')
  weather_data$conditions[grepl('(Freezing|Ice|Hail)', weather_data$conditions)] <- as.factor('Ice_Pellets')
  weather_data$conditions[grepl('(Mist|Haze)', weather_data$conditions)] <- as.factor('Fog')
  weather_data$conditions[grepl('(Drizzle)', weather_data$conditions)] <- as.factor('Drizzle')
  weather_data$conditions[grepl('(g|y) Snow', weather_data$conditions)] <- as.factor('Snow')
  # I could do more merging, for example Overcast, Partly Cloudy, and Mostly Cloudy
  
  # remove hours with multiple measurements
  weather_data <- weather_data[!duplicated(weather_data[c("start_date", 'start_hour')]),]
  
  weather_data <- data.table(weather_data)
  
  # delete the conditions that no longer are present
  conditions_sum <- table(weather_data$conditions) < 10
  weather_data <- weather_data[!conditions %in% names(conditions_sum[conditions_sum])]
  weather_data$conditions <- factor(weather_data$conditions) # delete unused factors
  
  # delete some more extraneous columns
  weather_data$wind_dir_deg <- NULL
  weather_data$events <- NULL
  weather_data$dew_pt <- NULL # this is colinear with temp and humidity
  
  return(weather_data)
}

## Calculate whether a station gains or loses bikes over time
calc_trips_per_station <- function(trips, stations) {
  # calculate the total number of trips starting at each stations
  trip_starts <- as.data.frame(table(trips$strt_statn))
  names(trip_starts) <- c('id', 'total_starts')
  stations <- merge(stations, trip_starts, by='id')
  
  # calculate number of trips ending at each stations
  trip_ends <- as.data.frame(table(trips$end_statn) )# this will be slightly different length that trip_starts, depending on reality
  names(trip_ends) <- c('id', 'total_ends')
  stations <- merge(stations, trip_ends, by='id')
  
  # add the net_trips
  stations <- transform(stations, net_trips = total_ends - total_starts )
  stations <- transform(stations, norm_net_trips = 2 * net_trips / (total_starts + total_ends) )
  
  loginfo('Verifying the normalized net trips per station is a normal distribution, first plotting the histogram of values:')
  qplot(norm_net_trips, data = stations, geom = 'histogram') + theme(axis.text.x = element_text(size = 14),
                                                                     axis.text.y = element_text(size = 16),
                                                                     axis.title.x = element_text(size = 18),
                                                                     axis.title.y = element_text(size = 18))
  loginfo(paste("It looks ok! Unfortunately, it does not pass the Shapiro-Wilkes test, with a p-value of:",
                as.character(shapiro.test(stations$norm_net_trips)[2]), "Let's plot the quantile-quantile distribution."))

  qqnorm(stations$norm_net_trips)
  qqline(stations$norm_net_trips)
  loginfo(paste("The data is normal in the middle, but has fat tails. I tried a few normalization techniques,",
                "but couldn't find a good one, so I'm going to stick to non-normal data."))
  
  
  return(stations)
}

##  Does elevation influence net_trips?
elevation_analysis <- function(stations, n_boot = 1000) {

  # now fit it with linear model
  # plot the net trips by elevation
  elev_plot <- plot(stations$elev,  stations$norm_net_trips, xlab = 'Station elevation (m)', ylab ='Normalized net trips')
  elev_fit <- lm(norm_net_trips ~ elev, data = stations)
  abline(elev_fit)
  print(summary(elev_fit))
  loginfo(paste('The model produced a fit with a significant p-value, but the R-squared was only',
                as.character(summary(elev_fit)$r.squared), 'The slope was:',
                as.character(elev_fit$coefficients[2])))
  
  loginfo(paste('Residuals fail the Shapiro-Wilkes test, with p =', as.character(shapiro.test(elev_fit$residuals)[2]),
                '. Plotting residuals of the elevation fit.'))
  par(mfrow=c(2,2))
  plot(elev_fit)
  par(mfrow=c(1,1))
  
  # the elevation parameter definitely has high leverage points; are they overly influencing the regression?
  # this bootstrap code was adapted from Quick-R
  boot_func <- function(formula, data, indices) {
    data_sample <- data[indices, ] # select a sample
    fit <- lm(formula, data = data_sample)
    return(coef(fit))
  }
  boot_elev <- boot(data = stations, statistic = boot_func,
                    R=n_boot, formula= norm_net_trips ~ elev)
  
  plot(boot_elev, index = 2)
  print(boot.ci(boot_elev, type = 'basic', index = 2))
  loginfo('The 95% confidence interval of the bootstrapped slope does not overlap with 0!')
}

## Does when a station is active influence whether people take one-way trips?
time_of_day_analysis <- function(trips, stations) {
  # get the average time of day that trips start at each station
  circ_start_times = trips[, circadian.mean(as.numeric(start_hour)), by='strt_statn']
  # rename columns for the merge, then merge by station
  setnames(circ_start_times, c('strt_statn', 'V1'), c('id', 'circ_start_time') )
  stations <- merge(stations, circ_start_times, by = 'id')
  
  # now fit it!
  plot(stations$circ_start_time, stations$norm_net_trips, xlab = 'Mean trip start time (hour)', ylab = 'Norm net trips')
  time_model <- lm(norm_net_trips ~ circ_start_time, data = stations)
  print(summary(time_model))
  abline(time_model)
}

# merge the weather data frame with the trips data frame
merge_weather_trips <- function(trips, weather, allx_flag = TRUE) {
  # get sum of all trips every hour
  trips_per_hour <- trips[, list(sum_trips = .N), by = c( 'start_hour', 'start_date')] # total trips per hour
  
  # merge weather and trips into single data.frame
  # merge and keep all.x, since we want to know when no one rode because of the weather
  trips_per_hour <- merge(weather, trips_per_hour, all.x=allx_flag,  by = c('start_hour', 'start_date'))
  trips_per_hour$sum_trips[is.na(trips_per_hour$sum_trips)] <- 0 # sometimes no one takes a trip!
  trips_per_hour$conditions <- factor(trips_per_hour$conditions) # lost some ice and snow conditions since no one rode then
  trips_per_hour$start_hour <- as.factor(trips_per_hour$start_hour)
  
  return(trips_per_hour)
}

if(interactive()) {
  # set up logging
  basicConfig()
  addHandler(writeToFile, file="Michael Patterson - hubway.log", level='DEBUG') 
  
  ## Load all the data, and munge it
  loginfo('Loading bike trips, then cleaning')
  trips <- load_trips()
  trips <- clean_datetimes_trips(trips)
  loginfo('Loading bike stations, and merging with trips')
  stations <- load_stations()
  stations <- calc_trips_per_station(trips, stations)
  
  loginfo('Loading weather data')
  weather <- load_weather()
  
  # Does elevation influence net_trips?
  loginfo('Analyzing influence of elevation on trips')
  elevation_analysis(stations)
  
  # Does time of day influence net trips?
  loginfo('Analyzing influence of time of day on trips')
  time_of_day_analysis(trips, stations)
  
  ## How does weather influence bike rides?
  # first, let's see distribution of trips per hour, since that is dependent variable
  loginfo('Merging weather and trip data')
  trips_per_hour_all <- merge_weather_trips(trips, weather, allx_flag = TRUE)
  
  # let's explore by plotting average trips per hour in select conditions
  loginfo('Plotting mean trips / hour under various conditions')
  trips_by_condition <- trips_per_hour_all[, list(mean_trips = mean(sum_trips)), by =c('start_hour', 'conditions')]
  conditions_per_hour_plot <- qplot(start_hour, mean_trips, data = trips_by_condition[conditions %in% c('Clear', 'Overcast','Rain', 'Light Snow')],
                                    color = conditions) + geom_point(size=5) + xlab('Hour') + ylab('Mean trips / hour')
  conditions_per_hour_plot + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 16),
                                   axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
  
  loginfo('Plotting distribution of dependent variable, the sum of trips / hour')
  par(mfrow=c(2,2))
  hist(trips_per_hour_all$sum_trips, 100, xlab='Trips / hour (raw)', main = '')
  hist(sqrt(trips_per_hour_all$sum_trips), 100, xlab='Square root trips / hour', main = '')
  hist(log(trips_per_hour_all$sum_trips), 100, xlab = 'log trips / hour', main = '')
  hist(log(trips_per_hour_all$sum_trips+1), 100)
  loginfo('No matter how we normalize, the 0 point is dominating!')
  
  # since I don't have tools to handle biased 0 distributions, let's take out zero times
  loginfo('Merging weather and trip data, ignoring 0 trip hours')
  trips_per_hour <- merge_weather_trips(trips, weather, allx_flag = FALSE)
  
  loginfo('Plotting distribution of dependent variable, the sum of trips / hour')
  par(mfrow=c(2,2))
  hist(trips_per_hour$sum_trips, 100, xlab='Trips / hour (raw)', main = '')
  hist(sqrt(trips_per_hour$sum_trips), 100, xlab='Square root trips / hour', main = '')
  hist(log(trips_per_hour$sum_trips), 100, xlab = 'log trips / hour', main = '')
  
  loginfo('Plotting qq of different normalizations')
  par(mfrow=c(2,2))
  qqnorm(log(trips_per_hour$sum_trips), main = 'Log normalized')
  qqline(log(trips_per_hour$sum_trips))
  qqnorm(sqrt(trips_per_hour$sum_trips), main = 'Sqrt normalized')
  qqline(sqrt(trips_per_hour$sum_trips))
  qqnorm((trips_per_hour$sum_trips) ^ (1/3), main = 'Cube root normalized')
  qqline((trips_per_hour$sum_trips) ^ (1/3))
  par(mfrow=c(1,1))
  loginfo("While the tails obviously don't fit, let's use the square root normaliztion")
  
  
  ## Now let's model it:
  loginfo('Now, we fit the trips per hour over all the variables')
  naive_conditions_model <- lm(sqrt(sum_trips) ~ . - start_date, data = trips_per_hour )
  print(summary(naive_conditions_model))
  loginfo(paste('There are a lot of significant variables! Overall, this model has an AIC of ',
                as.character(AIC(naive_conditions_model))))
  
  # that's a lot of significant variables! Many of them are probably correlated, like wind direction and speed
  # or temperature and conditions. let's see if we can winnow them down with a lasso regression
  # first need to normalize
  scale_cols <- c('temp', 'humidity', 'pressure', 'visibility', 'wind_speed', 'gust_speed', 'precipitation')
  scale_trips <- trips_per_hour[, (scale_cols) := lapply(.SD, scale), .SDcols = scale_cols]
  glm_matrix <- model.matrix(~ . - sum_trips - start_date, scale_trips)
 
  loginfo("Running a cross-validated lasso regression on the variables")
  cv_lasso_model <- cv.glmnet(glm_matrix, sqrt(trips_per_hour$sum_trips), alpha = 1, family ='gaussian')
  plot(cv_lasso_model)
  print(coef(cv_lasso_model, s = 'lambda.1se'))
  
  loginfo(paste("The lasso shows that visibility and wind speed don't  matter, but still thinks wind direction does, ",
                "especially for Eastern wind. Let's see if Eastern wind is correlated with precipitation."))
  print(trips_per_hour[,list(mean_precip = mean(precipitation)), by = 'wind_dir'])
  loginfo(paste("Easterly wind generally has more precipitation. Since most of the wind directions weren't,",
                "significant in the lasso, and it's collinear with precipitation, let's drop wind direction from the fit."))
  
  # choose the coefficients of the best model  
  lassoed_conditions_model <- lm(sqrt(sum_trips) ~ . -visibility-wind_dir-wind_speed-gust_speed- start_date, data = trips_per_hour )
  print(summary(lassoed_conditions_model))
  loginfo(paste("Oh no, the new model without wind direction is worse! The new AIC is:",
                as.character(AIC(lassoed_conditions_model))))
  loginfo(paste("Trying a model with the wind direction, but less the other lassoed out variables."))
  wind_conditions_model <- glm(sqrt(sum_trips) ~ . -visibility-wind_speed- gust_speed-start_date, family = 'gaussian', data = trips_per_hour )
  print(summary(wind_conditions_model))
  loginfo(paste("With the wind direction, the AIC is back to where it was originally, with an AIC of:",
                as.character(AIC(wind_conditions_model))))
  
  ## What influences ride duration?
  # finaqqnorm(sample(log(trips$duration), 10000 ))
  loginfo('Plotting trip duration histogram. It"s not normal!')
  hist(trips$duration / 60, 100, xlab = 'Trip duration (min)', ylab = 'Count', main = 'Histogram of trip duration')
  
  # since the data isn't normal, we can't use a t-test, so I will use a KS-test; this may be way too sensitive
  casual_duration <- trips$duration[trips$subsc_type == 'Casual'] / 60
  member_duration <- trips$duration[trips$subsc_type == 'Registered'] / 60
  loginfo('Plotting distribution of casual vs registered trips:')
  plot(ecdf(casual_duration), main = 'Cumulative distribution of trip lengths', xlab = 'Trip duration (min)',
       col='blue', lwd = 2)
  lines(ecdf(member_duration), lwd = 2)
  legend(600, 0.5, c('Casual', 'Member'), lwd=3, col=c('blue', 'black'))
  loginfo(paste('Casual trips are obviously longer. To test this, since it"t not normal,',
                'I use a KS-test, and get a p-value of:',
                as.character(ks.test(casual_duration, member_duration)[2]) ))
  
  # finally, I want to see if rider age influences trip duration; but since trip duration isn't normal, I have to normalize it
  loginfo('Plotting normalization quantile plots for trip duration')
  par(mfrow=c(2,2))
  qqnorm(log(sample(member_duration, 5000)), main = 'Log normalized')
  qqline(log(sample(member_duration, 5000)))
  qqnorm(sqrt(sample(member_duration, 5000)), main = 'Square root normalized')
  qqline(sqrt(sample(member_duration, 5000)))
  qqnorm((sample(member_duration, 5000))^(1/3), main = 'Cube root normalized')
  qqline((sample(member_duration, 5000))^(1/3))
  qqnorm((sample(member_duration, 5000))^(1/10), main = '10th root normalized')
  qqline((sample(member_duration, 5000))^(1/10))
  par(mfrow=c(1,1))
  
  # I choose the log normalization, and fit
  member_trips <- trips[trips$subsc_type == 'Registered']
  sample_trips <- member_trips[sample(nrow(member_trips), 1000)]
  plot( sample_trips$birth_date, log(sample_trips$duration), xlab = 'Birth date', ylab = 'Log(trip duration (sec))')
  
  age_model <- lm(log(duration) ~ birth_date, member_trips)
  abline(age_model)
  loginfo(paste('The model yielded a significant fit, but the R-squared was', as.character(summary(age_model)$r.squared),
          ', which means it accounts for basically 0 of the variance.'))
  
  ## How I would start doing graph analysis
  graph_input <- trips[,.(strt_statn, end_statn)]
  trips_graph <- graph.data.frame(graph_input)
  trips_graph <- simplify(trips_graph, remove.loops = FALSE, edge.attr.comb= 'sum')
}