#These functions are from Rezende et al. 2020 Science https://doi.org/10.5061/dryad.stqjq2c1r
 

#Function to estimate thermal tolerance landscape from static assays
# General procedure
# Step 1: Calculate CTmax and z from TDT curve
# Step 2: Calculate average log10 time and Ta (mean x and y for interpolation purposes)
# Step 3: Interpolating survival probabilities to make them comparable across treatments 
# Step 4: Overlap all survival curves into a single one by shifting each curve to mean x and y employing z
# Step 5: Build expected survival curve with mean x and y pooling all data
# Step 6: Expand expected curve to multiple Ta (with 0.1ºC difference for predictive purposes)

tolerance.landscape.np <- function(ta,time){
  
  data <- data.frame(ta,time)
  data <- data[order(data$ta,data$time),]
  
  # Step 1: Calculate CTmax and z from TDT curve
  ta <- as.numeric(levels(as.factor(data$ta)))			
  model <- lm(log10(data$time) ~ data$ta); summary(model)
  ctmax <- -coef(model)[1]/coef(model)[2]
  z <- -1/coef(model)[2]
  
  # Step 2: Calculate average log10 time and Ta (mean x and y for interpolation purposes)
  time.mn <- mean(log10(data$time))
  ta.mn <- mean(data$ta)
  
  # Step 3: Interpolating survival probabilities to make them comparable across treatments 
  time.interpol <- matrix(,1001,length(ta))
  for(i in 1:length(ta)){	
    time <- c(0,sort(data$time[data$ta==ta[i]]))
    p <- seq(0,100,length.out = length(time))
    time.interpol[,i] <- approx(p,time,n = 1001)$y}			
  
  # Step 4: Overlap all survival curves into a single one by shifting each curve to mean x and y employing z
  # Step 5: Build expected survival curve with median survival time for each survival probability
  shift <- (10^((ta - ta.mn)/z))
  time.interpol.shift <- t(t(time.interpol)*shift)[-1,]
  surv.pred <- 10^apply(log10(time.interpol.shift),1,median) 	
  
  # Step 6: Expand predicted survival curves to measured Ta (matrix m arranged from lower to higher ta)
  # Step 7: Obtain predicted values comparable to each empirical measurement
  m <- surv.pred*matrix ((10^((ta.mn - rep(ta, each = 1000))/z)), nrow = 1000)
  out <-0
  for(i in 1:length(ta)){
    time <- c(0,data$time[data$ta==ta[i]]); p <- seq(0,100,length.out = length(time))
    out <- c(out,approx(seq(0,100,length.out = 1000),m[,i],xout=p[-1])$y)}
  data$time.pred <- out[-1]
  colnames(m) <- paste("time.at",ta,sep=".")
  m <- cbind(surv.prob=seq(1,0.001,-0.001),m)
  
  for(i in 1:length(ta)){
    time <- c(0,sort(data$time[data$ta==ta[i]])); p <- seq(100,0,length.out = length(time))
    time <- c(0,sort(data$time.pred[data$ta==ta[i]]))}
  rsq <- round(summary(lm(log10(data$time) ~ log10(data$time.pred)))$r.square,3)
  list(ctmax = as.numeric(ctmax), z = as.numeric(z), ta.mn = ta.mn,  S = data.frame(surv=seq(0.999,0,-0.001),time=surv.pred),
       time.obs.pred=cbind(data$time,data$time.pred), rsq = rsq)}			
	
	

	# Function to estimate survival probability from tolerance landscapes and environmental temperature data


	
	dynamic.landscape.np.tc <- function(ta, tolerance.landscape) {
	  surv <- tolerance.landscape$S[, 2]
	  ta.mn <- tolerance.landscape$ta.mn
	  z <- tolerance.landscape$z
	  tc<-tolerance.landscape$tc
	  shift <- 10^((ta.mn - ta) / z)
	  time.rel <- 0
	  alive <- 100
	  
	  for (i in 1:length(ta)) {
	    if (ta[i] <= tc) {
	      # Set survival to 100% for temperatures lower than 26 degrees
	      alive <- c(alive, 100)
	      time.rel <- c(time.rel, 1)
	    } else {
	      if (alive[length(na.omit(alive))] > 0) {
	        alive <- try(
	          c(alive, approx(c(0, shift[i] * surv), seq(100, 0, length.out = length(c(0, surv))), xout = time.rel[i] + 1)$y),
	          silent = TRUE
	        )
	        time.rel <- try(
	          c(time.rel, approx(seq(100, 0, length.out = length(c(0, surv))), c(0, shift[i + 1] * surv), xout = alive[i + 1])$y),
	          silent = TRUE
	        )
	      } else {
	        alive <- 0
	      }
	    }
	  }
	  
	  out <- data.frame(cbind(ta = ta[1:(length(alive) - 1)], time = (1:length(ta))[1:(length(alive) - 1)], alive = alive[1:(length(alive) - 1)]))
	}
	

# Function to predict time based on dynamic curve to compare against observed data



###########################


	rezende_min<-function(tempseries,TDT_object,tc){
	  #first, make sure your temp time series are minute durations.  
	  time.min <- tempseries$time	
	  ta.min <- tempseries$temp
	  
	  #next, create a tolerance landscape
	  tl_object<-tolerance.landscape.np(TDT_object$temp,TDT_object$time)
	  tl_object$tc<-tc
	  
	  #repeat day number 1440 times (minutes in a day) 
	  day <- rep(1:(length(ta.min)/1440),each=1440)
	  daily.data <- vector("list",length(ta.min)/1440)
	  # Vectorized calculation of daily.data using lapply
	  daily.data <- lapply(split(ta.min, day), function(subset_ta) {
	    dynamic.landscape.np.tc(subset_ta, tl_object)
	  })
	  
	  # Convert the list to a data frame
	  daily.data <- do.call(rbind, daily.data)
	  
	  # Replace NAs with 0
	  daily.data[is.na(daily.data)] <- 0
	  
	  
	  #dataframe of mortality
	  daily.df<-data.frame("hours"=time.min,"mort"=100-daily.data[,3])
	  
	  #cumulative mortality
	  cum.mort<-data.frame("hours"=24*(1:max(day)),
	                       "mort"=100*(10^(cumsum(log10(daily.data[1440*(1:max(day)),3]/100)))),
	                       "type"="cum.mort")
	  
	  #repeats hourly cum.mort measure 1440 times for each minute of each day
	  cum.mort2<-cum.mort%>%slice(rep(1:n(), 1440))%>%arrange((hours))
	  #now that it is sorted and arranged by hour, give the cum.mort2 dataframe the sliced hour values
	  cum.mort2$hours<-daily.df$hours
	  daily.df$type<-"daily.mort"
	  mort.df<-rbind(daily.df,cum.mort2)
	  #add date column for easy plotting. Note we need to add 59 minutes to give us the lost hour between 11pM and midnight
	  mort.df<-mort.df%>%
	    mutate("date"=rep(seq(from=min(tempseries$time,na.rm=T),
	                          to=max(tempseries$time,na.rm=T), 
	                          by = "1 min"),times=2))
	  return(mort.df)
	  beep()
	}
	
	