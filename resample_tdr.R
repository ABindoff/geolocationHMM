library(dplyr)
downsample <- function(dat, rate = 5*60){
  dat$j <- cut(dat$Date, breaks = seq(dat$Date[1], dat$Date[nrow(dat)], rate))
  k <- group_by(dat, j)
  k <- summarise(k, max(Light), mean(Depth), mean(Temp))
  colnames(k) <- c("Date","Light", "Depth", "Temp")
  k$Date <- as.POSIXct(k$Date, tz = "GMT") + rate
  k$Date[nrow(k)] <- k$Date[nrow(k)-1]+rate
  return(k)
} 


tag <- "TDR-MK10-08A0023-PTT-86374.tab"
d.lig <- read_delim(tag, delim = " ", skip = 4, 
                    col_names = c("obs", "Date", "Depth", "Temp", "Light"))
head(d.lig, n= 30)
d.lig$Date <- as.POSIXct(d.lig$Date, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
d.lig <- downsample(d.lig, rate = 2*60)
d.lig$Date <- d.lig$Date - 11*60*60  # tag clock was set to local time, not UTC
write_csv(d.lig, "tdr86374.csv")


gps <- "09-86373.GPS.raw.track.txt"
gdat <- read_tsv(gps) %>% select(Name, Day, Time, Longitude, Latitude, Longitude)
gdat$Day <- as.POSIXct(strptime(gdat$Day, "%d-%B-%Y"), tz = "UTC", "%Y-%m-%d")
gdat <- na.omit(gdat)
write_csv(gdat, "gps86373.csv")
