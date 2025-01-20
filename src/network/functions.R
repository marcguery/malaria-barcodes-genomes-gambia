##############SEASONS##############

c_trans <- function(a, b, breaks = b$breaks, format = b$format, domain = b$domain) {
  a <- as.trans(a)
  b <- as.trans(b)
  
  name <- paste(a$name, b$name, sep = "-")
  
  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))
  
  trans_new(name, trans, inverse = inv, breaks = breaks, format=format, domain = domain)
  
}
rev_date <- c_trans("reverse", "time")
ori_date <- c_trans("identity", "time")

#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- c(2014, 2015, 2016, 2017)
wetseason <- as.Date(as.character(as.Date(c(paste0(years[1], "-01-01"),
                                            unlist(lapply(years, function(x){
                                              paste(x, wetseason, "15", sep="-")
                                            })),
                                            paste0(years[length(years)], "-12-31")), "%Y-%m-%d")), "%Y-%m-%d")

seasons <- data.frame(wetseason[-length(wetseason)], wetseason[-1])
seasons$type <- rep(c("dry", "wet"), length.out=nrow(seasons))
seasons$cycle <- rep(c(-1,-1), length.out=nrow(seasons))
seasons$cycle <- ifelse(rep(seasons$type[1]=="wet", length.out=nrow(seasons)),
                        seasons$cycle+floor((0:(nrow(seasons)-1))/2)+1,
                        seasons$cycle+floor((1:nrow(seasons))/2)+1)
colnames(seasons) <- c("date1", "date2", "type", "cycle")
mindate <- as.Date("2014-12-01")
maxdate <- as.Date("2017-05-01")
#Time zone warning but it is ok
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]

seasons.combination <- expand.grid(paste(seasons$date1, seasons$date2),
                                   paste(seasons$date1, seasons$date2), stringsAsFactors = F)
seasons.combination$date1start <- as.POSIXct(as.character(as.Date(sub("\\s\\S+", "", seasons.combination[,1]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination$date1end <- as.POSIXct(as.character(as.Date(sub("\\S+\\s", "", seasons.combination[,1]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination$date2start <- as.POSIXct(as.character(as.Date(sub("\\s\\S+", "", seasons.combination[,2]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination$date2end <- as.POSIXct(as.character(as.Date(sub("\\S+\\s", "", seasons.combination[,2]), "%Y-%m-%d")), "%Y-%m-%d")
seasons.combination <- seasons.combination[,-c(1,2)]

seasons.combination <- merge(seasons.combination, seasons,
                             by.x=c("date1start", "date1end"), by.y=c("date1", "date2"))
colnames(seasons.combination)[c((ncol(seasons.combination)-1):ncol(seasons.combination))] <- c("season1", "cycle1")
seasons.combination <- merge(seasons.combination, seasons, 
                             by.x=c("date2start", "date2end"), by.y=c("date1", "date2"))
colnames(seasons.combination)[c((ncol(seasons.combination)-1):ncol(seasons.combination))] <- c("season2", "cycle2")
############################
