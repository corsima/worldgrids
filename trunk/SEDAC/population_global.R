## Human population trends;

library(RCurl)
# verify the certificate:
curl <- getCurlHandle()
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)


## Download the population data:
## https://www.google.com/fusiontables/data?docid=1poMNMl-lMBt2RXqvHOcUMs3npjf2qBUmXxjZWLQ
## download the fusion table:
sql = URLencode("select * from 1poMNMl-lMBt2RXqvHOcUMs3npjf2qBUmXxjZWLQ")
cat(getURL(paste("https://www.google.com/fusiontables/api/query/?sql", sql, sep="=")), file="pop.csv") 
pop <- read.csv("pop.csv")
str(pop)

## fit a GLM:
pop.glm <- glm(mean ~ year, family=gaussian(link=log), data = pop)
summary(pop.glm)
x <- -5000:2100

## goodness of fit:
par(mar=c(4.5,4.5,.5,.5))
plot(type="l", x=pop$year, predict(pop.glm, pop), lwd=2, lty=2, col="grey", xlab="Year", ylab="Total population in M", ylim=c(-22,10), xlim=c(0,2020))
points(x=pop$year, y=pop.glm$linear.predictors, pch=21)
dev.off()
par(mar=c(4.5,4.5,.5,.5))
plot(type="l", x=x, predict(pop.glm, data.frame(year=x), type = "response"), lwd=2, lty=2, col="grey", xlab="Year", ylab="Total population in M", xlim=c(0,2100))
points(x=pop$year, y=pop$mean, pch=21)
## population by 2100:
predict(pop.glm, data.frame(year=2100), type = "response")

## end of script;
