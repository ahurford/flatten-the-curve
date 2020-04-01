
NLData <- data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv', fill=TRUE)[Province == 'NL']

NLData <-NLData[-9,]
New.Dates = data.frame(date = NLData$Date, stringsAsFactors = FALSE)
# This %>% is from the tidyr package
New.Dates = New.Dates %>% separate(date, sep="-", into = c("year", "month", "day"))
New.Dates$day = as.numeric(New.Dates$day)
New.Dates$month = as.numeric(New.Dates$month)
New.Dates$year = as.numeric(New.Dates$year)
Days.Since = julian(New.Dates$month,New.Dates$day, New.Dates$year,  c(month = 3,day= 16,year = 2020))
# Replace "NA" with 0 for comfirmed and presumptive cases.
NLData$presumptive_positive[is.na(NLData$presumptive_positive)] <- 0
NLData$confirmed_positive[is.na(NLData$confirmed_positive)] <- 0
TotalCases = NLData$presumptive_positive+NLData$confirmed_positive
NewCases = TotalCases-c(0,head(TotalCases,-1))
DaysBetween = c(1,tail(Days.Since,-1)-head(Days.Since,-1))
NewCasesPerDay = NewCases/DaysBetween
# index when first over 100 cases
i = min(which(TotalCases>=100))
# index of last obs.
i1 = length(Days.Since)

plot(Days.Since[i:i1], NewCasesPerDay[i:i1]/TotalCases[i:i1-1], typ = "l", ylab = "Change in cases per case", xlab = "Days sine March 16")