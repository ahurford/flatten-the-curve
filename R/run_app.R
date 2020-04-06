# Run app
dataNL <- data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv',
										fill = TRUE)[Province == 'NL']
