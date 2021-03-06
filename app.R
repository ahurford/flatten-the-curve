### The math behind flatten the curve ===
# by Amy Hurford and Alec Robitaille (Memorial University)

### Packages ----
library(deSolve)
library(shiny)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(data.table)
library(curl)
library(bbmle)
library(tidyr)
library(chron)

### To deploy ----
#rsconnect::deployApp("/Users/amyhurford/Desktop/flatten-the-curve")
# in R Console

### Theme ----
theme_set(theme_light() +
            theme(
              axis.text = element_text(size = 14),
              legend.text = element_text(size = 14),
              legend.spacing.y = unit(0.1, 'mm')
            ))

cols <- c('No changes implemented' = '#a6cee3',
          'With social distancing' = '#b2df8a',
          'Hospital capacity' = 'grey')

areaAlpha <- 0.6

### Functions ----
source('R/SIR.R')
source('R/SIHR.R')
source('R/SEIR.R')


### Global variables ----
# NL
dataNL <- data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID19_Canada.csv',
										fill = TRUE)[Province == 'NL']
dataNL[, Date := as.IDate(Date)]
dataNL[, jul := julian(Date)]
dataNL[, daysSince := jul - min(jul)]

dataNL[, positive := sum(presumptive_positive, confirmed_positive, na.rm = TRUE),
			 by = seq.int(nrow(dataNL))]
dataNL[, casesPerDay := positive - shift(positive)]

dataNL[, testsdaily := total_testing - shift(total_testing)]

today <- dataNL[, max(daysSince)]

# Parameters are taken from Bolker & Dushoff model
gamma <- 1 / 13
chi <- 0.03
v <- gamma * chi / (1 - chi)
H <- 15
c <- 0.4
a <- 0.5
I0 <- 0.005
S0 <- 1 - I0
# Parameters from Hill 2020
a1 <- 0.2
gamma1 <- 0.133
v1 <- gamma1 * chi / (1 - chi)

mintime <- 0
maxtime <- 250

# For data.frame to print
R_0 <- round(a * c / (v + gamma), 1)
DT <- round(log(2) / (a * c - v - gamma), 1)


## SIHR
# Yang: Lancet Respiratory Medicine - Clinical course and outcomes
# Survial of non-survivors 1-2 weeks. Median ICU to death 7 days
# 61.5% died before 28 days.
# 52/201 with pneumonia included.
# Assume 20% infections are severe.
rho <- 1 / 7
# 0.62 = vH/(vH + rho)
# <=> 0.62*(1/7) = vH*(1-0.62)
vH <- 0.62 * (1 / 7) / (1 - 0.62)
# 0.2*(52/201) = sigma/(v + gamma + sigma)
# 0.2*(52/201)*(v + gamma) = sigma*(1 - 0.2*(52/201))
sigma <- 0.2 * (52 / 201) * (0.00238 + 1 / 13) / (1 - 0.2 * (52 / 201))
pop.size <- 519716

### Server ----
server <- function(input, output) {

  dataSIR <- reactive({

      parms <- c(a = a, m1 = input$m1/100, c = c, gamma = gamma,
                 v = v, H = H)
      out <- ode(y = c(Sx = S0, Ix = I0, Fx = 0, Cx = 0, S = S0, I = I0,
      								 FS = 0, C = 0),
      					 times = seq(mintime, maxtime, 1), SIR, parms)
    })

  output$SIR <- renderPlot({

    df <- data.table(dataSIR())

    # Plot percent population infected
    areaAlpha <- 0.6
    g1 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Ix * 100, fill = 'No changes implemented'), show.legend = TRUE, alpha = areaAlpha - 0.2) +
      geom_area(aes(y = I * 100, fill = 'With social distancing'), show.legend = TRUE, alpha = areaAlpha) +
      geom_hline(aes(yintercept = H, color = 'Hospital capacity'), show.legend = TRUE, alpha = 0.8, size = 2.5) +
      labs(x = NULL, y = NULL, title = "Percent of the population currently infected", fill = NULL, color = NULL) +
      guides(fill = guide_legend(override.aes = list(linetype = 0),
                                 nrow = 1),
             color = guide_legend(override.aes = list(fill = 'white')))

    # Plot cumulative fatalities
    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Fx * 100, fill = 'No changes implemented'), alpha = areaAlpha - 0.2) +
      geom_area(aes(y = FS * 100, fill = 'With social distancing'), alpha = areaAlpha) +
      labs(x = "time (days)", y = NULL, title = "Cumulative fatalities (% of the population)", fill = NULL) +
      guides(fill = FALSE, color = FALSE)

    # Combine plots and table with patchwork
    (g1 /
      g2 &
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)),
                           labels = function(x) paste0(x, "%")) &
        scale_x_continuous(expand = expand_scale(mult = c(0, 0))) &
        scale_fill_manual(values = cols) &
        scale_color_manual(values = cols)) /
      guide_area() +
      plot_layout(guides = 'collect', heights = c(5, 5, 3))
  })

  output$SIRtab <- renderTable({
    out <- dataSIR()

    R_2 <- round((1 - input$m1/100) * R_0, 1)
    DT_2 <- max(round(log(2) / ((1 - input$m1/100) * a * c - v - gamma), 1), 0)
    fat <- round(100 * out[length(out[, 1]), 4], 1)
    fat_2 <- round(100 * out[length(out[, 1]), 8], 1)

    data.frame(
      " " = c("no distancing", "with distancing"),
      "doubling time" = c(DT, DT_2),
      "R0" = c(R_0, R_2),
      "fatalities" = c(fat, fat_2),
      check.names = FALSE
    )
  })

  dataSIHR <- reactive({
    parms <- c(a = a, m2 = input$m2/100, c = c, gamma = gamma,
               v = v, H2 = input$H2, rho = rho, vH = vH, sigma = sigma)
    out <- ode(y = c(Sx = S0, Ix = I0, Hx1 = 0, Hx = 0, Ux = 0, Cx = 0, S = S0,
    								 I = I0, H1 = 0, H0 = 0, U = 0, C = 0, Hcum = 0, Hcumx = 0),
    					 times = seq(mintime, maxtime, 1),
    					 SIHR, parms)
  })

  output$SIHR <- renderPlot({
    df <- data.table(dataSIHR())

    g1 <- ggplot(df, aes(x = time)) +
    	geom_area(aes(y = Ix * 100, fill = 'No changes implemented'), alpha = areaAlpha - 0.2) +
    	geom_area(aes(y = I * 100, fill = 'With social distancing'), alpha = areaAlpha) +
    	#geom_hline(aes(yintercept = input$H2/100), alpha = 0.2, size = 3) +
    	labs(x = NULL, y = NULL, fill = NULL, title = "Infected (% of population)") +
      guides(fill = FALSE)

    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = 100*Hx1, fill = 'No changes implemented'), show.legend = TRUE, alpha = areaAlpha - 0.2) +
      geom_area(aes(y = 100*H1, fill = 'With social distancing'), show.legend = TRUE, alpha = areaAlpha) +
    	geom_hline(aes(yintercept = input$H2, color = 'Hospital capacity'), show.legend = TRUE, alpha = 0.9, size = 3) +
      labs(x = NULL, fill = NULL, color = NULL, title = "Requiring critical care (% of population)", y = NULL) +
      guides(fill = guide_legend(override.aes = list(linetype = 0),
                                 nrow = 1),
             color = guide_legend(override.aes = list(fill = 'white')))

    g3 <- ggplot(df, aes(x = time)) +
    	geom_area(aes(y = Ux* 100, fill = 'No changes implemented'), alpha = areaAlpha - 0.2) +
    	geom_area(aes(y = U * 100, fill = 'With social distancing'), alpha = areaAlpha) +
    	#geom_hline(aes(yintercept = input$H2/100), alpha = 0.2, size = 3) +
    	labs(x = "time (days)", fill = NULL, y = NULL, title = "Cumulative unmet need (% of population)") +
      guides(fill = FALSE)

    # Combine plots and table with patchwork
    (g1 /
        g2 /
    			g3 &
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)),
                           labels = function(x) paste0(x, "%")) &
        scale_x_continuous(expand = expand_scale(mult = c(0, 0))) &
        scale_fill_manual(values = cols) &
        scale_color_manual(values = cols)) /
      guide_area() +
      plot_layout(guides = 'collect', heights = c(5, 5, 5, 3))
  })

  output$SIHRtab <- renderTable({
    df <- data.table(dataSIHR())

    # For SIHRtab
    final.unmet.x = round(last(df$Ux)*100,2)
    final.unmet = round(last(df$U)*100,2)
    final.hosp.x = round(last(df$Hcumx)*100,2)
    final.cases.x = round(last(df$Cx)*100,0)
    final.cases = round(last(df$C)*100,0)
    final.hosp = round(last(df$Hcum)*100,2)

    data.frame(
        " " = c("no distancing", "with distancing"),
        "Final unmet need (%)" = c(final.unmet.x, final.unmet),
        "Final critical care need (%):" = c(final.hosp.x, final.hosp),
        "Final infected (%)" = c(final.cases.x, final.cases),
        check.names = FALSE
      )
  })

  dataSEIR <- reactive({
  	parms <- c(beta = input$R0*a1 * (gamma1 + v1) / (gamma1 + v1 + a1), gamma1 = gamma1,
  						 v1 = v1, a1 = a1, today = today,
  						 beta1 = input$R01 * a1 * (gamma1 + v1) / (gamma1 + v1 + a1))
  	I0 = 50 / pop.size
  	E0 = 2 * I0
  	out <- ode(y = c(S = 1 - I0, E = E0, I = I0, C = I0), times = seq(11, maxtime, .5), SEIR, parms)

  })
  dataSEIRnull <- reactive({
  	parmsnull <- c(beta = input$R0*a1*(gamma1+v1)/(gamma1 + v1+a1), gamma1 = gamma1,
  						 v1 = v1, a1=a1)
  	I0 = 50/pop.size
  	E0 = 2*I0
  	outnull <- ode(y = c(S = 1 - I0, E = E0, I = I0, C = I0), times = seq(11, maxtime, .5), SEIRnull, parmsnull)
  })

  output$scrapeTab <- renderTable({
  	dataNL[order(-Date), .(Date = as.character(Date),
  												 positive,
  												 negative)]
  })

  output$scrapePlot <- renderPlot({
    df <- data.frame(dataSEIR())
    dfnull <- data.frame(dataSEIRnull())

    gnew <- ggplot(dataNL) +
    	geom_point(aes(daysSince, casesPerDay)) +
    	labs(x = "", y = "New cases") +
    	#geom_line(aes(time, (C - shift(C)) * pop.size, color = 'With social distancing'), show.legend = TRUE, data = df) +
    	#geom_line(aes(time, (C - shift(C)) * pop.size, color = 'No changes implemented'), show.legend = TRUE, data = dfnull) +
    	coord_cartesian(xlim = c(0, max(dataNL$daysSince, na.rm = TRUE) + 10),  ylim = c(0, max(dataNL$casesPerDay, na.rm = TRUE) + 10)) +
    	geom_vline(aes(xintercept = 11), color = 'grey', alpha = 0.9)

    gcumu <- ggplot(dataNL) +
    	geom_line(aes(daysSince, positive / 1000), size = 1.5) +
    	labs(x = "", y = "Cumulative cases (in thousands)") +
    	geom_line(aes(time, C * pop.size / 1000, color = 'With social distancing'), show.legend = TRUE, data = df) +
    	geom_line(aes(time, C * pop.size / 1000, color = 'No changes implemented'), show.legend = TRUE, data = dfnull)

    gtests <- ggplot(dataNL) +
    	geom_point(aes(Date, testsdaily)) +
    	labs(x = NULL, y = "Daily tests reported")

    (gcumu /
    		gnew /
    		gtests &
    		scale_fill_manual(values = cols) &
    		scale_color_manual(values = cols) &
    		guides(color = FALSE, fill = FALSE)) /
    	plot_layout(heights = c(5, 5, 5))

    #plot(tail(df$time, -1), c(diff(df$C)*pop.size), typ="l", ylab = "new cases",las=1, xlab = "days since first case on March 16", lwd=4, col='#b2df8a', xlim = c(0, max(Days.Since)+10),ylim = c(0,35))
  #ylim =c(0,max(diff(df$C)*pop.size, diff(dfnull$C)*pop.size)

    #points(tail(NLData$presumptive_positive+NLData$confirmed_positive, -1) - head(NLData$presumptive_positive+NLData$confirmed_positive, -1))


    # t <-seq(11,100, .1)
    # lambda <- 1
    # plot(t, lambda*t+log(102*exp(-lambda*11)), col = "dodgerblue", lwd=4, ylab = "log(cumulative cases)", xlab = "days since first case", typ="l", ylim = c(0,4*max(log(NLData$confirmed_positive))), xlim = c(0, max(Days.Since)+10))
    # points(Days.Since,log(NLData$presumptive_positive+NLData$confirmed_positive), pch = 16)
   #  for (j in names(NL)) set(NL, which(is.na(NL[[j]])), j, 0)
   #
   #  # Plot cases in NL
   #  cols <- c('Presumptive Positive' = '#881a58',
   #            'Confirmed Positive' = '#0e288e')
   #  (ggplot(NL, aes(x = Date, group = 1)) +
   #    #geom_line(aes(y = presumptive_positive + confirmed_positive, color = '#881a58'), size = 1.5) +
   #    geom_point(aes(y = presumptive_positive + confirmed_positive, color = '#881a58'), size = 4) +
   #    #geom_line(aes(y = confirmed_positive, color = 'Confirmed Positive'), size = 1.5,
   #              #show.legend = TRUE) +
   #    #geom_point(aes(y = confirmed_positive, color = 'Confirmed Positive'), size = 2) +
   #    labs(x = NULL, y = NULL, title = "Cases in NL (presumptive + confirmed)", color  = NULL) +
   #    #scale_color_manual(values = cols) +
   #    scale_y_continuous())  /
   #    guide_area() +
   #    plot_layout(guides = 'collect', heights = c(5, 3))#expand = expand_scale(mult = c(0, 0.1)))
    })

  output$SEIRtab <- renderTable({
  	out <- data.frame(dataSEIR())
  	i1 = min(which(out$time>=30))
  	finalC <- round(last(out$C)*100,digits=0)
  	month1<- round(out$C[i1]*100,digits=0)
  	i2 = min(which(out$time>=61))
  	month2<-round(out$C[i2]*100,digits=0)
  	i3 = min(which(out$time>=92))
  	month3<-round(out$C[i3]*100,digits=0)
  	idown = min(which(diff(out$I)<=0))
  	ipeak = min(which(diff(out$I)<=0))
  	peakweek = round(out$time[ipeak]/7,0)
  	newcases = diff(out$C)*pop.size
  	peakcases = round(newcases[ipeak],0)
		#### Null stats
  	outnull <-data.frame(dataSEIRnull())
  	finalCnull <- round(last(outnull$C)*100,digits=0)
  	month1null<- round(outnull$C[i1]*100,digits=0)
  	month2null<-round(outnull$C[i2]*100,digits=0)
  	month3null<-round(outnull$C[i3]*100,digits=0)
  	idown = min(which(diff(outnull$I)<=0))
  	ipeak = min(which(diff(outnull$I)<=0))
  	peakweeknull = round(outnull$time[ipeak]/7,0)
  	newcases = diff(outnull$C)*pop.size
  	peakcasesnull = round(newcases[ipeak],0)
  	data.frame(
  		" " = c("no change", "with change"),
  		"1 month (%)" = as.character(c(month1null, month1)),
  		"2 months (%)" = as.character(c(month2null, month2)),
  		"3 months (%)" = as.character(c(month3null, month3)),
  		"end (%)" = as.character(c(finalCnull,finalC)),
  		"week of peak" = as.character(c(peakweeknull,peakweek)),
  		"new cases at peak"= as.character(c(peakcasesnull,peakcases)),
  		check.names = FALSE)

  })

}


### UI ----
ui <- fluidPage(title = "The math behind flatten the curve",
  # Title row
  fluidRow(column(
    6,
    h1("The math behind flatten the curve"),
    p("by Amy Hurford and Alec Robitaille (Memorial University)")),
          column(
    6,
    tags$br(),
    p("Source code: ", tags$a(href = "https://github.com/ahurford/flatten-the-curve", icon('github', "fa-2x"))),
    p("Anyone interested in contributing should contact ahurford-at-mun-dot-ca."))),

  # Tabsets
  tabsetPanel(

    # Left column
    tabPanel("Social distancing",
             column(5,
                    # Text in sidebar
                    p("Have you heard the remark:"),
                    p(tags$b(" 'We'll never know the effect that social distancing has had;
                      we'll never know how many lives were saved' ")),
                    p("While we can never know with certainty,
                      we can get some idea using epidemic models. 'Flatten the curve'
                      argues that effective social distancing (green curves) will lessen the maximum number of infected people during an
                      epidemic, so that hospital resources are not overwhelmed (grey line).
                      On the right, we show that 'flatten the curve' arises from a mathematical model that describes the dynamics of an
                      epidemic, specifically ",
                      tags$a(href = "https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model", "the SIR equations.")),
                    p("The 'flatten the curve' graphic is not simply a drawing of an idea;
                       rather it can arise based on disease characteristics and the interactions between
                       susceptible and infected people. The SIR equations have been used for decades,
                       and today these equations remain a good start for understanding how the
                       number of infected people changes over time. The graphs on the right
                       are the computational output due to solving the SIR equations.
                       To make these graphs, we needed to define disease
                       characterisitcs such as the duration of infectivity (assumed to be 13 days),
                       and the percentage of infections that lead to fatalities (assumed to be 3%)."),
                    p("Not all the 'flatten the curve' graphs that have appeared in the media arise from SIR or related epidemic models",
                      tags$a(href = "http://ms.mcmaster.ca/~bolker/misc/peak_I_simple.html", "(Bolker and Dushoff 2020)"), "but none-the-less, as we have shown in the
                       graphs on the right, the 'flatten the curve' idea is consistent with
                       the epidemic models commonly found in textbooks.")),
             column(7,

           # Output SIR plot and help text below:
           plotOutput("SIR"),

           # Slider input: social distancing
           sliderInput("m1", "social distancing: 0%: no efforts - 100%: complete isolation",
                       min = 0, max = 100, step = 1, value = 20,
                       width = '100%'),

           tableOutput("SIRtab"),

           helpText("Cumulative fatalities does not account for an increased death rate when health resourses are exceeded."),
           helpText("Doubling time: Early on in the epidemic, the days for the number of infected people to double."),
           helpText("R0: Early on in the epidemic, the average number of people subsequently infected by an infected person."),
           helpText("Fatalities: The percentage of the population that has died from COVID-19 after 250 days, however this
                    does not consider an increased death rate when the hospital capacity is exceeded. Epidemic models, such as SIR, suggest
                    that, even aside from preventing exceeding the hospital capacity, a smaller percentage of the population
                    will die from COVID-19 under social distancing."),
           helpText("The parameterization for this SIR model was taken from Bolker and Dushoff (2020)."))),
    tabPanel("Your questions",
                      column(5,
                             # Text in sidebar
                             p("Here we answer some of your questions we received by email."),
                             p(tags$b("Q1. How can we estimate the hospital capacity?")),
                             p(tags$a(href = "https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model", "The SIR equations"), "do not distinguish between infected individuals
                               that require hospitalization and those who don't. Since hospital capacity is
                               a fundamental component of 'flatten the curve', I addressed this question by using
                               an 'SIHR' model to distinguish between infected people that require
                               critical care (H) and those who do not (I)."),
                             p("The top panel on the right shows the percentage of the population that is infected at a
                               given time, while the middle panel shows the percentage of the population that
                              requires critical care. Note that both graphs have nearly the same shape, but that
                              the percentage of people requiring critical care at any time is much smaller (<0.3% of the population).
                               "),
                             p("In terms of what you can 'flatten', as the average citizen it is unlikely that you will require
                               critical care, so you can flatten the infections curve (the top panel). However, a consequence
                              of flattening the infections curve is to also flatten the critical care curve (the middle panel).
                               This is central to the idea of flattening the curve: we implement social distancing not necessarily to protect ourselves, but to protect the most vulernable."),
                             p("The hospital capacity line is relevant to the percentage of the population
                               requiring critical care; the default value shows hosptial capacity as 0.2% of the population, i.e,
                               if more than 2 out of every 1000 people require critical care, then hospital capacity is exceeded.
                               Admittedly, this sounds like a very well-funded healthcare system, so please see the 'More models' tab to understand that the numbers generated by this app should not be taken too literally."),
                             p(tags$b("Q2. What does social distancing equal to 20% mean?")),
                             p("For COVID-19, a contact is defined as coming within 1-2m of another person.
                               Therefore, if at your baseline level (social distancing =0%) you contacted 50 people each day,
                               social distancing at 20% means that you reduce your contacts to 80% of your
                              baseline level (=40 people per day). Please remember that you
                               are not supposed to take these numbers too literally: many of the graphs here suggest that
                               social distancing at 20% will prevent hospital overload: these models are too simple for
                               their results to be understood that precisely (see the 'More models' tab)."),
                             p("My answer above is not very good. What about becoming infected via touching a contaminated
                              surface? What about touching a contaminated surface 1 day after the contamination? And what about coming
                              into contact with someone else for 1 minute vs. 1 hour? These details should matter, but classically
                              this component of epidemiological models has been difficult to pin down.")
                      ),

                      column(7,
                             plotOutput("SIHR"),
                             tags$br(),



                             # Slider input: social distancing
                             sliderInput("m2", "social distancing: 0%: no efforts - 100%: complete isolation",
                                         min = 0, max = 100, step = 1, value = 20,
                                         width = '100%'),
                             sliderInput("H2", "Hospital capacity (% of population):",
                                         min = 0, max = 0.3, step = .01, value = 0.2,
                                         width = '100%'),

                             tableOutput("SIHRtab"),

                             helpText("Cumulative unmet need is the percentage of the population that had an unmet need because they required
                                      critical care, when the capacity to provide critical care was exceeded. For the default parameters no
                                      green curve appears because the hospital capacity is never exceeded with social distancing at 20%."),
                             helpText("Final unmet need (%): after 250 days, the percentage of the population that had an unmet need for critical care."),
                             helpText("Final critical care need (%): after 250 days, the percentage of the population that required critical care"),
                             helpText("Final infected (%): after 250 days, the percentage of the population was infected. Note these
                                      numbers are much too high. Please see the 'More models' tab for a disclaimer.")
    )),
    # More models tab
    tabPanel("More models",
             column(10,
                    p(""),

                    p("The SIR model is a great starting point for describing the dynamics of an epidemic. Yet, you are probably wondering:"),

                    p("- 'How much social distancing is enough?'"),

                    p("- 'How long will the social distancing need to last?' and"),

                    p("- 'How many infections will there be next week?'"),

                    p("The exact numbers that arise from an SIR model shouldn't be taken too literally. The SIR model makes
                      important points on a more general level, i.e., that social distancing can have a large effect and prevent
                      overwhelming hospital resources."),

                    p("Many of the world's top experts are working on answering your questions above, but these answers require
                      more than a simple SIR model. For example, the SIR model fails to consider people who are not showing symptoms, but can infected others. For one approach to improving the consistency of the epidemic model with characteristics
                      of COVID-19, see",
                      tags$a(href = "https://alhill.shinyapps.io/COVID19seir/", "Hill (2020).")))),
    									# this link looks good: http://gabgoh.github.io/COVID/index.html
    # Newfoundland tab
    tabPanel("Newfoundland & Labrador",
             column(12,
                    p(),
										p("On this tab you can consider the computational output of a", tags$a(href="https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology", "SEIR epidemic model"), "relative to the recorded number of COVID-19 cases in Newfoundland and Labrador (black dots). However, since the SEIR model assumes
											only community spread of infections, we compare only to the reported cases beginning on March 27, when cumulative cases first exceeded 100 in Newfoundland and Labrador. This app allows you to consider different future scenarios (green curve) regarding
											better or worse infection rates than in the past. Please note that this is a simple epidemic model, and predictions will give a general picture: the exact numerical values are unlikely to come to pass. See the right for a discussion of the model's limitations.")
										),
             column(7,
             			 # Slider input: social distancing
             			 sliderInput("R0", "R0 past ---------- shift the slider to achieve agreement of the curve with the reported cases (black dots)",
             			 						min = 1, max =3, step = .01, value = 2,
             			 						width = '100%'),
             			 sliderInput("R01", "R0 future ---------- shift the slider to consider a different future scenario",
             			 						min = 1, max = 3, step = .01, value = 1.5,
             			 						width = '100%'),
             			 helpText("R0 is the average number of people that are infected per one infected person, when most other people are susceptible.
             			 				 R0 past: considers March 27 (day 11 since the first case) to today. R0 future: considers from today, forward."),
                    plotOutput("scrapePlot", width = "100%"),
             			 helpText("Blue curve (no change): assumes that 'R0 past' continues into the future. Green curve (with change): assumes that after today, R0 changes to 'R0 future'."),
             			 # comma needs to be inserted above
             			 tableOutput("SEIRtab"),
             			 helpText("1 month (%): percentage of the NL population that have been infected 1 month from the first case (2 month, 3 month, and end (%) are definited similarly)."),
									 helpText("Week of peak: the week when the peak number of cases occurs. The week of the first case is week 0."),
									 helpText("New cases at peak: the number of cases reported on the day when the epidemic reaches its peak. Note that the default slider settings give
									 				 values that seem much too high.")


                    ),
             column(5,
             			 h4("Model limitations"),
             			 p("The model does not differentiate based on where people live: R0 averages across the entire provincial population.
             			 	However, the rate that St. John's residents contact others will be much different than this same rate in rural Newfoundland and Labrador."),
             			 p("The model does not consider asymptomatic infections."),
             			 p("The model parameters correspond to a 5 day period, where a person is infectious, but has not yet developed symptoms.
     A person assumed to be symptomatic and infectious for 7.5 days. Of infected people, it is assumed that 3% do not survive.
             			 These values are based on", tags$a(href = "https://alhill.shinyapps.io/COVID19seir/", "Hill (2020).")),
             			 helpText("The data below are filtered from a dataset by Dr. Michael Li", tags$a(href = "https://github.com/wzmli/COVID19-Canada/blob/master/README.md", "(here)."),
             			 				 "These data are compiled by recording information from provincal health websites around 11.30pm NST.
             			 				 "),
                    tableOutput("scrapeTab")


             )
)))

### Run app ----
shinyApp(ui = ui, server = server)

