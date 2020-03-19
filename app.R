### The math behind flatten the curve ===
# by Amy Hurford and Alec Robitaille (Memorial University)

### Packages ----
library(deSolve)
library(shiny)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(data.table)

### To deploy ----
#rsconnect::deployApp("/Users/amyhurford/Desktop/flatten-the-curve")
# in R Console

### Theme ----
theme_set(theme_light())

### Functions ----
SIR <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dSx <- -a * c * Sx * Ix
    dIx <- a * c * Sx * Ix - gamma * Ix - v * Ix
    dFx <- v * Ix
    dCx <- a*c*Sx*Ix
    dS <- -a * (1 - m1) * c * S * I
    dI <-  a * (1 - m1) * c * S * I - gamma * I - v * I
    dF <- v * I
    dC <- a*c*(1-m1)*S*I

    list(c(Sx = dSx, Ix = dIx, Fx = dFx, Cx = dCx, S = dS,
           I = dI, Fs = dF, C=dC))
  })
}

SIR2 <- function(t, y, p) {
  with(as.list(c(y, p)), {
  	# The x subscript indicates "no social distancing".
    dSx <- -a * c * Sx * Ix
    dIx <- a * c * Sx * Ix - gamma * Ix - sigma*Ix - v*Ix
    # Hospitalized without a capacity cap for reference
    dHx1 <- sigma*Ix - vH*Hx1 - rho*Hx1
    dCx <- a * c * Sx * Ix
    dHcumx <- sigma*Ix
    if(Hx < H2/100){
    	# Hospitalized with a capacity cap
      dHx <- sigma*Ix - vH*Hx - rho*Hx
      # Unmet need
      dUx <- 0
    } else {
      dHx <- sigma*H2/100 - vH*Hx - rho*Hx
      # Cumulative unmet need
      dUx <- sigma*(Ix-H2/100)
    }

    # with social distancing
    dS <- -a * (1 - m2) * c * S * I
    dI <-  a * (1 - m2) * c * S * I - gamma * I - v * I - sigma*I
    dHcum <- sigma*I
    # without a cap on hospital resources
    dH1 <- sigma*I - vH*H1 - rho*H1
    dC <- a * (1 - m2) * c * S * I

    if(H0 < H2/100){
      dH0 <- sigma*I - vH*H0 - rho*H0
      dU <- 0
    } else {
      dH0 <- sigma*H2/100 - vH*H0 - rho*H0
      dU <- sigma*(I-H2/100)
    }
    list(c(Sx = dSx, Ix = dIx, Hx1 = dHx1, Hx=dHx, Ux=dUx,Cx=dCx, S = dS, I = dI, H1 = dH1, H0 = dH0, U=dU,C=dC, Hcum=dHcum, Hcumx=dHcumx))
  })
}

### Server ----
server <- function(input, output) {
  output$SIR <- renderPlot({
    # Parameters are taken from Bolker & Dushoff model
    gamma <- 1/13
    chi <- 0.03
    v <- gamma * chi / (1 - chi)
    H <- 15
    c <- 0.4
    a <- 0.5
    parms <- c(a = a, m1 = input$m1, c = c, gamma = gamma,
               v = v, H = H)
    I0 <- 0.005
    S0 <- 1 - I0
    mintime <- 0
    maxtime <- 250
    out <- ode(y = c(Sx = S0, Ix = I0, Fx = 0, Cx=0, S = S0, I = I0, FS = 0, C=0), times = seq(mintime, maxtime, 1), SIR, parms)
    df <- data.table(out)

    # Plot percent population infected
    areaAlpha <- 0.6
    g1 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Ix * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = I * 100), fill = '#b2df8a', alpha = areaAlpha) +
      geom_hline(aes(yintercept = H), alpha = 0.2, size = 3) +
      labs(x = NULL, y = NULL, title = "Percent of the population currently infected")

    # Plot cumulative fatalities
    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Fx * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = FS * 100), fill = '#b2df8a', alpha = areaAlpha) +
      labs(x = "time (days)", y = NULL, title = "Cumulative fatalities (% of the population)")

    # Generate data.frame to print
    R_0 <- round(a * c / (v + gamma), 1)
    R_2 <- round((1 - input$m1) * R_0, 1)
    DT <- round(log(2) / (a * c - v - gamma), 1)
    DT_2 <- max(round(log(2) / ((1 - input$m1) * a * c - v - gamma), 1), 0)
    fat <- round(100 * out[length(out[, 1]), 4], 1)
    fat_2 <- round(100 * out[length(out[, 1]), 7], 1)
    toprint <-
      data.frame(
        " " = c("no distancing", "with distancing"),
        "doubling time" = c(DT, DT_2),
        "R0" = c(R_0, R_2),
        "fatalities" = c(fat, fat_2),
        check.names = FALSE
    )

    # Combine plots and table with patchwork
    (g1 /
      g2 &
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)),
                         labels = function(x) paste0(x, "%")) &
        scale_x_continuous(expand = expand_scale(mult = c(0, 0)))) /
      tableGrob(toprint, rows = NULL, theme = ttheme_minimal())
  })

  output$AreaUnder <- renderPlot({
    # Parameters are taken from Bolker & Dushoff model
    gamma <- 1/13
    chi <- 0.03
    v <- gamma * chi / (1 - chi)
    c <- 0.4
    a <- 0.5
    # Yang: Lancet Respiratory Medicine - Clinical course and outcomes
    # Survial of non-survivors 1-2 weeks. Median ICU to death 7 days
    # 61.5% died before 28 days.
    # 52/201 with pneumonia included.
    # Assume 20% infections are severe.
    rho <- 1/7
    # 0.62 = vH/(vH + rho)
    # <=> 0.62*(1/7) = vH*(1-0.62)
    vH <- 0.62*(1/7)/(1-0.62)
    # 0.2*(52/201) = sigma/(v + gamma + sigma)
    # 0.2*(52/201)*(v + gamma) = sigma*(1 - 0.2*(52/201))
    sigma <- 0.2*(52/201)*(0.00238 + 1/13)/(1 - 0.2*(52/201))

    parms <- c(a = a, m2 =input$m2, c = c, gamma = gamma,
               v = v, H2 = input$H2, rho = rho, vH = vH, sigma = sigma)
    I0 <- 0.005
    S0 <- 1 - I0
    mintime <- 0
    maxtime <- 250
    out <- ode(y = c(Sx = S0, Ix = I0, Hx1 = 0, Hx=0, Ux=0, Cx=0, S = S0, I = I0, H1 = 0, H0 = 0, U = 0, C=0, Hcum=0, Hcumx=0), times = seq(mintime, maxtime, 1), SIR2, parms)
    df <- data.table(out)
		print(tail(df))

    areaAlpha <- 0.6

    g1 <- ggplot(df, aes(x = time)) +
    	geom_area(aes(y = Ix * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
    	geom_area(aes(y = I * 100), fill = '#b2df8a', alpha = areaAlpha) +
    	#geom_hline(aes(yintercept = input$H2/100), alpha = 0.2, size = 3) +
    	labs(x = NULL, y = NULL, title = "Infected (percentage of population)")

    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = 100*Hx), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = 100*H0), fill = '#b2df8a', alpha = areaAlpha) +
      geom_line(aes(y=100*Hx1), alpha = 1, size = 3, col = '#a6cee3') +
    	geom_line(aes(y=100*H1), alpha = 1, size = 3, col = '#b2df8a') +
    	geom_hline(aes(yintercept = input$H2), alpha = 0.2, size = 3) +
      labs(x = NULL, title = "Requiring critical care (percentage of population)", y= NULL)

    g3 <- ggplot(df, aes(x = time)) +
    	geom_area(aes(y = Ux* 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
    	geom_area(aes(y = U * 100), fill = '#b2df8a', alpha = areaAlpha) +
    	#geom_hline(aes(yintercept = input$H2/100), alpha = 0.2, size = 3) +
    	labs(x = "time (days)", y = NULL, title = "Cumulative unmet need (percentage of population)")

    final.unmet.x = round(last(df$Ux)*100,2)
    final.unmet = round(last(df$U)*100,2)
    final.hosp.x = round(last(df$Hcumx)*100,2)
    final.cases.x = round(last(df$Cx)*100,0)
    final.cases = round(last(df$C)*100,0)
    final.hosp = round(last(df$Hcum)*100,2)
    print(tail(df))

    toprint <-
      data.frame(
        " " = c("no distancing", "with distancing"),
        "Final unmet need (%)" = c(final.unmet.x, final.unmet),
        "Final critical care need (%):" =c(final.hosp.x, final.hosp),
        "Final cases (%)" = c(final.cases.x, final.cases),
        check.names = FALSE
      )

    # Combine plots and table with patchwork
    (g1 /
        g2 /
    			g3
      &
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)),
                           labels = function(x) paste0(x, "%")) &
        scale_x_continuous(expand = expand_scale(mult = c(0, 0)))) /
      tableGrob(toprint, rows = NULL, theme = ttheme_minimal())

  })

  output$selected_var <- renderText({
    "You have selected this"
  })

  output$scrapeTab <- renderTable({
    invalidateLater(24 * 60 * 60 * 1000)
    data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv')[Province == 'NL']
  })

  output$scrapePlot <- renderPlot({
    invalidateLater(24 * 60 * 60 * 1000)
    NL <- data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv')[Province == 'NL']

    # Plot cases in NL
    ggplot(NL, aes(x = Date, group = 1)) +
      geom_area(aes(y = log(presumptive_positive)), color = 'grey', alpha = 0.5) +
      geom_area(aes(y = log(confirmed_positive)), color = 'black', alpha = 0.6) +
      # geom_hline(aes(yintercept = H), alpha = 0.2, size = 3) +
      labs(x = "date", y = NULL, title = "log(cases in NL)") +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)),
                         labels = function(x) paste0(x, "%"))
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
                    p(""),

                    # Slider input: social distancing
                    sliderInput("m1", "social distancing:",
                                min = 0, max = 1, step = 0.01, value = .2,
                                width = '100%'),
                    helpText("0: no efforts to enact social distancing"),
                    helpText("1: fully effective isolation"),

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

           helpText("Blue curve: no changes implemented; Green curve: with social distancing; Grey line: capacity of the health care system."),
           helpText("Cumulative fatalities does not account for an increased death rate when health resourses are exceeded."),
           helpText("Doubling time: Early on in the epidemic, the time for the number of infected people to double."),
           helpText("R0: Early on in the epidemic, the average number of people subsequently infected by an infected person."),
           helpText("Fatalities: The percentage of the population that has died from COVID-19 after 250 days, however this
                    does not consider an increased death rate when health resources are overwhelmed. Epidemic models, such as SIR, suggest
                    that, even aside from preventing overwhelming the health care system, a smaller percentage of the population
                    will die from COVID-19 under social distancing."),
           helpText("The parameterization for this SIR model was taken from Bolker and Dushoff (2020)."))),
    # Area under the curve
    tabPanel("Area under the curve",
             # Left column
             tabPanel("Social distancing",
                      column(5,
                             p(""),

                             # Slider input: social distancing
                             sliderInput("m2", "social distancing:",
                                         min = 0, max = 1, step = 0.01, value = .2,
                                         width = '100%'),
                             helpText("0: no efforts to enact social distancing"),
                             helpText("1: fully effective isolation"),
                             sliderInput("H2", "Hospital capacity (%):",
                                         min = 0, max = 0.3, step = .01, value = 0.2,
                                         width = '100%'),

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
             plotOutput("AreaUnder"),
             textOutput("selected_var")
    ))),
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
    # Newfoundland tab
    tabPanel("Newfoundland",
             column(10,
                    p(""),
                    p("We aim to make some Newfoundland-specific graphs, but this work
                    is currently in progress")
                    ),

             plotOutput("scrapePlot"),
             tags$br(),
             tags$br(),
             tags$br(),
             tableOutput("scrapeTab"))
))

### Run app ----
shinyApp(ui = ui, server = server)

