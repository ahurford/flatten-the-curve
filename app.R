library(deSolve)
library(shiny)
library(ggplot2)
library(patchwork)
library(gridExtra)

# To deploy
#rsconnect::deployApp("/Users/amyhurford/Desktop/flatten-the-curve")
# in R Console

theme_set(theme_light())

SIR <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dSx <- -a * c * Sx * Ix
    dIx <- a * c * Sx * Ix - gamma * Ix - v * Ix
    dFx <- v * Ix
    dS <- -a * (1 - m1) * c * S * I
    dI <-  a * (1 - m1) * c * S * I - gamma * I - v * I
    dF <- v * I

    list(c(Sx = dSx, Ix = dIx, Fx = dFx, S = dS,
           I = dI, Fs = dF))
  })
}

server <- function(input, output) {
  output$SIR <- renderPlot({
    # Parameters are taken from Bolker & Dushoff model
    gamma <- 1/13
    chi <- 0.03
    v <- gamma * chi / (1 - chi)
    H <- 15
    c <-0.4
    a<-0.5
    parms <- c(a = a, m1 = input$m1, c = c, gamma = gamma,
               v = v, H = H)
    I0 <- 0.005
    S0 <- 1 - I0
    mintime <- 0
    maxtime <- 150
    out <- ode(y = c(Sx = S0, Ix = I0, Fx = 0, S = S0, I = I0, FS = 0), times = seq(mintime, maxtime, 1), SIR, parms)
    df <- data.frame(out)

    areaAlpha <- 0.6
    g1 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Ix * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = I * 100), fill = '#b2df8a', alpha = areaAlpha) +
      geom_hline(aes(yintercept = H), alpha = 0.2, size=3)+
      labs(x = NULL, y = NULL, title = "Percent of population infected")

    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Fx * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = FS * 100), fill = '#b2df8a', alpha = areaAlpha) +
      labs(x = "time (days)", y = NULL, title = "Cumulative fatalities")
    
    label<-""
    R_0 <- round(a*c/(v+gamma),2)
    R_2 <- round((1-input$m1)*R_0,2)
    DT <- round(log(2)/(a*c - v - gamma),2)
    DT_2 <- round(log(2)/((1-input$m1)*a*c - v - gamma),2)
    fat <- 2
    fat_2 <- 2
    toprint <- data.frame("no distancing"=label,"doubling time" = DT,"R0"=R_0,"fatalities" = fat, "with distancing"=label, "R0"=R_2, "doubling time" = DT_2, "fatalities"=fat_2)


    (g1 /
      g2 &
      scale_y_continuous(expand = expand_scale(c(0, 0.1)), labels = function(x) paste0(x,"%")) &
        scale_x_continuous(expand = c(0, 0)) ) /
      tableGrob(toprint, rows = NULL, theme = ttheme_minimal())
  })
} # End server function

# Define UI for app that draws a histogram ----
ui <- fluidPage(title = "The math behind flatten the curve",
  fluidRow(
    column(6,
           h1("The math behind flatten the curve"),
           p("by Amy Hurford, Alec Robitaille, and Joseph Baafi (Memorial University)"),
           p("Anyone interested in contributing should contact ahurford-at-mun-dot-ca.")),
    column(6,

           )
    ),
  tabsetPanel(
    tabPanel("Social distancing",
             column(5,
                    p(""),
                    
                    # Input: Slider for the number of bins ----
                    sliderInput("m1", "social distancing:",
                                min = 0, max = 1, step = 0.01, value = .2,
                                width = '100%'),helpText("0: no efforts to enact social distance"),helpText("1: fully effective isolation"),
                    p(""),
                    p("Have you heard the remark:"),
                    p(tags$b("We'll never know the effect that social distancing has had;
                      we'll never know how many lives were saved")
                    ),
                    p("We can never know this with absolute certainty,
                      but we can get some idea using mathematical models. "),
           p("Below we show that the 'flatten the curve' graphic arises from a mathematical model:"),

p(tags$a(href = "https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model", "The SIR equations")),
        p("The 'flatten the curve' graphic is not simply a drawing of an idea,
        rather, it is based on epidemiological characteristics such the duration of infectivity,
        and the disease mortality rate. The SIR equations have been used for decades,
and today, these equations remain a good start for understanding the time course of epidemics. The graphics below
are the computational output due to solving the SIR equation.
        "),

),
           column(7,

           # Output:
           plotOutput("SIR"),
           helpText("Cumulative fatalities does not account for increased death rate when health resourses are exceeded")

           )
    ),
    tabPanel("More models",

             column(6,
                    p(""),
                    
                    p("We aim to make some Newfoundland-specific graphs, but this work
                      is currently in progress")
                    )
           
       ),
  tabPanel("Newfoundland",


           column(6,
                  p(""),

                  p("We aim to make some Newfoundland-specific graphs, but this work
                    is currently in progress")
                  )

           )
))

shinyApp(ui = ui, server = server)

