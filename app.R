library(deSolve)
library(shiny)
library(ggplot2)
library(patchwork)
library(gridExtra)

theme_set(theme_light())


SIR <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dSx <- -a * c * Sx * Ix
    dIx <- a * c * Sx * Ix - gamma * Ix - v * Ix
    dFx <- v * Ix
    if (Ix > H / 100) {
      dHx <- a * c * Sx * Ix
    } else {
      dHx <- 0
    }

    dS <- -a * (1 - m1) * c * S * I
    dI <-  a * (1 - m1) * c * S * I - gamma * I - v * I
    dF <- v * I
    if (I > H / 100) {
      dH <- a * (1 - m1) * c * S * I
    } else {
      dH <- 0
    }
    list(c(Sx = dSx, Ix = dIx, Fx = dFx, Hx = dHx, S = dS,
           I = dI, Fs = dF, HS = dH))
  })
}

# TODO: is there a better name for this function?
SIR2 <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dSx <- -(1 - m1) * (1 - m2) * a * c * Sx * Ix
    dIx <- (1 - m1) * (1 - m2) * a * c * Sx * Ix - gamma * Ix - v * Ix
    dS <- -a * c * S * I
    dI <-  a * c * S * I - gamma * I - v * I
    list(c(
      Sx = dSx,
      Ix = dIx,
      S = dS,
      I = dI
    ))
  })
}


server <- function(input, output) {
  output$SIR <- renderPlot({
    gamma = 1
    chi = 0.03
    v = gamma * chi / (1 - chi)
    H = 40
    c = 7
    a = 1
    parms <- c(a = a, m1 = input$m1, c = c, gamma = gamma,
               v = v, H = H)
    H = 40
    I0 = 0.01
    S0 = 1 - I0
    mintime <- 0
    maxtime <- 6
    out <- ode(y = c(Sx = S0, Ix = I0, Fx = 0, Hx = 0, S = S0, I = I0, FS = 0, HS = 0), times = seq(mintime, maxtime, .1), SIR, parms)
    df <- data.frame(out)

    areaAlpha <- 0.6
    g1 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Ix * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = I * 100), fill = '#b2df8a', alpha = areaAlpha) +
      geom_hline(aes(yintercept = H), alpha = 0.2) +
      labs(x = NULL, y = NULL, title = "Percent of population infected")

    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = Fx * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = FS * 100), fill = '#b2df8a', alpha = areaAlpha) +
      labs(x = "time (months)", y = NULL, title = "Cumulative fatalities")

    # xTODO (Alec #2): I would also like to print out R_0 1, R_2
    # doubling time 1, doubling time 2, and final size....
    # TODO (AH): update this placeholder
    R_0 <- 10
    R_2 <- 12
    toprint <- data.frame(R_0, R_2)


    (g1 /
      g2 &
      scale_y_continuous(expand = expand_scale(c(0, 0.1)), labels = function(x) paste0(x,"%")) &
        scale_x_continuous(expand = c(0, 0)) ) /
      tableGrob(toprint, rows = NULL, theme = ttheme_minimal())
  })

  output$more <- renderPlot({
    gamma = 0.5
    v = gamma * input$chi / (1 - input$chi / 100) / 100
    parms <- c(a = input$a, m1 = input$m1, c = input$c, m2 = input$m2, gamma = gamma, v =
                 v)
    H = input$H
    I0 = 0.01
    S0 = 1 - I0
    out2 <- ode(y = c(Sx = S0, Ix = I0, S = S0, I = I0), times = seq(0, 12, .1), SIR2, parms)
    df2 <- data.frame(out2)

    areaAlpha <- 0.6
    ggplot(df2, aes(x = time)) +
      geom_area(aes(y = Ix * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = I * 100), fill = '#b2df8a', alpha = areaAlpha) +
      geom_hline(aes(yintercept = H), alpha = 0.2) +
      labs(x = NULL, y = NULL, title = "Percent of population infected") + scale_y_continuous(expand = expand_scale(c(0, 0.1)), labels = function(x) paste0(x,"%")) +
      scale_x_continuous(expand = c(0, 0))
  })
} # End server function

# Define UI for app that draws a histogram ----
ui <- fluidPage(title = "The math behind flatten the curve",
  fluidRow(
    column(6,
           h1("The math behind flatten the curve"),
           p("by Amy Hurford, Alec Robitialle, and Joseph Baafi (Memorial University)"),
           p("Anyone interested in contributing should contact ahurford-at-mun-dot-ca.")),
    column(6,
             p("Have you heard the remark:"),
             p("\"We'll never know the effect that social distancing has had;
                we'll never know how many lives were saved\""
             ),
             p("We can never know this with absolute certainty,
                but we can get some idea using mathematical models. ")
           )
    ),
  tabsetPanel(
    tabPanel("Social distancing",
             column(4,
           p("Below we show that the 'flatten the curve' graphic arises from a mathematical model:"),

p(tags$a(href = "https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model", "The SIR equations")),
        p("The 'flatten the curve' graphic is not simply a drawing of an idea,
        rather, rather it is based on epidemiological characteristics such the duration of infectivity,
        and the disease mortality rate. The SIR equations have been used for decades,
and today, these equations remain a good start for understanding the time course of epidemics. The graphics below
are the computational output due to solving the SIR equation.
        "),

           # Input: Slider for the number of bins ----
           sliderInput("m1", "social distancing:",
                       min = 0, max = 1, step = 0.01, value = .2,
                       width = '100%'),helpText("0: no efforts to enact social distance"),helpText("1: fully effective isolation"), ),
           column(8,

           # Output:
           plotOutput("SIR"),
           helpText("Cumulative fatalities does not account for increased death rate when health resourses are exceeded")

           )
    ),
    tabPanel("More models",

           # TODO (AH): update this p with "to the left" etc
           column(4,
             p("Above we showed that the flatten the curve graph arises from a well-established epidemiological
        model. However, the shape of the curves depend characteristics of the disease. Below we let you choose
        the characteristics of the disease.
        "),
             p("Similar resources: https://alhill.shinyapps.io/COVID19seir/"),
             # TODO (AH): where does this go?
             p("In my CBC St. John's Morning Show talk, I discussed exponential growth"),

           # Input: Slider for the number of bins ----
           sliderInput("H", "hospital capacity:", min = 0, max = 1, value = 0.4),
           sliderInput("a", "hygiene:", min = 0, max = 1, step = 0.01, value = .5),
           sliderInput("m1", "hygiene improvement factor:", min = 0, max = 1, step = 0.01, value = .2),
           sliderInput("c", "contact rate:", min = 0, max = 10, step = 0.1, value = 10),
           sliderInput("m2", "social distancing improvement factor:", min = 0, max = 1, step = 0.01, value = .2),
           sliderInput("chi", "case fatality (%):", min = 0, max = 10, step = 0.1, value = 3)),
           column(8,
                  plotOutput("more")
                  # TODO (AH): include second plot here
       )
  ),
  tabPanel("Newfoundland",

           column(4,

                  # TODO: sidebar content

                  ),

           column(8,

                  # TODO: main content
                  )

           )
))


# Alec #3: I would like a layout with: a sidebar and a plot, and then another sidebar and a plot

server <- function(input, output) {
  output$SIR <- renderPlot({
    gamma = 0.5
    v = gamma*input$chi/(1-input$chi/100)/100
    parms <- c(a=input$a,m1=input$m1,c=input$c, m2=input$m2, gamma=gamma, v=v)
    H=input$H
    I0 = 0.01
    S0 = 1-I0
    out <- ode(y = c(Sx=S0, Ix=I0, S=S0, I=I0), times=seq(0, 12, .1), SIR, parms)
    #matplot.0D(out)
    plot(out[,1], out[,3], typ ="l", ylim = c(0,S0), ylab = "prop. of population infected", xlab = "time (months)")
    lines(out[,1], out[,5], lty=2, col = "red")
    lines(c(0,12), c(H,H))
    text(12,0.8, "Final size: print")
  })
} # End server function
#
# # Define UI for app that draws a histogram ----
# ui <- fluidPage(
#
#   # App title ----
#   titlePanel("The math behind flatten the curve"),
#
#   sidebarLayout(
#
#     # Sidebar panel for inputs ----
#     sidebarPanel(
#
#       # Input: Slider for the number of bins ----
#       sliderInput("H", "hospital capacity:", min = 0, max = 1, value = 0.4),
#       sliderInput("a", "hygiene:", min = 0, max = 1,step=0.01, value = .5),
#       sliderInput("m1", "hygiene improvement factor:", min = 0, max = 1,step=0.01, value = .2),
#       sliderInput("c", "contact rate:", min = 0, max = 10,step=0.1, value = 10),
#       sliderInput("m2", "social distancing improvement factor:", min = 0, max = 1,step=0.01, value = .2),
#       sliderInput("chi", "case fatality (%):", min = 0, max = 10, step = 0.1, value = 3)
#     ),
#
#     # Main panel for displaying outputs ----
#     mainPanel(
#
#       p("Below we show that the 'flatten the curve' graphic can be
#         produced from mathematical equations that explain the time course of
#         an epidemic: an SIR model (link to wikipedia).
#         "),
#
#       # Output:
#       plotOutput("SIR")
#   )# end main panel
# )# end layout
# )# end ui

shinyApp(ui=ui, server=server)

