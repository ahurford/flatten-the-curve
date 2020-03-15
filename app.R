library("deSolve")
library("shiny")

SIR <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dSx <- -(1-m1)*a*c*Sx*Ix
    dIx <- (1-m1)*a*c*Sx*Ix - gamma*Ix - v*Ix
    dS <- -a*c*S*I 
    dI <-  a*c*S*I - gamma*I - v*I
    list(c(Sx=dSx, Ix=dIx, S=dS, I=dI))
  })
}

server <- function(input, output) {
  output$SIR <- renderPlot({
    gamma = 0.5
    chi = 0.03
    v = gamma*chi/(1-chi)
    H = 0.4
    c=5
    a=1  
    parms <- c(a=a,m1=input$m1,c=c, gamma=gamma, v=v)
    H=0.4
    I0 = 0.01
    S0 = 1-I0
    out <- ode(y = c(Sx=S0, Ix=I0, S=S0, I=I0), times=seq(0, 12, .1), SIR, parms)
    #matplot.0D(out)
    plot(out[,1], out[,3], typ ="l", ylim = c(0,S0), ylab = "prop. of population infected", xlab = "time (months)")
    lines(out[,1], out[,5], lty=2, col = "red")
    lines(c(0,12), c(H,H))
    text(12,0.8, "Final size: print")
    # Alec #1: can you make the plot a bit prettier.
    # Alec #2: I would also like to print out R_0 1, R_2,
    # doubling time 1, doubling time 2, and final size. I can put the formulas for
    # these in later, if you can make some dummy outputs appear.
  })
} # End server function

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("The math behind flatten the curve"),
  
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput("m1", "social distancing (0=none ---> 1=complete isolation):", min = 0, max = 1,step=0.01, value = .2)
      # Alec do you know how to make a slider note?
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      p("Have you heard the remark:"),
      
      p("We'll never know the effect that social distancing has had;
        we'll never know how many lives were saved?"
        ),
      
      p("But, while we can never know for sure, we can get some idea using mathematical
        models.
        "),
      
      p("Below we show that the 'flatten the curve' graphic arises from a mathematical model: the SIR equations (link to wikipedia).
        The lines in the 'flatten the curve' graphic are not simply a drawing of an idea;
        rather, they are based on epidemiological characteristics such the duration of infectivity,
        and the disease mortality rate.
        "),
      
      # Output:
      plotOutput("SIR"),
      
      p("Above we showed that the flatten the curve graph arises from a well-established epidemiological
        model. However, the shape of the curves depend characteristics of the disease. Below we let you choose
        the characteristics of the disease.
        
        "),
      
      p("In my CBC talk, I discussed exponential growth"),
      
      p("This app was made by XXX, and anyone interested in contributing should contact ahurford-at-mun-dot-ca")
      
      )# end main panel
    )# end layout
  )# end ui


# Alec #3: I would like a layout with: a sidebar and a plot, and then another sidebar and a plot
# This code below is the start of what I would like for the SECOND sidebar and plot.

# SIR <- function(t, y, p) {
#   with(as.list(c(y, p)), {
#     dSx <- -(1-m1)*(1-m2)*a*c*Sx*Ix
#     dIx <- (1-m1)*(1-m2)*a*c*Sx*Ix - gamma*Ix - v*Ix
#     dS <- -a*c*S*I 
#     dI <-  a*c*S*I - gamma*I - v*I
#     list(c(Sx=dSx, Ix=dIx, S=dS, I=dI))
#   })
# }
# 
# server <- function(input, output) {
#   output$SIR <- renderPlot({
#     gamma = 0.5
#     v = gamma*input$chi/(1-input$chi/100)/100
#     parms <- c(a=input$a,m1=input$m1,c=input$c, m2=input$m2, gamma=gamma, v=v)
#     H=input$H
#     I0 = 0.01
#     S0 = 1-I0
#     out <- ode(y = c(Sx=S0, Ix=I0, S=S0, I=I0), times=seq(0, 12, .1), SIR, parms)
#     #matplot.0D(out)
#     plot(out[,1], out[,3], typ ="l", ylim = c(0,S0), ylab = "prop. of population infected", xlab = "time (months)")
#     lines(out[,1], out[,5], lty=2, col = "red")
#     lines(c(0,12), c(H,H))
#     text(12,0.8, "Final size: print")
#   })
# } # End server function
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

