library("deSolve")
library("shiny")

SIR <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dSx <- -betax*Sx*Ix
    dIx <- betax*Sx*Ix - gamma*Ix
    dS <- -beta*S*I
    dI <- beta*S*I - gamma*I
    list(c(Sx=dSx, Ix=dIx, S=dS, I=dI))
  })
}

server <- function(input, output) {
  output$SIR <- renderPlot({
    parms <- c(beta=input$beta,betax=input$betax,gamma=0.05)
    N0=input$N0
    H=input$H
    I0 = 1/N0/1000
    S0 = 1-I0
    out <- ode(y = c(Sx=S0, Ix=I0, S=S0, I=I0), times=seq(0, 12, .1), SIR, parms)
    #matplot.0D(out)
    plot(out[,1], out[,3], typ ="l", ylim = c(0,S0), ylab = "prop. of population infected", xlab = "time (months)")
    lines(out[,1], out[,5], lty=2, col = "red")
    lines(c(0,12), c(H,H))
  })
}

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # App title ----
  titlePanel("Hello Shiny!"),
  
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput("H", "hospital capacity:", min = 0, max = 1, value = 0.4),
      numericInput("N0", "population size (thousands):", value = 120, step = 100),
      sliderInput("beta", "beta:", min = 0, max = 5,step=0.1, value = 5),
      sliderInput("betax", "betax:", min = 0, max = 5,step=0.1, value = 5), width = 4
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output:
      plotOutput("SIR")
  )
)
)

shinyApp(ui=ui, server=server)

