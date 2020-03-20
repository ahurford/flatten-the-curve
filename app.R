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

### To deploy ----
#rsconnect::deployApp("/Users/amyhurford/Desktop/flatten-the-curve")
# in R Console

### Theme ----
theme_set(theme_light() +
            theme(
              axis.text = element_text(size = 14),
              legend.text = element_text(size = 14)
            ))

cols <- c('No changes implemented' = '#a6cee3',
          'With social distancing' = '#b2df8a',
          'Hospital capacity' = 'grey')

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

SIHR <- function(t, y, p) {
  with(as.list(c(y, p)), {
  	# The x subscript indicates "no social distancing".
    dSx <- -a * c * Sx * Ix
    dIx <- a * c * Sx * Ix - gamma * Ix - sigma*Ix - v*Ix
    # Hospitalized without a capacity cap for reference
    dHx1 <- sigma*Ix - vH*Hx1 - rho*Hx1
    dCx <- a * c * Sx * Ix
    dHcumx <- sigma*Ix
# Much better to code the below as min/max. Poor computation times
    # with if esle.
    	# Hospitalized with a capacity cap
    # with cap on hospital admissions
    dHx <- sigma*min(Ix,H2) - vH*H0 - rho*H0
    # cumulative unmet need
    dUx <- sigma*max(0,(Ix-H2))

    # with social distancing
    dS <- -a * (1 - m2) * c * S * I
    dI <-  a * (1 - m2) * c * S * I - gamma * I - v * I - sigma*I
    dHcum <- sigma*I
    # without a cap on hospital resources
    dH1 <- sigma*I - vH*H1 - rho*H1
    # culmulative cases
    dC <- a * (1 - m2) * c * S * I
    # with cap on hospital admissions
    dH0 <- sigma*min(I,H2) - vH*H0 - rho*H0
    # cumulative unmet need
    dU <- sigma*max(0,(I-H2))
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
    parms <- c(a = a, m1 = input$m1/100, c = c, gamma = gamma,
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
        scale_color_manual(values = cols) &
        guides(fill = guide_legend(override.aes = list(linetype = 0),
                                   nrow = 1),
               color = guide_legend(override.aes = list(fill = 'white')))) /
      guide_area() +
      plot_layout(guides = 'collect', heights = c(5, 5, 2))
  })

  output$SIRtab <- renderTable({
    # Generate data.frame to print
    R_0 <- round(a * c / (v + gamma), 1)
    R_2 <- round((1 - input$m1/100) * R_0, 1)
    DT <- round(log(2) / (a * c - v - gamma), 1)
    DT_2 <- max(round(log(2) / ((1 - input$m1) * a * c - v - gamma), 1), 0)
    fat <- round(100 * out[length(out[, 1]), 4], 1)
    fat_2 <- round(100 * out[length(out[, 1]), 7], 1)
    data.frame(
        " " = c("no distancing", "with distancing"),
        "doubling time" = c(DT, DT_2),
        "R0" = c(R_0, R_2),
        "fatalities" = c(fat, fat_2),
        check.names = FALSE
      )
  })

  output$SIHR <- renderPlot({
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

    parms <- c(a = a, m2 =input$m2/100, c = c, gamma = gamma,
               v = v, H2 = input$H2, rho = rho, vH = vH, sigma = sigma)
    I0 <- 0.005
    S0 <- 1 - I0
    mintime <- 0
    maxtime <- 250
    out <- ode(y = c(Sx = S0, Ix = I0, Hx1 = 0, Hx=0, Ux=0, Cx=0, S = S0, I = I0, H1 = 0, H0 = 0, U = 0, C=0, Hcum=0, Hcumx=0), times = seq(mintime, maxtime, 1), SIHR, parms)
    df <- data.table(out)

    areaAlpha <- 0.6

    g1 <- ggplot(df, aes(x = time)) +
    	geom_area(aes(y = Ix * 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
    	geom_area(aes(y = I * 100), fill = '#b2df8a', alpha = areaAlpha) +
    	#geom_hline(aes(yintercept = input$H2/100), alpha = 0.2, size = 3) +
    	labs(x = NULL, y = NULL, title = "Infected (% of population)")

    g2 <- ggplot(df, aes(x = time)) +
      geom_area(aes(y = 100*Hx1), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
      geom_area(aes(y = 100*H1), fill = '#b2df8a', alpha = areaAlpha) +
    	geom_hline(aes(yintercept = input$H2), alpha = 0.2, size = 3) +
      labs(x = NULL, title = "Requiring critical care (% of population)", y= NULL)

    g3 <- ggplot(df, aes(x = time)) +
    	geom_area(aes(y = Ux* 100), fill = '#a6cee3', alpha = areaAlpha - 0.2) +
    	geom_area(aes(y = U * 100), fill = '#b2df8a', alpha = areaAlpha) +
    	#geom_hline(aes(yintercept = input$H2/100), alpha = 0.2, size = 3) +
    	labs(x = "time (days)", y = NULL, title = "Cumulative unmet need (% of population)")

    final.unmet.x = round(last(df$Ux)*100,2)
    final.unmet = round(last(df$U)*100,2)
    final.hosp.x = round(last(df$Hcumx)*100,2)
    final.cases.x = round(last(df$Cx)*100,0)
    final.cases = round(last(df$C)*100,0)
    final.hosp = round(last(df$Hcum)*100,2)

    toprint <-
      data.frame(
        " " = c("no distancing", "with distancing"),
        "Final unmet need (%)" = c(final.unmet.x, final.unmet),
        "Final critical care need (%):" = c(final.hosp.x, final.hosp),
        "Final infected (%)" = c(final.cases.x, final.cases),
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


  output$scrapeTab <- renderTable({
    invalidateLater(24 * 60 * 60 * 1000)
    data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv')[
      Province == 'NL'][, .(Province, Date, negative, presumptive_negative, presumptive_positive, confirmed_positive)]
  })

  output$scrapePlot <- renderPlot({
    invalidateLater(24 * 60 * 60 * 1000)
    NL <- data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv')[Province == 'NL']

    for (j in names(NL)) set(NL, which(is.na(NL[[j]])), j, 0)

    # Plot cases in NL
    cols <- c('Presumptive Positive' = '#881a58',
              'Confirmed Positive' = '#0e288e')
    ggplot(NL, aes(x = Date, group = 1)) +
      geom_line(aes(y = presumptive_positive, color = 'Presumptive Positive'),
                show.legend = TRUE) +
      geom_point(aes(y = presumptive_positive), color = '#881a58') +
      geom_line(aes(y = confirmed_positive, color = 'Confirmed Positive'),
                show.legend = TRUE) +
      geom_point(aes(y = confirmed_positive), color = '#0e288e') +
      labs(x = NULL, y = NULL, title = "Cases in NL", color  = NULL) +
      scale_color_manual(values = cols) +
      scale_y_continuous(#expand = expand_scale(mult = c(0, 0.1)),
                         limits = c(0, max(c(NL$presumptive_positive, NL$confirmed_positive), na.rm = TRUE) + 2))
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

    # tabPanel('Introduction',
             # )

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
           helpText("Doubling time: Early on in the epidemic, the time for the number of infected people to double."),
           helpText("R0: Early on in the epidemic, the average number of people subsequently infected by an infected person."),
           helpText("Fatalities: The percentage of the population that has died from COVID-19 after 250 days, however this
                    does not consider an increased death rate when the hospital capacity is exceeded. Epidemic models, such as SIR, suggest
                    that, even aside from preventing exceeding the hospital capacity, a smaller percentage of the population
                    will die from COVID-19 under social distancing."),
           helpText("The parameterization for this SIR model was taken from Bolker and Dushoff (2020)."))),
    tabPanel("Your questions",
             tabPanel("Your questions",
                      column(5,
                             p(""),

                             # Slider input: social distancing
                             sliderInput("m2", "social distancing (%):",
                                         min = 0, max = 100, step = 1, value = 20,
                                         width = '100%'),
                             helpText("The green curve shows the effect of social distancing"),
                             helpText("0%: no efforts to enact social distancing"),
                             helpText("100%: complete isolation"),
                             sliderInput("H2", "Hospital capacity (% of population):",
                                         min = 0, max = 0.3, step = .01, value = 0.2,
                                         width = '100%'),
                             helpText("Move the slider and the grey line will change"),

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
                              the percentage of people requiring critical care at anytime is much smaller (<0.3% of the population).
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
                             helpText("Blue curve: no social distancing; Green curve: with social distancing; Grey line: hospital capacity."),
                             helpText("Cumulative unmet need is the percentage of the population that had an unmet need because they required
                                      critical care, when the capacity to provide critical care was exceeded. For the default parameters no
                                      green curve appears because the hospital capacity is never exceeded with social distancing at 20%."),
                             helpText("Final unmet need (%): after 250 days, the percentage of the population that had an unmet need for critical care."),
                             helpText("Final critical care need (%): after 250 days, the percentage of the population that required critical care"),
                             helpText("Final infected (%): after 250 days, the percentage of the population was infected. Note these
                                      numbers are much too high. Please see the 'More models' tab for a disclaimer."),
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
    									# this link looks good: http://gabgoh.github.io/COVID/index.html
    # Newfoundland tab
    tabPanel("Newfoundland",
             column(12,
                    p(""),
                    p("We aim to make some Newfoundland-specific graphs and analysis, but this work
                    is currently in progress")
                    ),

             plotOutput("scrapePlot", width = "55%"),
             tags$br(),
             tags$br(),
             tags$br(),
             tableOutput("scrapeTab"))
))

### Run app ----
shinyApp(ui = ui, server = server)

