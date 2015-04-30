library(shiny)
library(shinyBS)

shinyUI(pageWithSidebar(
    headerPanel("App for choosing optimal FP2"),
    sidebarPanel(
        numericInput("p", "protein production rate p (molecules per minute)", 0.5, min=0, max=Inf),
        bsTooltip("p", "The number of proteins that are produced per minute"),
        numericInput("sigma", "additive noise sigma", 1, min=0, max=Inf),
        bsTooltip("sigma", "Additive noise ill be drawn from N(0,sigma^2)"),
        numericInput("nRealizations", "Number of simulation realizations", 100, min=1, max=1000),
        bsTooltip("nRealizations", "Run more simulations to get smoother profiles"),
        numericInput("tA", "Protein half-life A (minutes)", 30, min=0, max=1000),
        bsTooltip("tA", "Lower bound of protein half-lives in the experiment"),
        numericInput("tB", "Protein half-life B (minutes)", 180, min=0, max=1000),
        bsTooltip("tB", "Upper bound of protein half-lives in the experiment"),
        numericInput("t1", "Fluorophore 1 maturation time (minutes)", 5,
                     min=0, max=1000),
        bsTooltip("t1", "Maturation time of the faster-maturing fluorophore (FP1 in paper). Convert to maturation rate using transformation log(2)/T, where T is the maturation time."),
        sliderInput("t2range", "Fluorophore 2 maturation time range (minutes)",
                    min=0, max=1000, value=c(20, 200)),
        bsTooltip("t2range", "Maturation time range of the slow-maturing fluorophore. Choose a range which covers the maturation times of the fluorophores that you are considering using."),
        numericInput("E", "FRET coefficient (Fluorophore 1 -> 2 transfer)", 
                    min=0, max=1, value=0, step=0.1),
        bsTooltip("E", "The proportion of FP1 fluorescence that is transferred to FP2."),
        numericInput("f", "Prefactor f", min=0, value=1),
        bsTooltip("f", "Equal to f2/f1, where f1 and f2 are prefactors for FP1 and FP2 respectively that relate fluorescence intensity to the number of fluorescent molecules. The calculated of optimal FP2 maturation time is unaffected by this number.")
    ),

    mainPanel(
        helpText("This app accompanies the paper 'TimerQuant: A modelling approach to tandem fluorescent timer design and data interpretation for measuring protein turnover in embryos' (2015). The purpose is to guide users on how to choose the FP2 (slow-maturing fluorophore) so that they find the best tradeoff between signal and noise."),
        plotOutput("timerSignal", height="800px", width="500px")
    )
))

