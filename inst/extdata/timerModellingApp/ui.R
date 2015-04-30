library(shiny)
library(shinyBS)

shinyUI(pageWithSidebar(

    headerPanel("Timer modelling app"),
    sidebarPanel(
        selectInput("prodModel", "Production Model:", list("constant", "linear", "step"), selected="step"),
        numericInput("p0", "Protein production rate (molecules per minute)", 15, min=0, max=10000),
        numericInput("dp", "Total change production (molecules per minute)", 10, min=-10000, max=10000),
        bsTooltip("dp", title="Set to be positive for production increases and negative for production decreases."),
        numericInput("startProd", "Start time production change (minutes):", 100, min=0),
        bsTooltip("startProd", title="When the change in production starts to take effect."),
        numericInput("endProd", "End time production change (minutes):", 110, min=0),
        bsTooltip("endProd", title="When the change in production finishes. For bursts select step model and set this to be close to the start time of production changes."),
        selectInput("degModel", "Degradation Model:", list("constant", "linear", "step"), selected="constant"),
        numericInput("k0", "Degradation rate (1/minutes)", 0.01, min=log(2)/(1/60), max=log(2)/(7*24*60)),
        bsTooltip("k0", "Convert to half-life using the transformation log(2)/k, where k is the degradation rate."),
        numericInput("dk", "Total change degradation (1/minutes)", -0.005, min=-50, max=50),
        bsTooltip("dk", "Set to be positive for degradation rate increases (removes more proteins) and negative for degradation rate decreases (removes less proteins)"),
        numericInput("startDeg", "Start time degradation change (minutes):", 0, min=0),
        bsTooltip("startDeg", title="When the change in degradation rate starts to take effect."),
        numericInput("endDeg", "End time degradation change (minutes):", 600, min=0),
        bsTooltip("endDeg", title="When the change in degradation rate finishes. For bursts select step model and set this to be close to the start time of degradation changes."),
        numericInput("t1", "Fluorophore 1 maturation time (minutes)", 5, min=0, max=1000),
        numericInput("t2", "Fluorophore 2 maturation time (minutes)", 60, min=0, max=1000),
        numericInput("Tf", "Maximum time (minutes)", min=0, max=7*24*60, value=10*60, step=60),
        bsTooltip("Tf", "Simulation end time. Increase to see how readouts develop at future times.")
    ),

    mainPanel(
        helpText("This app accompanies the paper 'TimerQuant: A modelling approach to tandem fluorescent timer design and data interpretation for measuring protein turnover in embryos' (2015). The purpose is to guide users on how to interpret timer readouts when production and degradation rates are changing over time."),
        plotOutput("modelPlot"),
        bsPopover("modelPlot", "Production and degradation rate profiles", "These plots show production and degradation behaviour as a function of time. Use the side panel to input different parameters and manipulate the shape of the profiles.", placement="left"),
        plotOutput("solutionPlot"),
        bsPopover("solutionPlot", "Model solutions", "Timer readouts are shown as a function of time. FP1 (fast-maturing fluorophore) profiles are shown in green, FP2 (slow-maturing fluorophore) in red, and the FP2/FP1 ratio in blue. Use both green and blue profiles together when comparing different models.", placement="left")
    )
))

