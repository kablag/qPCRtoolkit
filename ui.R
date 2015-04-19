library(shiny)
library(ggvis)

shinyUI(fluidPage(
  includeCSS("www/css/styles.css"),
  includeScript("www/js/responsiveTable.js"),
  # Application title
  titlePanel("qPCRtoolkit"),
  
  sidebarLayout(    
    # Sidebar with a slider input
    sidebarPanel(
      
      conditionalPanel(
        condition = "input.mainTabsetPanel == 'experiments'",
        fileInput("rdmlFile",
                  "Add File...")),
      #       shinyFilesButton("rdmlFile", 
      #                        "Добавить файл...", 
      #                        'Please select a file', multiple = TRUE)
      #       uiOutput("useFilesListUi"),
      conditionalPanel(
        condition = "input.mainTabsetPanel == 'absQuant'",
        uiOutput("limitCyclePanel"),
        uiOutput("backgroundPanel"),
        checkboxInput("useBackgroundSubstruction",
                      "Apply Background Substraction",
                      TRUE),
        uiOutput("thresholdPanel"),
        
        radioButtons("showPlot", "Show:",
                     c("Amplification" = "ampCurve",
                       "Calibration" = "calib"),
                     selected = "ampCurve"),
        uiOutput("filterTargetsUi"),
        downloadButton("downloadCq", "Export .xlsx")
      ),
      conditionalPanel(
        condition = "input.mainTabsetPanel == 'relQuant'",
        wellPanel(
          uiOutput("refGenesUi"),
          uiOutput("runCalibUi"),
          uiOutput("studyCalibUi")),
        wellPanel(
          uiOutput("showRelGenesUi"),
          uiOutput("showRelSamplesUi"),
          uiOutput("showRelConditionsUi"),
          selectInput("splitOpts", "Split by",
                      choices = c("Sample" = "sample",
                                  "Target" = "target",
                                  "Condition" = "conditions"),
                      selected = c("Sample" = "sample",
                                   "Target" = "target",
                                   "Condition" = "conditions"),
                      multiple = TRUE),
          radioButtons("displayRatio", "Show Ratio:",
                       c("Relative" = "ratio_mean",
                         "Normalized" = "normalized.ratio_mean",
                         "Scaled" = "scaled.ratio_mean"),
                       selected = "ratio_mean"),
          radioButtons("curveTypeRel", "Curve Type:",
                       c("Linear" = "linear",
                         "Log" = "log"),
                       selected = "linear"),
          checkboxInput("showErrRel",
                        "Show Errors", value = FALSE)),
        downloadButton("downloadRel", "Export .xlsx")
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        id = "mainTabsetPanel",
        tabPanel("Experiments",
                 value = "experiments",
                 uiOutput("experimentsChecksUi")),
        tabPanel("Absolute Quantification",
                 value = "absQuant",
                 conditionalPanel(
                   condition = "input.showPlot == 'ampCurve'",
                   ggvisOutput("curvesPlot"),
                   fluidRow(
                     column(6,
                            uiOutput("plateTblExpSelectorUi"),
                            uiOutput("plateTbl")),
                     
                     column(3,
                            tags$div(class = "verticalLine",
                                     radioButtons("colorBy", "Color by:",
                                                  c("Type" = "sample.type",
                                                    "Target" = "target"),
                                                  selected = "target"))),
                     column(3,
                            radioButtons("curveType", "Curve Type:",
                                         c("Linear" = "linear",
                                           "Log" = "log"),
                                         selected = "linear")))
                 ),
                 conditionalPanel(
                   condition = "input.showPlot == 'calib'",
                   ggvisOutput("calibrationPlot"),
                   htmlOutput("calibText")
                 ),
                 dataTableOutput("cqConcTbl")
        ), 
        tabPanel("Relative Quantification",
                 value = "relQuant",
                 plotOutput("relativeQuantPlot"),
                 dataTableOutput("relQuantTbl")
        )
      )
    )
  )
))