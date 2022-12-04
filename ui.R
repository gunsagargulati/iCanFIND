data <- readRDS("./www/data.rds")
seer <- data$seer
mut_prediction <- data$mut_prediction
cna_type <- data$cna_type
fusion_type <- data$fusion_type
years <- data$years
cancers <- data$cancers
seer <- data$seer
hugo_symbol <- data$hugo_symbol
source("www/functions.R")

ui <- bootstrapPage(
  navbarPage(theme = shinythemes::shinytheme("flatly"), collapsible = TRUE,
             title = div(img(src = "iCanFindlogo.png",
                   height = 60,
                   style = "margin-top: -14px;
                               padding-right:10px;
                               padding-bottom:10px"
                 )),
                 
             windowTitle = "iCanFIND",
             
             tabPanel("Home",
                      div(class="outer",
                          tags$head(includeCSS("styles.css")),
                          
                          tabPanel("Region plots",
                                   
                                   sidebarLayout(
                                     sidebarPanel(
                                       
                                       span(tags$i(h6("Please select from the menus below to generate bar plots on the right.")), style="color:#045a8d"),

                                       shiny::sliderInput("year_select",
                                                   "Year:",
                                                   sep= "",
                                                   min = min(years),
                                                   max = max(years),
                                                   value= 2022),
                                       shinyWidgets::pickerInput("cancer_select", "Cancer subtype:",   
                                                   choices = cancers, 
                                                   options = list(`actions-box` = TRUE, `none-selected-text` = "Please make a selection!"),
                                                   selected = "BLCA",
                                                   multiple = TRUE),
                                       shinyWidgets:: pickerInput("gene_select", "Gene:",   
                                                   choices = hugo_symbol,
                                                   options = list(`live-search`=TRUE, `none-selected-text` = "Please make a selection!"),
                                                   selected = "TSC1",
                                                   multiple = TRUE), 
                                       shinyWidgets::pickerInput(inputId = "mutselect", 
                                                   label = "Pick subgroup", 
                                                   choices = c("None", "Inactivating","Likely inactivating", "Select my own mutations"),
                                                   options = list(`none-selected-text` = "Please make a selection!"),
                                                   selected = "Inactivating",
                                                   multiple = FALSE),
                                       shiny::conditionalPanel(condition = "input.mutselect == 'Select my own mutations'",
                                                        uiOutput("aminoacids")),
                                                        uiOutput("cnawidget"),
                                                        uiOutput("fusionwidget")),
                                   mainPanel(
                                     tabsetPanel(type = "tabs",
                                     tabPanel("Frequency",
                                     fluidPage(plotly::plotlyOutput("plot_seer1"),
                                               downloadButton('resultsDL', label = "Download frequency"))),
                                     tabPanel("Incidence",
                                    fluidPage(plotly::plotlyOutput("plot_seer2"),
                                              downloadButton('resultsSEERDL', label = "Download incidence")))),
                                     
                                       ))))),

             tabPanel("Data",
                      div(class="outer",
                          tags$head(includeCSS("styles.css"))),
                      fluidPage(
                        titlePanel("Download data", windowTitle = "Data"),
                        mainPanel(
                          p(
                            "The following files contain processed genomic data stored as .rds (RDS = R data single files), which can be accessed via R. The datasets include samples from TCGA (n = 9,809) and the AACR GENIE database (e.g. MSK, n = 15,333 and DFCI, n = 5,720) with paired mutation, copy number, and structural variant data. "
                          ),
                          htmltools::HTML('<br>'),
                          downloadButton('tcgaDL', label = "Download TCGA data"),
                          htmltools::HTML('<br><br>'),
                          downloadButton('mskpDL', label = "Download MSK data"),
                          htmltools::HTML('<br><br>'),
                          downloadButton('dfciDL', label = "Download DFCI data"),
                          htmltools::HTML('<br>'),
                          htmltools::HTML('<hr style="color: black;">'),
                          p(
                            "The following file includes demographic data from the 30,862 samples above. "
                          ),
                          htmltools::HTML('<br>'),
                          downloadButton('demoDL', label = "Download demographic data"),
                          htmltools::HTML('<br>'),
                          htmltools::HTML('<hr style="color: black;">'),
                          p(
                            "The following file include processed data from SEER, predicting annual U.S. incidence of 31 cancer subtypes from 2020-2045. "
                          ),
                          htmltools::HTML('<br>'),
                          downloadButton('seerDL', label = "Download SEER data")
                        ))),
             tabPanel("Contact Us",
                      div(class="outer",
                          tags$head(includeCSS("styles.css"))),
                      fluidPage(
                        titlePanel("Contact us", windowTitle = "Contact us"),
                        mainPanel(
                          htmltools::HTML('<br>'),
                          htmltools::HTML('<hr style="color: black;">'),
                          p(
                            "If you have any questions or suggestions, please email us directly at gsgulati@bwh.harvard.edu. "
                          )
             
             
  )))))



