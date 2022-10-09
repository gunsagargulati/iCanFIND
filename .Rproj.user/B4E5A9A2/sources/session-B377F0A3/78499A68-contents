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
                                                   min = min(years),
                                                   max = max(years),
                                                   value= min(years)),
                                       shinyWidgets::pickerInput("cancer_select", "Cancer subtype:",   
                                                   choices = cancers, 
                                                   options = list(`actions-box` = TRUE, `none-selected-text` = "Please make a selection!"),
                                                   selected = cancers[1],
                                                   multiple = TRUE),
                                       shinyWidgets:: pickerInput("gene_select", "Gene:",   
                                                   choices = hugo_symbol,
                                                   options = list(`live-search`=TRUE, `none-selected-text` = "Please make a selection!"),
                                                   selected = hugo_symbol[1],
                                                   multiple = FALSE), 
                                       shinyWidgets::pickerInput(inputId = "mutselect", 
                                                   label = "Pick subgroup", 
                                                   choices = c("Inactivating","Likely inactivating", "Select my own mutations"),
                                                   options = list(`none-selected-text` = "Please make a selection!"),
                                                   selected = "Inactivating",
                                                   multiple = FALSE),
                                       shiny::conditionalPanel(condition = "input.mutselect == 'Select my own mutations'",
                                                        uiOutput("aminoacids")),
                                       shinyWidgets::pickerInput(
                                       inputId = "cnaselect",
                                       label = "Copy number:",
                                       choices = c("Any", "None", "Amplification", "Deletion"),
                                       selected = "Any",
                                       options = list(`none-selected-text` = "Please make a selection!"),
                                       multiple = FALSE),
                                       shinyWidgets::pickerInput(
                                       inputId = "fusionselect",
                                       label = "Fusion:",
                                       selected = "Any",
                                       choices = c("Any", "None", "Likely Loss-of-function", "Loss-of-function", "Likely Gain-of-function", "Gain-of-function"),
                                       options = list(`none-selected-text` = "Please make a selection!"),
                                       multiple = FALSE)),
 
                                   mainPanel(
                                     fluidPage(plotly::plotlyOutput("plot_seer1"),
                                                plotly::plotlyOutput("plot_seer2")),
                                       ))))),
             
             tabPanel("Methods",
                      div(class="outer",
                          tags$head(includeCSS("styles.css"))),
                      fluidPage("Under construction")),

             tabPanel("Data",
                      div(class="outer",
                          tags$head(includeCSS("styles.css"))),
                      fluidPage("Under construction")),
             tabPanel("Contact Us",
                      div(class="outer",
                          tags$head(includeCSS("styles.css"))),
                      fluidPage("Under construction"))
             
             
  ))



