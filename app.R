#app.R

#Libraries
library(shinythemes)
library(shinyWidgets)
library(data.table)
library(lubridate) 
library(plotly)
library(shiny)
library(ggplot2)
library(dplyr)
library(htmltools)


source("ui.R")
source("server.R")

runApp(ui = ui, server = server)
