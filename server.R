library(shiny)
library(ggplot2)

server <- function(input, output, session) {
source("www/functions.R")
#data
  data <- readRDS("./www/data.rds")
  seer <- data$seer
  mut_prediction <- data$mut_prediction
  cna_type <- data$cna_type
  fusion_type <- data$fusion_type
  years <- data$years
  cancers <- data$cancers
  seer <- data$seer
  hugo_symbol <- data$hugo_symbol
  
  withProgress(message = 'Loading packages, functions, and datasets', value = 0, {
    incProgress(1/4)
    
  tcga <- readRDS("www/TCGA_variants2.rds")
  incProgress(1/4)
  
  dfci <- readRDS("www/DFCI_variants2.rds")
  incProgress(1/4)
  
  mskp <- readRDS("www/MSKP_variants2.rds")
  incProgress(1/4)
  
  geneMuts <- readRDS("www/geneMuts.rds")
  
  Sys.sleep(12)
  })
  
mutselect <- reactive({input$mutselect})
aa_select <- reactive({input$aa_select})
gene_select <- reactive({input$gene_select})
fusionselect <- reactive({input$fusionselect})
cnaselect <- reactive({input$cnaselect})
cancer_select <- reactive({input$cancer_select})
year_select <- reactive({input$year_select})


output$aminoacids <- shiny::renderUI({
  aa <- as.character(unlist(geneMuts[gene_select()]))
  aa <- aa[which(nchar(aa) > 0)]
  aa <- sort(aa)
  shinyWidgets::pickerInput(inputId = "aa_select",
              label = "Select mutation:",
              choices = aa,
              selected = NULL,
              options = list(`actions-box` = TRUE, `none-selected-text` = "Select mutations", `live-search` = TRUE, title = "Select mutations"),
              multiple = TRUE)
})



  
tcga_mut <- reactive({mut_value(tcga[[1]], gene_select(), mutselect(), aa_select())})
tcga_cna <- reactive({cna_value(tcga[[2]], gene_select(), cnaselect())})
tcga_fusion <- reactive({fusion_value(tcga[[3]], gene_select(), fusionselect())})

dfci_mut <- reactive({mut_value(dfci[[1]], gene_select(), mutselect(), aa_select())})
dfci_cna <- reactive({cna_value(dfci[[2]], gene_select(), cnaselect())})
dfci_fusion <- reactive({fusion_value(dfci[[3]], gene_select(), fusionselect())})

mskp_mut <- reactive({mut_value(mskp[[1]], gene_select(), mutselect(), aa_select())})
mskp_cna <- reactive({cna_value(mskp[[2]], gene_select(), cnaselect())})
mskp_fusion <- reactive({fusion_value(mskp[[3]], gene_select(), fusionselect())})



tcga_results <- reactive({aggregate(cbind(tcga_mut(),tcga_cna(), tcga_fusion()), list(tcga[[4]]), mean, na.rm = T)})
dfci_results <- reactive({aggregate(cbind(dfci_mut(),dfci_cna(), dfci_fusion()), list(dfci[[4]]), mean, na.rm = T)})
mskp_results <- reactive({aggregate(cbind(mskp_mut(),mskp_cna(), mskp_fusion()), list(mskp[[4]]), mean, na.rm = T)})

results <- reactive({
  results_mean <- data.frame((as.matrix(tcga_results()[,-1]) + as.matrix(dfci_results()[,-1]) + as.matrix(mskp_results()[,-1]))/3)
  rownames(results_mean) <- tcga_results()[,1]
  results_mean <- data.frame(results_mean)[cancer_select(),]
  colnames(results_mean) <- c("mut", "cna", "fusion")
  results_mean <- results_mean[order(-apply(results_mean, 1, mean, na.rm = T)),]
  results_mean <- reshape2::melt(cbind.data.frame(Tumor = rownames(results_mean), results_mean))
  results_mean <- cbind.data.frame(Tumor = results_mean$Tumor, Alteration = results_mean$variable, 
                              mean = results_mean$value*100)
  return(results_mean)
})


results_seer <- reactive({
results_mean <- data.frame((as.matrix(tcga_results()[,-1]) + as.matrix(dfci_results()[,-1]) + as.matrix(mskp_results()[,-1]))/3)
rownames(results_mean) <- tcga_results()[,1]
results_mean <- data.frame(results_mean)[cancer_select(),]
colnames(results_mean)<- c("mut", "cna", "fusion")
results_mean <- results_mean[order(-apply(results_mean, 1, mean, na.rm = T)),]
seer_df <- cbind.data.frame(tumor = unlist(cancer_select()), 
                            year = unlist(year_select()), 
                            value = as.numeric(unlist(data.frame(seer[which(seer[,2] %in% cancer_select()),which(colnames(seer) %in% year_select())]))))
seer_df <- seer_df[order(-apply(results_mean, 1, mean, na.rm = T)),]
seer_mean <- results_mean*seer_df$value
seer_mean <- seer_mean2 <- seer_mean[order(-apply(seer_mean, 1, mean, na.rm = T)),]
seer_mean <- reshape2::melt(cbind.data.frame(Tumor = rownames(seer_mean), seer_mean))
results_seer <- cbind.data.frame(Tumor = seer_mean$Tumor, Alteration = seer_mean$variable, 
                            mean = seer_mean$value)
})

palette <- reactive({
palette <- c("#F3D671", "#DD5E46", "#3EADA3")
names(palette) <- c("mut", "cna", "fusion")
return(palette)
})

resultsOrder <- reactive({
  results_mean <- data.frame((as.matrix(tcga_results()[,-1]) + as.matrix(dfci_results()[,-1]) + as.matrix(mskp_results()[,-1]))/3)
  rownames(results_mean) <- tcga_results()[,1]
  results_mean <- data.frame(results_mean)[cancer_select(),]
  colnames(results_mean)<- c("mut", "cna", "fusion")
  results_mean <- results_mean[order(-apply(results_mean, 1, mean, na.rm = T)),]
  return(rownames(results_mean))
})

resultsSeerOrder <- reactive({
  results_mean <- data.frame((as.matrix(tcga_results()[,-1]) + as.matrix(dfci_results()[,-1]) + as.matrix(mskp_results()[,-1]))/3)
  rownames(results_mean) <- tcga_results()[,1]
  results_mean <- data.frame(results_mean)[cancer_select(),]
  colnames(results_mean) <- c("mut", "cna", "fusion")
  results_mean <- results_mean[order(-apply(results_mean, 1, mean, na.rm = T)),]
  seer_df <- cbind.data.frame(tumor = unlist(cancer_select()), 
                              year = unlist(year_select()), 
                              value = as.numeric(unlist(data.frame(seer[which(seer[,2] %in% cancer_select()),which(colnames(seer) %in% year_select())]))))
  seer_df <- seer_df[order(-apply(results_mean, 1, mean, na.rm = T)),]
  seer_mean <- results_mean*seer_df$value
  seer_mean <- seer_mean2 <- seer_mean[order(-apply(seer_mean, 1, mean, na.rm = T)),]
  return(rownames(seer_mean))
})

# PLOTS

output$plot_seer1 <- plotly::renderPlotly({
p =  ggplot2::ggplot(results(), aes(fill=factor(Alteration, levels = rev(c("mut", "cna", "fusion"))), 
                      y=mean, x= factor(Tumor, levels = resultsOrder()),
                      text = paste(
                        "Frequency:", round(mean, 2))))  
p = p + ggplot2::geom_bar(position="stack", stat="identity", color = "black")
p = p + ggplot2::scale_fill_manual(values = palette(),
                               guide = ggplot2::guide_legend(frame.colour = "black",
                                                             override.aes = list(size = 0.5)))
    #  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3)+
p = p + labs(x = "Tumor subtype", y = "Frequency (%)")
p = p + ggpubr::theme_pubr()
p <- p + ggplot2::theme(
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 21, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.y = ggplot2::element_text(size=16),
      axis.text.x = ggplot2::element_text(size=16, angle = 50, vjust = 1, hjust = 1),
      legend.position="right",
      plot.margin = ggplot2::unit(c(0.5,1,0.5,1), "cm")
    )
p = plotly::ggplotly(p, tooltip = "text")
p
  })

output$plot_seer2 <- plotly::renderPlotly({
p = ggplot2::ggplot(results_seer(), aes(fill=factor(Alteration, levels = rev(c("mut", "cna", "fusion"))), 
                           y=mean, x= factor(Tumor, levels = resultsSeerOrder()),
           text = paste("U.S. incidence in", year_select(),":", round(mean, 2)))) 
p = p + ggplot2::geom_bar(position="stack", stat="identity", color = "black")
p = p + ggplot2::scale_fill_manual(values = palette(),
                               guide = ggplot2::guide_legend(frame.colour = "black",
                                                             override.aes = list(size = 0.5)))
    #  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3)+
p = p + ggplot2::labs(x = "Tumor subtype", y = paste("U.S. incidence in", year_select()))
p = p + ggpubr::theme_pubr()
p = p + ggplot2::theme(
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 21, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.y = ggplot2::element_text(size=16),
      axis.text.x =ggplot2::element_text(size=16, angle = 50, vjust = 1, hjust = 1),
      legend.position="right",
      plot.margin = ggplot2::unit(c(0.5,1,0.5,1), "cm")
    )
p = plotly::ggplotly(p, tooltip = "text")
p
})


}


                        