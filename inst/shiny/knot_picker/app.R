library(shiny)
library(shinythemes)
library(bslib)
library(markdown)
library(dplyr)
library(ggplot2)
library(grid)


ui <- fixedPage(

  theme = shinytheme("flatly"),
  br(),
  titlePanel(
    HTML(markdownToHTML(fragment.only = TRUE,
                        text = c(
                          "# **KNOT PICKER for SENNA**")))),
  br(),
  hr(),
  br(),
  fixedRow(
    class = "well",
    column(width = 12,
           HTML(markdownToHTML(fragment.only = TRUE,
                               text = c(
                                 "> **Welcome to** ***knot picker*** **!**")))),
    column(width = 12,
           textInput(inputId = "dotsize", 
                     label = "DOT SIZE FOR CLICKABLE POINTS (default=1):", 
                     value = "1")),
    column(width = 6,
           numericInput(inputId = "pah", 
                        label = "TISSUE MAP ALPHA:", 
                        value = 0.7, min = 0, max = 1, step = 0.1)),
    column(width = 6,
           numericInput(inputId = "kah", 
                        label = "CANVAS ALPHA:", 
                        value = 0.7, min = 0, max = 1, step = 0.1)),
    column(width = 6,
           input_switch(
             id = "image",
             label = "Background Image",
             value = FALSE
           ))),
  br(),
  hr(),
  br(),
  fixedRow(

    column(width = 12,
           column(width = 6,
                  HTML(markdownToHTML(fragment.only = TRUE,
                                      text = c(
                                        "## Tissue Map"
                                      ))),
                  plotOutput("palette", height = 400,
                             brush = brushOpts(
                               id = "segmentation",
                               resetOnNew = FALSE)
                  )
           ),

           column(width = 6,
                  HTML(markdownToHTML(fragment.only = TRUE,
                                      text = c(
                                        "## Canvas for Click"
                                      ))),
                  plotOutput("picker", height = 400,
                             click = clickOpts(
                               id = "pick"
                             )
                  )
           )
    )
  ),

  tags$hr(),


  # Table info

  fixedRow(
    column(width = 12,
           fluidRow(
             column(width = 8,
                    HTML(markdownToHTML(fragment.only = TRUE,
                                        text = c(
                                          "### Selected Points (Marked in black)"
                                        )))),
             column(width = 2,
                    actionButton("erase", label = "Erase")),
             column(width = 2,
                    downloadButton("export", "Export")),
           ),
           br(),
           verbatimTextOutput("knots")),
  ),


  fixedPanel(width = 400, top = 10, right = 10,
             height = "auto",
             draggable = TRUE,
             wellPanel(
               HTML(markdownToHTML(fragment.only = TRUE,
                                   text = c(
                                     "**Current points (Marked in 'X')**

※ This panel is draggable."))),
               verbatimTextOutput("point_info"),
               actionButton("add", label = "Add")),
             style = "opacity : 0.9"
  ),
  br(),
  br(),
  br(),
  br(),
  br(),
  br()
)



server <- function(input, output) {
  
  bg_active <- reactive({
    input$image
  })
  
  ds <- reactive({
    if(bg_active()) min(app_image$idim) * as.numeric(input$dotsize)
    else as.numeric(input$dotsize)
    })
  

  # Masked range

  ranges <- reactiveValues(x1 = NULL, x2 = NULL)
  observe({
    bsh <- input$segmentation
    if(!is.null(bsh)) {
      ranges$x1 <- c(bsh$xmin, bsh$xmax)
      ranges$x2 <- c(bsh$ymin, bsh$ymax)
    } else {
      if(bg_active()){
        ranges$x1 <- c(0, app_image$idim[1])
        ranges$x2 <- c(0, app_image$idim[2])
      }
      else{
        ranges$x1 <- c(0, 1)
        ranges$x2 <- c(0, 1) 
      }
    }
  })



  # Plot rendering
  
  output$palette <- renderPlot({
    pp <- ggplot()
    
    if(bg_active()){
      pp <- pp +
        annotation_custom(app_image$bgr, 
                          xmin = 0, 
                          xmax = min(app_image$idim), 
                          ymin = 0, 
                          ymax = min(app_image$idim)) +
        lims(x = c(0, ymax = min(app_image$idim)), 
             y = c(0, ymax = min(app_image$idim))) +
        coord_fixed(ratio = 1,
                    expand = FALSE)
      
      dsc <- app_image$grid
      kdf <- knotdf2$df
    } else {
      dsc <- app_dat
      kdf <- knotdf$df
    }
    
    if(nrow(clk$df) != 0){
      if(nrow(knotdf$df) != 0){
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor") {
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 2,
                         data = clk$df) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else if(class(dsc$CLS) == "character") {
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_identity() +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 2,
                         data = clk$df) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else{
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_gradientn(colors = colorset) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 2,
                         data = clk$df) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         data = kdf) +
              theme_void() +
              theme(legend.position = "left") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          }
          
        } else{
          pp <- pp +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$pah,
                       size = 0.8,
                       data = dsc) +
            geom_point(aes(X1, X2),
                       shape = 4,
                       data = clk$df) +
            geom_point(aes(X1, X2),
                       col = 'black',
                       data = kdf) +
            theme_void() +
            coord_fixed(ratio = 1,
                        expand = TRUE)
        }
        
      } else{
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor") {
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 2,
                         data = clk$df) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else if(class(dsc$CLS) == "character") {
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_identity() +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 2,
                         data = clk$df) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else{
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_gradientn(colors = colorset) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 2,
                         data = clk$df) +
              theme_void() +
              theme(legend.position = "left") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          }
          
        } else{
          pp <- pp +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$pah,
                       size = 0.8,
                       data = dsc) +
            geom_point(aes(X1, X2),
                       shape = 4,
                       data = clk$df) +
            theme_void() +
            coord_fixed(ratio = 1,
                        expand = TRUE)
        }
      }
    } else {
      if(nrow(kdf) != 0){
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor") {
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else if(class(dsc$CLS) == "character"){
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_identity() +
              geom_point(aes(X1, X2),
                         col = 'black',
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else{
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_gradientn(colors = colorset) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         data = kdf) +
              theme_void() +
              theme(legend.position = "left") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          }
          
        } else{
          pp <- pp +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$pah,
                       size = 0.8,
                       data = dsc) +
            geom_point(aes(X1, X2),
                       col = 'black',
                       data = kdf) +
            theme_void() +
            coord_fixed(ratio = 1,
                        expand = TRUE)
        }
      } else{
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor"){
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else if(class(dsc$CLS) == "character"){
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_identity() +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          } else{
            pp <- pp +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$pah,
                         size = 0.8,
                         data = dsc) +
              scale_color_gradientn(colors = colorset) +
              theme_void() +
              theme(legend.position = "left") +
              coord_fixed(ratio = 1,
                          expand = TRUE)
          }
        } else{
          pp <- pp +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$pah,
                       size = 0.8,
                       data = dsc) +
            theme_void() +
            coord_fixed(ratio = 1,
                        expand = TRUE)
        }
      }
    }
    
    return(pp)
  })
  
  
  output$picker <- renderPlot({
    pk <- ggplot()
    
    if(bg_active()){
      pk <- pk +
        annotation_custom(app_image$bgr, 
                          xmin = 0, 
                          xmax = min(app_image$idim), 
                          ymin = 0, 
                          ymax = min(app_image$idim)) +
        lims(x = c(0, ymax = min(app_image$idim)), 
             y = c(0, ymax = min(app_image$idim))) +
        coord_fixed(ratio = 1,
                    expand = FALSE)
      
      dsc <- app_image$grid
      kdf <- knotdf2$df
    } else {
      dsc <- app_dat
      kdf <- knotdf$df
    }
    
    if(nrow(clk$df) != 0){
      if(nrow(knotdf$df) != 0){
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                                 max(ranges$x2) - min(ranges$x2)),
                         data = clk$df) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else if(class(dsc$CLS) == "character"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              scale_color_identity() +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                                 max(ranges$x2) - min(ranges$x2)),
                         data = clk$df) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else{
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              
              scale_color_gradientn(colors = colorset) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                                 max(ranges$x2) - min(ranges$x2)),
                         data = clk$df) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          }
        } else{
          pk <- pk +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$kah,
                       size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                         max(ranges$x2) - min(ranges$x2)),
                       data = dsc) +
            geom_point(aes(X1, X2),
                       shape = 4,
                       size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                               max(ranges$x2) - min(ranges$x2)),
                       data = clk$df) +
            geom_point(aes(X1, X2),
                       col = 'black',
                       size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                         max(ranges$x2) - min(ranges$x2)),
                       data = kdf) +
            theme_void() +
            coord_fixed(ratio = 1,
                        xlim = ranges$x1,
                        ylim = ranges$x2,
                        expand = FALSE)
        }
      } else {
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor") {
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                                 max(ranges$x2) - min(ranges$x2)),
                         data = clk$df) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else if(class(dsc$CLS) == "character"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              scale_color_identity() +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                                 max(ranges$x2) - min(ranges$x2)),
                         data = clk$df) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else{
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              scale_color_gradientn(colors = colorset) +
              geom_point(aes(X1, X2),
                         shape = 4,
                         size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                                 max(ranges$x2) - min(ranges$x2)),
                         data = clk$df) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          }
        } else{
          pk <- pk +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$kah,
                       size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                         max(ranges$x2) - min(ranges$x2)),
                       data = dsc) +
            geom_point(aes(X1, X2),
                       shape = 4,
                       size = 1.5 * ds() / max(max(ranges$x1) - min(ranges$x1),
                                               max(ranges$x2) - min(ranges$x2)),
                       data = clk$df) +
            theme_void() +
            coord_fixed(ratio = 1,
                        xlim = ranges$x1,
                        ylim = ranges$x2,
                        expand = FALSE)
        }
      }
    } else {
      if(nrow(kdf) != 0){
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else if(class(dsc$CLS) == "character"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              scale_color_identity() +
              geom_point(aes(X1, X2),
                         col = 'black',
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else{
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc) +
              scale_color_gradientn(colors = colorset) +
              geom_point(aes(X1, X2),
                         col = 'black',
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = kdf) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          }
        } else{
          pk <- pk +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$kah,
                       size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                         max(ranges$x2) - min(ranges$x2)),
                       data = dsc) +
            geom_point(aes(X1, X2),
                       col = 'black',
                       size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                         max(ranges$x2) - min(ranges$x2)),
                       data = kdf) +
            theme_void() +
            coord_fixed(ratio = 1,
                        xlim = ranges$x1,
                        ylim = ranges$x2,
                        expand = FALSE)
        }
      } else{
        if("CLS" %in% colnames(dsc)){
          if(class(dsc$CLS) == "factor"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc)  +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else if(class(dsc$CLS) == "character"){
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc)  +
              scale_color_identity() +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          } else{
            pk <- pk +
              geom_point(aes(X1, X2, colour = CLS),
                         alpha = input$kah,
                         size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                           max(ranges$x2) - min(ranges$x2)),
                         data = dsc)  +
              scale_color_gradientn(colors = colorset) +
              theme_void() +
              theme(legend.position = "none") +
              coord_fixed(ratio = 1,
                          xlim = ranges$x1,
                          ylim = ranges$x2,
                          expand = FALSE)
          }
        } else{
          pk <- pk +
            geom_point(aes(X1, X2),
                       col = "grey",
                       alpha = input$kah,
                       size = ds() / max(max(ranges$x1) - min(ranges$x1),
                                         max(ranges$x2) - min(ranges$x2)),
                       data = dsc) +
            theme_void() +
            coord_fixed(ratio = 1,
                        xlim = ranges$x1,
                        ylim = ranges$x2,
                        expand = FALSE)
        }
      }
    }
    
    return(pk)
  })




  # Click logic

  clk <- reactiveValues(
    df = data.frame(
      X1 = NULL,
      X2 = NULL
    )
  )

  observeEvent(input$pick, {
    req(input$pick)
    if (bg_active()) {
      dsc <- app_image$grid
    } else {
      dsc <- app_dat
    }
    
    clk$df <- nearPoints(df = dsc,
                         coordinfo = input$pick,
                         xvar = "X1", 
                         yvar = "X2", 
                         addDist = FALSE)
  })


  output$point_info <- renderPrint({
    as.data.frame(clk$df)
  })





  # knots data logic

  knotdf <- reactiveValues(
    df = data.frame(
      X1 = NULL,
      X2 = NULL
    )
  )
  
  knotdf2 <- reactiveValues(
    df = data.frame(
      X1 = NULL,
      X2 = NULL
    )
  )

  # Add
  observeEvent(input$add, {
    if(bg_active()){
      knotdf2$df <- BiocGenerics::rbind(knotdf2$df,
                                        clk$df)
      knotdf$df <- app_dat[rownames(knotdf2$df),,drop = FALSE]
    }
    
    else{
      knotdf$df <- BiocGenerics::rbind(knotdf$df,
                                       clk$df) 
    }
  })

  output$knots <- renderPrint(knotdf$df)


  # Erase

  observeEvent(input$erase, {
    knotdf$df <- knotdf$df[-nrow(knotdf$df),]
    
    if(bg_active()) {
      knotdf2$df <- knotdf2$df[-nrow(knotdf2$df),]
    }
  })

  output$knots <- renderPrint({
    as.data.frame(knotdf$df)
  })


  # Export csv

  output$export <- downloadHandler(
    filename = function(){
      default <- "knots.csv"
    },
    content = function(file) {
      write.csv(knotdf$df, file, row.names = FALSE)
    }
  )
}





# Run the application
shinyApp(ui = ui, server = server)
