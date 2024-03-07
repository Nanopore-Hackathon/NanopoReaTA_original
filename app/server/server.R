#### server.R

# ______________________________________________________________________________
# DESCRIPTION OF BUTTONS

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
names_colorblind_palette <- c("light blue", "light red", "mustard", "green", "violett", "magenta", 
                              "turquoise", "dirty yellow", "dark pink", "merlot", "sky", "light grey")

# ______________________________________________________________________________
# SERVER ---------------------------------------------------------
server <- function(input, output, session) {
  ## 0. Server settings ####
  cores = 16
  myProcess <- NULL
  
  
  ### 0.2 Cores for DEA ####
  # Set number of cores of default cores used for R process (from input in DEA tab)
  observe({register(MulticoreParam(max(c(4,round(cores * 0.6)))))})
  
  
  ## RAM in use (currently not embedded into ui)
  output$memuse.out <- renderText({
    invalidateLater(5000, session)
    as.character(Sys.meminfo()[[2]])
  })
  
  ## INPUT FILE PATH
  metadata_path = "metadata.tsv"
  bed_file = "hg38_GENCODE_V42_Basic.bed"
  gtf_file = "converted_gtf.csv"
  
  ## READ INPUT FILES 
  metadata <- data.table::fread(metadata_path, data.table = F)
  row.names(metadata) = metadata$Samples
  gtf_df <- data.table::fread(gtf_file, data.table = T, sep = ",", header = T)
  gtf_df$V1 <- NULL
  gtf_df <- unique(gtf_df)
  gtf_df <- as.data.frame(gtf_df)
  
  ## 1. Preprocessing #######
  
  # Not really needed
  alignment = reactiveVal()
  alignment(0)
  
  ### Buttons Actions ####
  
  observeEvent(input$run_dea_B, {
    updateTabsetPanel(session, "dea.box",
                      selected = "dea.tab")
  })

  observeEvent(input$run_dte_B, {
    updateTabsetPanel(session, "dea.box",
                      selected = "dte.tab")
  })

  observeEvent(input$run_dt_preprocess, {
    updateTabsetPanel(session, "dea.box",
                      selected = "deu.tab")
  })
  

  
  output$metadata.tables.out <- renderUI({
    metadata_table_box(input, output, session)
  })
  
  output$design_column_in <- renderUI({
    if (is.null(metadata)) return()
    selectInput(inputId = 'design_column', 
                label = 'Select condition column', 
                choices = colnames(metadata), selected = colnames(metadata)[2])
  })
  
  output$feature_column_A <- renderUI({
    if (is.null(metadata)) return()
    x = unique(metadata[,input$design_column])
    print(x[1])
    selectInput(inputId = 'feature_A', 
                label = 'Select condition A', selected = x[1],
                choices = unique(metadata[,input$design_column]))
  })
  
  output$feature_column_B <- renderUI({
    if (is.null(metadata)) return()
    x = unique(metadata[,input$design_column])
    print(x[2])
    selectInput(inputId = 'feature_B', 
                label = 'Select condition B', 
                selected = x[2],
                choices = x[!(x %in% input$feature_A)])
  })
  

  output$color_feature_A <- renderUI({
    if (is.null(metadata)) return()
    choices_df = data.frame(
      names = names_colorblind_palette, 
      hex = safe_colorblind_palette
    )
    selectInput("color_A", 
                label = "Select color A", 
                choices = c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
                            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"),
                selected = "#88CCEE")
  })
  
  output$color_feature_B <- renderUI({
    if (is.null(metadata)) return()
    choices_df = data.frame(
      names = names_colorblind_palette, 
      hex = safe_colorblind_palette
    )
    selectInput(inputId = "color_B", 
            label = "Select color B", 
            choices = c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
                            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"),
            selected = "#DDCC77")
  })

  condi_cols = reactive({
    req(metadata)
    condi_cols = c(input$color_A, input$color_B)
    names(condi_cols) = c(input$feature_A, input$feature_B)
    condi_cols
  })
  
  output$design.matrix.out <- renderUI({
    design_matrix_box(input, output, session, F)
  })

  
  output$metadata.out <- renderDT({
    if (is.null(metadata)){
      h5("No metadata")
    } else {
      metadata #< reactive metadata value
    }
  })
  

  ## 1.1 Output from Preprocessing ####
  
  ### Load counts file ####
  dds_file = reactiveValues(df = NULL, df_trans = NULL)
  rld_file = reactiveValues(df = NULL, df_trans = NULL)
  

  observeEvent({input$run_dea_B + input$submit_gene_selection}, {
    print(paste0("Current tab: ", input$tabs))

    print("###################################")
    print("##     Load/Update Counts file   ##")
    print("###################################")
    print(paste0(Sys.time(), ": Starting count file loading/updating"))

    counts_file = "gene.RData"

    print(sprintf("Count file does exists?: %s", file.exists(counts_file)))
    if (!file.exists(counts_file)) {
      showModal(modalDialog(size = "l",
                            title = "No gene count file",
                            "No gene count file found! Wait until the preprocessing pipeline has generated a gene count file! This can take a few minutes",
                            easyClose = TRUE,
                            footer = NULL
      ))
    } else {
      load(counts_file, verbose = T) # loads dds and rld object
      dds_file$df = dds
      rld_file$df = rld
      
    } 
  }, 
  priority = 2
  )
  
  
  observeEvent({input$run_dt_preprocess + input$run_dte_B}, {
    req(metadata)
    print(paste0("Current tab: ", input$tabs))
    print("###################################")
    print("##     Load/Update Transcript Counts file   ##")
    print("###################################")
    print(paste0(Sys.time(), ": Starting transcript count file loading/updating"))
    
    counts_file = "transcript.RData"
    
    print(sprintf("Salmon Count file does exists?: %s", file.exists(counts_file)))
    if (file.exists(counts_file)) {
      load(counts_file, verbose = T)
      dds_file$df_trans = dds
      rld_file$df_trans = rld
      
      print("Transcript count file loaded successfully!")

    }
  }, 
  priority = 2
  )
  
  
  ### LOAD DTU RDATA SCRIPT ####
  
  dtu_file = reactiveValues(d = NULL, dxd = NULL)
  
  observeEvent({input$run_dt_preprocess + input$run_dte_B}, {
    req(metadata)
    print(paste0("Current tab: ", input$tabs))
    print("###################################")
    print("##     Load/Update Transcript Usage file   ##")
    print("###################################")
    print(paste0(Sys.time(), ": Starting transcript usage file loading/updating"))
    
    DTU_file = "transcript_usage.RData"
    
    print(sprintf("Trancsript Usage file does exists?: %s", file.exists(DTU_file)))
    if (file.exists(DTU_file)) {
      load(DTU_file, verbose = T)
      dtu_file$d = d_object
      dtu_file$dxd = dxd_object
      
      print("Transcript usage file loaded successfully!")
      
    }
  }, 
  priority = 2
  )


  ## 3. Run Overview ####
  
  ### Read length distribution ####
  # store inputs from reactiveFileReader later
  overview_variables <- reactiveValues(
    readLengthTable = NULL
  )

  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$tabs
               },
               handlerExpr = {
                 if (input$tabs == "run_overview"){
                 overview_variables$readLengthTable <- reactivePoll(10000, session,
                                                                       checkFunc = function(){
                                                                        inputSizeDir = "ReadLengthFolder"
                                                                        if (dir.exists(inputSizeDir)){
                                                                          filesList <- list.files(inputSizeDir, full.names = TRUE)
                                                                          file.info(filesList[1])$mtime[1]
                                                                        } else {
                                                                          print("Directory does not exists!")
                                                                        }
                                                                       },
                                                                       valueFunc = function(){
                                                                        inputSizeDir = "ReadLengthFolder"
                                                                        if (dir.exists(inputSizeDir)){
                                                                          filesList <- list.files(inputSizeDir, full.names = TRUE)
                                                                          if (length(filesList) == 0) return()

                                                                          names(filesList) = gsub("_read_lengths_pass.txt", "", basename(filesList))

                                                                            # Load files and add pass/fail column
                                                                          readLengths <- lapply(filesList, read.table, header = T)
                                                                          readLengths
                                                                        } else {
                                                                          print("Directory does not exists!")
                                                                          return()
                                                                        }
                                                                      }
                                                                    )
                 }
               }
  )

  readLengthTab <- reactive({
    
    req(overview_variables$readLengthTable())
    print("### Overview Length reactive")
    print(length(overview_variables$readLengthTable()))
    if (is.null(overview_variables$readLengthTable())) {
      return()
    } else {
      readLengths = overview_variables$readLengthTable()
      # Filter out lists of no reads
      readLengths = readLengths[unlist(lapply(readLengths, function(x) nrow(x) > 0))]
      # Create data.frame of list
      readLengths_df = data.table::rbindlist(readLengths, idcol = "Sample")
      upper_bound <- quantile(readLengths_df$Length, 0.99)
      readLengths_df_filt = readLengths_df[readLengths_df$Length < upper_bound,]
      readLengths_df_filt$Sample = gsub("_.*", "", readLengths_df_filt$Sample)
      readLengths_df_filt

    }
  })

  output$samplewise_read_length.out <- renderPlot({
    req(readLengthTab())
    print("Samplewise read length plot")
    print(length(readLengthTab()))
    if(is.null(readLengthTab())) return()
    g = createLengthPlots(readLengthTab(), metadata, input$design_column, c(input$feature_A, input$feature_B), c(input$color_A, input$color_B))
    g[[1]]
  }, bg = "transparent")

  output$groupwise_read_length.out <- renderPlot({
    req(readLengthTab())
    if(is.null(readLengthTab())) return()
    g = createLengthPlots(readLengthTab(), metadata, input$design_column, c(input$feature_A, input$feature_B), c(input$color_A, input$color_B))
    g[[2]]
  }, bg = "transparent")

  # downloads

  output$samplewise_read_length.down <- downloadHandler(
    filename = function(){
      paste0("Read_length_distribution_samplewise_", Sys.Date(), ".", input$samplewise_read_length.type)
    },
    content = function(file){
      if(is.null(readLengthTab())) return()
      plot = samplewise_read_length.download(readLengthTab(), metadata, input$design_column, c(input$feature_A, input$feature_B))
      ggsave(file, plot = plot, device = input$inner_var_sample.type, height = 5, width = 14, units = "in")
    })
  output$groupwise_read_length.down <- downloadHandler(
    filename = function(){
      paste0("Read_length_distribution_groupwise_", Sys.Date(), ".", input$groupwise_read_length.type)
    },
    content = function(file){
      if(is.null(readLengthTab())) return()
      plot = groupwise_read_length.download(readLengthTab(), metadata, input$design_column, c(input$feature_A, input$feature_B), c(input$color_A, input$color_B))
      ggsave(file, plot = plot, device = input$groupwise_read_length.type, height = 5, width = 14, units = "in")
    })

  output$process_time.down <- downloadHandler(
    filename = function(){
      paste0("Process_time_tools_over_time_", Sys.Date(), ".", input$process_time.type)
    },
    content = function(file){
      req(process_time$df())
      if (!is.null(process_time$df())){
        print(head(process_time$df()))
        x = process_time$df()
        # replace _ with whitespace in tool names for a cleaner legend
        x = x %>% mutate(Tool = gsub("_", " ", Tool))
        # get the x-tick interval based on the number of iterations in the data
        max_it <- max(x$Iteration)
        xticks <- 1:max_it
        if (max_it <= 10) { # every tick
          xticks <- xticks
        } else if (max_it <= 50) { # every 5th tick
          xticks <- xticks[xticks %% 5 == 0]
        } else if (max_it <= 100) { # every 10th tick
          xticks <- xticks[xticks %% 10 == 0]
        } else if (max_it <= 1000) { # every 50th tick
          xticks <- xticks[xticks %% 50 == 0]
        } else if (max_it <= 10000) { # every 100th tick
          xticks <- xticks[xticks %% 100 == 0]
        } else { # every 1000th tick
          xticks <- xticks[xticks %% 1000 == 0]
        }
        plot = ggplot(x, aes(x = Iteration, y = Time, fill = Tool)) +
          geom_bar(stat = "identity") +
          scale_x_continuous(breaks = xticks) + # label each x tick individually + # label each x tick individually
          ylab("Time [s]") + # add [s] to y-axis label
          theme_bw() +
          scale_fill_manual("Preprocessing steps", values = c("#CC6677", "#d8c66c", "#117733", "#88CCEE", "#AA4499", "#24011e","#092a36")) +
          theme(
          panel.grid.minor = element_blank(), # remove minor grid lines
          panel.grid.major.x = element_blank(),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(angle = 45, hjust = 1, size = 17),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
          axis.title = element_text(size = 23))
        ggsave(file, plot = plot, device = input$process_time.type, height = 10, width = 10, units = "in")
      }
  })
  
  ### Gene Expression Variability ####
  
  # store inputs from reactiveFileReader later
  overview_variables_infer <- reactiveValues(
    inner_var_sample.table = NULL, 
    expGenes_counted.table = NULL
  )
  
  observeEvent(input$tabs,{
    if (input$tabs == "run_overview"){
      
      overview_variables_infer$inner_var_sample.table <- reactivePoll(10000, session,
                                     checkFunc = function(){
                                       inner_var="inner_variability_per_sample.csv"
                                      if (file.exists(inner_var)){   
                                        file.info(inner_var)$mtime[1]
                                       } else {
                                         print("Inner variability file does not exists!")
                                       }
                                     },
                                     valueFunc = function(){
                                       inner_var="inner_variability_per_sample.csv"
                                       
                                       if (file.exists(inner_var)) {
                                        read.table(inner_var, sep="\t", header=T)
                                       } else {
                                        NULL
                                       }
                                     })
      overview_variables_infer$expGenes_counted.table <- reactivePoll(10000, session,
                                                                checkFunc = function(){
                                                                  genes_counted="exp_genes_counted_per_sample.csv"
                                                                  if (file.exists(genes_counted)){   
                                                                        file.info(genes_counted)$mtime[1]
                                                                      
                                                                  } else {
                                                                    print("Counts file does not exists!")
                                                                  }
                                                                },
                                                                valueFunc = function(){
                                                                  
                                                                  genes_counted="exp_genes_counted_per_sample.csv"
                                                                  
                                                                  if (file.exists(genes_counted)) {
                                                                    read.table(genes_counted, sep="\t", header=T)
                                                                  } else {
                                                                    NULL
                                                                  }
                                                                })
    }
  },
  priority = -1
  )
  

  output$inner_var_sample.out <- renderPlot({
    if (is.null(overview_variables_infer$inner_var_sample.table())){ 
      inner_var_plot_per_sample(data.frame())
    }
    else{
    inner_var_plot_per_sample(overview_variables_infer$inner_var_sample.table())
    }
  }, bg = "transparent")
      
  output$inner_var_condition.out <- renderPlot({
    req(metadata)
    if (is.null(overview_variables_infer$inner_var_sample.table())){
        inner_var_plot_per_condition(data.frame(), metadata, c(input$color_A, input$color_B))
    } 
    else{
    inner_var_plot_per_condition(overview_variables_infer$inner_var_sample.table(), metadata, c(input$color_A, input$color_B))
    }
  }, bg = "transparent")
  
  output$total_genes_sample.out <- renderPlot({
    req(overview_variables_infer$expGenes_counted.table())

    if (is.null(overview_variables_infer$expGenes_counted.table())){
        total_genes_counted_plot_per_sample(data.frame())
      }
    else{
      total_genes_counted_plot_per_sample(overview_variables_infer$expGenes_counted.table())
    }
    }, bg = "transparent")
      
  output$total_genes_condition.out <- renderPlot({
    req(metadata)
    if(is.null(overview_variables_infer$expGenes_counted.table())){
      total_genes_counted_plot_per_condition(data.frame(), metadata, c(input$color_A, input$color_B))
    }
    else{
      total_genes_counted_plot_per_condition(overview_variables_infer$expGenes_counted.table(), metadata, c(input$color_A, input$color_B))
    }
  
  }, bg = "transparent")
    
  # downloads 
  

  output$inner_var_sample.down <- downloadHandler(
   filename = function(){
    paste0("Expression_variability_over_time_samplewise_", Sys.Date(), ".", input$inner_var_sample.type)
    }, 
   content = function(file){
    if(is.null(overview_variables_infer$inner_var_sample.table())) return()
    plot = inner_var_plot_per_sample.download(overview_variables_infer$inner_var_sample.table())
    ggsave(file, plot = plot, device = input$inner_var_sample.type, height = 5, width = 14, units = "in")
   })



  output$inner_var_condition.down <- downloadHandler(
    filename = function(){
      paste0("Expression_variability_over_time_conditionwise_", Sys.Date(), ".", input$inner_var_condition.type)
    }, 
    content = function(file){
      if(is.null(overview_variables_infer$inner_var_sample.table()) | is.null(metadata)) return()
      plot = inner_var_plot_per_condition.download(overview_variables_infer$inner_var_sample.table(), metadata, c(input$color_A, input$color_B))
      ggsave(file, plot = plot, device = input$inner_var_condition.type, height = 5, width = 14, units = "in")
    })

  

  output$total_genes_sample.down <- downloadHandler(
    filename = function(){
      paste0("Gene_expression_counts_over_time_samplewise_", Sys.Date(), ".", input$total_genes_sample.type)
    }, 
    content = function(file){
      if(is.null(overview_variables_infer$expGenes_counted.table())) return()
      plot = total_genes_counted_plot_per_sample.download(overview_variables_infer$expGenes_counted.table())
      ggsave(file, plot = plot, device = input$total_genes_sample.type, height = 5, width = 14, units = "in")
    })


 
  output$total_genes_condition.down <- downloadHandler(
    filename = function(){
      paste0("Gene_expression_counts_over_time_conditionwise_", Sys.Date(), ".", input$total_genes_condition.type)
    }, 
    content = function(file){
      if(is.null(overview_variables_infer$expGenes_counted.table()) | is.null(metadata)) return()
      plot = total_genes_counted_plot_per_condition.download(overview_variables_infer$expGenes_counted.table(), metadata, c(input$color_A, input$color_B))
      ggsave(file, plot = plot, device = input$total_genes_condition.type, height = 5, width = 14, units = "in")
    })

  
  
  
  
  ## 2. Quality Control ####
  ### Buttons ####
  table_of_normCounts <- reactiveValues(df_norm = NULL)
  
  output$submit_gene_selection.out <- renderUI({
    if (!is.null(gtf_df)){
         actionButton("submit_gene_selection", "Submit genes",
                      class = "btn btn-primary",
                      style='font-size:200%; color: white;
                             background-color: #e76f51; 
                             border-radius: 5px')
      } else {
        actionButton("submit_gene_selection_disabled", "Submit genes",
                     class = "btn btn-primary",
                     style='font-size:200%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
    }
  })
  
  observeEvent({input$run_dea_B + input$submit_gene_selection}, {
    print(paste0("Norm current tab: ", input$tabs))
    req(dds_file$df)
    print(!is.null(dds_file$df))
    if (is.null(dds_file$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    
    print(!is.null(metadata))
    if (!is.null(dds_file$df) & !is.null(metadata)){
      if (is.null(input$design_column)){
        print("input$design_column is null")
        return()
      }
      print(input$design_column)
      table_of_normCounts$df_norm <- counts(dds_file$df, normalize = T)
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      table_of_normCounts$df_norm <- NULL
      updateTabItems(session, "tabs", "run_overview")
      updateTabsetPanel(session, "main_quality_control",
                      selected = "read_length_dist_panel")
    }
  }, 
  priority = 3)
  
   
  output$counts_plots_box <- renderUI({
    input$submit_gene_selection
    if (is.null(input$submit_gene_selection)) return()
    if (input$submit_gene_selection <= 0){
      box(title = "Counts plots", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = T, 
          width = 12,
          fluidRow(
            column(12, radioButtons("disp", "Display",
                                    choices = c("Boxplot" = "boxplot",
                                                "Dotplot" = "dotplot",
                                                "Violin plot" = "violinplot"),
                                    selected = "boxplot", inline = T)),
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_tea')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("tea_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("teaPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1")),
            column(12, br(), p(em(fig_des$gene_counts)))
          ))
    } else {
      box(title = "Counts plots", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = F, 
          width = 12,
          fluidRow(
            column(12, radioButtons("disp", "Display",
                                    choices = c("Boxplot" = "boxplot",
                                                "Dotplot" = "dotplot",
                                                "Violin plot" = "violinplot"),
                                    selected = "boxplot", inline = T)),
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_tea')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("tea_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("teaPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1")),
            column(12, br(), p(em(fig_des$gene_counts)))
          ))
    }
  })
  
  
  genes.list <- reactive({
    gtf_genes = unique(gtf_df[, c("gene_id", "gene_name")])
    
    print(input$table_of_genes_df_rows_selected)
    gtf_genes[input$table_of_genes_df_rows_selected, ]
  })
  tea_res = reactiveValues(df_res = NULL)
  observeEvent({input$submit_gene_selection}, {
    print(paste0("TEA current tab: ", input$tabs))

    cat(input$submit_gene_selection)
    if (is.null(input$submit_gene_selection)) return()
    if (!(input$tabs == "qualitycontrol")) return()
    req(dds_file$df)
    if (is.null(dds_file$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (is.null(genes.list())){
      print("genes list is null")
      return()
    }
    req(genes.list())
    
    print(head(metadata))
    print(paste0(Sys.time(), ": Starting TEA"))
    print(genes.list()[,1])
    print(condi_cols())
    print(head(counts(dds_file$df)))
    print(metadata)
    print(condi_cols())
    tea_res$df_res <- TEA(counts = counts(dds_file$df),
        norm_counts = counts(dds_file$df, normalize = T), 
        genes.list = genes.list(), 
        metadata = metadata,
        pvalue = 0.05, 
        output.dir = "TEA", condi_cols())
  }, priority = 2)
  
  output$table_of_genes_df <- DT::renderDataTable({
    config = gtf_df
    
    config_out = config
    config_out = unique(config_out)
    DT::datatable(
      unique(config_out[, c("gene_id", "gene_name")]),
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
      rownames = T)
  }, server = FALSE)
  
  output$selectedGenes <- DT::renderDataTable({
    req(genes.list())
    config = genes.list()
    config_out = config
    DT::datatable(
      config_out, 
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 't',
                     fixedColumns = TRUE),
      rownames = T)
    
  }, server = FALSE)
  
    output$teaPlot <- renderPlot({
    req(tea_res$df_res)
    if (input$disp == "boxplot"){
      return(tea_res$df_res[["Boxplot"]])
    }
    if (input$disp == "dotplot"){
      return(tea_res$df_res[["Dotplot"]])
    }
    if (input$disp == "violinplot"){
      return(tea_res$df_res[["Violinplot"]])
    }
  }, bg="transparent")
  
  # Save button for TEA plot
  observe({
    req(tea_res$df_res)
    if (input$disp == "boxplot"){
      output$down_tea = downloadHandler(
        filename = function() {
          paste0("Gene_counts_boxplot_", Sys.Date(), ".", input$tea_down.type)
        },
        content = function(file) {
          plot <- tea_res$df_res[["Boxplot.down"]]
          ggsave(file, plot, device = input$tea_down.type, bg="white", height = 8, width = 20, units = "in")
        }
      )
    }
    if (input$disp == "dotplot"){
      output$down_tea = downloadHandler(
        filename = function() {
          paste0("Gene_counts_dotplot_", Sys.Date(), ".", input$tea_down.type)
          
        },
        content = function(file) {
          plot <- tea_res$df_res[["Dotplot.down"]]
          ggsave(file, plot, device = input$tea_down.type, bg="white", height = 8, width = 20, units = "in")
        }
      )
    }
    if (input$disp == "violinplot"){
      output$down_tea = downloadHandler(
        filename = function() {
          
          paste0("Gene_counts_violinplot_", Sys.Date(), ".", input$tea_down.type)
          
        },
        content = function(file) {
          plot <- tea_res$df_res[["Violinplot.down"]]
          ggsave(file, plot, device = input$tea_down.type, bg="white", height = 8, width = 20, units = "in")
        }
      )
    }
    
  })
  ## Gene Body Coverage ####

  output$table_of_genes_df_gC <- DT::renderDataTable({
    config = gtf_df
    
    config_out = config
    config_out = unique(config_out)
    DT::datatable(
      unique(config_out[, c("gene_id", "gene_name")]),
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 200,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
                     selection = list(mode = "single", target = "row"),
      rownames = T)
  },server = FALSE)

   genes.list_gC <- reactive({
    gtf_genes = unique(gtf_df)
    print(input$table_of_genes_df_gC_rows_selected)
    gtf_genes[input$table_of_genes_df_gC_rows_selected, ]
  })

  output$selectedGene_gC <- DT::renderDataTable({
    req(genes.list_gC())
    config = genes.list_gC()
    config_out = config
    DT::datatable(
      config_out, 
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 200,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
      selection = list(mode = "single"
                      , target = "row"),
      rownames = T)
    
  }, server = FALSE)

  output$submit_gene_selection_gC.out <- renderUI({
    if (!is.null(gtf_df)){
         actionButton("submit_gene_selection_gC", "Submit gene",
                      class = "btn btn-primary",
                      style='font-size:200%; color: white;
                             background-color: #e76f51; 
                             border-radius: 5px')
      } else {
        actionButton("submit_gene_selection_gC_disabled", "Submit gene",
                     class = "btn btn-primary",
                     style='font-size:200%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
    }
  })
  
  ### Start genebodycoverage pipeline ####
  gB_results = reactiveVal(NULL)
  observeEvent(input$submit_gene_selection_gC, {
    req(input$submit_gene_selection_gC)
    if (!file.exists("g_percentiles.json")) return()
    if (!is.null(input$submit_gene_selection_gC)){

      # Save config file
      currentwd = dirname(getwd())
      print(currentwd)
      print(genes.list_gC())
      bamFileList = "bam_genome_merged/"
      print(bamFileList)
      converted_gtf = "converted_gtf.csv"
      print(converted_gtf)
      myProcess <<- processx::process$new(command = "server/bash_scripts/run_genebodycoverage.sh", 
                                          args = c("server/python_scripts/get_geneBody_coverage.py", bamFileList, genes.list_gC()$gene_id, converted_gtf, "./"), echo_cmd = T,stderr = "error_GC.log", stdout = "index_GC.log")
      myProcess$wait()
      df = "samples.geneBodyCoverage.txt"
      if (!file.exists(df)) print("something went wrong")
      geneBodyCov = read.table(df, header = T)
      print(geneBodyCov)
      geneBodyCov.plots = geneBodyCov.plot(geneBodyCov, genes.list_gC(), metadata, condi_cols())
      gB_results(geneBodyCov.plots)
      print("GB is done!")
    }
  }, priority = 2)
  
  output$geneBodyCoveragePlot <- renderUI({
    input$submit_gene_selection_gC
    if (is.null(input$submit_gene_selection_gC)) return()
    if (input$submit_gene_selection_gC <= 0){
      box(title = "Gene Body Coverage", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = T, 
          width = 12,
          fluidRow(
            column(12, radioButtons("disp", "Display",
                                    choices = c("Barplot" = "barplot",
                                                "Dotplot" = "dotplot",
                                                "Line plot" = "lineplot"),
                                    selected = "lineplot", inline = T)),
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_geneBod')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("geneBod_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("geneBodyPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1"))
          ))
    } else {
      box(title = "Gene Body Coverage", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = F, 
          width = 12,
          fluidRow(
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_geneBod')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("geneBod_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("geneBodyPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1"))
          ))
    }
  })

  output$geneBodyPlot <- renderPlot({   
    req(gB_results())
    ggarrange(gB_results()[[1]], gB_results()[[2]], ncol = 2)
  }, bg="transparent")
  
  # Save button for Coverage plot
  output$down_geneBod = downloadHandler(
    filename = function() {
        paste0("GeneBodyCoverage_",genes.list_gC()$gene_name, "_", Sys.Date(), ".", input$geneBod_down.type)
  
    },
    content = function(file) {
      g1 = gB_results()[[1]] + theme_bw() +
          theme(
        legend.title = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"), # remove rotation from x tick labels
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"))

      g2 = gB_results()[[2]] + theme_bw() +
          theme(
        legend.title = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"), # remove rotation from x tick labels
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"))

      g_out = ggarrange(g1, g2, ncol = 2)
     ggsave(file, plot = g_out, device = input$geneBod_down.type, bg="white", height = 10, width = 25, units = "in")
    }
  )


  # Infer experiment #####
  infer_plots <- reactive({
    req(dds_file$df)
    if (is.null(dds_file$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (!is.null(dds_file$df) & !is.null(metadata)){
      counts_file = counts(dds_file$df)
      if (is.null(input$design_column) | is.null(input$design_column) | is.null(input$design_column)){
        return()
      }
      p1 = total_genes_counted_plot_per_condition(counts_file, metadata)
      p2 = total_genes_counted_plot_per_sample(counts_file)
      list(p1, p2)
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      return()
      
    }
  })
  
  
  output$gene_counts_develop <- renderPlot({
    if (is.null(dds_file$df) | is.null(metadata)) return()
    req(dds_file$df)
    req(metadata)
    infer_plots()[[1]]
    
  }, bg="transparent")
  output$gene_counts_develop2 <- renderPlot({
    if (is.null(dds_file$df)) return()
    req(dds_file$df)
    infer_plots()[[2]]
  }, bg="transparent")
  
  
  
  ## 3. DEA ####
  ### Buttons & Actions ####
  #Confirm button for DEA settings
  output$submit_dge <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dea_B" , label = "Gene Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:120%;
                                                        background-color: #e76f51;")
    } else {
      actionButton(inputId = "run_dea_B_disabled" , label = "Gene Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:120%;
                                                        background-color: gray;")
    }
  })


  output$submit_dte <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dte_B" , label = "Transcript Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:120%;
                                                        background-color: #e76f51;")
    } else {
      actionButton(inputId = "run_dte_B_disabled" , label = "Transcript Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:120%;
                                                        background-color: gray;")
    }
  })

  
  output$submit_dt_preprocess <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dt_preprocess" , label = "Transcript Usage Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:120%;
                                                        background-color: #e76f51;")
    } else {
      actionButton(inputId = "run_dt_preprocess_disabled" , label = "Transcript Usage Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:120%;
                                                        background-color: gray;")
      
    }
  })

  
  output$submit_dtu <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dtu_B" , label = "Submit gene selection", class = "btn btn-primary",
                      style='font-size:200%; color: white;
                             background-color: #e76f51; 
                             border-radius: 5px')
    } else {
      actionButton(inputId = "run_dtu_B_disabled" , label = "Submit gene selection", class = "btn btn-primary",
                     style='font-size:200%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
    }
  })
  
  observeEvent(input$run_dea_B, {
    updateTabsetPanel(session, "dea.box",
                      selected = "dea.tab")
  }, priority = 10)

  observeEvent(input$run_dte_B, {
    updateTabsetPanel(session, "dea.box",
                      selected = "dte.tab")
  }, priority = 10)
  
  
  observeEvent(input$run_dt_preprocess, {
    updateTabsetPanel(session, "dea.box",
                      selected = "deu.tab")
  }, priority = 10)

  
  # Dialog is only run if table_of_genes is loaded
  observeEvent(input$run_dea_B,{
    showModal(modalDialog(size = "l",
                          title = "Differential Expression Analysis started",
                          "The DEA has started!",
                          easyClose = TRUE,
                          footer = NULL
    ))
    
    if (input$run_dea_B == 1){
      print(head(gtf_df))
    }
  })

  observeEvent(input$run_dte_B,{
    showModal(modalDialog(size = "l",
                          title = "Differential Expression Analysis started",
                          "The DTE has started!",
                          easyClose = TRUE,
                          footer = NULL
    ))
    
    if (input$run_dte_B == 1){
      print(head(gtf_df))
    }
  })
  ### Select Colors####
    
  output$select_color_heat <- renderUI({
    fluidRow(
      box(
      column(12, div(style = "font-size: 13px", selectInput(inputId = "main_color_heat",label = "Select Colorpalette", choices = c("RdBu", "PiYG", "PRGn", "PuOr", "Blue-Yellow",
                                                                                                                  "Teal", "Sunset", "Viridis")))),
      width = 12, solidHeader = T, collapsible = T, collapsed = T, status = "primary", title = div("Customize Colors", style = "font-size: 16px")
    )
    )
  })

  outputOptions(output, "select_color_heat", suspendWhenHidden = FALSE)

  output$select_color_dte_heat <- renderUI({
    fluidRow(
      box(
      column(12, div(style = "font-size: 13px", selectInput(inputId = "main_color_dte_heat",label =  "Select Colorpalette", choices = c("RdBu", "PiYG", "PRGn", "PuOr", "Blue-Yellow",
                                                                                                                  "Teal", "Sunset", "Viridis")))),
      width = 12, solidHeader = T, collapsible = T, collapsed = T, status = "primary", title = div("Customize Colors", style = "font-size: 16px")
    )
    )
  })

  outputOptions(output, "select_color_dte_heat", suspendWhenHidden = FALSE)

  ######## DTU ANALYSIS #########
  preProcTrans <- reactiveValues(pre_list = NULL)
    observeEvent({input$run_dt_preprocess}, {
    req(dtu_file)
      
    print(input$tabs)
    print(input$dea.box)
    if (input$run_dt_preprocess > 0){
      preProcTrans$pre_list <- DRIM_seq_prep(dtu_file$d,
                                    condition_col = input$design_column, 
                                    first.level = input$feature_A, 
                                    ref.level = input$feature_B,
                                    cores = max(c(4,round(cores * 0.6)))
      )
      print("DRIMSEQPREP DONE!!!")
      print(length(preProcTrans$pre_list))
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      updateTabItems(session, "tabs", "run_overview")
      updateTabsetPanel(session, "main_quality_control",
                      selected = "read_length_dist_panel")
    }
    }, priority  = 0)
  

  #############DEXSeq Analysis##########################
  
  DTU_general_run <- reactiveValues(df_res_dte = NULL)
  
  observeEvent({input$run_dt_preprocess}, {
    req(dtu_file$dxd)
    if (is.null(dtu_file$dxd)) return()
    print(names(dtu_file$dxd))
 
    DTU_general_run$df_res_dte <- DTU_general(dtu_file$dxd,
                                  condition_col = input$design_column, 
                                  first.level = input$feature_A, 
                                  ref.level = input$feature_B, 
                                  gtf_table = gtf_df,
                                  cores = max(c(4,round(cores * 0.6))),
                                  pvalue_input = as.numeric(input$pvalue))    
  })
    
    
  
  ##################################################
  
  
  ###########DRIMSeq Analysis##########################  
  DTU_run <- reactiveValues(df_res_dtu = NULL)

  observeEvent(input$run_dtu_B, {
    print(input$tabs)
    print(preProcTrans$pre_list)
    req(preProcTrans$pre_list)
    
    if (is.null(preProcTrans$pre_list)) return()

    s = input$dtu_tab_rows_selected
    print("Input DTU rows")
    print(s)
    if (is.null(gtf_df$gene_id[input$dtu_tab_rows_selected])){
      print("Please select gene of interest in datatable")
    } else{

      DTU_run$df_res_dtu <- DTU_special(preProcTrans$pre_list,
                                        condition_col = input$design_column, 
                                        first.level = input$feature_A, 
                                        ref.level = input$feature_B,
                                        goi_id = gtf_df$gene_id[s],
                                        gtf_tab = gtf_df,
                                        cores = max(c(4,round(cores * 0.6))),
                                        pvalue_input = as.numeric(input$pvalue)
      )
    }
  })
  
  
  
  
  
  #################################################
  dea_res_preprocess <- reactiveValues(df_res = NULL)
  dte_res_preprocess <- reactiveValues(df_res = NULL)
  
  
  observeEvent({input$run_dea_B}, {
    req(dds_file$df)
    req(rld_file$df)
  
    if (is.null(dds_file$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (!(is.null(dds_file$df)) & !is.null(metadata)){
        print(input$pvalue)
        dea_res_preprocess$df_res <- run_preprocessing_dea(dds = dds_file$df,
                              rld = rld_file$df, 
                              condition.col = input$design_column,
                              first.level = input$feature_A, 
                              ref.level = input$feature_B, 
                              pvalue = as.numeric(input$pvalue), 
                              gtf_file = gtf_df
                              )
          save_rds(dea_res_preprocess$df_res, ".")
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      updateTabItems(session, "tabs", "run_overview")
      updateTabsetPanel(session, "main_quality_control",
                      selected = "read_length_dist_panel")
    }

  })
  
  dea_res <- reactive({
    req(dea_res_preprocess$df_res)    
    run_dea(condition.col = input$design_column, 
            first.level = input$feature_A, 
            ref.level = input$feature_B, 
            rld = rld_file$df, 
            dds = dds_file$df, 
            res_df = dea_res_preprocess$df_res$res_df,
            gtf_file = gtf_df, 
            condi_cols = condi_cols()
    )
  })
  
  
  ### Download buttons ####
  # Save button for PCA plot
  output$down_pca = downloadHandler(
    
    filename = function() {
        paste0("pca_", Sys.Date(), ".", input$down_pca.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)
      plot <- createPCA(rld_file$df, 
              first.level = input$feature_A, 
              ref.level = input$feature_B, 
              condi_cols()) + theme_bw() +
              theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_pca.type, width = 10, height = 10, bg = "white")
    }
  )
  
  # Save button for volcano plot
  output$down_vol = downloadHandler(
    
    filename = function() {
      paste0("DEA_volcano_plot_", Sys.Date(), ".", input$down_vol.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)
      plot <- createVolcano(dea_res_preprocess$df_res$res_df, condi_cols()) + theme_bw() +
              theme(
                legend.title = element_blank(),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_vol.type, bg = "white", width = 10, height = 7)
    }
  )
  
  # Save button for Sam2sam plot
  output$down_s2s = downloadHandler(
    
    filename = function() {
      paste0("DEA_sam2sam_", Sys.Date(), ".", input$down_s2s.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)

      if (input$down_s2s.type == "pdf"){
        pdf(file = file, width = 10, height = 10)
      } 
      if (input$down_s2s.type == "png"){
        png(file = file, width = 10, height = 20, units = "cm", res = 300)
      }
      draw(createSam2Sam(rld_file$df)[[2]])
      dev.off()
    }
  )
  
  # Save button for heatmap
  output$down_heat = downloadHandler(
    
    filename = function() {
      paste0("DEA_expression_heatmap_", Sys.Date(), ".", input$down_heat.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)
      color_plot <- condi_cols()
      plot <- createHeatmap(dds_file$df,
                           dds_file$rld, 
                           condi_col = color_plot, 
                           genes = dea_res_preprocess$df_res$res_df[1:20, ], 
                           main_color = input$main_color_heat)
      if (input$down_heat.type == "pdf"){
        pdf(file = file,width= 20, height = 20)
      } 
      if (input$down_heat.type == "png"){
        png(file = file,width= 20, height = 20, units = "cm", res = 300)
      }
      print(plot$heat.down)
      dev.off()
    }
  )

  ## Download button for DEX Volcano
  output$down_vol_dex = downloadHandler(
    
    filename = function() {
      paste0("DEX_volcano_plot_", Sys.Date(), ".", input$down_vol_dex.type)
    },
    content = function(file) {
      req(DTU_general_run$df_res_dte$volcano_plot)
      plot = DTU_general_run$df_res_dte$volcano_plot + 
        theme_bw() +
        theme(legend.title = element_blank(),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_vol_dex.type, bg = "white", width = 10, height = 7)
    }
  )


  ###############################################
    ## Differential transcript usage

  observeEvent({input$run_dte_B}, {
    req(dds_file$df_trans)
    print(dds_file$df_trans)
    if (is.null(dds_file$df_trans)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (!(is.null(dds_file$df_trans)) & !is.null(metadata)){
        cat(unlist(metadata))
        dte_res_preprocess$df_res <- run_preprocessing_dte(dds_file$df_trans, rld_file$df_trans,
                              condition.col = input$design_column,
                              first.level = input$feature_A, 
                              ref.level = input$feature_B, 
                              pvalue = as.numeric(input$pvalue), 
                              gtf_file = gtf_df
                              )
          save_rds_dte(dte_res_preprocess$df_res, "./")
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      updateTabItems(session, "tabs", "run_overview")
      updateTabsetPanel(session, "main_quality_control",
                      selected = "read_length_dist_panel")
    }

  })
  
  dte_res <- reactive({
    req(dte_res_preprocess$df_res)    
    run_dte(metadata = dte_res_preprocess$df_res$metadata, 
            counts = dte_res_preprocess$df_res$counts,
            condition.col = input$design_column, 
            first.level = input$feature_A, 
            ref.level = input$feature_B, 
            rld = dte_res_preprocess$df_res$rld, 
            dds = dte_res_preprocess$df_res$dds, 
            res_df = gtf_df, 
            condi_cols = condi_cols()
    )
  })
  
  output$text <- renderText({
    req(dte_res())
    if (length(dte_res()) == 0){
      "No file loaded!"
    } else {
      paste0(Sys.time(), " ", paste0(names(dte_res()), collapse = ","))
      
    }
  })
  
  ### Download buttons ####
  # Save button for PCA plot
  output$down_pca_dte = downloadHandler(
    
    filename = function() {
        paste0("pca_dte_", Sys.Date(), ".", input$down_pca_dte.type)
    },
    content = function(file) {
      req(dte_res_preprocess$df_res)
      plot <- createPCA_DTE(dte_res_preprocess$df_res$rld, 
              first.level = input$feature_A, 
              ref.level = input$feature_B, 
              condi_cols()) + theme_bw() +
              theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_pca_dte.type, width = 10, height = 10, bg = "white")
    }
  )
  
  # Save button for volcano plot
  output$down_vol_dte = downloadHandler(
    
    filename = function() {
      paste0("DTE_volcano_plot_", Sys.Date(), ".", input$down_vol_dte_type)
    },
    content = function(file) {
      req(dte_res_preprocess$df_res)
      plot <- createVolcano_DTE(dte_res_preprocess$df_res$res_df, condi_cols()) + theme_bw() +
              theme(
                legend.title = element_blank(),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_vol_dte_type, bg = "white", width = 10, height = 7)
    }
  )
  
  # Save button for Sam2sam plot
  output$down_s2s_dte = downloadHandler(
    
    filename = function() {
      paste0("DTE_sam2sam_", Sys.Date(), ".", input$down_s2s_dte_type)
    },
    content = function(file) {
      req(dte_res_preprocess$df_res)

      if (input$down_s2s_dte_type == "pdf"){
        pdf(file = file, width = 10, height = 10)
      } 
      if (input$down_s2s_dte_type == "png"){
        png(file = file, width = 10, height = 20, units = "cm", res = 300)
      }
      draw(createSam2Sam_DTE(dte_res_preprocess$df_res$rld)[[2]])
      dev.off()
    }
  )
  
  # Save button for heatmap
  output$down_heat_dte = downloadHandler(
    
    filename = function() {
      paste0("DTE_expression_heatmap_", Sys.Date(), ".", input$down_heat_dte_type)
    },
    content = function(file) {
      req(dte_res_preprocess$df_res)
      color_plot <- condi_cols()
      plot_dte <- createHeatmap_DTE(dte_res_preprocess$df_res$dds,
                           dte_res_preprocess$df_res$rld, 
                           condi_col = color_plot, 
                           transcripts = dte_res_preprocess$df_res$res_df[1:20, ], 
                           main_color = input$main_color_dte_heat)
      if (input$down_heat_dte_type == "pdf"){
        pdf(file = file,width= 20, height = 20)
      } 
      if (input$down_heat_dte_type == "png"){
        png(file = file,width= 20, height = 20, units = "cm", res = 300)
      }
      print(plot_dte$heat_dte.down)
      dev.off()
    }
  )

  ## Download button for DEX Volcano
  output$down_vol_dex = downloadHandler(
    
    filename = function() {
      paste0("DEX_volcano_plot_", Sys.Date(), ".", input$down_vol_dex.type)
    },
    content = function(file) {
      req(DTU_general_run$df_res_dte$volcano_plot)
      plot = DTU_general_run$df_res_dte$volcano_plot + 
        theme_bw() +
        theme(legend.title = element_blank(),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_vol_dex.type, bg = "white", width = 10, height = 7)
    }
  )

    ## Download button for DEX Volcano
  output$down_dtu_boxplot = downloadHandler(
    
    filename = function() {
      req(input$dtu_tab_rows_selected)
      paste0("DTU_", gtf_df$gene_id[input$dtu_tab_rows_selected], "_", Sys.Date(), ".", input$down_dtu_boxplot.type)
    },
    content = function(file) {
      req(DTU_run$df_res_dtu$bp)
      plot = DTU_run$df_res_dtu$bp + 
        theme_bw() +
        theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 8),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_dtu_boxplot.type, bg = "white", width = 10, height = 7)
    }
  )
  
  ## Sequencing run table ####
  seqRunInfos = reactiveValues(df = NULL) 
  
  observeEvent(input$tabs,ignoreNULL = T,{
    if (input$tabs == "run_overview"){
              seqRunInfos$df <- reactivePoll(10000, session,
                checkFunc = function(){
                  mapStats="mapping_stats.txt"
                  if (file.exists(mapStats)){   
                    file.info(mapStats)$mtime[1]
                  } else {
                    print("Not processed yet")
                  }
                },
                valueFunc = function(){
                  
                  df_out = data.frame("Samples" = metadata$Samples, "mapped_reads" = 0, "gene_counts" = 0, "transcript_counts" = 0)
                  #print(df_out)

                  # Num Q reads and mapped reads
                  if (file.exists("mapping_stats.txt")){
                    print("Mapping statistics exists.")
                    mapStats = read.table("mapping_stats.txt", header = T, sep = "\t")
                    print(head(mapStats))
                    colnames(mapStats) = c("Samples","mapped_reads", "X")
                    mapStats = mapStats[match(df_out$Samples, mapStats$Samples),]
                    df_out$mapped_reads = mapStats$mapped_reads
                  }

                  ## Num gene counts 
                  if (!is.null(dds_file$df)){
                    counts = counts(dds_file$df)
                    print(head(counts))
                    print(colSums(counts))
                    df_out$gene_counts =  as.integer(colSums(counts)[c(df_out$Samples)])
                  } else {
                    print("Counts file not present")
                  }
                  
                  ## Num gene counts 
                  if (!is.null(dds_file$df_trans)){
                    print(head(dds_file$df_trans))
                    counts_trans = counts(dds_file$df_trans)
                    print(colSums(counts_trans))
                    print(as.integer(colSums(counts_trans)[c(df_out$Samples)]))
                    df_out$transcript_counts =  as.integer(colSums(counts_trans)[c(df_out$Samples)])
                  } else {
                    print("Transcript counts file not present")
                  }
                  df_out
                })
               }
    },
               priority = -1
  )


  ### Tables and Plots ####
  # Show table of seq. run infos
  output$table_seq_run <- renderDT({

    #print(seqRunInfos$df())
    DT::datatable(
      seqRunInfos$df(),
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  },server = FALSE)


  process_time = reactiveValues(df = NULL)
  
  observeEvent(input$tabs,ignoreNULL = T,{
    if (input$tabs == "run_overview"){
              process_time$df <- reactivePoll(10000, session,
                checkFunc = function(){
                  if (file.exists("processing_time_table.csv")){
                    file.info("processing_time_table.csv")$mtime[1]
                  }
                },
                valueFunc = function(){
                  if (file.exists("processing_time_table.csv")){
                    read.table("processing_time_table.csv", header = T, sep = ",")
                  } else {
                    NULL
                  }
                })
               }},
               priority = -1
  )
  ## Plot Process time per Iteration
   output$process_time.out <- renderPlot({
    req(process_time$df())
    if (!is.null(process_time$df())){
      print(head(process_time$df()))
      x = process_time$df()
      # replace _ with whitespace in tool names for a cleaner legend
      x = x %>% mutate(Tool = gsub("_", " ", Tool))
      
      # get the x-tick interval based on the number of iterations in the data
      max_it <- max(x$Iteration)
      xticks <- 1:max_it
      if (max_it <= 10) { # every tick
        xticks <- xticks
      } else if (max_it <= 50) { # every 5th tick
        xticks <- xticks[xticks %% 5 == 0]
      } else if (max_it <= 100) { # every 10th tick
        xticks <- xticks[xticks %% 10 == 0]
      } else if (max_it <= 1000) { # every 50th tick
        xticks <- xticks[xticks %% 50 == 0]
      } else if (max_it <= 10000) { # every 100th tick
        xticks <- xticks[xticks %% 100 == 0]
      } else { # every 1000th tick
        xticks <- xticks[xticks %% 1000 == 0]
      }

      ggplot(x, aes(x = Iteration, y = Time, fill = Tool)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = xticks) + # label each x tick individually
        ylab("Time [s]") + # add [s] to y-axis label
        theme_bw() +
         scale_fill_manual("Preprocessing steps", values = c("#CC6677", "#d8c66c", "#117733", "#88CCEE", "#AA4499", "#24011e","#092a36")) +
        theme(
        panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.title = element_text("", size = 20, color = "white"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 20, color = "white"),
        axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
        axis.title = element_text(size = 23, color = "white"))

    }
  }, bg="transparent")
  
  
  
  # Show table of DEGs 
  output$degs_tab <- renderDT({
    req(dea_res_preprocess$df_res)
    res = dea_res_preprocess$df_res$res_df
    degs = res[res$padj < as.numeric(input$pvalue),]
    colnames(degs) = c("gene name","gene base mean","log 2 fold change","lfcSE", "stat", "p-value", "p-adjusted","significance","gene ID")
    DT::datatable(
      degs,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  }, server = FALSE)

  output$dtes_tab <- renderDT({
    req(dte_res_preprocess$df_res)
    res = dte_res_preprocess$df_res$res_df
    dtes = res[res$padj < as.numeric(input$pvalue),]
    colnames(dtes) = c("transcript name","transcript base mean","log 2 fold change","lfcSE", "stat", "p-value", "p-adjusted","significance","gene ID")
    DT::datatable(
      dtes,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  }, server = FALSE)
  

  output$drim_tab <- renderDT({
    req(preProcTrans$pre_list)
    res = data.frame(preProcTrans$pre_list$res_df)
    res = na.omit(res)
    res = res[res$adj_pvalue < as.numeric(input$pvalue),]
    gtf_table <- gtf_df
    res$gene_name = gtf_table[match(res$feature_id, gtf_table$transcript_id), "gene_name"]
    res$transcript_name = gtf_table[match(res$feature_id, gtf_table$transcript_id), "transcript_name"]
    dtus = data.frame(gene_name = as.character(res$gene_name), transcript_name = as.character(res$transcript_name),gene_id = as.character(res$gene_id),feature_id = as.character(res$feature_id),lr = as.numeric(res$lr), pvalue = as.numeric(res$pvalue),adj_pvalue = as.numeric(res$adj_pvalue))
    colnames(dtus) = c("gene name","transcript name","gene ID","feature ID", "lr","pvalue", "p-adjusted")
    DT::datatable(
      dtus,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  }, server = FALSE)



  # Show table of DTUs
  output$dex_tab <- renderDT({
    req(DTU_general_run$df_res_dte$dxr_df)
    res = DTU_general_run$df_res_dte$dxr_df
    dtus = res[res$padj < as.numeric(input$pvalue),]
    dtus = data.frame(gene_name = as.character(dtus$gene_name), transcript_name = as.character(dtus$transcript_name),group_id = as.character(dtus$groupID), feature_id = as.character(dtus$featureID),exon_base_mean = as.numeric(as.character(dtus$exonBaseMean)),p_adjusted = as.numeric(as.character(dtus$padj)), dispersion = as.numeric(as.character(dtus$dispersion)))
    colnames(dtus) = c("gene name","transcript name","gene ID","feature ID", "exon base mean", "p-adjusted", "dispersion")
    DT::datatable(
      dtus,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  }, server = FALSE)
  
  output$dtu_tab <- DT::renderDataTable({
    config = gtf_df
    
    config_out = config
    DT::datatable(
      config_out,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     pageLength = 5,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
      selection = list(mode = "single"
                       , target = "row"),
      rownames = T)
  }, server = FALSE)
  
  output$selectedGene_dtu <- DT::renderDataTable({
    s = gtf_df[input$dtu_tab_rows_selected,]
    config = s
    config_out = config
    DT::datatable(
      config_out, 
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 't',
                     fixedColumns = TRUE),
      rownames = T)
    
  }, server = FALSE)
  
  # Section One ---------------------------------
  
  # Plot volcano plot stored in RDS object
  output$volcano_plot <- renderPlot({
    req(dea_res_preprocess$df_res)
    createVolcano(dea_res_preprocess$df_res$res_df, condi_cols())
  }, bg="transparent")
  
  
  output$pca_plot <- renderPlot({
    print("PCA is made")
    req(dea_res_preprocess$df_res)
    createPCA(dea_res_preprocess$df_res$rld, 
              first.level = input$feature_A, 
              ref.level = input$feature_B, 
              condi_cols())   
  }, bg="transparent")
  
  output$s2s_plot <- renderPlot(bg = "transparent",{
    print("Heatmap")

    req(dea_res_preprocess$df_res)
    createSam2Sam(dea_res_preprocess$df_res$rld)[[1]]
    
  })
  
  output$heat_plot <- renderPlot(bg = "transparent", {
    print("Heatmap")
    req(dea_res_preprocess$df_res)
    plot = createHeatmap(dea_res_preprocess$df_res$dds,
                  dea_res_preprocess$df_res$rld, 
                  condi_col = condi_cols(), 
                  genes = dea_res_preprocess$df_res$res_df[1:20, ], 
                  main_color = input$main_color_heat)
    plot$heat
  })

  ### Output for DTE
  output$volcano_plot_dte <- renderPlot({
    req(dte_res_preprocess$df_res)
    createVolcano_DTE(dte_res_preprocess$df_res$res_df, condi_cols())
  }, bg="transparent")
  
  
  output$pca_plot_dte <- renderPlot({
    print("PCA is made")
    req(dte_res_preprocess$df_res)
    createPCA_DTE(dte_res_preprocess$df_res$rld, 
              first.level = input$feature_A, 
              ref.level = input$feature_B, 
              condi_cols())   
  }, bg="transparent")
  
  output$s2s_plot_dte <- renderPlot(bg = "transparent",{
    print("Heatmap")
    
    req(dte_res_preprocess$df_res)
    createSam2Sam_DTE(dte_res_preprocess$df_res$rld)[[1]]
    
  })
  
  output$heat_plot_dte <- renderPlot(bg = "transparent", {
    print("Heatmap")
    req(dte_res_preprocess$df_res)
    plot_dte = createHeatmap_DTE(dte_res_preprocess$df_res$dds,
                  dte_res_preprocess$df_res$rld, 
                  condi_col = condi_cols(), 
                  transcripts = dte_res_preprocess$df_res$res_df[1:20, ], 
                  main_color = input$main_color_dte_heat)
    plot_dte$heat_dte
  })
  
  output$volcano_dex_output <- renderPlot(bg = "transparent",{
    print("DEX Volcano plot")
    req(DTU_general_run$df_res_dte$volcano_plot)
    DTU_general_run$df_res_dte$volcano_plot
  })
  
  output$dtu_boxplot <- renderPlot(bg = "transparent",{
    print("DTU Boxplot")
    req(DTU_run$df_res_dtu)
    DTU_run$df_res_dtu$bp
  })
    
  output$dtu_boxplot_UI <- renderUI({
      if (!is.null(DTU_run$df_res_dtu)) {
        if (is.null(DTU_run$df_res_dtu$bp)) {
          h3(p(em(strong("No transcripts found!"))))
        } else {
          column(12, 
          div(style = "display: inline-block; vertical-align: top", downloadButton("down_dtu_boxplot")),
          div(style = "display: inline-block; vertical-align: top; width: 80px", selectInput("down_dtu_boxplot.type", NULL, choices = c("png", "pdf"))),
          plotOutput("dtu_boxplot") %>% withSpinner(color = "#0dc5c1"),
          br(), p(em(fig_des$single_gene_deu)))
        }
      }
  })

  output$dea_results.out <- renderUI({
    if (is.null(input$run_dea_B)) return()
    if(input$run_dea_B > 0){
      fluidRow(
          column(12, box(
            title = "Differential expressed genes",
            width = 12,
            status = "primary",
            solidHeader = T,
            collapsible = T,
            collapsed = F,
            DT::dataTableOutput("degs_tab") %>% withSpinner(color = "#0dc5c1")
          )),
          column(
            12,
            tabBox(
              width = 12, title = "Plots",
              id = "tabset.plots",
              tabPanel(
                title = "PCA", value = "pca.tab",
                fluidRow(column(
                  12,
                  div(style = "display:inline-block; vertical-align: top", downloadButton("down_pca")),
                  div(style = "display:inline-block; vertical-align: top; width: 80px", selectInput("down_pca.type", NULL, choices = c("png", "pdf"))),
                  plotOutput("pca_plot") %>% withSpinner(color = "#0dc5c1")
                ),
                column(12, br(), p(em(fig_des$pca_genes))))
              ),
              tabPanel(
                title = "Volcano", value = "volcano.tab",
                fluidRow(column(
                  12,
                  div(style = "display: inline-block; vertical-align: top", downloadButton("down_vol")),
                  div(style = "display: inline-block; vertical-align: top; width: 80px", selectInput("down_vol.type", NULL, choices = c("png", "pdf"))),
                  plotOutput("volcano_plot") %>% withSpinner(color = "#0dc5c1")
                ),
                column(12, br(), p(em(fig_des$volcano_genes))))
              ),
              tabPanel(
                title = "Sample2Sample", value = "s2s.tab",
                fluidRow(column(
                  12,
                  div(style = "display:inline-block; vertical-align: top", downloadButton("down_s2s")),
                  div(style = "display:inline-block; vertical-align: top; width: 80px", selectInput("down_s2s.type", NULL, choices = c("png", "pdf"))),
                  plotOutput("s2s_plot") %>% withSpinner(color = "#0dc5c1")
                ),
                column(12, br(), p(em(fig_des$sample2sample_genes))))
              ),
              tabPanel(
                title = "Heatmap", value = "heatmap.tab",
                fluidRow(column(
                  12,
                  div(style = "display:inline-block; vertical-align:top", downloadButton("down_heat")),
                  div(style = "display:inline-block; vertical-align:top; width:80px", selectInput("down_heat.type", NULL, choices = c("png", "pdf"))),
                  uiOutput("select_color_heat"),
                  plotOutput("heat_plot") %>% withSpinner(color = "#0dc5c1")
                ),
                column(12, br(), p(em(fig_des$heatmap_genes))))
              )
            ) # tabBox Plots
          )
        ) # fluidRowcolumn
    }
  })


  output$dte_results.out <- renderUI({
  if (is.null(input$run_dte_B)) return()
  if(input$run_dte_B > 0){
    fluidRow(
        column(12, box(
          title = "Differential expressed transcripts",
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = T,
          collapsed = F,
          DT::dataTableOutput("dtes_tab") %>% withSpinner(color = "#0dc5c1")
        )),
        column(
          12,
          tabBox(
            width = 12, title = "Plots",
            # The id lets us use input$tabset1 on the server to find the current tab
            tabPanel(
              title = "PCA", value = "pca_dte.tab",
              fluidRow(column(
                12,
                div(style = "display:inline-block; vertical-align: top", downloadButton("down_pca_dte")),
                div(style = "display:inline-block; vertical-align: top; width: 80px", selectInput("down_pca_dte.type", NULL, choices = c("png", "pdf"))),
                plotOutput("pca_plot_dte") %>% withSpinner(color = "#0dc5c1")
              ),
              column(12, br(), p(em(fig_des$pca_transcripts))))
            ),
            tabPanel(
              title = "Volcano", value = "volcano_dte.tab",
              fluidRow(column(
                12,
                div(style = "display: inline-block; vertical-align: top", downloadButton("down_vol_dte")),
                div(style = "display: inline-block; vertical-align: top; width: 80px", selectInput("down_vol_dte_type", NULL, choices = c("png", "pdf"))),
                plotOutput("volcano_plot_dte") %>% withSpinner(color = "#0dc5c1")
              ),
              column(12, br(), p(em(fig_des$volcano_transcripts))))
            ),
            tabPanel(
              title = "Sample2Sample", value = "s2s_dte.tab",
              fluidRow(column(
                12,
                div(style = "display:inline-block; vertical-align: top", downloadButton("down_s2s_dte")),
                div(style = "display:inline-block; vertical-align: top; width: 80px", selectInput("down_s2s_dte_type", NULL, choices = c("png", "pdf"))),
                plotOutput("s2s_plot_dte") %>% withSpinner(color = "#0dc5c1")
              ),
              column(12, br(), p(em(fig_des$sample2sample_transcrips))))
            ),
            tabPanel(
              title = "Heatmap", value = "heatmap_dte.tab",
              fluidRow(column(
                12,
                div(style = "display:inline-block; vertical-align:top", downloadButton("down_heat_dte")),
                div(style = "display:inline-block; vertical-align:top; width:80px", selectInput("down_heat_dte_type", NULL, choices = c("png", "pdf"))),
                uiOutput("select_color_dte_heat"),
                plotOutput("heat_plot_dte") %>% withSpinner(color = "#0dc5c1")
              ),
              column(12, br(), p(em(fig_des$heatmap_transcripts))))
            )
          ) # tabBox Plots
        )
      ) # fluidRowcolumn
    }
  })

 

  
  # Kill all processes when app is closed
  onStop(function(){
    cat(sprintf("Session %s was closed\n", session$token))
    if(!is.null(myProcess)){
      if(myProcess$is_alive()){
        myProcess$kill()
        unlink("work", recursive = TRUE) 
      }
    }
    
  })
  
}
