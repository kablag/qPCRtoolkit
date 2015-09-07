library(shiny)
library(plyr)
library(dplyr)
library(tidyr)
library(rlist)

library(ggvis)
library(ggplot2)
library(scales)
library(chipPCR)
library(xlsx)

# if(installed.packages(.Library) %>%
#    as.data.frame %>%
#    filter(Package == "RDML",
#           Version == "0.8-4") %>% nrow() %>% "=="(0)) {
#   library(devtools)
#   install_github("kablag/RDML")
# }

library(RDML)

source('df2html.R')
source("relativeCalc.R")

smart.round <- function(x) {
  if(is.na(x)) return("")
  if(x == 0) return("0")
  if(x < 1 || x > 1000)
    return(format(x, digits = 3, scientific = TRUE))
  sprintf("%.2f", x)
}

smart.round.col <-
  function(col) {
    if(is.numeric(col))
      sapply(col, smart.round)
    else col
  }

mylog10_trans <- function(from=0) 
{
  base = 10
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}

shinyServer(function(input, output, session) {
  
  values <- reactiveValues(experiments = list(),
                           merged = NULL,
                           tbl = NULL,
                           tblCq = NULL,
                           calib = NULL,
                           calcRelative = FALSE)
  
  
  
  
  
  observe({
    if(is.null(input$rdmlFile))    
      return(NULL)    
    cat("adding experiment\n")
    path <- input$rdmlFile$datapath
    isolate({
      values$experiments <- 
        c(values$experiments,
          RDML$new(path,
                   conditions.sep = "@"))
      values$mergedRDML <- NULL
      values$tbl <- NULL
      values$tblCq <- NULL
      values$calib <- NULL
      if(!is.null(substrCalced))
        substrCalced <- FALSE
    })
  })
  
  observe({
    if (input$useExampleFileBtn == 0 )
      return(NULL)
    cat("adding example experiment\n")
    isolate({
      values$experiments <-
        c(values$experiments,
          RDML$new("www/data/Demo_Rel Quant Dual Color.lc96p",
                   conditions.sep = "@"))
      values$mergedRDML <- NULL
      values$tbl <- NULL
      values$tblCq <- NULL
      values$calib <- NULL
      if(!is.null(substrCalced))
        substrCalced <- FALSE
    })
  })
  
  output$experimentsChecksUi <- renderUI({
    if(length(values$experiments) == 0)
      return(NULL)
    cat("creating experimentsChecksUi\n")
    isolate({
      values$expIds <- c()
      values$expTexts <- c()
      for(i in 1:length(values$experiments)){
        values$expIds <- 
          c(values$expIds,
            values$experiments[[i]]$experiment[[1]]$id$id)
        values$expTexts <- 
          c(values$expTexts,                 
            sprintf("%s",
                    values$experiments[[i]]$dateMade)
          )
      }
      names(values$expIds) <- values$expTexts
      names(values$expTexts) <- values$expIds
      checkboxGroupInput(
        "experimentsChecks", 
        "Experiments:",
        values$expIds,
        selected = values$expIds)
    })
  })
  
  observe({  
    if(is.null(input$experimentsChecks) ||
       length(input$experimentsChecks) == 0) {
      values$mergedRDML <- NULL
      return(NULL)
    }
    cat("creating mergedRDML\n")
    # print(input$experimentsChecks)
    isolate({
      values$mergedRDML <- MergeRDMLs(
        # values$experiments[[which(values$expIds == input$experimentsChecks)]]
        values$experiments
      )
        # RDML$new()
#       l_ply(input$experimentsChecks,
#             function(experimentName) {
#               values$mergedRDML$Merge(
                # values$experiments[[which(values$expIds == experimentName)]])
            # })
      #       values$calcRelative <- {
      #         types <- laply(values$mergedRDML$target,
      #                        function(target) target$type)
      #         if("ref" %in% types && "toi" %in% types)
      #           TRUE
      #         else FALSE
      #       }
    })
  })
  
  output$filterTargetsUi <- renderUI({
    if(is.null(values$tbl))
      return(NULL)
    cat("creating filterTargetsUi\n")
    unique(values$tbl$target) %>%      
      selectInput("filterTargets", "Targets",
                  choices = .,
                  selected = .,
                  multiple = TRUE)
    
  })
  
  output$plateTblExpSelectorUi <- renderUI({
    if((input$experimentsChecks %>% length) == 0)
      return(NULL)
    cat("creating plateTblExpSelectorUi\n")
    selectInput("plateTblExpSelector", "Experiment",
                choices =  values$expIds[values$expTexts[input$experimentsChecks]])
  })
  
  observe({
    if(is.null(values$mergedRDML))
      return(NULL)
    cat("creating mergedTbl\n")
    values$tbl <-
      values$mergedRDML$AsTable(
        name.pattern = paste(
          react$Position(run$pcrFormat),
          react$sample$id,
          private$.sample[[react$sample$id]]$type$value,
          data$tar$id,
          sep = "_"),
        minCyc = min(data$adp$fpoints[, "cyc"]),
        maxCyc = max(data$adp$fpoints[, "cyc"]),
        conditions = {
          if(id[[1]]$publisher == "Roche Diagnostics") {
            cond <- list.filter(sample[[react$sample$id]]$annotation,
                           .$property == sprintf("Roche_condition_at_%s",
                                               react$id$id))
            if (length(cond) == 0)
                   NA
            else
              cond[[1]]$value
          }
          else NA
        }, 
        quantity = 
          if(id[[1]]$publisher == "Roche Diagnostics") {
            quantity <- list.filter(sample[[react$sample$id]]$annotation,
                               .$property == sprintf("Roche_quantity_at_%s_%s",
                                                   target[[data$tar$id]]$dyeId$id,
                                                   react$id$id))
            if (length(quantity) == 0)
                   NA
            else
              quantity[[1]]$value
          }
        else
          sample[[react$sample$id]]$quantity$value
      )
  })
  
  output$limitCyclePanel <- renderUI({
    if(is.null(values$tbl) ||
       is.null(values$mergedRDML))
      return(NULL)
    cat("creating limitCyclePanel\n")    
    
    c(min(values$tbl$minCyc),
      max(values$tbl$maxCyc)) %>%        
      sliderInput("limitCycleFromTo",
                  "Limit Cycles",
                  value = .,
                  min = .[1],
                  max = .[2],
                  step = 1)
    
  })
  
  output$plateTbl<- renderUI({
    if(is.null(input$plateTblExpSelector))
      return(NULL)
    cat("creating plateTbl\n")
    cellClass <- aaply(1:96, 1,
                       function(i) {
                         values$tbl %>% 
                           filter(
                             exp.id == input$plateTblExpSelector,
                             react.id == i) %>%
                           select(sample.type) %>% 
                           .[1, 1]
                       }) %>% 
      matrix(8, 12, byrow = TRUE,
             dimnames = list(LETTERS[1:8], 1:12))
    
    cellText <- cellClass
    cellText[is.na(cellText)] <- "    "
    cellText <- as.data.frame(cellText) %>% 
      cbind(LETTERS[1:8], .)
    colnames(cellText) <- c("", 1:12)
    
    cellClass <- cbind(rep(NA, 8), cellClass)
    
    HTML(df2html(
      cellText,
      class = "tbl selCell",
      id = "plateTbl",
      cellClass = cellClass)
    )
  })
  
  output$backgroundPanel <- renderUI({
    if(is.null(input$limitCycleFromTo))
      return(NULL)
    cat("creating backgroundPanel\n")
    sliderInput("backgroundFromTo",
                "Background Cycles",
                value = c(3, 10),
                min = input$limitCycleFromTo[1],
                max = input$limitCycleFromTo[2],
                step = 0.1)
    
  })
  #   
  substrCalced <- reactive({    
    if(is.null(input$backgroundFromTo) ||
       is.null(values$mergedRDML) ||
       input$limitCycleFromTo[2] < input$backgroundFromTo[2] ||
       input$limitCycleFromTo[1] > input$backgroundFromTo[1])
      return(NULL)
    cat("rcalc substr\n")
    fdata <- values$mergedRDML$GetFData(values$tbl,
                                        limits = input$limitCycleFromTo) 
    substr <-       
      cbind(cyc = fdata[, 1],
            apply(fdata[, -1], 2,
                  function(x) CPP(fdata[, 1],
                                  x,
                                  bg.range =                                          
                                    input$backgroundFromTo,
                                  smoother = FALSE,
                                  method.norm = "none")$y))
    tbl <- values$tbl
    tbl$exp.id <- sprintf("%s_substr", tbl$exp.id)
    values$mergedRDML$SetFData(substr, tbl)
    values$tbl2 <- 
      values$mergedRDML$AsTable(
        name.pattern = paste(
          react$Position(run$pcrFormat),
          react$sample$id,
          private$.sample[[react$sample$id]]$type$value,
          data$tar$id,
          sep = "_"),
        minCyc = min(data$adp$fpoints[, "cyc"]),
        maxCyc = max(data$adp$fpoints[, "cyc"]),
        conditions = {
          if(id[[1]]$publisher == "Roche Diagnostics") {
            cond <- list.filter(sample[[react$sample$id]]$annotation,
                                .$property == sprintf("Roche_condition_at_%s",
                                                      react$id$id))
            if (length(cond) == 0)
              NA
            else
              cond[[1]]$value
          }
          else NA
        }, 
        quantity = 
          if(id[[1]]$publisher == "Roche Diagnostics") {
            quantity <- list.filter(sample[[react$sample$id]]$annotation,
                                    .$property == sprintf("Roche_quantity_at_%s_%s",
                                                          target[[data$tar$id]]$dyeId$id,
                                                          react$id$id))
            if (length(quantity) == 0)
              NA
            else
              quantity[[1]]$value
          }
        else
          sample[[react$sample$id]]$quantity$value
      )
    TRUE    
  })
  #  
  output$thresholdPanel <- renderUI({
    if( is.null(substrCalced()) ||
        substrCalced() == FALSE)
      return(NULL)
    if(input$useBackgroundSubstruction == TRUE)
      data.type = ".*_substr$"
    else
      data.type = ".*"
    
    
    cat("create threshold panel\n")
    
    tryCatch(
      lapply(unique(values$tbl$target),
             function(trgt){
               
               ftbl <-  filter(values$tbl2,
                               grepl(data.type, exp.id),
                               target == trgt)
               fluo <-          
                 values$mergedRDML$GetFData(ftbl,
                   limits = input$limitCycleFromTo)[, -1]
               minfluo <- min(fluo) %>% smart.round %>% as.numeric
               maxfluo <- max(fluo) %>% smart.round %>% as.numeric
               sliderInput(paste("thresholdValue", trgt, sep = "_"), 
                           sprintf("Threshold %s",
                                   trgt), 
                           min = minfluo,
                           max = maxfluo,
                           value = maxfluo * 0.1,
                           step = (maxfluo - minfluo) / 100,
                           sep = NA)
             }
      ),
      error = function(e) return(NULL)
      
    )
    
  })  
  #  
  observe({
    if( is.null(substrCalced()) ||
        substrCalced() == FALSE ||
        is.null(input[[sprintf("thresholdValue_%s", 
                               unique(values$tbl$target)[1])]]))
      return(NULL)
    cat("calc Cq\n")
    isolate({
      values$tblCq <- NULL
      values$calib <- NULL
      if(input$useBackgroundSubstruction == TRUE)
        data.type = ".*_substr$"
      else
        data.type = ".*"
      fdata <- values$mergedRDML$GetFData(filter(values$tbl2,
                                                 grepl(data.type, exp.id)),
                                          limits = input$limitCycleFromTo) 
      
      # calc Cq
      
      tryCatch({
        values$tblCq <- 
          values$tbl2 %>%
          group_by(fdata.name) %>%
          filter(grepl(data.type, exp.id)) %>% 
          mutate(Cq = as.numeric(tryCatch(
            th.cyc(fdata[, 1],
                   fdata[, fdata.name],
                   # auto = TRUE
                   r = input[[sprintf("thresholdValue_%s", 
                                      target)]]
            )@.Data[1],
            error = function(e) NA))
          ) %>% 
          group_by(sample, target, quantity, conditions)  %>%
          mutate(Cq_Mean = mean(Cq, na.rm=TRUE),
                 Cq_SD = sd(Cq, na.rm=TRUE),
                 Cq_Var = var(Cq, na.rm=TRUE) %>% 
                   ifelse(is.na(.),
                          0,
                          .))},
        warning = function(w) NULL,
        error = function(e) NULL
      )
      cat("calc calib\n")
      if(!("std" %in% values$tblCq$sample.type))
        return(NULL)
      targets <- unique(values$tblCq$target)
      values$calib <- 
        lapply(targets,
               function(trgt) {        
                 data <- filter(values$tblCq, 
                                target == trgt,
                                sample.type == "std")
                 if(nrow(data) == 0) return(NA)
                 q.list <- list()        
                 for(n in unique(data$quantity)){
                   q.list[[as.character(n)]] <- data$Cq[which(data$quantity == n)]
                 }
                 lmax <- max(sapply(q.list, length))
                 q.table <-c()
                 for(i in 1:length(q.list)) {
                   l <- length(q.list[[i]])
                   if(l < lmax) 
                     q.table <- rbind(q.table,
                                      c(as.numeric(names(q.list)[i]),
                                        q.list[[i]],
                                        rep(NA, lmax-l)))
                   else 
                     q.table <- rbind(q.table,
                                      c(as.numeric(names(q.list)[i]),
                                        q.list[[i]]))
                 }
                 tryCatch({
                   eff <- effcalc(q.table[, 1], q.table[, 2:ncol(q.table)]) 
                   list(eff = eff,
                        slope = eff@regression$coefficients[["res[, 1]"]],
                        intercept = eff@regression$coefficients[["(Intercept)"]],
                        efficiency = 10 ^ (-1 / eff@regression$coefficients[["res[, 1]"]]))
                   
                 },
                 error = function(e) {
                   # print(e)
                   NULL
                 }
                 )
               }
        )
      names(values$calib) <- targets
      cat("calc conc\n")
      if(is.null(values$tblCq))
        return(NULL)
      values$tblCq <-
        values$tblCq %>%
        group_by(fdata.name) %>%
        mutate(quantity = 
                 as.numeric({
                   if(sample.type == "std") {
                     quantity
                   } else if (is.na(Cq) ||
                              is.null(values$calib[[target]]) ||
                              is.na(values$calib[[target]])) {
                     NA
                   } else {
                     10 ^ ((Cq - values$calib[[target]]$intercept)/values$calib[[target]]$slope)
                   }
                 }),
               log10quantity = log10(quantity)
        ) %>% 
        group_by(sample, target)  %>%     
        mutate(quantity_mean = 
                 as.numeric({
                   if("std" %in% sample.type)
                     quantity
                   else
                     mean(quantity, na.rm=TRUE)}),
               quantity_SD = as.numeric({
                 if("std" %in% sample.type)
                   0
                 # TODO if na then na
                 else
                   sd(quantity, na.rm=TRUE)})
        )
    })
  })
  
  
  observe({ 
    if(is.null(substrCalced()) ||
       substrCalced() == FALSE || 
       length(input$filterTargets) == 0 ||
       is.null(input[[sprintf("thresholdValue_%s", 
                              input$filterTargets[[1]])]]))
      return(NULL)
    cat("plot curves\n")
    if(input$useBackgroundSubstruction == TRUE)
      data.type = ".*_substr$"
    else
      data.type = ".*"
    
    if(is.null(input$plateTbl)) {
      filteredData <- filter(values$tbl2,
                             exp.id == {
                               if(input$useBackgroundSubstruction == TRUE)
                                 paste0(input$plateTblExpSelector, "_substr")
                               else
                                 input$plateTblExpSelector
                             },
                             target %in% input$filterTargets)
    } else {
      react.ids <- 
        ((input$plateTbl[1:(length(input$plateTbl)/2)] - 1 ) * 12 + 
           input$plateTbl[(length(input$plateTbl)/2 + 1):length(input$plateTbl)]) %>% 
        as.character
      filteredData <- 
        filter(values$tbl2,
               exp.id == {
                 if(input$useBackgroundSubstruction == TRUE)
                   paste0(input$plateTblExpSelector, "_substr")
                 else
                   input$plateTblExpSelector
                 },
               target %in% input$filterTargets,
               react.id %in% react.ids)
    }
    if(nrow(filteredData) == 0)
      return(NULL)
    fluo <- 
      values$mergedRDML$GetFData( 
        filteredData,
        limits = input$limitCycleFromTo,
        long.table = TRUE)
    fluo$color <- fluo[[input$colorBy]]
    threshold <- input$filterTargets %>%      
      ldply(function (trgt) 
        data.frame(target = rep(trgt, 2),
                   thresholdValY = rep(input[[sprintf("thresholdValue_%s", 
                                                      trgt)]], 2),
                   thresholdValX = input$limitCycleFromTo)
      )
    fluo$gr <-as.character(fluo$fdata.name)
    if(input$curveType == "log"){
      fluo$fluo <- log(fluo$fluor)
      gr <- 0
      fd.name <- ""
      for(irow in 1:nrow(fluo)) {
        if(fluo$fdata.name[irow] != fd.name) {
          fd.name <- fluo$fdata.name[irow]
          gr <- gr + 1
        }
        if(is.nan(fluo$fluor[irow]))
          gr <- gr + 1
        fluo[irow, "gr"] <- as.character(gr)
      }
      threshold$thresholdValY <- log(threshold$thresholdValY)
    }
    tryCatch({
      minfluo <- round(min(fluo$fluor, na.rm = TRUE), digits = 2)
      maxfluo <- round(max(fluo$fluor, na.rm = TRUE), digits = 2)
    },
    warning = function(w) return(NULL)      
    )
    background <- data.frame(
      x1 = input$backgroundFromTo[1],
      x2 = input$backgroundFromTo[2],
      y1 = max(fluo$fluor, na.rm = TRUE),
      y2 = min(fluo$fluor, na.rm = TRUE))
    fluo %>% 
      filter(!is.nan(fluor)) %>%      
      group_by(as.factor(gr)) %>%      
      ggvis(~cyc, ~fluor) %>%
      layer_paths(stroke.hover = ~color,
                  stroke = ~color,
                  strokeWidth.hover := 4,
                  strokeWidth := 2) %>%      
      layer_paths(data = group_by(threshold, target),
                  x = ~thresholdValX,
                  y = ~thresholdValY,
                  stroke = ~target) %>%      
      add_tooltip(function(data){
        if("as.factor(gr)" %in% names(data)) {
          fdata <- filter(fluo, gr == data[["as.factor(gr)"]])[1,]        
          paste(fdata$position,
                fdata$sample,
                fdata$target,
                fdata$sample.type,
                sep = "<br>")
        }
        else {
          data[["target"]]
        }
      }, "hover") %>%
      layer_rects(data = background,
                  ~x1, ~y1, x2 = ~x2, y2 = ~y2,
                  fillOpacity := 0.1,
                  strokeOpacity := 0.1) %>%
      add_axis("y", title = "RFU",
               title_offset = 50) %>% 
      add_axis("x", title = "Cycles") %>% 
      add_legend(scales = "stroke", title = "") %>% 
      bind_shiny("curvesPlot", "curvesPlot_ui")
    if(values$calib %>% length == 0)
      return(NULL)
    cat("plot effsd\n")
    effsd <- values$calib %>%
      ldply(function(calib) {
        if(!is.na(calib) && !is.null(calib$eff))
          data.frame(calib$eff@.Data)
      },
      .id = "target") %>%
      filter(target %in% input$filterTargets) %>%        
      group_by(target, Concentration)
    minmax <-
      adply(input$filterTargets, 1,
            function(trgt) {
              if(!is.na(values$calib[[trgt]]) && !is.null(values$calib[[trgt]]$intercept)) {
                cqvals <- filter(values$tblCq,
                                 target == trgt,
                                 !is.na(Cq))$Cq
                data.frame(target = rep(trgt,2),
                           Cq  = c(
                             min(cqvals) - (max(cqvals) - min(cqvals)) * .05,
                             max(cqvals) + (max(cqvals) - min(cqvals)) * .05),
                           quant = c(
                             ((min(cqvals) - (max(cqvals) - min(cqvals)) * .05 - values$calib[[trgt]]$intercept)/values$calib[[trgt]]$slope),
                             ((max(cqvals) + (max(cqvals) - min(cqvals)) * .05 - values$calib[[trgt]]$intercept)/values$calib[[trgt]]$slope)
                           )
                )
              }
            }) %>%
      group_by(target)
    all_values <- function(x) {
      if(is.null(x)) return(NULL)
      paste0(names(x), ": ", format(x), collapse = "<br />")
    }
    if(is.null(input$selectedRows)) {
      filteredData <- values$tblCq
    } else {      
      l <- length(input$selectedRows)
      fdata.names <- paste(
        input$selectedRows[seq(1, l, by = 10)],
        input$selectedRows[seq(2, l, by = 10)],
        input$selectedRows[seq(4, l, by = 10)],
        input$selectedRows[seq(3, l, by = 10)],
        sep = "_")
      filteredData <- 
        filter(values$tblCq,
               fdata.name %in% fdata.names)
    }
    effsd %>%     
      ggvis() %>%
      layer_points(x = ~Concentration, 
                   y = ~Location..Mean.,
                   fill = ~target) %>%
      layer_rects(width := 1,
                  x = ~Concentration,
                  y = ~Location..Mean. + Deviation..SD.,
                  y2 = ~Location..Mean. - Deviation..SD.,
                  fill = ~target,
                  stroke = ~target) %>%
      layer_paths(data = minmax,
                  x = ~quant,
                  y = ~Cq,
                  stroke = ~target,
                  strokeWidth := 2) %>%
      layer_points(data = filter(filteredData,
                                 target %in% input$filterTargets,
                                 !is.na(Cq),
                                 sample.type != "std"),
                   x = ~log10quantity,
                   y = ~Cq,
                   fill = ~target,
                   size := 125,
                   shape := "cross") %>%
      add_tooltip(all_values, "hover") %>%
      add_axis("x", title = "Log Quantity") %>% 
      hide_legend(scales = c("shape", "stroke")) %>% 
      add_legend(scales = c("fill"), title = " ") %>% 
      bind_shiny("calibrationPlot", "calibrationPlot_ui")
  })
  
  output$calibText <- renderUI({
    if(values$calib %>% length == 0)
      return(NULL)
    compact <- function(x) Filter(Negate(is.null), x)
    sapply(names(values$calib), function(trgt) {
      if(is.null(values$calib[[trgt]]$slope))
        return(NULL)
      sprintf("%s - Slope %.2f; Intersect %.2f; R<sup>2</sup> %.4f; Efficiency %.2f",
              trgt,
              values$calib[[trgt]]$slope,
              values$calib[[trgt]]$intercept,
              values$calib[[trgt]]$eff@correlation.test$estimate  %>% abs,
              values$calib[[trgt]]$efficiency)
    }) %>% compact %>% paste(collapse = "<br/>") %>% HTML
  })
  
  output$cqConcTbl<-
    renderDataTable({
      if(is.null(values$tblCq)
         || is.null(values$mergedRDML))
        return(NULL)
      cat("view tbl\n")
      isolate({
        values$cqTblOut <- values$tblCq[c("position",
                                          "sample", 
                                          "target",
                                          "target.dyeId",
                                          "sample.type",
                                          if(!all(is.na(values$tblCq$conditions)))
                                            c("conditions"),
                                          if(!all(is.na(values$tblCq$quantity)))
                                            c("quantity",
                                              "quantity_mean",
                                              "quantity_SD"),
                                          "Cq",
                                          "Cq_Mean",
                                          "Cq_SD"
        )]
        
        colnames(values$cqTblOut) <- c("Pos.",
                                       "Sample", 
                                       "Target",
                                       "Dye",
                                       "Type",
                                       if(!all(is.na(values$tblCq$conditions)))
                                         c("Condition"),
                                       if(!all(is.na(values$tblCq$quantity)))
                                         c("Quantity",
                                           "Quantity Mean",
                                           "Quantity SD"),
                                       "Cq",
                                       "Cq Mean",
                                       "Cq SD") 
        tbl <- values$cqTblOut %>% mutate_each(funs(smart.round.col))
        if(tail(colnames(tbl), 1) != "Cq SD") 
          tbl <- tbl[, -ncol(tbl)]
        tbl[is.na(tbl)] <- ""
        tbl
      })
    })
  
  # Relative quantification -------------------------------------------------
  
  output$refGenesUi <- renderUI({
    if (is.null(values$tblCq))
      return(NULL)
    cat("creating refGenesUi\n")
    refGenes <- c()
    for (target in values$mergedRDML$target) {
      if (target$type$value == "ref")
        refGenes <- c(refGenes, target$id$id)
    }
    values$tblCq %>% 
      filter(sample.type %in% c("unkn", "pos")) %>%
      ungroup %>% 
      select(target) %>% unlist %>%  
      unique %>%      
      selectInput("refGenes", "Reference Genes",
                  choices = .,
                  selected = refGenes,
                  multiple = TRUE)
    
  })
  
  output$runCalibUi <- renderUI({
    if(is.null(values$tblCq))
      return(NULL)
    cat("creating runCalibUi\n")
    values$tblCq %>% 
      filter(sample.type %in% c("unkn", "pos")) %>%
      ungroup %>% 
      select(sample) %>% unlist %>%  
      unique %>% c(NA) %>% 
      selectInput("runCalib", "Run Calibrator",
                  choices = .,
                  selected = NA,
                  multiple = FALSE)
    
  })
  
  output$studyCalibUi <- renderUI({
    if(is.null(values$tblCq))
      return(NULL)
    cat("creating studyCalibUi\n")
    unique(values$tblCq$conditions) %>%      
      selectInput("studyCalib", "Study Calibrator",
                  choices = .,
                  selected = NA,
                  multiple = FALSE)
    
  })
  
  output$showRelGenesUi <- renderUI({
    if(is.null(values$tblCq))
      return(NULL)
    cat("creating showRelGenesUi\n")
    unique(values$tblCq$target) %>%
      setdiff(input$refGenes) %>% 
      selectInput("showRelGenes", "Show Targets",
                  choices = .,
                  selected = .,
                  multiple = TRUE)
    
  })
  
  output$showRelSamplesUi <- renderUI({
    if(is.null(values$tblCq))
      return(NULL)
    cat("creating showRelSamplesUi\n")
    values$tblCq %>% 
      filter(sample.type %in% c("unkn", "pos")) %>% 
      .$sample %>% unique %>%
      selectInput("showRelSamples", "Show Samples",
                  choices = .,
                  selected = .,
                  multiple = TRUE)
    
  })
  
  output$showRelConditionsUi <- renderUI({
    if(is.null(values$tblCq))
      return(NULL)
    cat("creating showRelConditionsUi\n")
    unique(values$tblCq$conditions) %>%
      selectInput("showRelConditions", sprintf("Show Conditions"),
                  choices = .,
                  selected = .,
                  multiple = TRUE)
    
  })
  
  
  
  output$relQuantTbl <- renderDataTable({
    if(is.null(values$tblCq))
      return(NULL)
    cat("Calc relative\n")
    
    values$tblRelative <- calcRelative(values$tblCq,
                                       ref.genes = input$refGenes,
                                       run.calib = input$runCalib,
                                       study.calib.cond = input$studyCalib,
                                       calib = values$calib)
    if(is.null(values$tblRelative)) return(NULL)
    isolate({
      values$tblRelativeOut <- 
        values$tblRelative[, c("position",
                               "sample",
                               "conditions",
                               "target",
                               "sample.type",
                               "ratio",
                               "ratio.error",
                               "normalized.ratio",
                               "normalized.ratio.error",
                               "scaled.ratio",
                               "scaled.ratio.error")]
      colnames(values$tblRelativeOut) <- c("Pos.",
                                           "Sample", 
                                           "Condition", 
                                           "Target",
                                           "Type",
                                           "Ratio",
                                           "Ratio Error",
                                           "Norm. Ratio",
                                           "Norm. Ratio Error",
                                           "Scaled Ratio",
                                           "Scaled Ratio Error") 
      tbl <- values$tblRelativeOut %>% mutate_each(funs(smart.round.col))
      if(tail(colnames(tbl), 1) != "Scaled Ratio Error") 
        tbl <- tbl[, -ncol(tbl)]
      tbl[is.na(tbl)] <- ""
      tbl
    })
  })
  
  
  
  output$relativeQuantPlot <- renderPlot({
    if(is.null(values$tblRelative) ||
       length(input$splitOpts) != 3)
      return(NULL)
    #     split.opts <- c("Образец" = "sample",
    #                     "Мишень" = "target",
    #                     "Услови\u044F" = "conditions")
    
    splitNames <- c(sample = "Sample",
                    target = "Target",
                    conditions = "Condition")
    values$tblRelative %>%
      unite_("split", input$splitOpts[1:2], 
             sep = "\n", remove = FALSE) %>%
      mutate(conditions = 
               ifelse(is.na(conditions),
                      "NA",
                      conditions),
             conditions = as.factor(conditions),
             split = factor(split)) %>% 
      filter(target %in% input$showRelGenes,
             sample %in% input$showRelSamples,
             sample.type %in% c("unkn", "pos"),
             conditions %in% input$showRelConditions,
             !(sample == input$runCalib &
                 input$displayRatio %in% c("normalized.ratio_mean",
                                           "scaled.ratio_mean")),
             !(conditions == input$studyCalib &
                 input$displayRatio == "scaled.ratio_mean")) %>% 
      ggplot(aes_string(x = "split", y = input$displayRatio)) + 
      aes_string(fill = sprintf("as.factor(%s)" , input$splitOpts[3])) +
      geom_bar(stat = "identity", position = "dodge") +
      { if(input$showErrRel)
        geom_errorbar(aes_string(ymin = 
                                   sprintf("%s - %s_mean",
                                           input$displayRatio,
                                           gsub("_mean", ".error",
                                                input$displayRatio)),
                                 ymax = 
                                   sprintf("%s + %s_mean",
                                           input$displayRatio,
                                           gsub("_mean", ".error",
                                                input$displayRatio))),
                      position=position_dodge(width=0.9),
                      width=0.25,
                      size = 1) 
      } +
      { if(input$curveTypeRel == "log") 
        scale_y_continuous(
          trans = mylog10_trans( 
            from = values$tblRelative[[input$displayRatio]] %>%
              min(na.rm = TRUE) %>% log10 %>% "-"(.5))) } +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("Ratio") + 
      xlab(sprintf("%s, %s",
                   splitNames[input$splitOpts[1]],
                   splitNames[input$splitOpts[2]])) +
      guides(fill = guide_legend(title = splitNames[input$splitOpts[3]]))
  })
  
  # downloads ---------------------------------------------------------------
  
  output$downloadCq <- downloadHandler(
    filename = function() { paste0(values$expTexts[input$experimentsChecks[1]], '.xlsx') },
    content = function(file) { write.xlsx(values$cqTblOut, file) }
  )
  
  output$downloadRel <- downloadHandler(
    filename = function() { paste0(values$expTexts[input$experimentsChecks[1]], '.xlsx') },
    content = function(file) { write.xlsx(values$tblRelativeOut, file) }
  )
  
})



