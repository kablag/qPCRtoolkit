calcRelative <- function(cq.data,
                         ref.genes,
                         run.calib,
                         study.calib.cond,
                         calib) {
  if(is.null(ref.genes)) return(NULL)
  eff <- function(ref.gene) {
    if(!is.null(calib[[ref.gene]]$efficiency) &&
       !is.na(calib[[ref.gene]]$efficiency))
      calib[[ref.gene]]$efficiency
    else 2
  }
  
  CqRef <- function(smpl,
                    cond,
                    gene) {
    cq.data %>% filter(sample == smpl,
                       if(is.na(cond)) is.na(conditions)
                       else conditions == cond,
                       sample.type %in% c("unkn", "pos"),
                       target == gene) %>%  .[1, "Cq_Mean"]
  }
  
  CqVar <- function(smpl,
                    cond,
                    gene) {
    cq.data %>% filter(sample == smpl,
                       if(is.na(cond)) is.na(conditions)
                       else conditions == cond,
                       sample.type %in% c("unkn", "pos"),
                       target == gene) %>%  .[1, "Cq_Var"]
  }
  
  ratio.numerator <- function(smpl,
                              cond) {
    sapply(ref.genes, 
           function(ref.gene) 
             eff(ref.gene) ^ CqRef(smpl,
                                   cond,
                                   ref.gene)) %>%
    prod %>% '^'(1 / length(ref.genes))
  }
  
  norm.denominator <- function(smpl,
                              cond) {
    sapply(ref.genes, 
           function(ref.gene) 
             eff(ref.gene) ^ (CqRef(run.calib,
                                   NA,
                                   ref.gene) - 
             CqRef(smpl, cond, ref.gene))) %>%
      prod %>% '^'(1 / length(ref.genes))
  }
  
  scaled.denominator <- function(smpl,
                                 cond) {
    sapply(ref.genes, 
           function(ref.gene) 
             eff(ref.gene) ^ (CqRef(smpl,
                                    study.calib.cond,
                                    ref.gene) - 
                                CqRef(smpl, cond, ref.gene))) %>%
             prod %>% '^'(1 / length(ref.genes))
  }
  
  ratio.er.calc <- function(smpl,
                            cond,
                            trgt,
                            ratio) {
    
    sapply(c(ref.genes, trgt),
           function(gene) {
             (((ratio / 1) * log(eff(gene))) ^ 2) * CqVar(smpl,
                                                          cond,
                                                          gene)
           }) %>% unlist %>% sum %>% sqrt
  }
  
  cq.data <- cq.data %>%
    group_by(fdata.name) %>%
    mutate(
      ratio = {
        ifelse(target %in% ref.genes || !(sample.type %in% c("unkn", "pos")),
               NA,
               ratio.numerator(sample, conditions) / 
                 eff(target) ^ Cq) %>% as.numeric
      },
      ratio.error = {
        ifelse(is.na(ratio),
               NA,
               ratio.er.calc(sample, conditions, target, ratio)) %>%
                 as.numeric
      },
      normalized.ratio = {
        ifelse(run.calib %>% is.null ||
                 target %in% ref.genes ||
                 #sample.type == "unkn",
                 !(sample.type %in% c("unkn", "pos")),
               NA,
               eff(target) ^ (CqRef(run.calib, NA, target) - Cq) /
                 norm.denominator(sample, conditions)) %>% 
          as.numeric
      },
      normalized.ratio.error = {
        ifelse(is.na(normalized.ratio),
               NA,
               ratio.er.calc(sample, conditions, target, 
                             normalized.ratio)) %>%
          as.numeric
      },
      scaled.ratio = {
        ifelse(study.calib.cond %>% is.null ||
                 is.na(conditions) ||
                 target %in% ref.genes ||
                 #sample.type == "unkn",
                 !(sample.type %in% c("unkn", "pos")),
               NA,
               eff(target) ^ (CqRef(sample, study.calib.cond, target) - Cq) /
                 scaled.denominator(sample, conditions)) %>% 
          as.numeric
      },
      scaled.ratio.error = {
        ifelse(is.na(scaled.ratio),
               NA,
               ratio.er.calc(sample, conditions, target, 
                             scaled.ratio)) %>%
          as.numeric
      }
    ) %>%
#     mutate(nre = {
#       ifelse(run.calib %>% is.null || is.na(normalized.ratio),
#              NA,
#              ratio.er.calc(sample, conditions, target, 
#                            ratio)) %>%
#         as.numeric
#     }) %>% 
    group_by(sample,
             conditions,
             target) %>%
    mutate(ratio_mean = 
             as.numeric(
                 mean(ratio, na.rm=TRUE)),
           normalized.ratio_mean = 
             as.numeric(
               mean(normalized.ratio, na.rm=TRUE)),
           scaled.ratio_mean = 
             as.numeric(
               mean(scaled.ratio, na.rm=TRUE)),
           ratio.error_mean = 
             as.numeric(
               mean(ratio.error, na.rm=TRUE)),
           normalized.ratio.error_mean = 
             as.numeric(
               mean(normalized.ratio.error, na.rm=TRUE)),
           scaled.ratio.error_mean = 
             as.numeric(
               mean(scaled.ratio.error, na.rm=TRUE))
           )
    
  return(cq.data)
}