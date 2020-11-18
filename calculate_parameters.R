
library(janitor)


n_boot = 10

get_sample_sizes <- function(dim.students, dim.questions) {
  
  #dim.students <- 525
  #dim.questions <- 939
  
  #dim.students <- 525
  #dim.questions <- 939
  
  student_question_ratios = c( 0.3, 0.4, 0.5, 1,2,4, 6, 8, 10,  15, 20, 25, 30, 100) #0.1, 0.2,   # Student/Questions
  
  sample_sizes <- round(dim.questions * student_question_ratios)
  allow_sizes <- sample_sizes <= dim.students & sample_sizes > 5
  
  sample_sizes <- sample_sizes[allow_sizes]
  student_question_ratio_all <- student_question_ratios[allow_sizes]
  
  sample_sizes <- c(sample_sizes, dim.students)
  #student_question_ratio_all <- c(student_question_ratio_all, round(dim.students/dim.questions,1))
  
  #student_question_ratios <- factor(c(as.character(student_question_ratio_all)),
  #                                  # So that missing levels are also included
  #                                  levels = union(as.character(student_question_ratios), round(dim.students/dim.questions,1)), 
  #                                  ordered = TRUE)
  
  return(sample_sizes)
}

### Bootstrapping (without replacement) 
get_parameter_estimates <- function(X.p,Q_reduced, group = "Group", Q_names) {
  
  df.cdm <- CDM::din(X.p, Q_reduced, progress=FALSE)
  
  fit <- IRT.modelfit(df.cdm)
  
  df.p2.t1 <- tibble("value" = df.cdm$slip$est, 
                     "key" = Q_names) %>% 
    spread(key = "key", value = "value") %>%
    select(!!!Q_names) %>% 
    mutate(parameter = "Slip", group = group, type = "item") %>% 
    gather(key = "entities", value = "value", -group, -parameter, -type)
  
  
  
  
  df.p2.t2 <- tibble("value" = df.cdm$guess$est, 
                     "key" = Q_names) %>% 
    spread(key = "key", value = "value") %>% 
    select(!!!Q_names) %>% 
    mutate(parameter = "Guess", group =group, type = "item")  %>% 
    gather(key = "entities", value = "value", -group, -parameter, -type)
  
  
  #Item Parameters
  df.p2.t <- rbind(df.p2.t1, df.p2.t2) 
  
  # Fit Parameters
  df.p2.fit1 <- data.frame(model_fit = fit$modelfit.stat['SRMSR',] 
  ) %>% mutate(parameter = "SRMSR", type = "fit_val", group = group) %>% 
    gather(key = "entities", value = "value", -group, -parameter, -type)
  
  df.p2.fit2 <- data.frame(
    model_fit = fit$modelfit.test[1,2]
    
  ) %>% mutate(parameter = "CHI_SQR", type = "fit_val", group = group) %>% 
    gather(key = "entities", value = "value", -group, -parameter, -type)
  
  
  df.p2.fit3 <- data.frame( model_fit = fit$modelfit.test[1,3] 
  ) %>% mutate(parameter = "CHI_SQR", type = "p_val", group = group) %>% 
    gather(key = "entities", value = "value", -group, -parameter, -type)
  
  
  
  df.p2.fit <- rbind(df.p2.fit1, df.p2.fit2, df.p2.fit3) 
  
  
  
  #"sqr" = round(dim(df.X.p1)[1] / dim(df.X.p1)[2] ,1))
  
  
  # Combining all results together
  df.sim <- rbind(df.p2.t, df.p2.fit)
  
  return (df.sim)
  
}

# This function would partition data into 2 sub-populations and calculation DINA for each sub-population
cdm_fn <- function(df, Q, sample_size = -1, i = 1:dim(df)[1], do_print = FALSE) {
  
  #print(dim(df))
  #browser()
  #df.t <- X
  #Q_reduced = Q
  
  X <- df[i,] #Generating bootstrapped data by using bootstrapped index provided by boot function
  X.s <- head(X, sample_size) # Taking only 1:sample_size data to create smaller sample from bootstrapped X.
  
  X.p1.s <- X.s %>% head(round(dim(X.s)[1]/2)) # Partitition sample into two sub-populations.
  X.p2.s <- X.s %>% tail(round(dim(X.s)[1]/2))
  
  X.p1.s <- X.p1.s %>% remove_empty(.,which = "cols")
  X.p2.s <- X.p2.s %>% remove_empty(.,which = "cols")
  
  # This code helps to remove the questions that has not been answered by both the groups.
  Q_dist <-  X.p1.s %>% gather() %>% distinct(key) %>% inner_join(
   ( X.p2.s %>% gather()%>% distinct(key)), by = "key"
  ) %>% distinct(key)
  
  Q_reduced <- Q %>% mutate(key = colnames(X)) %>% semi_join(Q_dist) %>% select(-key)
  Q_names <- Q %>% mutate(key = colnames(X)) %>% semi_join(Q_dist) %>% select(key)
  Q_names <- Q_names$key
  
  
  df.X.p1 <- X.p1.s %>% mutate(id = row_number()) %>% gather(key = "key", value = "value", -id) %>% 
    semi_join(Q_dist, by = "key") %>% spread(key = "key", value = "value") %>% select(-id) # To select items available in Q matrix, because some students in smaller sample may never have answered 
  
  df.X.p2 <- X.p2.s %>% mutate(id = row_number()) %>% gather(key = "key", value = "value", -id) %>% 
    semi_join(Q_dist, by = "key") %>% spread(key = "key", value = "value") %>% select(-id)
  
  if(do_print) {
    print(paste("sample p1:", paste0(dim(df.X.p1),collapse = ","), 
                "; sample p2:", paste0(dim(df.X.p2), collapse = ","), 
                "Q Reduced:", paste0(dim(Q_reduced), collapse = ",")))
    
  }  
  ####################################################################
  # Calculation Attempts per Question
  
  df.attempts <- df.X.p1 %>% 
    
    summarise_all(funs(sum((!is.na(.))))) %>% 
    gather(key="entities", value = "attempts") %>% 
    mutate(group = "Partition 1") %>% 
    
    union (
    
      df.X.p2 %>% 
        summarise_all(funs(sum((!is.na(.))))) %>% 
        gather(key="entities", value = "attempts") %>% 
        mutate(group = "Partition 2")   
  )
    
  

  ####################################################################
  # Estimating parameters for Partition 1
  
  df.p1 <- get_parameter_estimates(df.X.p1, Q_reduced, "Partition 1", Q_names)
  
  ####################################################################
  # Estimating parameters for Partition 2
 
  df.p2 <- get_parameter_estimates(df.X.p2, Q_reduced, "Partition 2", Q_names)
  
  df.data.sim <- rbind(df.p1, df.p2) %>% left_join(df.attempts)
  
  # Returning item and model fit statistics for both partitions.
  return(df.data.sim)

  
} 



# Need to repeat this block with difference student size in X. 
#sample_size = 0.1

get_boot_index <- function(X) {
  
  
  
  X.bt <- boot(data = X , statistic = function(X, i) return(i), R = n_boot, stype = "i") # R needs to be 2 atleast for below code to work correctly
  
  return (X.bt$t %>% as.data.frame())
}


# This function calculate non-parameter bootstrapped statistics for a given sample size
get_mean_sample_error <- function(X, Q, sample_size) {
  
    print( c("Starting ", n_boot, " bootstrap for sample size:", sample_size) )
  
    # Bootstrap indexes are created from the whole data. However, DINA parameters are estimated only for the given sample size
    X.bt <- get_boot_index(X) %>% as.data.frame()
    X.index <- X.bt[1,]  %>% gather() 
    X.index <- X.index$value
    
    
    df.data.sim = cdm_fn(X, Q, sample_size = sample_size, X.index, do_print = TRUE) %>% mutate(sim_no = 1)

    #browser()
    for (i_val in 2: nrow(X.bt)) {
    
      X.index <- X.bt[i_val,]  %>% gather() 
      X.index <- X.index$value
      
      
      df.data.sim.t <- cdm_fn(X, Q, sample_size = sample_size, X.index)  %>% mutate(sim_no = i_val)
      
      df.data.sim = df.data.sim %>% bind_rows(df.data.sim.t)
    
    } # Gathering data for all sample sizes

    
  #browser()
  df.data.sim.item <- df.data.sim %>% filter(type == "item") %>% 
    select(-type) %>% rename(questions = entities, item_parameters = value) %>% 
    mutate(sample_size = sample_size)
  
  # Calculating Stats like Mean and Sampling Error for each Item Parameter in both subpopulations in a given sample
  df.data <- df.data.sim.item %>%
    
    # e.g. group_by ("Partition 1", "Slip", Question 1"). this function is called for single sample size value only
    group_by(group, parameter, questions, sample_size) %>% 
    summarise( 
      total_count = n(),   # Simulation Count
      #na_count = sum(is.na(item_parameters)),
      
      sampling_mean = mean(item_parameters), 
      sampling_error = sqrt(sum((item_parameters - sampling_mean)^2)/(n() - 1)), 
      sampling_variance = sampling_error^2,

      attempts = mean(attempts) # To average out attempts across simulations.
      
    ) %>% 
    
    
    ungroup()  
  
  
  # Calculating Stats like Mean and Sampling Error for Overall Item Parameters in both subpopulations in a given sample.
  #df.data.agg <- NA
  df.data.agg <- df.data %>% 
    
    group_by(group, parameter) %>% 
    summarise(sampling_error_mean = sqrt(sum(sampling_variance)/(n())),
              avg_attempts = round(mean(attempts))
              ) %>%
    
    ungroup() %>% 
    mutate(sample_size = sample_size) #dim(df.X.p1)[1])
  
  #browser()
  # Stats for Fit Parameters
  df.data.fit <- df.data.sim %>% filter(entities == "model_fit") %>% select(-entities) %>% rename(fit = type, fit_value = value) %>%
    group_by(group, parameter, fit) %>%
    summarise( 
     
      sampling_mean = mean(fit_value, na.rm = TRUE), 
      sampling_error = sqrt(sum((fit_value - sampling_mean)^2, na.rm = TRUE)/(n() - 1)), 
     
    ) %>% mutate(sample_size = sample_size) %>% inner_join(  
      df.data %>% distinct(group, attempts)
      )
  
  return(list(df.data,df.data.agg, df.data.fit, df.data.sim.item))
  
}
