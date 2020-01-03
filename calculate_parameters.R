
get_sample_sizes <- function(dim.students, dim.questions) {
  
  #dim.students <- 525
  #dim.questions <- 939
  
  student_question_ratios = c(0.1, 0.2, 0.3, 0.4, 0.5, 1,2,5, 10, 15)
  
  sample_sizes <- round(dim.questions * student_question_ratios)
  allow_sizes <- sample_sizes <= dim.students & sample_sizes > 5
  
  sample_sizes <- sample_sizes[allow_sizes]
  student_question_ratio_all <- student_question_ratios[allow_sizes]
  
  sample_sizes <- c(sample_sizes, dim.students)
  student_question_ratio_all <- c(student_question_ratio_all, round(dim.students/dim.questions,1))
  
  student_question_ratios <- factor(c(as.character(student_question_ratio_all)),
                                    # So that missing levels are also included
                                    levels = union(as.character(student_question_ratios), round(dim.students/dim.questions,1)), 
                                    ordered = TRUE)
  
  return(list(ratio = student_question_ratios, sample_size = sample_sizes))
}


total_cdm_fn <- function(df, i = 1:dim(df)[1]) {
  
  #print(dim(df))
  #browser()
  #df.t <- X
  df.t <- df[i,]
  
  df.t.s1 <- df.t #%>% select(E1:E28) 
  
  ###########################
  
  df.cdm <- CDM::din(df.t.s1, Q, progress=FALSE)
  
  
  df.t1 <- tibble("value" = df.cdm$slip$est, "key" = paste0("s1_slip_E", seq(1,length(df.cdm$slip$est),1))) %>% 
    spread(key = "key", value = "value") %>% 
    select(!!!paste0("s1_slip_E", seq(1,length(df.cdm$slip$est),1)))
  
  
  df.t2 <- tibble("value" = df.cdm$guess$est, "key" = paste0("s1_guess_E", seq(1,length(df.cdm$guess$est),1))) %>% 
    spread(key = "key", value = "value") %>% 
    select(!!!paste0("s1_guess_E", seq(1,length(df.cdm$guess$est),1)))
  
  
  df.t <- cbind(df.t1, df.t2)
  
  return(as.matrix(df.t))
  
} 



# Need to repeat this block with difference student size in X. 
sample_size = 0.1


get_boots_1 <- function(X.p1, X.p2, sample_size) {
  df.X.p1 <- X.p1 %>% union(X.p2) %>% sample_n(sample_size)
  
  print(paste0("All of X:", dim(df.X.p1)[1]))
  
  print("Starting Boot 1")
  X.bt.1 <- boot(data = df.X.p1 , statistic = total_cdm_fn, R = 100, stype = "i") # R needs to be 2 atleast for below code to work correctly
  
  #print("Starting Boot 2")
  #X.bt.2 <- boot(data = df.X.p2 , statistic = total_cdm_fn, R = 10, stype = "i")
  
  question_size = dim(X)[2]
  
  df.s1_slip <- X.bt.1$t[,1:question_size] %>%  as_tibble() %>% mutate(parameter = "Slip", group = "All Data")
  df.s1_guess <- X.bt.1$t[,(question_size+1):ncol(X.bt.1$t)] %>% as_tibble() %>% mutate(parameter = "Guess", group = "All Data")
  
  
  #df.s2_slip <- X.bt.2$t[,1:question_size] %>%  as_tibble() %>% mutate(parameter = "Slip", group = "Partition 2")
  #df.s2_guess <- X.bt.2$t[,(question_size+1):ncol(X.bt.1$t)] %>% as_tibble() %>% mutate(parameter = "Guess", group = "Partition 2")
  
  
  df.data.sim <- df.s1_slip %>% bind_rows(df.s1_guess) #%>% bind_rows(df.s2_slip) %>% bind_rows(df.s2_guess) 
  
}

get_boots_2 <- function(X.p1, X.p2, sample_size) {
  df.X.p1 <- X.p1 %>% sample_n(sample_size)
  df.X.p2 <- X.p2 %>% sample_n(sample_size)
  
  print(paste0("sample p1:", dim(df.X.p1)[1], "; sample p2:", dim(df.X.p2)[1]))
  
  print("Starting Boot 1")
  X.bt.1 <- boot(data = df.X.p1 , statistic = total_cdm_fn, R = 100, stype = "i") # R needs to be 2 atleast for below code to work correctly
  
  print("Starting Boot 2")
  X.bt.2 <- boot(data = df.X.p2 , statistic = total_cdm_fn, R = 100, stype = "i")
  
  question_size = dim(X)[2]
  
  df.s1_slip <- X.bt.1$t[,1:question_size] %>%  as_tibble() %>% mutate(parameter = "Slip", group = "Partition 1")
  df.s1_guess <- X.bt.1$t[,(question_size+1):ncol(X.bt.1$t)] %>% as_tibble() %>% mutate(parameter = "Guess", group = "Partition 1")
  
  
  df.s2_slip <- X.bt.2$t[,1:question_size] %>%  as_tibble() %>% mutate(parameter = "Slip", group = "Partition 2")
  df.s2_guess <- X.bt.2$t[,(question_size+1):ncol(X.bt.1$t)] %>% as_tibble() %>% mutate(parameter = "Guess", group = "Partition 2")
  
  
  df.data.sim <- df.s1_slip %>% bind_rows(df.s1_guess) %>% bind_rows(df.s2_slip) %>% bind_rows(df.s2_guess) 
  
}

get_mean_sample_error <- function(X.p1, X.p2, sample_size, sqr, n_boot) {
  
  if(n_boot == 1) {
    df.data.sim <- get_boots_1(X.p1, X.p2, sample_size)
  } else {
    df.data.sim <- get_boots_2(X.p1, X.p2, sample_size)
  }
  
  
  col_names <- colnames(df.data.sim %>% select(-group, -parameter) ) 
  
  df.data <- df.data.sim %>% gather(key = "questions", value="item_parameters", -group, -parameter) %>% 
    
    mutate(questions = factor(questions, 
                              labels = colnames(X),
                              levels = col_names, # All unique values of questions header return by boot function v1:v_
                              ordered = TRUE)) %>% 
    
    group_by(group, parameter, questions) %>% 
    summarise( 
      total_count = n(),   
      na_count = sum(is.na(item_parameters)),
      sampling_mean = mean(item_parameters), 
      sampling_error = sqrt(sum((item_parameters - sampling_mean)^2)/(n() - 1)), 
      sampling_variance = sampling_error^2
    ) %>% 
    
    ungroup() %>% 
    mutate(sample_size = sample_size, ratio = sqr) #dim(df.X.p1)[1])
  
  df.data.agg <- df.data.sim %>% gather(key = "questions", value="item_parameters", -group, -parameter) %>% 
    
    mutate(questions = factor(questions, 
                              labels = colnames(X),
                              levels = col_names, # All unique values of questions header return by boot function v1:v_
                              ordered = TRUE)) %>% 
    
    group_by(group, parameter, questions) %>% 
    summarise( sampling_mean = mean(item_parameters), 
               sampling_error = sqrt(sum((item_parameters - sampling_mean)^2)/(n() - 1)), 
               sampling_variance = sampling_error^2
    ) %>% 
    
    group_by(group, parameter) %>% 
    summarise(sampling_error_mean = sqrt(sum(sampling_variance)/(n()-1))) %>%
    
    ungroup() %>% 
    mutate(sample_size = sample_size, ratio = sqr) #dim(df.X.p1)[1])
  
  
  return(list(df.data,df.data.agg))
  
}
