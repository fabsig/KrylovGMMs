library(xtable)
library(ggplot2)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data_sets <- c("cars", "chicago_building_permits",
               "instEval", "amazon_employee_access", "KDDCup09_upselling")

###Table: data info#############################################################

#import results
info_data <- data.frame()
for(d in 1:length(data_sets)){
  my_data <- readRDS(paste0("./", data_sets[d],".rds"))
  if(d==1){
    info_data <- my_data$data_info
  }
  else{
    info_data <- rbind(info_data, my_data$data_info)
  }
}

print(xtable(info_data))

###Table: runtimes##############################################################

#import results
results <- data.frame()
c <- c(1,2,4)
for(d in 1:length(data_sets)){
  my_data <- readRDS(paste0("./", data_sets[d],".rds"))
  new_results <-  my_data$all_results[,c]
  new_results <- cbind(ds_name = c(data_sets[d], rep(NA, nrow(new_results)-1)), new_results)
  if(d==1){
    results <- new_results
  }
  else{
    results <- rbind(results, new_results)
  }
}

print(xtable(results))

###Plots########################################################################

#import results
results <- data.frame()
c <- c(1,2,4)
for(d in 1:length(data_sets)){
  my_data <- readRDS(paste0("./", data_sets[d],".rds"))
  new_results <-  my_data$all_results[,c]
  new_results$ds_name <- data_sets[d]
  if(d==1){
    results <- new_results
  }
  else{
    results <- rbind(results, new_results)
  }
}

results$method <-  as.factor(results$method)
results$ds_name <- as.factor(results$ds_name)

#Plot time for each data set
ggplot(results, aes(x=ds_name, y=time_estimation, group = method)) +
  geom_point(aes(color=method, shape=method), size=2.5) + 
  geom_line(aes(color=method), linewidth=1) + theme_bw() + 
  scale_y_log10(breaks=c(1,3,10,30,100,300,1000,3000,8000)) + theme(legend.position = "top", legend.title=element_blank()) + 
  xlab("Data set") + ylab("Time (s)") + scale_color_brewer(type = "qual", palette=6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Plot average relative difference to the fastest model for each model
results <- results[!(results$ds_name %in% c("KDDCup09_upselling")),]  #remove the data sets KDDCup09_upselling
results <- results %>% group_by(ds_name) %>% mutate(min_time = min(time_estimation), 
                                                    rel_time = (time_estimation - min_time)/min_time)
mean_rel_diff <- results %>% group_by(method) %>% summarize(mean_rel_time = mean(rel_time))

ggplot(mean_rel_diff, aes(x=method, y=mean_rel_time)) +
  geom_point(aes(color=method, shape=method), size=2.5) + theme_bw() + 
  theme(legend.position="none") + 
  scale_y_continuous(breaks=c(0,seq(10,100,by=10))) +
  xlab("") + ylab("Average relative difference") + scale_color_brewer(type = "qual", palette=6)
