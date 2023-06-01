# Recreating the grid plot
rm(list=ls())
library(tidyverse)
grid_knots <- 1:20

plotting <- data.frame(n_rep = NA,
           rmse = NA,
           crps = NA,
           pi = NA,
           nIknots = NA,
           intercept = NA)

plotting <- plotting[-1,]

for(i in 1:length(grid_knots)){

     # Merging the data
     load_obj <- readRDS(paste0("~/Documents/adding_new_folder_file/backupsplines/sbart2/R/results/faithful/faithful_spatial_N_100_tuning_n_tree_10_nIknots_",grid_knots[i],".Rds"))
     plotting <- plotting %>% rbind( cbind(load_obj$metrics_df, nIknots = grid_knots[i],intercept = "BASE" ))

     # Merging the data
     # load_obj <- readRDS(paste0("~/Documents/adding_new_folder_file/backupsplines/sbart2/R/results/N_100_tuning_n_tree_10_nIknots_",grid_knots[i],".Rds"))
     # plotting <- plotting %>% rbind( cbind(load_obj$metrics_df %>% filter(model == "SBART"), nIknots = grid_knots[i], intercept = "BS" ))

}

plotting %>% mutate(nIknots = as.factor(nIknots)) %>%
ggplot()+
     geom_boxplot(mapping = aes(x = nIknots, y = rmse, col = model))+
     ggtitle("cars")+
     theme_bw()


# plotting %>% mutate(nIknots = as.factor(nIknots)) %>%
#      group_by(nIknots) %>% summarise(m = mean(rmse)) %>%
#      ggplot()+
#      geom_boxplot(mapping = aes(x = nIknots, y = m, col = intercept))

plotting %>% group_by(nIknots) %>% summarise( m = mean(rmse)) %>% arrange(m) %>% head(10)
plotting %>% group_by(nIknots) %>% summarise( m = median(rmse)) %>% arrange(m) %>% head(10)
