setwd("C:/Users/jwahnzavalet/20220406_biotin_labeled_quant_mzML_analysis")

library(stringr)
library(ggplot2)
library(tibble)
library(magrittr)
library(plyr)
library(purrr)
library(dplyr)
library(ggpubr)
library(ggrepel)

##### Reading data #####################################################################################
data = read.csv("data_clean.csv", header = T, row.names = 1)
#data = data[,-1]

conditions = c("G13D_39", "G13D_40", "G13D_BOTH", "G13D_CTRL",
               "PAR_39", "PAR_40", "PAR_BOTH", "PAR_CTRL",
               "WT_39", "WT_40", "WT_BOTH", "WT_CTRL")

##### Functions ####################################################################################
#Selects proteins in common between condition and the wt_ctrl
inter_cond = function(cond1, cond2){
  idx = which(data[,cond1] != 0) # & data[,cond1]<data[,cond2]
  intersec = data.frame(data[idx,cond1], data[idx, cond2])
  rownames(intersec) = rownames(data[idx,])
  colnames(intersec) = c(cond1, cond2)
  
  return(intersec)
}

#Quantile normalization of the data
quantile_normalization <- function(df){
  
  # Find rank of values in each column
  df_rank <- map_df(df,rank,ties.method="average")
  # Sort observations in each column from lowest to highest 
  df_sorted <- map_df(df,sort)
  # Find row mean on sorted columns
  df_mean <- rowMeans(df_sorted)
  
  # Function for substiting mean values according to rank 
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  # Replace value in each column with mean according to rank 
  df_final <- map_df(df_rank,index_to_mean, my_mean=df_mean)
  
  return(df_final)
}

#T-test
t_test = function(dt, cond1, cond2){
  control = cond1
  condition = cond2
  #Subset control group and convert to numeric
  x = dt[,control] %>% unlist %>% as.numeric()
  
  #Subset condition group and convert to numeric
  y = dt[,condition] %>% unlist %>% as.numeric()
  
  result = try(t.test(x,y), silent = T)
  
  
  if (is(result, "try-error")){
    return(tibble(p_val = 1))}
  
  else
    p_vals = tibble(p_val = result$p.value)
  return(p_vals*nrow(dt))
}

Fc_2 = function(df, cond1, cond2){
  control = cond1
  condition = cond2
  dt = df %>% 
    select(-c(p_val)) %>%   #select columns
    log2()                  #log data
  
  data = bind_cols(dt, df[,"p_val"])
  
  
  mean_control = rowMeans(data[,control])
  mean_condition = rowMeans(data[,condition]) 
  
  log_fc = mean_control - mean_condition
  log_pval = -1*log10(data$p_val)
  
  data = data %>%
    mutate(log_Fc = log_fc, log_pval = log_pval)
  
  return(data)
}

#Preping data
volcano_prep = function(cond1, cond2){
  
  #keeping only proteins present in both
  df = inter_cond(cond1, cond2)

  genes = rownames(df)
  
  #Quantile Normalization
  df = quantile_normalization(df)
  #Hypothesis testing
  df =  plyr::adply(df,.margins = 1, .fun = t_test,
                    cond1, cond2) %>% as_tibble()
  #Fold change
  df = Fc_2(df, cond1, cond2)
  df$Genes = genes
  
  #if ("-Inf" %in% df$log_Fc | "Inf" %in% df$log_Fc){
  #  df = df[-which(df$log_Fc == "-Inf" |df$log_Fc == "Inf"),]
  #}
  
  return(df)
}

#Graphing volcanos
volcano_plot = function(kras, database = "global"){
  conditions = c()
  
  
  for (fix in c("_39", "_40", "_BOTH")){
    conditions = c(conditions, paste(kras, fix, sep=""))
  }
  
  for (cond in conditions){
    if (database == "glut"){
      var = paste("glut_", kras, "_CTRLvs", cond, sep="")
    }
    else{
      var = paste(kras, "_CTRLvs", cond, sep="")
    }
    
    assign(cond,get(var) %>%
             # Add a threhold for significant observations
             mutate(threshold = if_else(log_Fc >= 2 & log_pval >= 1.3 |
                                          log_Fc <= -2 & log_pval >= 1.3,"A", "B")) %>%
             # Plot with points coloured according to the threshold
             ggplot(aes(log_Fc,log_pval, colour = threshold)) +
             geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
             # Add dotted lines to indicate the threshold, semi-transparent
             geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
             geom_vline(xintercept = 2, linetype = 2, alpha = 0.5) +
             geom_vline(xintercept = -2, linetype = 2, alpha = 0.5) +
             # Set the colour of the points
             scale_colour_manual(values = c("A"= "deeppink1", "B"= "black")) +
             xlab("log2 fold change") + ylab("-log10 p-value") + # Relabel the axes
             theme_minimal() + # Set the theme
             theme(legend.position="none") + # Hide the legend
             xlim(-5.0,5.0) +
             ggtitle(var))
    
    
  }
  
  
  ggarrange(get(conditions[1]), get(conditions[2]), get(conditions[3]) + rremove("x.text"), 
            labels = c("A", "B", "C"),
            ncol = 2, nrow = 2)
  
  
}

extract_prots_VP = function(kras, Fc="down", file_name){
  
  vars = c()
  file.create(file_name)
  
  for (fix in c("_39", "_40", "_BOTH")){
    vars = c(vars, paste(kras, fix, sep=""))
  }
  
  prots = c()
  for (cond in vars){
    var = paste(kras, "_CTRLvs", cond, sep="")
    cat(var,"\n", file=file_name,append=T, sep="\n")
    
    var = get(var)
    
    if (Fc == "up"){
      result = var$Genes[which(var$log_Fc > 2)]
      prots = c(prots,result)
      cat(result,"\n", file=file_name, append=T, sep="\n")
    }
    
    if (Fc == "down"){
      result = var$Genes[which(var$log_Fc < (-2))]
      prots = c(prots,result)
      cat(result,"\n", file=file_name, append=T, sep="\n")
    }
    
    
  }
  print("Proteins in common")
  print(table(prots)[as.logical(table(prots)==3)])
}

##### WT_CTRL vs everything #####################################################################################
### Creates data frames containing the intensity values of wt_ctrl and condition
for (cond in conditions){
  var = paste("WTCTRL_", cond, sep="")
  assign(var, volcano_prep("WT_CTRL",cond))
}


##### DMSO CONTROL vs everything #####################################################################################
### Creates data frames containing the intensity values of wt_ctrl and condition
for (cond in conditions){
  var = paste("DMSO_", "cond", sep="")
  assign(var, volcano_prep("CONTROL",conditions))
}

volcano_plot()

##### CTRL OF EACH KRAS vs conditions #####################################################################################
###### WT
for (cond in c("WT_39", "WT_40", "WT_BOTH")){
  assign(paste("WT_CTRLvs", cond, sep=""), volcano_prep("WT_CTRL", cond))
}

volcano_plot("WT")

extract_prots_VP("WT", "down", "extractedProts_WT_down.txt")
extract_prots_VP("WT", "up", "extractedProts_WT_up.txt")

###### G13D
for (cond in c("G13D_39", "G13D_40", "G13D_BOTH")){
  assign(paste("G13D_CTRLvs", cond, sep=""), volcano_prep("G13D_CTRL", cond))
}

volcano_plot("G13D")
extract_prots_VP("G13D", "down", "extractedProts_G13D_down.txt")
extract_prots_VP("G13D", "up", "extractedProts_G13D_up.txt")

###### PAR
for (cond in c("PAR_39", "PAR_40", "PAR_BOTH")){
  assign(paste("PAR_CTRLvs", cond, sep=""), volcano_prep("PAR_CTRL", cond))

}

volcano_plot("PAR")
extract_prots_VP("PAR", "down", "extractedProts_PAR_down.txt")
extract_prots_VP("PAR", "up", "extractedProts_PAR_up.txt")




#################### other plots ######################################################################

barplot(as.matrix(data[which(rownames(data) == "PRDX5_HUMAN"),]), las=2, cex.names=0.8, col = "grey", main = "PRDX5 Intensity")
data[which(rownames(data) == "LANC1_HUMAN"),]


####### Number of proteins detected
nbr_prots = matrix(data = 1:12, nrow = 1, ncol=12)
colnames(nbr_prots) = conditions

for (cond in conditions){
  nbr_prots[1,cond] = length(which(data[,cond] != 0))
}

barplot(nbr_prots, cex.names = 0.8, las=2, main = "Number of proteins detected per condition", 
        ylim=c(0,1000), col="grey")

###### Intensities per condition
plot(sort(data[com_prots_uni,"G13D_39"]),col="blue", type="l", 
     main = "Intensities of common proteins between G13D cell lines", ylab = "Intensity",
     ylim = c(0,8e+08), xlim = c(100,222))
lines(1:222, sort(data[com_prots_uni,"G13D_40"]), col="red")
lines(1:222, sort(data[com_prots_uni,"G13D_BOTH"]), col="green")
lines(1:222, sort(data[com_prots_uni,"G13D_CTRL"]), col="black")
legend(x = "topleft", legend = conditions[1:4], lty = 1, col = c("blue", "red", "green", "black"))

plot(sort(data[com_prots_uni,"WT_39"]),col="blue", type="l",
     main = "Intensities of common proteins between WT cell lines", ylab = "Intensity",
     ylim = c(0,7e+08), xlim = c(100,222))
lines(1:222, sort(data[com_prots_uni,"WT_40"]), col="red")
lines(1:222, sort(data[com_prots_uni,"WT_BOTH"]), col="green")
lines(1:222, sort(data[com_prots_uni,"WT_CTRL"]), col="black")
legend(x = "topleft", legend = conditions[9:12], lty = 1, col = c("blue", "red", "green", "black"))

plot(sort(data[com_prots_uni,"PAR_39"]),col="blue", type="l", 
     main = "Intensities of common proteins between PAR cell lines", ylab = "Intensity",
     ylim = c(0,3e+08), xlim = c(100,222))
lines(1:222, sort(data[com_prots_uni,"PAR_40"]), col="red")
lines(1:222, sort(data[com_prots_uni,"PAR_BOTH"]), col="green")
lines(1:222, sort(data[com_prots_uni,"PAR_CTRL"]), col="black")
legend(x = "topleft", legend = conditions[5:8], lty = 1, col = c("blue", "red", "green", "black"))
