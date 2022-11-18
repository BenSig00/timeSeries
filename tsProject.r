
library(igraph)
library(rmatio)
library(depmixS4)
library(ggplot2)
library(glue)

#Function creating formula for depmix
make_formula <- function(n_states, data_name = "celldata"){
  formula <- list()
  for(i in 1:n_states){
    formula[[i]] <- as.formula(paste(data_name,glue("[,{i}]~1")))
  }
  return(formula)
}
#Function creating list of distribution for depmix
make_fam_list <- function(n_states){
  fam_list <- list()
  for(i in 1:n_states){
    fam_list[[i]] <- poisson()
  }
  return(fam_list)
}

#Function visualizing the result
depmix_and_angle_data <- function(model, angle_data){
  
  #model: fitted depmix model
  #angle data: vector of head angles. Same dimension as model (states)
  
  states = model@posterior$state
  n_states = length(unique(states))
  t_end = length(states)
  
  if(t_end != length(angle_data)){print("Length of states does not match observed data.")}
  
  #Transition matrix
  P = t(matrix(getpars(model)[(n_states +1): (n_states + n_states^2 )], nrow = n_states))
  
  t = c(1:t_end) / t_end + 0.2
  df = data.frame(angle_data,t)
  
  #Polar plot
  
  gg = ggplot(df, aes(x = angle_data, y = t)) + labs(color="State") +
    geom_point(aes(color = as.factor(states))) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + xlab(" ") + ylab(" ") + coord_polar()
  
  diag(P) = 0
  g = graph_from_adjacency_matrix(P, mode = "directed", weighted = T)
  plot(g, layout = layout_with_fr(g), edge.width = .1, edge.arrow.size = 0.001, vertex.label = " ")
  
  return(gg)
}

# Load data
data = read.mat("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\Mouse_data_and_python\\Mouse28-140313_BS0150_HMMready.mat") #Change!
celldata = data.frame(t(data$celldata))
had = data$resampledAwakeHeadAngleData

#Reduce data
num_zero = rep(0, ncol(celldata))
for(i in 1:ncol(celldata)){  num_zero[i] = sum(celldata[,i] == 0) }
celldata_reduced <- celldata[, (num_zero<13000)]
thr = 100
col_sum = colSums(celldata)
celldata_r <- celldata[,(col_sum>thr)]
ncol(celldata_r)
n_neur <- ncol(celldata_r)

# Train model
states = 16
formula <- make_formula(n_neur, data_name = "celldata_r")
fam_list <- make_fam_list(n_neur)
mod <- depmix(formula, data = celldata_r, nstates = states, family = fam_list)
fm <- fit(mod, verbose = T)

# Visualize model compared do observed data
gg = depmix_and_angle_data(fm, had)
gg

#-----------------------------------------------

tm10 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\transition_matrix_10.csv")
states10 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\df_10.csv")

tm14 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\transition_matrix_14.csv")
states14 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\df_14.csv")

tm28 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\transition_matrix_28.csv")
states28 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\df_28.csv")

states16 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\df_16.csv")
tm16 = read.csv("C:\\Users\\oy019\\OneDrive\\Documents\\SKOLE\\Time Series\\Project\\transition_matrix_16.csv")

transmat = as.matrix(tm28)
P = transmat
diag(P) = 0
g = graph_from_adjacency_matrix(P, mode = "directed", weighted = T)
plot(g, layout = layout_with_fr(g), edge.width = .1, edge.arrow.size = 0.001, vertex.label = " ")


#----------------------------------------------
states = states28
#Polar plot

gg = ggplot(states, aes(x = angle_data, y = t)) + labs(color="State") +
  geom_point(aes(color = as.factor(states))) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlab(" ") + ylab(" ") + coord_polar()
gg
