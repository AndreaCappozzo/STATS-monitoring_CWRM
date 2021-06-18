library(flexCWM)
library(here)
library(tidyverse)
library(patchwork)
library(parallel)
library(foreach)
library(plotly)
library(patchwork)
library(scatterplot3d)
library(GGally)

source(file = "code/00-sim_CWM.R")
# Source trim.CWM
source(here("code/robustCWM.R"))

######################################################################
# 
# penalty to be employed within tBIC
# 
######################################################################
# df must be 4 columns: the first three the variables and the fourth is pop aka the classification

ggplot_solutions_CWRM <- function(df) {
  dens_ggpairs <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_density(..., alpha = 0.5)
  }
  
  points_ggpairs <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(..., alpha = 0.5)
  }
  
  legend_plot <-
    ggplot(data = mutate(df, pop = factor(ifelse(pop == 0, "out", pop))), aes(x =
                                                                                X1, y = y, col = pop)) +
    geom_point() +
    scale_color_manual(
      values = c(
        "out" = "black",
        "1" = palette_plotly[1],
        "2" = palette_plotly[2],
        "3" = palette_plotly[3],
        "4" = palette_plotly[4],
        "5" = palette_plotly[5]
      )
    ) + theme_bw() +
    labs(color = "") +
    theme(legend.position = "bottom")
  
  gg_pairs_plot <- ggpairs(
    data = df[, -4],
    aes(col = factor(df$pop), pch=ifelse(df$pop==0, 4,19)),
    upper = list(continuous = points_ggpairs),
    lower = list(continuous = points_ggpairs),
    diag = list(continuous = dens_ggpairs)
    # legend = grab_legend(legend_plot)
  ) +
    scale_color_manual(
      values = c(
        "0" = "black",
        "1" = palette_plotly[1],
        "2" = palette_plotly[2],
        "3" = palette_plotly[3],
        "4" = palette_plotly[4],
        "5" = palette_plotly[5]
      )
    ) +
    scale_fill_manual(
      values = c(
        "0" = "black",
        "1" = palette_plotly[1],
        "2" = palette_plotly[2],
        "3" = palette_plotly[3],
        "4" = palette_plotly[4],
        "5" = palette_plotly[5]
      )
    ) +
    theme_bw() +
    scale_shape_identity()+
    theme(legend.position = "bottom")
  
  gg_pairs_plot[1, 1] <-
    df %>%
    mutate(pop = factor(ifelse(pop == 0, "out", pop))) %>%
    filter(pop != "out") %>%
    ggplot(aes(x = get(names(df[1])), fill = pop)) +
    geom_density(alpha = .3, color = NA) +
    scale_fill_manual(
      values = c(
        "0" = "black",
        "1" = palette_plotly[1],
        "2" = palette_plotly[2],
        "3" = palette_plotly[3],
        "4" = palette_plotly[4],
        "5" = palette_plotly[5]
      )
    ) +
    theme_bw()
  gg_pairs_plot[3, 3] <-
    df %>%
    mutate(pop = factor(ifelse(pop == 0, "out", pop))) %>%
    filter(pop != "out") %>%
    ggplot(aes(x = get(names(df[3])), fill = pop)) +
    geom_density(alpha = .3, color = NA) +
    scale_fill_manual(
      values = c(
        "0" = "black",
        "1" = palette_plotly[1],
        "2" = palette_plotly[2],
        "3" = palette_plotly[3],
        "4" = palette_plotly[4],
        "5" = palette_plotly[5]
      )
    ) +
    theme_bw()
  
  gg_pairs_plot[2, 2] <-
    df %>%
    mutate(pop = factor(ifelse(pop == 0, "out", pop))) %>%
    filter(pop != "out") %>%
    ggplot(aes(x = get(names(df[2])), fill = pop)) +
    geom_density(alpha = .3, color = NA) +
    scale_fill_manual(
      values = c(
        "0" = "black",
        "1" = palette_plotly[1],
        "2" = palette_plotly[2],
        "3" = palette_plotly[3],
        "4" = palette_plotly[4],
        "5" = palette_plotly[5]
      )
    ) +
    theme_bw()
  gg_pairs_plot
}

pnlt_term<-function(c1,c2,k,p){
  fixed_term<-(k-1) +k*p +k*(p+1)
  #(k-1) mixture weights, k*p means in X,k*(p+1) beta coeff for regression bgx+b0g,
  var_term <- 1 + ( (k*p-1) + k*p*(p-1)/2 ) * (1-1/c1)
  # this is the part relative to modelling X
  # 1+ (k*p-1)  1 free eigenvalue and k*p-1 constrained eigenvalues,   
  # k*p*(p-1)/2 rotation matrices for Sigma_g,
  # (1-1/c1) take into account constrained eigenvalues 
  var_term <- var_term + 1 + (k-1) * (1-1/c2)
  # this is the part relative to modelling Y|X
  # 1+ (k-1)*(1-1/c2) one free σ^2g and k-1 constrained σ^2g  
  return(fixed_term+var_term)
}



###########################################################################
# Model validation
# based on the trimmed likelihood
# OR
# on the 
# to validate the obtained clustering
# Ingrassia punzo 2019
# Cluster validation for Mixtures of 
# Regressions via the tot sum of squares
############################################################################
var_dec<-function(a,Y,X){
  # in a$czise we have  sumCols(z_ij)
  rob_mean_Y<-mean(Y*rowSums(a$z_ij))*length(Y)/length(which(rowSums(a$z_ij)==1))
  #  iter$center[k,] = (t(a$z_ij[,k]) %*% Y) / iter$csize[k]
  means_y_g = (t(a$z_ij) %*% Y) / a$csize
  #  means_y_g <- a$b[,1] + a$center * a$b[,2]
  K<-length(a$cw)
  sc_mu_j_da_mu<-as.matrix(means_y_g-rep(rob_mean_Y,K),dim=c(1,K))
  BSS<- sum((sc_mu_j_da_mu)^2 * a$csize)
  #  X1<-cbind(rep(1,length(X)),X)
  #  explained<-X1 %*% t(a$b)- rep(1,length(X)) %*% t(means_y_g)
  #  explained^2
  one=matrix(1,nrow=length(Y),ncol=1)
  X1=cbind(one,X)
  #\mu(x_i;\hat(beta)_g) in formula 30
  # beta_0g+beta_1g'x
  forESS<-one%*%t(means_y_g)-X1%*%t(a$b)
  EWSS<-sum(a$z_ij*(forESS^2))
  #primo modo per calcolare gli errori
  onek<-rep(1,K)
  err=Y%*%t(onek)-X1%*%t(a$b)
  RWSS<-sum(a$z_ij*(err^2))
  # secondo modo per calcolare gli errori
  residual<-Y %*% t(rep(1,K))-X1 %*% t(a$b)
  # sono uguali gli errori=scarti residui?
  # sum(abs(residual-err))
  # cbind(Y,X1 %*% t(a$b),a$z_ij)[1:5,]
  scartirob<-(Y-rob_mean_Y)^2*rowSums(a$z_ij)
  # cbind((Y-mean(Y))^2,rowSums(a$z_ij))
  TSS<-sum(scartirob)
  devtot<-var(Y)*(length(Y)-1)
  # devtot
  # TSS
  # BSS  # the soft between sum of squares 
  # # (variability of Y explained by the latent group variable)
  # # BSS can be seen as a separation measure along the Y axis
  # EWSS  # the soft within-group sum of squares explained 
  # # by the model (thanks to the covariates)
  # RWSS  # the soft residual within-group sum of squares can 
  # # WSS=EWSS+RWSS can be seen as a compactness measure
  # BSS+EWSS+RWSS
  # TSS-(BSS+EWSS+RWSS)
  return(c(BSS,EWSS,RWSS))
}

outliers_adder <- function(n_out){
  p <- ncol(df)
  out <- matrix(0, nrow = n_out, ncol = p)
  iter <- 1
  
  while(iter<=n_out){
    
    out_cand <-
      sapply(1:p, function(dim)
        runif(1, min = out_range[1, dim], max = out_range[2, dim]))
    
    is_out <- dCWM_KL(yX = out_cand,
                      tau = tau,
                      beta = beta,
                      mu = mu,
                      Sigma = Sigma,
                      sigma = sigma2,
                      log = TRUE
    ) < min_log_d
    
    if(is_out){
      out[iter,] <- out_cand
      iter <- iter+1
      
    }
  }
  out
}

#  Parameters setting ----------------------------------------------------

p_x <- 2
K <- 4
n <- 980
tau <- c(.2,.4,.35,.05)
n_g <- tau*n

beta = matrix(c(10, 3,4,20,6,7,20,6,7,40,-6,-7), nrow = p_x+1, ncol = K)
mu = matrix(c(2,2,2,-2,1,-1,3,3), nrow = p_x, ncol = K)
Sigma = array(c(1,0,0,1,2,.5,.5,1,.5,0,0,.5,3,.5,.5,2), dim = c(p_x, p_x, K),)

eigen_Sigma <- apply(Sigma,3,function(x)eigen(x,only.values = TRUE)$values)
c_X_true <- max(eigen_Sigma)/min(eigen_Sigma)

# sigma = c(1:3,1)
sigma2 = c(1:3*5,1)
c_y_true <- max(sigma2)/min(sigma2)

# Data sim

set.seed(44)
df <- rCWM(n = n,tau = tau,beta = beta,mu = mu,Sigma = Sigma,sigma = sigma2, ng_fixed = TRUE)
class <- attributes(df)$clust
pairs(df, col=class)

min_log_d <- min(dCWM(
  y = df$y,
  X = as.matrix(df[, 2:3]),
  tau = tau,
  beta = beta,
  mu = mu,
  Sigma = Sigma,
  sigma = sigma2,log = TRUE
))

out_range <- apply(df, 2, range)
n_outliers <- 20
outliers <- outliers_adder(n_out = n_outliers)

df_out <- rbind(as.matrix(df), outliers)
class_out <- c(class,rep(0,n_outliers))

pairs(df_out)
pairs(df_out, col=class_out+1)

plotyy_df <- data.frame(cbind(df_out, pop=class_out)) %>% 
  as_tibble() 

palette_plotly <- RColorBrewer::brewer.pal(n = 8, name = "Set2")

# With plotly
scatter_3d <- plot_ly(
  data = plotyy_df,
  x =  ~ X1,
  y =  ~ X2,
  z =  ~ y,
  color =  ~ pop,
  colors = c("black", palette_plotly[1], palette_plotly[2],palette_plotly[3],palette_plotly[4]),
  size=.2
) 

scatter_3d
# # with 3dplot
# 
scatterplot3d(
  x=plotyy_df$X1,
  y=plotyy_df$X2,
  z=plotyy_df$y,
  # plotyy_df[,-4],
  pch = 16,
  grid = TRUE,
  box = FALSE,
  color = rep(c(
    palette_plotly[1],
    palette_plotly[2],
    palette_plotly[3],
    palette_plotly[4],
    "black"
  ), table(plotyy_df$pop)[c(2:5,1)]),
  angle=30
)
# 
# scatter_3d

# Pairs plot

dens_ggpairs <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_density(..., alpha = 0.5)
}

points_ggpairs <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_point(..., alpha = 0.5) 
}

legend_plot <-
  ggplot(data = mutate(plotyy_df, pop = factor(ifelse(pop == 0, "out", pop))), aes(x =X1, y = y, col = pop)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "out" = "black",
      "1" = palette_plotly[1],
      "2" = palette_plotly[2],
      "3" = palette_plotly[3],
      "4" = palette_plotly[4]
    )
  ) + theme_bw() +
  labs(color="") +
  theme(legend.position = "bottom")

gg_pairs_plot <- ggpairs(data = plotyy_df[,-4],aes(col=factor(plotyy_df$pop)),
                         upper = list(continuous=points_ggpairs), 
                         lower = list(continuous=points_ggpairs),diag = list(continuous=dens_ggpairs),
                         legend = grab_legend(legend_plot)) +
  scale_color_manual(
    values = c(
      "0" = "black",
      "1" = palette_plotly[1],
      "2" = palette_plotly[2],
      "3" = palette_plotly[3],
      "4"=palette_plotly[4]
    )
  ) +
  scale_fill_manual(
    values = c(
      "0" = "black",
      "1" = palette_plotly[1],
      "2" = palette_plotly[2],
      "3" = palette_plotly[3],
      "4"=palette_plotly[4]
    )
  ) +
  theme_bw() +
  theme(legend.position = "bottom") 

gg_pairs_plot[1, 1] <-
  plotyy_df %>%
  mutate(pop = factor(ifelse(pop == 0, "out", pop))) %>%
  filter(pop != "out") %>%
  ggplot(aes(x = get(names(plotyy_df[1])), fill = pop)) +
  geom_density(alpha = .3, color = NA) +
  scale_fill_manual(
    values = c(
      "0" = "black",
      "1" = palette_plotly[1],
      "2" = palette_plotly[2],
      "3" = palette_plotly[3],
      "4" = palette_plotly[4]
    )
  ) +
  theme_bw() 
gg_pairs_plot[3, 3] <-
  plotyy_df %>%
  mutate(pop = factor(ifelse(pop == 0, "out", pop))) %>%
  filter(pop != "out") %>%
  ggplot(aes(x = get(names(plotyy_df[3])), fill = pop)) +
  geom_density(alpha = .3, color = NA) +
  scale_fill_manual(
    values = c(
      "0" = "black",
      "1" = palette_plotly[1],
      "2" = palette_plotly[2],
      "3" = palette_plotly[3],
      "4" = palette_plotly[4]
    )
  ) +
  theme_bw() 

gg_pairs_plot[2, 2] <-
  plotyy_df %>%
  mutate(pop = factor(ifelse(pop == 0, "out", pop))) %>%
  filter(pop != "out") %>%
  ggplot(aes(x = get(names(plotyy_df[2])), fill = pop)) +
  geom_density(alpha = .3, color = NA) +
  scale_fill_manual(
    values = c(
      "0" = "black",
      "1" = palette_plotly[1],
      "2" = palette_plotly[2],
      "3" = palette_plotly[3],
      "4" = palette_plotly[4]
    )
  ) +
  theme_bw() 

gg_pairs_plot

ggsave(
  filename = "output/sim_pairs.eps",
  device = cairo_ps,
  width = 15.3*.4,
  height = 10.5*.4
)

c1<-c2<-c(1,2,4,16,32,64)
k=1:5

N <- nrow(df_out)

alpha <- seq(0,0.10,.002)
true_alpha_no_out <- n_outliers/N
true_alpha_no_forth_group <- (n_outliers+n_g[4])/N

parameters_grid <- expand.grid(num_groups=k,restr_factor_X=c1,restr_factor_Y=c2,trim_level=alpha)

# NOT RUN: it takes a while. Uncomment (and prepare to wait) to generate the file STATS_monitoring_CWRM_STRONG_overlapping.Rds
# clustvarsel::startParallel(10)
# i <- NULL
# 
# out <- foreach(
#   i = 1:nrow(parameters_grid)
# ) %dopar% {
#   
#   trim.cwm(X = df_out[,2:3],Y = df_out[,1],
#            K = parameters_grid[i,"num_groups"],
#            niter = 100,
#            Ksteps = 50,
#            maxfactx = parameters_grid[i,"restr_factor_X"],
#            maxfacty = parameters_grid[i,"restr_factor_Y"],
#            zero.tol = 1e-16,
#            zero.tol1 = 1e-16,
#            zero.tol2 = 1e-90,
#            trace = 0,
#            alpha=parameters_grid[i,"trim_level"],
#            init = 3
#   )
# }
# 
# saveRDS(out,file = "code/results/STATS_monitoring_CWRM_STRONG_overlapping.Rds")

# *.Rds files not tracked in the repo
models <- readRDS(file = "code/results/STATS_monitoring_CWRM_STRONG_overlapping.Rds")

model_results <- tibble(parameters_grid) %>%
  mutate(model = models) %>%
  rowwise() %>%
  mutate(
    TBIC = -2 * model$obj +
      pnlt_term(restr_factor_X, restr_factor_Y, num_groups, p = p_x) *
      log(N * (1 - trim_level)),
    decomp = list(var_dec(
      a = model, Y = df_out[, 1], X = df_out[, 2:3]
    )),
    BSS = decomp[1],
    EWSS = decomp[2],
    RWSS = decomp[3],
    R = EWSS / (EWSS + RWSS)
  ) %>%
  select(-decomp) %>%
  ungroup()


# First step: detecting the best alpha/alphas -----------------------------

best_alpha_wise <- model_results %>% 
  group_by(trim_level) %>% 
  filter(TBIC==min(TBIC)) %>% 
  ungroup() %>% 
  arrange(trim_level)

best_alpha_wise %>% 
  print(n=Inf)

# Model with correct level of trimming

correct_trim_level_model <- best_alpha_wise %>% 
  filter(trim_level==.02) %>% 
  pull(model) %>% 
  pluck(1)

table(correct_trim_level_model$assig, class_out)
# Highest trim level for finding deepest obs
most_robust_model <- best_alpha_wise$model[[nrow(best_alpha_wise)]]

# Manual relabeling of the most robust modeling parameters
most_robust_relabeler <- c(3,2,1)

most_robust_model$cw=most_robust_model$cw[most_robust_relabeler]
most_robust_model$b=most_robust_model$b[most_robust_relabeler,]
most_robust_model$center=most_robust_model$center[most_robust_relabeler,]
most_robust_model$sigma=most_robust_model$sigma[,,most_robust_relabeler]
most_robust_model$v=most_robust_model$v[most_robust_relabeler]

table(most_robust_model$assig, class_out)

# I relabel it according to table(most_robust_model$assig, class_out)



pairs(df_out, col=most_robust_model$assig+1)
table(most_robust_model$assig, class_out)

G_most_robust_modeling <- nrow(most_robust_model$center)

y_representative_fake_obs <- sapply(1:G_most_robust_modeling,function(g)
  most_robust_model$b[g,] %*% c(1,most_robust_model$center[g,]))
group_representative_fake_obs <- cbind(y_representative_fake_obs, most_robust_model$center)


z_CWM <-
  estep_CWM(
    y = group_representative_fake_obs[,1],
    X = group_representative_fake_obs[,-1,drop=FALSE],
    tau = most_robust_model$cw,
    beta = t(most_robust_model$b),
    mu = t(most_robust_model$center),
    Sigma = most_robust_model$sigma,
    sigma = most_robust_model$v
  )

map_CWM <- mclust::map(z_CWM)

n_models <- nrow(best_alpha_wise)

relabeler_list <- vector(mode = "list", length = n_models)
relabeler_list[[n_models]] <- 1:length(map_CWM)

relabeler_list[[n_models]] <- most_robust_relabeler

labels_relabeled_df <- map_dfc(.x = 1:n_models, .f = ~identity(best_alpha_wise$model[[.x]]$assig))

labels_to_be_changed <- labels_relabeled_df[[n_models]]

for(k in 1:length(relabeler_list[[n_models]])){
  labels_to_be_changed[labels_relabeled_df[[n_models]]==k] <- which(relabeler_list[[n_models]]==k)
}

labels_relabeled_df[[n_models]] <-   labels_to_be_changed

table(labels_relabeled_df[[n_models]], class_out)

# View(labels_relabeled_df)

for(n_mod in (n_models-1):1){
  
  cluster_relabeler_output <- cluster_relabeler_via_CWM_density(
    X = df_out,
    group_representative_fake_obs = group_representative_fake_obs,
    old_partition = labels_relabeled_df[[n_mod+1]],
    new_partition = labels_relabeled_df[[n_mod]],
    new_model = best_alpha_wise$model[[n_mod]],
    old_model =best_alpha_wise$model[[n_mod+1]]
  )
  # Save relabeler
  relabeler_list[[n_mod]] <- cluster_relabeler_output$relabeler
  
  # Relabel model n_mod
  labels_to_be_changed <- labels_relabeled_df[[n_mod]]
  for(k in 1:length(relabeler_list[[n_mod]])){
    labels_to_be_changed[labels_relabeled_df[[n_mod]]==k] <- which(relabeler_list[[n_mod]]==k)
  }
  labels_relabeled_df[[n_mod]] <-   labels_to_be_changed
  
  # Update groupwise fake obs
  group_representative_fake_obs <- cluster_relabeler_output$group_representative_fake_obs
}

# View(labels_relabeled_df)

# Plotting with relabeling ------------------------------------------------

theme_set(theme_bw())

# Mixing proportions

# gg_mix <- best_alpha_wise %>%
#   hoist(.col=model, "cw") %>% 
#   mutate(cw_relabeled= map(1:nrow(best_alpha_wise), ~cw[[.x]][relabeler_list[[.x]]])) %>% 
#   unnest_wider(cw_relabeled, names_sep = ".") %>%
#   select(trim_level, starts_with("cw_relabeled")) %>%
#   pivot_longer(-trim_level) %>%
#   mutate(name=str_extract(name,".$"), title="pi") %>% 
#   ggplot(aes(x=trim_level,y=value,fill=name)) +
#   geom_col() +
#   facet_wrap(~title, labeller = label_parsed) +
#   scale_x_continuous(breaks = alpha) +
#   labs(x="", y="", fill="G") +
#   scale_fill_brewer(palette = "Set2") +
#   theme(axis.text.x = element_blank(), legend.position = "bottom")
color_palette_set_2 <- RColorBrewer::brewer.pal(n=8,name = "Set2")
cluster_prop_colors <- c("0"="black",
                         "1"=color_palette_set_2[1],
                         "2"=color_palette_set_2[2],
                         "3"=color_palette_set_2[3],
                         "4"=color_palette_set_2[4],
                         "5"=color_palette_set_2[5],
                         "6"=color_palette_set_2[6],
                         "7"=color_palette_set_2[7],
                         "8"=color_palette_set_2[8])
# Cluster proportions including trimming level
gg_mix <- best_alpha_wise %>%
  hoist(.col=model, "csize") %>% 
  mutate(csize_relabeled= map(1:nrow(best_alpha_wise), ~csize[[.x]][relabeler_list[[.x]]]/N)) %>% 
  rowwise() %>% 
  mutate(csize_relabeled=list(c(trim_level,csize_relabeled))) %>% 
  ungroup()  %>% 
  unnest_wider(csize_relabeled, names_sep = ".") %>%
  select(trim_level, starts_with("csize_relabeled")) %>%
  pivot_longer(-trim_level) %>%
  mutate(name=factor(as.numeric(str_extract(name,".$"))-1), title="Groups proportion") %>% 
  ggplot(aes(x=trim_level,y=value,fill=name)) +
  geom_col() +
  facet_wrap(~title) +
  scale_x_continuous(breaks = alpha) +
  labs(x="", y="", fill="G") +
  scale_fill_manual(breaks = names(cluster_prop_colors[-1]), values = cluster_prop_colors, 
  ) +
  theme(axis.text.x = element_blank(), legend.position = "bottom")
gg_mix
best_alpha_wise$model[[1]]$csize

# Beta and sigma coefficients


beta1_df <- map(1:nrow(best_alpha_wise),~best_alpha_wise$model[[.x]]$b[,2]) %>% 
  enframe() %>% 
  mutate(value_relabeled=map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>% 
  unnest_wider(col = value_relabeled, names_sep = "Beta1") %>% 
  select(-value)

beta2_df <- map(1:nrow(best_alpha_wise),~best_alpha_wise$model[[.x]]$b[,3]) %>% 
  enframe() %>% 
  mutate(value_relabeled=map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>% 
  unnest_wider(col = value_relabeled, names_sep = "Beta2") %>% 
  select(-value)

sigma_regr_df <- map(1:nrow(best_alpha_wise),~sqrt(best_alpha_wise$model[[.x]]$v)) %>% 
  enframe() %>% 
  mutate(value_relabeled=map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>% 
  unnest_wider(col = value_relabeled, names_sep = "sigma_regr") %>% 
  select(-value)

det_sigma_df <- map(1:nrow(best_alpha_wise),~apply(best_alpha_wise$model[[.x]]$sigma,3,function(s) det(s)^(1/p_x))) %>% 
  enframe() %>% 
  mutate(value_relabeled=map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>% 
  unnest_wider(col = value_relabeled, names_sep = "det_sigma") %>% 
  select(-value)

param_df <- best_alpha_wise %>% 
  select(trim_level) %>% 
  # bind_cols(inner_join(inner_join(beta1_df,beta2_df,"name"),sigma_regr_df,"name")) %>% 
  bind_cols(plyr::join_all(list(beta1_df,beta2_df, sigma_regr_df, det_sigma_df),by="name",type="inner")) %>% 
  select(-name) %>% 
  pivot_longer(-trim_level) %>% 
  mutate(param=case_when(str_detect(name,pattern = "Beta1")~"b[g]^1",
                         str_detect(name,pattern = "Beta2") ~ "b[g]^2",
                         str_detect(name,pattern = "sigma_regr")~"sigma[g]",
                         str_detect(name,pattern = "det_sigma")~"abs(Sigma[g])^{1/d}"),
         group=str_extract(string = name,pattern = ".$")) 


gg_param <- param_df %>% 
  mutate(param=fct_relevel(.f = param,"abs(Sigma[g])^{1/d}",after = Inf)) %>% 
  ggplot(aes(trim_level,value, group=group, col=group)) +
  geom_point() +
  geom_line(aes(lty=group)) +
  # facet_wrap(~param,ncol = 1,labeller = label_parsed,scales = "free") +
  facet_grid(rows = vars(param),labeller = label_parsed,scales = "free") +
  scale_x_continuous(breaks = alpha) +
  labs(x=quote(alpha), y="", col="G", lty="G") +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.text.y = element_text(angle = 0))

gg_param

gg_beta_1 <- param_df %>% 
  filter(param=="b[g]^1") %>% 
  ggplot(aes(trim_level,value, group=group, col=group)) +
  geom_point(show.legend = FALSE) +
  geom_line(aes(lty=group),show.legend = FALSE) +
  facet_wrap(~param,ncol = 1,labeller = label_parsed,scales = "free") +
  scale_x_continuous(breaks = alpha) +
  scale_color_brewer(palette = "Set2") +
  labs(x="", y="", col="G", lty="G") +
  theme(axis.text.x = element_blank())

gg_beta_2 <- param_df %>% 
  filter(param=="b[g]^2") %>% 
  ggplot(aes(trim_level,value, group=group, col=group)) +
  geom_point(show.legend = FALSE) +
  geom_line(aes(lty=group),show.legend = FALSE) +
  facet_wrap(~param, labeller = label_parsed) +
  scale_x_continuous(breaks = alpha) +
  scale_color_brewer(palette = "Set2") +
  labs(x="", y="", col="G", lty="G") +
  theme(axis.text.x = element_blank())

gg_sigma <- param_df %>% 
  filter(param=="sigma[g]") %>% 
  ggplot(aes(trim_level,value, group=group, col=group)) +
  geom_point(show.legend = FALSE) +
  geom_line(aes(lty=group),show.legend = FALSE) +
  facet_wrap(~param, labeller = label_parsed) +
  scale_x_continuous(breaks = alpha) +
  scale_color_brewer(palette = "Set2") +
  labs(x="", y="", col="G", lty="G") +
  theme(axis.text.x = element_blank())

gg_det_Sigma <- param_df %>% 
  filter(param=="abs(Sigma[g])^{1/d}") %>% 
  ggplot(aes(trim_level,value, group=group, col=group)) +
  geom_point(show.legend = FALSE) +
  geom_line(aes(lty=group),show.legend = FALSE) +
  facet_wrap(~param, labeller = label_parsed) +
  scale_x_continuous(breaks = alpha) +
  scale_color_brewer(palette = "Set2") +
  labs(x="", y="", col="G", lty="G") +
  theme(axis.text.x = element_blank())

# Normalized Three-Term Decomposition \cite{Ingrassia2020}

gg_decomp <- best_alpha_wise %>% 
  select(trim_level,BSS,EWSS,RWSS) %>% 
  mutate(TSS=BSS+EWSS+RWSS,
         across(-trim_level,.fns = ~.x/TSS)) %>% 
  select(-TSS) %>% 
  rename_with(.fn = ~paste0("N",.x),-trim_level) %>% 
  pivot_longer(-trim_level) %>% 
  mutate(title="Total Sum of Squares Decomposition") %>% 
  ggplot() +
  geom_col(aes(trim_level,value, fill=name))+
  scale_x_continuous(breaks = alpha) +
  facet_wrap(~title) +
  labs(x="", y="", fill="TSS Decomposition") +
  scale_fill_manual(values = c("NBSS"="darkblue", "NRWSS"="darkred", "NEWSS"="darkgreen")) +
  theme(axis.text.x = element_blank())

gg_ARI <-
  tibble(
    trim_level = alpha[-length(alpha)],
    ARI = map_dbl(
      1:(nrow(best_alpha_wise) - 1),
      ~ mclust::adjustedRandIndex(best_alpha_wise$model[[.x]]$assig, best_alpha_wise$model[[.x +1]]$assig)
    ), title="Adjusted Rand Index"
  ) %>% 
  ggplot(aes(trim_level, ARI)) +
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = alpha) +
  facet_wrap(~title) +
  labs(x=quote(alpha), y="")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

gg_combined <- ((gg_mix / gg_decomp) /
                  gg_beta_1/gg_beta_2/gg_sigma/gg_det_Sigma/
                  gg_ARI) &
  theme(legend.position = "bottom")


gg_combined + plot_layout(guides = "collect")

ggsave(
  filename = "output/step_1_sim_monitoring_WITH_TRIMMING.eps",
  device = cairo_ps,width = 210, height = 280,units = "mm")

# Second step A: detecting the best solutions conditioning on alpha =.02-----------------------------

# Fixing alpha
solution_1_alpha_fixed <- model_results %>% 
  filter(trim_level==.02) %>% 
  select(1:3,5:6) %>% 
  ungroup() %>% 
  mutate(rank_solution = rank(TBIC, ties.method = "first"))

threshold_ARI <- .7

output_second_step_LIGHTER <-
  solutions_ranking_identifier_LIGHTER(
    solutions_df = solution_1_alpha_fixed,
    number_of_optimal_sol = 4,
    epsilon_ARI = threshold_ARI
  )


output_second_step_LIGHTER$optimal_solutions

output_second_step_LIGHTER$df_for_plotting %>% 
  filter(num_groups!=1) %>% 
  mutate(to_fill=case_when((is_stable & is_best_interval)~TRUE,
                           is_stable~FALSE,
                           TRUE~NA),
         n_sol=factor(ifelse(is.na(to_fill),NA, n_sol)),
         num_groups=paste("G = ", num_groups)) %>% 
  ggplot(aes(x=restr_factor_X,y=restr_factor_Y))+
  geom_tile(aes(fill=n_sol,alpha=(to_fill)|!is.na(rank_best)), colour = "black", na.rm = TRUE, size = 0.5)+
  # geom_tile(aes(fill=n_sol),colour = "black", na.rm = TRUE, size = 0.5)+
  geom_text(aes(label=rank_best)) +
  labs(x=quote(c[X]),y=quote(c[y]), fill="Optimal solutions")+
  theme_bw() +
  guides(alpha=FALSE) +
  facet_wrap (~ num_groups,nrow = 1) +
  scale_fill_brewer(palette="Dark2",na.value = "white", breaks=1:(n_distinct(output_second_step_LIGHTER$df_for_plotting$n_sol)-1)) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(
  filename = "output/step_2_monitoring_sim.eps",
  device = cairo_ps,width = 13.7, height = 4.36)

table(output_second_step$optimal_solutions$model[[1]]$assig, class_out)


pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[4]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[5]]$assig+1)

table(output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[4]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[5]]$assig, class_out)


# Manual relabeling of the classes different sol


first_cl <- output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig
second_cl <- output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig
third_cl <- output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig

table(first_cl, class_out)


first_true_relabeler <- c(3,2,4,1)
first_cl_relabeled <- first_cl

for(k in 1:length(first_true_relabeler)){
  first_cl_relabeled[first_cl==k] <- first_true_relabeler[k]
}

table(first_cl_relabeled, class_out)

table(second_cl, class_out)
second_true_relabeler <- c(1,4,3,2)
second_cl_relabeled <- second_cl

for(k in 1:length(second_true_relabeler)){
  second_cl_relabeled[second_cl==k] <- second_true_relabeler[k]
}

table(second_cl_relabeled, class_out)

table(third_cl, class_out)

third_true_relabeler <- c(2,1,4)

third_cl_relabeled <- third_cl

for(k in 1:length(third_true_relabeler)){
  third_cl_relabeled[third_cl==k] <- third_true_relabeler[k]
}

table(third_cl_relabeled, class_out)

df_1 <- mutate(plotyy_df, pop=first_cl_relabeled)
df_2 <- mutate(plotyy_df, pop=second_cl_relabeled)
df_3 <- mutate(plotyy_df, pop=third_cl_relabeled)

ggplot_solutions_CWRM(df = df_1)+
  ggtitle(expression(paste("First optimal solution (", alpha==0.02,")")))

ggsave(
  filename = "output/first_opt_sol.eps",
  device = cairo_ps,width = 6.94, height = 6.44)

ggplot_solutions_CWRM(df_2)+
  ggtitle(expression(paste("Second optimal solution (", alpha==0.02,")")))

ggsave(
  filename = "output/second_opt_sol.eps",
  device = cairo_ps,width = 6.94, height = 6.44
)

ggplot_solutions_CWRM(df_3)+
  ggtitle(expression(paste("Third optimal solution (", alpha==0.02,")")))

ggsave(
  filename = "output/third_opt_sol.eps",
  device = cairo_ps,width = 6.94, height = 6.44)

# Second step B: detecting the best solutions conditioning on alpha =.08-----------------------------

# Fixing alpha
solution_1_alpha_fixed <- model_results %>% 
  filter(trim_level==.08) %>% 
  select(1:3,5:6) %>% 
  ungroup() %>% 
  mutate(rank_solution = rank(TBIC, ties.method = "first"))

threshold_ARI <- .7

output_second_step_LIGHTER <-
  solutions_ranking_identifier_LIGHTER(
    solutions_df = solution_1_alpha_fixed,
    number_of_optimal_sol = 4,
    epsilon_ARI = threshold_ARI
  )


output_second_step_LIGHTER$optimal_solutions
output_second_step_LIGHTER$df_for_plotting %>% 
  filter(num_groups!=1) %>% 
  mutate(to_fill=case_when((is_stable & is_best_interval)~TRUE,
                           is_stable~FALSE,
                           TRUE~NA),
         n_sol=factor(ifelse(is.na(to_fill),NA, n_sol)),
         num_groups=paste("G = ", num_groups)) %>% 
  ggplot(aes(x=restr_factor_X,y=restr_factor_Y))+
  geom_tile(aes(fill=n_sol,alpha=(to_fill)|!is.na(rank_best)), colour = "black", na.rm = TRUE, size = 0.5)+
  # geom_tile(aes(fill=n_sol),colour = "black", na.rm = TRUE, size = 0.5)+
  geom_text(aes(label=rank_best)) +
  labs(x=quote(c[X]),y=quote(c[y]), fill="Optimal solutions")+
  theme_bw() +
  guides(alpha=FALSE) +
  facet_wrap (~ num_groups,nrow = 1) +
  scale_fill_brewer(palette="Dark2",na.value = "white", breaks=1:(n_distinct(output_second_step_LIGHTER$df_for_plotting$n_sol)-1)) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggsave(
  filename = "output/step_2_B_monitoring_sim.eps",
  device = cairo_ps,width = 13.7, height = 4.36)

pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[4]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[5]]$assig+1)

table(output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[4]]$assig, class_out)


# Manual relabeling of the classes different sol

first_cl <- output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig
second_cl <- output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig
third_cl <- output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig

table(first_cl, class_out)

table(third_cl, class_out)

first_true_relabeler <- c(3,2,1)
first_cl_relabeled <- first_cl

for(k in 1:length(first_true_relabeler)){
  first_cl_relabeled[first_cl==k] <- first_true_relabeler[k]
}

table(first_cl_relabeled, class_out)

table(second_cl, class_out)

second_true_relabeler <- c(2,4,1,3)
second_cl_relabeled <- second_cl

for(k in 1:length(second_true_relabeler)){
  second_cl_relabeled[second_cl==k] <- second_true_relabeler[k]
}

table(second_cl_relabeled, class_out)

table(third_cl, class_out)
third_true_relabeler <- c(4,3,1,2)

third_cl_relabeled <- third_cl

for(k in 1:length(third_true_relabeler)){
  third_cl_relabeled[third_cl==k] <- third_true_relabeler[k]
}

table(third_cl_relabeled, class_out)

df_1 <- mutate(plotyy_df, pop=first_cl_relabeled)
df_2 <- mutate(plotyy_df, pop=second_cl_relabeled)
df_3 <- mutate(plotyy_df, pop=third_cl_relabeled)

ggplot_solutions_CWRM(df_1)+
  ggtitle(expression(paste("First optimal solution (", alpha==0.08,")")))

ggsave(
  filename = "output/first_opt_sol_B.eps",
  device = cairo_ps,width = 6.94, height = 6.44)

ggplot_solutions_CWRM(df = df_2)+
  ggtitle(expression(paste("Second optimal solution (", alpha==0.08,")")))

ggsave(
  filename = "output/second_opt_sol_B.eps",
  device = cairo_ps,width = 6.94, height = 6.44
)

ggplot_solutions_CWRM(df_3)+
  ggtitle(expression(paste("Third optimal solution (", alpha==0.08,")")))

ggsave(
  filename = "output/third_opt_sol_B.eps",
  device = cairo_ps,width = 6.94, height = 6.44)

# Second step C: detecting the best solutions conditioning on alpha =0-----------------------------

# Fixing alpha
solution_1_alpha_fixed <- model_results %>% 
  filter(trim_level==0) %>% 
  select(1:3,5:6) %>% 
  ungroup() %>% 
  mutate(rank_solution = rank(TBIC, ties.method = "first"))

threshold_ARI <- .7

output_second_step_LIGHTER <-
  solutions_ranking_identifier_LIGHTER(
    solutions_df = solution_1_alpha_fixed,
    number_of_optimal_sol = 4,
    epsilon_ARI = threshold_ARI
  )


output_second_step_LIGHTER$optimal_solutions
output_second_step_LIGHTER$df_for_plotting %>% 
  filter(num_groups!=2) %>%
  mutate(to_fill=case_when((is_stable & is_best_interval)~TRUE,
                           is_stable~FALSE,
                           TRUE~NA),
         n_sol=factor(ifelse(is.na(to_fill),NA, n_sol)),
         num_groups=paste("G = ", num_groups)) %>% 
  ggplot(aes(x=restr_factor_X,y=restr_factor_Y))+
  geom_tile(aes(fill=n_sol,alpha=(to_fill)|!is.na(rank_best)), colour = "black", na.rm = TRUE, size = 0.5)+
  # geom_tile(aes(fill=n_sol),colour = "black", na.rm = TRUE, size = 0.5)+
  geom_text(aes(label=rank_best)) +
  labs(x=quote(c[X]),y=quote(c[y]), fill="Optimal solutions")+
  theme_bw() +
  guides(alpha=FALSE) +
  facet_wrap (~ num_groups,nrow = 1) +
  scale_fill_brewer(palette="Dark2",na.value = "white", breaks=1:(n_distinct(output_second_step_LIGHTER$df_for_plotting$n_sol)-1)) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggsave(
  filename = "output/step_2_C_monitoring_sim.eps",
  device = cairo_ps,width = 13.7, height = 4.36)

pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig+1)
pairs(df_out,col=output_second_step_LIGHTER$optimal_solutions$model[[4]]$assig+1)

table(output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig, class_out)
table(output_second_step_LIGHTER$optimal_solutions$model[[4]]$assig, class_out)


# Manual relabeling of the classes different sol

first_cl <- output_second_step_LIGHTER$optimal_solutions$model[[1]]$assig
second_cl <- output_second_step_LIGHTER$optimal_solutions$model[[2]]$assig
third_cl <- output_second_step_LIGHTER$optimal_solutions$model[[3]]$assig

table(first_cl, class_out)

table(third_cl, class_out)

first_true_relabeler <- c(1,2,3,4,5)
first_cl_relabeled <- first_cl

for(k in 1:length(first_true_relabeler)){
  first_cl_relabeled[first_cl==k] <- first_true_relabeler[k]
}

table(first_cl_relabeled, class_out)
table(second_cl, class_out)
second_true_relabeler <- c(4,5,1,2)
second_cl_relabeled <- second_cl

for(k in 1:length(second_true_relabeler)){
  second_cl_relabeled[second_cl==k] <- second_true_relabeler[k]
}

table(second_cl_relabeled, class_out)

table(third_cl, class_out)
third_true_relabeler <- c(5,1,3,4,2)

third_cl_relabeled <- third_cl

for(k in 1:length(third_true_relabeler)){
  third_cl_relabeled[third_cl==k] <- third_true_relabeler[k]
}

table(third_cl_relabeled, class_out)

df_1 <- mutate(plotyy_df, pop=first_cl_relabeled)
df_2 <- mutate(plotyy_df, pop=second_cl_relabeled)
df_3 <- mutate(plotyy_df, pop=third_cl_relabeled)

ggplot_solutions_CWRM(df_1)+
  ggtitle(expression(paste("First optimal solution (", alpha==0,")")))

ggsave(
  filename = "output/first_opt_sol_C.eps",
  device = cairo_ps,width = 6.94, height = 6.44)

ggplot_solutions_CWRM(df = df_2)+
  ggtitle(expression(paste("Second optimal solution (", alpha==0,")")))

ggsave(
  filename = "output/second_opt_sol_C.eps",
  device = cairo_ps,width = 6.94, height = 6.44
)

ggplot_solutions_CWRM(df_3)+
  ggtitle(expression(paste("Third optimal solution (", alpha==0,")")))

ggsave(
  filename = "output/third_opt_sol_C.eps",
  device = cairo_ps,width = 6.94, height = 6.44)
