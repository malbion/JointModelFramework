# This functions spits out diagnostics for convergence
stan_diagnostic <- function(fit, 
                            results_folder,
                            ...) {
  
  # Print diagnostics
  check_hmc_diagnostics(fit)
  
  check_treedepth(fit)
  # if fit saturates the threshold, rerun with larger maximum tree depth and recheck saturation
  check_energy(fit) # check E-BFMI
  check_divergences(fit) # check divergence (validity concern)
  
  fit_summary <- summary(fit)$summary
  # get the range of rhat values
  print(c('rhat range: ', range(na.omit(fit_summary[ , 'Rhat']))))
  # and n_eff
  print(c('n_eff range: ', range(na.omit(fit_summary[ , 'n_eff']))))
  
  # print some diagnostic plots 
  png(paste0(results_folder, '/diag_rhat.png'), width=800, height=800)
  x <- stan_rhat(fit)  # distribution of Rhat
  print(x)
  dev.off()
  png(paste0(results_folder, '/diag_effSS_totalSS.png'), width=800, height=800)
  x <- stan_ess(fit)   # ratio of effective sample size to total sample size
  print(x)
  dev.off()
  png(paste0(results_folder, '/diag_MCse_postsd.png'), width=800, height=800)
  x <- stan_mcse(fit)  # ratio of Monte Carlo standard error to posterior standard deviation
  print(x)
  dev.off() 
  
  #Histograms of lp__ and accept_stat, and their joint distribution
  png(paste0(results_folder, '/diag_sample.png'), width=800, height=800)
  stan_diag(fit, 'sample')
  dev.off()
  # Violin plots showing the distributions of lp__ and accept_stat at each of the sampled step sizes (one per chain)
  png(paste0(results_folder, '/diag_stepsize.png'), width=800, height=800)
  stan_diag(fit, 'stepsize')
  dev.off()
  # Histogram of treedepth and violin plots showing the distributions of lp__ and accept_stat for each value of treedepth
  png(paste0(results_folder, '/diag_treedepth.png'), width=800, height=800)
  stan_diag(fit, 'treedepth')
  dev.off()
  # Violin plots showing the distributions of lp__ and accept_stat for iterations that encountered divergent transitions (divergent=1) and those that did not (divergent=0)
  png(paste0(results_folder, '/diag_divergence.png'), width=800, height=800)
  stan_diag(fit, 'divergence')
  dev.off()
  
  # NB: for further diagnostics, I can explore with
  # - stan_par(fit, 'specific parameter')
  # - stan_ac(fit, 'param') for auto-correlation
}



# This function checks the validity of the RE model

stan_model_check <- function(fit,
                             results_folder,
                             params = c('disp_dev', 'beta_i0', 'effect', 'response'),
                             ...) {
  
  fit_summary <- summary(fit)$summary
  
  # remove 'mu' from the parameter vector and do it after because it is too long
  params <- params[params != 'mu']
  
  sapply(params, function(x) { 
    
    # this is just to length the plots so they can show all parameters
    N <- length(grep(x, rownames(fit_summary)))
    # some exceptions: 
    if (x == 'beta_i0') {N <- 20}
    if (x == 're') {N <- 800}
    if (x == 'sigma') {N <- 10}
    
    # save traceplot of the posterior draws
    png(paste0(results_folder, '/tplot_', x, '.png'), width = 800, height = (N/5*100))
    tplot <- traceplot(fit, pars = x, inc_warmup = F, ncol = 5)
    print(tplot)
    dev.off()
    
    # plot posterior uncertainty intervals
    png(paste0(results_folder, '/postui_', x, '.png'), width = 600, height = (N/5*100))
    postint <- plot(fit, pars = x, show_density = T) 
    print(postint)
    dev.off()
    
    })
  
  
}


# This function plots the posterior predicted seed number vs the actual data

stan_post_pred_check <- function(post.draws,
                                 results_folder,
                                 stan.data,
                                 ...) {
  
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 
  
  # extract mu and phi
  mu <- post.draws$mu # matrix with nrow = draws and ncol = observations
  disp_dev <- post.draws$disp_dev
  phi <- (disp_dev^2)^(-1)
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i])  
    }
  }
  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                       max(density(stan.data$perform)$y)))
  
  # dev.new(noRStudioGD = T)
  png(paste0(results_folder, '/postpredch.png'), width=800, height=800)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ]), ylim = c(0, max.density), col = 'darkgrey',
                   ylab = 'Seed density',
                   main = 'Post. pred. check',
                   sub = '(grey = predicted, black = observed)') 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ]), col = 'darkgrey')
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data$perform), col = 'black', lwd = 2)  
  print(ppc.plot)
  dev.off()
  
}