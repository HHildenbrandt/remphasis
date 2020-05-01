### EMPHASIS functions
emphasis <- function(brts,
                     soc = 2,
                     model = "remphasisddd",
                     init_par = c(0.05, 0.5, 0.0),
                     lower_bound = c(0, 0, -Inf),  
                     upper_bound = c(Inf, Inf, Inf),
                     max_lambda = 500,
                     xtol = 0.001, 
                     tol = 0.01,
                     verbose = FALSE,
                     em_tol = 0.25,
                     sample_size_tol = 0.005,
                     return_trees = FALSE,
                     max_missing = 10000,  # maximum tree size
                     burnin_sample_size = 200,
                     pilot_sample_size = seq(100, 1000, by = 100),
                     burnin_iterations = 20,
                     num_threads = 0) {
  
  msg1 = paste("Initializing emphasis...")
  msg2 = paste("Age of the tree: ",max(brts))
  msg3 = paste("Number of speciations: ",length(brts))
  msg4 = paste("Diversification model to fit:",model)
  msg5 = "######################################"
  cat(msg1, msg2, msg3, msg4, msg5, sep = "\n")
  
  cat( "Performing Phase 1: burn-in", sep = "\n")
  mc = mcEM_step(brts = brts,
            pars = init_par,
            sample_size = burnin_sample_size,
            model = model,
            soc = soc,
            max_missing = max_missing,
            max_lambda = max_lambda,
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            xtol = xtol,
            num_threads = num_threads,
            return_trees = FALSE,
            verbose = FALSE,
            tol = em_tol,
            burnin = burnin_iterations)
  
  M = mc$mcem
  pars = c(mean(tail(M$par1,n = nrow(M)/2)),
           mean(tail(M$par2,n = nrow(M)/2)),
           mean(tail(M$par3,n = nrow(M)/2)),
           mean(tail(M$par4,n = nrow(M)/2)))
  
  cat("\n", msg5, sep = "\n")
  cat( "Phase 2: Assesing required MC sampling size \n")
  
  for (i in 1:length(pilot_sample_size)) {
    cat(paste("\n Sampling size: ", as.character(pilot_sample_size[i]),"\n"))
    mc = mcEM_step(brts = brts,
              pars = pars,
              sample_size = pilot_sample_size[i],
              model = model,
              soc = soc,
              max_missing = max_missing,
              max_lambda = max_lambda,
              lower_bound = lower_bound,
              upper_bound = upper_bound,
              xtol = xtol,
              num_threads = num_threads,
              return_trees = FALSE,
              verbose = FALSE,
              tol = em_tol,
              burnin = 10)
              
    ta = tail(mc$mcem,n = nrow(M)/2)
    pars = c(mean(ta$par1),mean(ta$par2),mean(ta$par3),mean(ta$par4))
    M = rbind(M,mc$mcem)
  }
  n.r = get_required_sampling_size(M[-(1:burnin_iterations), ], 
                                   tol = sample_size_tol)
  sample_size = max(pilot_sample_size + 2, n.r)
  n.r_old = -1
  j = 1
  while (n.r_old < n.r) {
    msg6 = paste0("Required sampling size: ",n.r)
    msg7 = paste0("Phase 3: Performing metaiteration: ",j)
    cat("\n", msg5, msg7, msg6, sep = "\n")
    mc = mcEM_step(brts = brts,
              pars = pars,
              sample_size = sample_size,
              model = model,
              soc = soc,
              max_missing = max_missing,
              max_lambda = max_lambda,
              lower_bound = lower_bound,
              upper_bound = upper_bound,
              xtol = xtol,
              num_threads = num_threads,
              return_trees = FALSE,
              verbose = FALSE,
              tol = em_tol,
              burnin = 2)
        
    M <- rbind(M, mc$mcem)
    n.r_old = n.r
    j = j + 1
    n.r = get_required_sampling_size(M[-(1:burnin_iterations), ], 
                                     tol = sample_size_tol)
    pars = as.numeric(colMeans(mc$mcem)[1:4])
    sample_size = n.r
  }
  
  cat(pars)
  return(list(pars = pars, MCEM = M))
}




mcEM_step <- function(brts, 
                      pars, 
                      sample_size = 10000, 
                      model = "DDD", 
                      soc = 2, 
                      max_missing = 10000,
                      max_lambda = 500,
                      lower_bound = c(0, 0, -Inf),
                      upper_bound = c(Inf, Inf ,Inf),
                      xtol = 0.001, # tolerance in M step
                      tol = 0.01,   #tolerance of mcEM step
                      burnin = 20,
                      num_threads = 0,
                      return_trees = FALSE, 
                      verbose = TRUE) {
  mcem = NULL
  sde = 10 
  i = 0
  times = NULL
  while (sde > tol) {
    i <- i + 1
    results = remphasis::e_and_m_step(brts,
                                  pars,
                                  sample_size,                       
                                  model,             
                                  soc,                           
                                  max_missing,               
                                  max_lambda,             
                                  lower_bound,  
                                  upper_bound,  
                                  xtol = 0.001,                     
                                  num_threads,
                                  return_trees,
                                  verbose)
    
    mcem <- rbind(mcem, data.frame(pars = results$estimates, 
                                   fhat = results$ll, 
                                   sample_size = sample_size))
    
    if (verbose) {
      print(paste("(mean of) loglikelihood estimation: ",mean(mcem$fhat)))
    }
    
    times <- c(times, mcem$time)
    time_p_it = mean(times)
    
    if (i > burnin) {
      mcem_est = mcem[floor(nrow(mcem)/2):nrow(mcem),]
      mcem_est = mcem_est[is.finite(mcem_est$fhat),]
      sde0 = sde
      sde = sd(mcem_est$fhat)/sqrt(nrow(mcem_est))
      mde = mean(mcem_est$fhat)
      msg = paste("Iteration:",i," SE of the loglikelihood: ",sde)
      cat("\r",msg) 
    } else {
      msg = paste("Remaining time (burn-in): ", 
                  round(time_p_it * (burnin - i), digits = 0),"sec")
      cat("\r",msg) 
    }
  }
  return(list(mcem = mcem))
}

get_required_sampling_size <- function(M, tol=.05){
  n <- M$sample_size
  f<-  M$fhat
  hlp<-MASS:::rlm(f~I(1/n),weights = n)
  ab<-coef(hlp)
  
  f.r<-ab[1]-tol
  n.r<-ceiling(ab[2]/(f.r-ab[1]))
  return(n.r)
}

e_and_m_step <- function(brts,
                         pars,
                         sample_size,                       
                         model,             
                         soc,                           
                         max_missing,               
                         max_lambda,             
                         lower_bound,  
                         upper_bound,  
                         xtol = 0.001,                     
                         num_threads,
                         return_trees,
                         verbose) 
 {
    return(remphasis::mcem_cpp(brts,
                               pars,
                               sample_size,
                               locate_plugin(model),
                               soc,                           
                               max_missing,               
                               max_lambda,             
                               lower_bound,  
                               upper_bound,  
                               xtol = 0.001,                     
                               num_threads,
                               return_trees,
                               verbose))
}


  



