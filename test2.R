library(remphasis)
library(remphasisrpd5c)

brts = 10:1
pars = c(0.1,0.8,-0.01,0)
model = "remphasisrpd5c"
it = remphasis::em_cpp(brts = brts,
                       init_pars = pars,
                       sample_size = 10000,
                       maxN = 100000,
                       plugin = locate_plugin(model),
                       soc = 2,
                       max_misssing = 1000000,
                       max_lambda = 100,
                       lower_bound = numeric(),
                       upper_bound = numeric(),
                       xtol_rel = 0.001,
                       num_threads=10,
                       copy_trees = FALSE)

#logf = sapply(it$trees,emphasis:::loglik.tree("rpd1"), pars=pars)
#logg = sapply(it$trees,emphasis:::sampling_prob, pars=pars,model="rpd1")
#log_weights = logf-logg
#prop_const = max(log_weights)
#log_weights_norm = log_weights - prop_const
#w = exp(log_weights_norm)
#log_fhat = log(mean(w)) + prop_const
show(it$fhat)
