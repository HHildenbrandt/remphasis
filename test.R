library(remphasis)
library(remphasisrpd1)
library(remphasisrpd5c)



brts_Megapodiidae = c(
	35.012472823,
	32.530356812,
	30.880632632,
	30.39947118,
	23.095866119,
	18.049627291,
	11.039829463,
	10.894029534,
	8.479030522,
	8.289813192,
	7.980299923,
	7.711562259,
	6.137002438,
	5.4937384316,
	4.191252593,
	3.078151366,
	3.026430002,
	2.456891288,
	1.836070379,
	1.262732134
)
pars = c(0.102054, 0.834852, -0.0361973)
sample_size = 100
so = locate_plugin("remphasisrpd5c")

cat('clade: Megapodiidae\n')
cat('model: rpd1 ', so, '\n')
cat('sample size: ', sample_size, '\n')
cat('initial pars:', pars, '\n')

#em <- em_cpp(brts_Megapodiidae, pars, sample_size, 10*sample_size, so, 2, 10000, 500, vector(), vector(), 0.001, 0, FALSE)
#show(em)

#e <- e_cpp(brts_Megapodiidae, pars, sample_size, 10*sample_size, so, 2, 10000, 500, vector(), vector(), 0.001, 0)
#m <- m_cpp(e, pars, so, vector(), vector(), 0.001, 0) 
#show(e)
#show(m)


# conditional closure
cond_closure = function(gam) {
  function(pars) {
    pr = as.numeric(mgcv:::predict.gam(gam,
                                       newdata = data.frame(mu=pars[1],
                                                            lambda=pars[2],
                                                            betaN=pars[3],
                                                            betaP=pars[4]),
                                       type = "response"))
    return(pr)
  }
}


e <- e_cpp(brts_Megapodiidae, pars, sample_size, 10*sample_size, so, 2, 10000, 500, vector(), vector(), 0.001, 0)
m <- m_cpp(e, pars, so, vector(), vector(), 0.001, 0, cond_closure(srv.gam)) 
#show(e)
#show(m)
