var documenterSearchIndex = {"docs":
[{"location":"#Replication-of-\"Public-Debt,-Interest-Rates,-and-Negative-Shocks\"-(Evans,-R.-2020)","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"This replication study was part of our evaluation for the course Numerical Methods at SciencesPo Paris in Spring 2021The functions used to replicate this paper are:","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"Modules = [Evans2020]","category":"page"},{"location":"#Evans2020.LN_pdf-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.LN_pdf","text":"LN_pdf\n\n#--------------------------------------------------------------------\n#This function gives the PDF of the lognormal distribution for xvals\n#given mu and sigma\n\n#(LN): f(x; mu, sigma) = (1 / (x * sigma * sqrt(2 * pi))) *\n#        exp((-1 / 2) * (((log(x) - mu) / sigma) ** 2))\n#        x in [0, infty), mu in (-infty, infty), sigma > 0\n#--------------------------------------------------------------------\n#INPUTS:\n#xvals = (N,) vector, data\n#mu    = scalar, mean of the ln(x)\n#sigma = scalar > 0, standard deviation of ln(x)\n\n#OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None\n\n#OBJECTS CREATED WITHIN FUNCTION:\n#pdf_vals        = (N,) vector, probability of each observation given\n#                  the parameter values\n\n#FILES CREATED BY THIS FUNCTION: None\n\n#RETURNS: pdf_vals\n#--------------------------------------------------------------------\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_1pr_MU_c2_pdf-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_1pr_MU_c2_pdf","text":"get_1pr_MU_c2_pdf\n\n#This function is the target for calculating the integral\n#(expectation): E[(1+r_{tp1})*(c_{2,t+1})**(-gamma)]. This function\n#returns the value of\n#(1 + r_{tp1})*((c_{2,t+1})**(-gamma)) * pdf(A|mu,sigma)\n#for a given value of A and k2tp1\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_C-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_C","text":"get_C\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_Eul_err-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_Eul_err","text":"get_Eul_err\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_ExpMU_c2_b-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_ExpMU_c2_b","text":"get_ExpMU_c2_b\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_ExpMU_c2tp1_k-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_ExpMU_c2tp1_k","text":"get_ExpMU_c2tp1_k\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_Hbar_err-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_Hbar_err","text":"`getHbarerr\n\nThis function is the error function that solves for the current\nperiod shock that sets w * n1 + x1 - c_min - K_min = Hbar. This is\nthe minimum shock that does not create default.\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_Ht-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_Ht","text":"get_Ht\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_MU_c2_pdf-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_MU_c2_pdf","text":"get_MU_c2_pdf\n\nThis function is the target for calculating the integral (expectation): E[(c{2,t+1})^(-gamma)].  This function returns the value of ((c{2,t+1})^(-gamma)) * pdf(A|mu,sigma) for a given value of A and k2tp1\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_MUc_CRRA-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_MUc_CRRA","text":"get_MUc_CRRA\n\n\n\nGenerate marginal utility(ies) of consumption with CRRA consumption utility and stitched function at lower bound such that the new hybrid function is defined over all consumption on the real line but the function has similar properties to the Inada condition. u'(c) = c ** (-sigma) if c >= epsilon      = g'(c) = 2 * b2 * c + b1 if c < epsilon\n\nsuch that g'(epsilon) = u'(epsilon) and g''(epsilon) = u''(epsilon)\n\nu(c) = (c ** (1 - sigma) - 1) / (1 - sigma) g(c) = b2 * (c ** 2) + b1 * c + b0 –––––––––––––––––––––––––––––––––– INPUTS: c  = scalar, individual consumption in a particular period gamma = scalar >= 1, coefficient of relative risk aversion for CRRA        utility function: (c**(1-gamma) - 1) / (1 - gamma) graph = boolean, =True if want plot of stitched marginal utility of         consumption function\n\nOTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None\n\nOBJECTS CREATED WITHIN FUNCTION: epsilon    = scalar > 0, positive value close to zero cs        = scalar, individual consumption cscnstr  = boolean, =True if cs < epsilon b1         = scalar, intercept value in linear marginal utility b2         = scalar, slope coefficient in linear marginal utility MUc       = scalar or (p,) vector, marginal utility of consumption              or vector of marginal utilities of consumption p          = integer >= 1, number of periods remaining in lifetime cveccnstr = (p,) boolean vector, =True for values of cvec < epsilon\n\nFILES CREATED BY THIS FUNCTION:     MUcstitched.png\n\nRETURNS: MU_c\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_Y-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_Y","text":"get_Y(k2t, zt, args)\n\n--------------------------------------------------------------------\nCalculate aggregate output\n--------------------------------------------------------------------\nINPUTS:\nKt   =\nLt   =\nzt   =\nargs =\n\nRETURNS: Yt\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_c1mgam-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_c1mgam","text":"get_c1mgam\n\n#--------------------------------------------------------------------\n#Generate marginal utility(ies) of consumption with CRRA consumption\n#utility and stitched function at lower bound such that the new\n#hybrid function is defined over all consumption on the real\n#line but the function has similar properties to the Inada condition.\n\n#f(c) = c ** (1-sigma) if c >= epsilon\n#g(c) = b2 * c + b1    if c < epsilon\n\n#    such that g(epsilon) = f(epsilon)\n#    and g'(epsilon) = f'(epsilon)\n\n#    f(c) = c ** (1 - sigma)\n#    g(c) = b2 * c + b1\n\n#    s.t. b2 = (1 - gamma) * (epsilon ** (-gamma))\n#         b1 = epsilon**(-gamma) - (1-gamma) * (epsilon ** (1-gamma))\n#--------------------------------------------------------------------\n#INPUTS:\n#c  = scalar, individual consumption in a particular period\n#gamma = scalar >= 1, coefficient of relative risk aversion for CRRA\n#        utility function: (c**(1-gamma) - 1) / (1 - gamma)\n#graph = boolean, =True if want plot of stitched marginal utility of\n#        consumption function\n\n#OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None\n\n#OBJECTS CREATED WITHIN FUNCTION:\n#epsilon    = scalar > 0, positive value close to zero\n#b1         = scalar, intercept value in linear marginal utility\n#b2         = scalar, slope coefficient in linear marginal utility\n#MU_c       = scalar or (p,) vector, marginal utility of consumption\n#             or vector of marginal utilities of consumption\n#p          = integer >= 1, number of periods remaining in lifetime\n#cvec_cnstr = (p,) boolean vector, =True for values of cvec < epsilon\n\n#FILES CREATED BY THIS FUNCTION:\n#    MU_c_stitched.png\n\n#RETURNS: f_c\n#--------------------------------------------------------------------\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_c2t-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_c2t","text":"get_c2t\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_c2tp1_1mgam_pdf-Tuple{Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_c2tp1_1mgam_pdf","text":"get_c2tp1_1mgam_pdf\n\nThis function is the target for calculating the integral (expectation): E[(c{2,t+1})^(1-gamma)].  This function returns the value of ((c{2,t+1})^(1-gamma)) * pdf(A|mu,sigma) for a given value of A and k2tp1\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_k2tp1-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_k2tp1","text":"get_k2tp1\n\n\n\nSolve for k2tp1 c1t + k2tp1 = wt * n1 - tau * w1 * n1 ––––––––––––––––––––––––––––––––––\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_r-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_r","text":"get_r\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_w-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_w","text":"get_w\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.get_zstar-Tuple{Any, Any, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.get_zstar","text":"get_zstar\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.sim_timepath-NTuple{28, Any}","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.sim_timepath","text":"simtimepath(     Hbar, beta, gamma, k20, sigma, x1, T, z0, zmin, rho, mu, nvec,     epsilon, alpha, delta, tau, cmin, Kmin, Amin, yrsinper,     Hind, risktypeind, riskvalind,     avgrtp1ind, avgrbartind, Sind, ztvec,     rand_seed)\n\nRuns the simulation\n\n\n\n\n\n","category":"method"},{"location":"#Evans2020.trunc_norm_draws","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Evans2020.trunc_norm_draws","text":"trunc_norm_draws(unifvals::Any, mu::Any, sigma::Any, cutlb::Any = nothing, cut_ub::Any = nothing)\n\n--------------------------------------------------------------------\nDraw (N x S) matrix of random draws from a truncated normal\ndistribution based on a normal distribution with mean mu and\nstandard deviation sigma and cutoffs (cut_lb, cut_ub). These draws\ncorrespond to an (N x S) matrix of randomly generated draws from a\nuniform distribution U(0,1).\n--------------------------------------------------------------------\nINPUTS:\nunif_vals = (N, S) matrix, (N,) vector, or scalar in (0,1), random\n            draws from uniform U(0,1) distribution\nmu        = scalar, mean of the nontruncated normal distribution\n            from which the truncated normal is derived\nsigma     = scalar > 0, standard deviation of the nontruncated\n            normal distribution from which the truncated normal is\n            derived\ncut_lb    = scalar or string, ='None' if no lower bound cutoff is\n            given, otherwise is scalar lower bound value of\n            distribution. Values below this cutoff have zero\n            probability\ncut_ub    = scalar or string, ='None' if no upper bound cutoff is\n            given, otherwise is scalar lower bound value of\n            distribution. Values below this cutoff have zero\n            probability\n\nOTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION:\n    scipy.stats.norm()\n\nOBJECTS CREATED WITHIN FUNCTION:\ncut_ub_cdf  = scalar in [0, 1], cdf of N(mu, sigma) at upper bound\n              cutoff of truncated normal distribution\ncut_lb_cdf  = scalar in [0, 1], cdf of N(mu, sigma) at lower bound\n              cutoff of truncated normal distribution\nunif2_vals  = (N, S) matrix, (N,) vector, or scalar in (0,1),\n              rescaled uniform derived from original.\ntnorm_draws = (N, S) matrix, (N,) vector, or scalar in (0,1),\n              values drawn from truncated normal PDF with base\n              normal distribution N(mu, sigma) and cutoffs\n              (cut_lb, cut_ub)\n\nFILES CREATED BY THIS FUNCTION: None\n\nRETURNS: tnorm_draws\n--------------------------------------------------------------------\n\n\n\n\n\n","category":"function"},{"location":"#Case-number-1-:-and-epsilon;-1-and-and-mu;-is-constant","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Case number 1 : &epsilon; = 1 and &mu; is constant","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Evans2020\n\njulia> ","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"end","category":"page"},{"location":"#Our-Replication-of-The-Cash-Paradox-(Jiang-and-Shao,-2019)","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"This replication study was part of our evaluation for the course Numerical Methods at SciencesPo Paris in Spring 2021","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"The functions used to replicate this paper are:","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"Modules = [CashParadox]","category":"page"},{"location":"#Replication-of-model-predictions:-Figures-5a-5d","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of model predictions: Figures 5a-5d","text":"","category":"section"},{"location":"#Figure-5a:-AUSTRALIA","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure 5a: AUSTRALIA","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> Fig5(2)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: fig5aus)","category":"page"},{"location":"#Figure-5b:-CANADA","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure 5b: CANADA","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> Fig5(0)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: fig5can)","category":"page"},{"location":"#Figure-5c:-UK","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure 5c:  UK","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> Fig5(3)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: fig5uk)","category":"page"},{"location":"#Figure-5d:-US","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure 5d: US","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> Fig5(1)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: fig5us)","category":"page"},{"location":"#Replication-of-model-predictions:-Figures-A2a-A2d","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of model predictions: Figures A2a-A2d","text":"","category":"section"},{"location":"#Figure-A2a:-AUSTRALIA","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure A2a: AUSTRALIA","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA2(2)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA2aus)","category":"page"},{"location":"#Figure-A2b:-CANADA","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure A2b: CANADA","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA2(0)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA2can)","category":"page"},{"location":"#Figure-A2c:-UK","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure A2c:  UK","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA2(3)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA2uk)","category":"page"},{"location":"#Figure-A2d:-US","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Figure A2d: US","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA2(1)","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA2us)","category":"page"},{"location":"#Replication-of-regime-changes-over-time:-Figure-A3","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of regime changes over time: Figure A3","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA3()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA3)","category":"page"},{"location":"#Replication-of-The-Value-of-ATM-withdrawals-over-CIC:-Figure-A4","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of The Value of ATM withdrawals over CIC: Figure A4","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA4()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA4)","category":"page"},{"location":"#Replication-of-Cash-receipts-from-circulation-in-the-Federal-Reserve-Banks:-Figure-A5","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of Cash receipts from circulation in the Federal Reserve Banks: Figure A5","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigA5()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figA5)","category":"page"},{"location":"#Replication-of-Different-measures-of-nominal-interest-rates:-Figure-D1","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of Different measures of nominal interest rates: Figure D1","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigD1()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figD1)","category":"page"},{"location":"#Replication-of-CIC/GDP-with-different-interest-rate-specifications:-Figure-D2","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of CIC/GDP with different interest rate specifications: Figure D2","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> FigD2()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figD2)","category":"page"},{"location":"#Replication-of-Table-1-and-Table-2:-Calibration-results-and-Cash-shares-relative-to-credit","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of Table 1 and Table 2: Calibration results and Cash shares relative to credit","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"Most of results in this table are very close to those from Jiang & Shao (2019). However, we encounter problems replicating the result for the NCF model with UK data.","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create a dataframe containing this results type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> Table1()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figG1)","category":"page"},{"location":"#Replication-of-Table-D.1:-Parameter-values-and-Table-D.2:-Model-performance-comparison","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of Table D.1: Parameter values and Table D.2: Model performance comparison","text":"","category":"section"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"In order to create a dataframe containing this results type","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"julia> using Cash Paradox\n\njulia> Table2()","category":"page"},{"location":"","page":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","title":"Replication of \"Public Debt, Interest Rates, and Negative Shocks\" (Evans, R. 2020)","text":"(Image: figD1D2)","category":"page"}]
}
