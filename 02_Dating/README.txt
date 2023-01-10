# DATING WITH MCMCTREE

# For the approximate likelihood method, run first using usedata=3, and then with usedata=2. See MCMCtree documentation.

#PARAMETERS/Priors

# Root age constraint

RootAge = '<7'  * safe constraint on root age, used if no fossil for root.
 # I use 700 Ma as a "safe" soft constraint, precambrian LCA of annelids, rotifers, mollusca.

# Prior for rate, follows a gamma distribution with parameters alpha and beta

rgene_gamma = 1 10    

	# mean rate = alpha/beta. I use here = 0.1 sub / unit time (here =100Ma), so prior for rate is 0.1 subs/100Ma based ~ on posterior rate
	# estimate in BEAST results of Kocot et al 2020 Monoplacophora Sci Rep = ~0.08 subs/100Ma
	# alpha = 2 is a "diffuse/uninformative prior", following MCMCTree manual. Variance of the distribution is alpha/beta^2 so increasing alpha&beta 
	# proportionally (i.e. keeping the mean constant) leads to lower variance = more constrained prior

	# alternative: run CODEML with few constraints and strict clock to have rough estimate of subs rate and use that as prior (PAML manual, page 43)
	# use @ or = to assign node ages in nwk file for this analysis. BECAREFUL, v4.9 doesn't work, use v4.8
	# Analysis in CODEML gave a mean subst for all 50 genes of ~0.09/100Ma which is rougly equal to BEAST estimates from Kocot. (RESULTS NOT INCLUDED IN SUPPLEMENTARY DATA)

 "baseml and codeml, clock dating.  When clock=1 and some internal
nodes are fixed at certain ages (by using the symbol '@X.XX' or '=X.XX'),
baseml and codeml rescale the estimated relative node ages (measured
in expected numbers of substitutions per site) into absolute node
ages.  In such cases there should be a block of output including
"Substitution rate is per time unit" and the absolute node ages below
"Detailed output identifying parameters". " https://github.com/abacus-gene/paml/blob/master/doc/pamlHistory.txt

# Prior for rate variation across branches, follows gamma distribution

sigma2_gamma = 1 10   * gammaDir prior for sigma^2 

	# I use 1 10 conservatively (= 0.1 mean rate variance). In Kocot et al 2020 BEAST results estimate a posterior or ~0.03 for the mean variance of the 
	# rates.


# priors on transition/transversions rate and among site variation

"Under HKY85 + >5, the transition/transversion rate ratio kappa is
assigned the gamma prior G(2, 0.5), with mean 4, while the shape parameter alpha for rate
variation among sites is assigned the prior G(2, 2), with mean 1. These two parameters are
reliably estimated in the data, so that the prior has little significance."  Yang's book


kappa_gamma = 6 2      * gamma prior for kappa
alpha_gamma = 2 2      * gamma prior for alpha

### NOTE ###
 # for gamma distribution rate = 1/scale

 # CALIBRATIONS

 NODE 106 'B(5.46,5.58)' following Kocot et al 2020 BEAST gamma distribution 2.5 tail values (Sup. Table 4). Diversification of Mollusks, first shell record
 NODE 134 'B(4.98,5.93)' following Kocot et al 2020 BEAST gamma distribution 2.5 tail values (Sup. Table 4). Split Cephalopoda/(Bivalvia.Gastro.Scapho)
 NODE 144 'B(5.27,5.42)' following Kocot et al 2020 BEAST gamma distribution 2.5 tail values (Sup. Table 4). Split Bivalvia/(Scapho.Gastro). Fordilla. Plectronoceras.
 NODE 145 'B(4.78,5.21)' following Kocot et al 2020 BEAST gamma distribution 2.5 tail values (Sup. Table 4). Diversification of Bivalvia.
 NODE 107 'B(4.78,5.21)' following Kocot et al 2020 BEAST gamma distribution 2.5 tail values (Sup. Table 4). Split Polyplacophora/Aplacophora. Orthriochiton.
 NODE 204 'B(3.47,3.85)' following Kocot et al 2020 BEAST gamma distribution 2.5 tail values (Sup. Table 4). Diversification of Scaphopoda. Pentalium.
 NODE 157 'L(4.75)'		 following Zapata et al 2014 lower bound Glyptarca serrata
 NODE 170 'L(4.18)'		 following Zapata et al 2014 lower bound Subtiloidea
 NODE 195 'B(1.30,1.52)' following Romero et al 2016 BEAST log-normal distribution 2.5 tail values. MRCA Ellobiidae. 
