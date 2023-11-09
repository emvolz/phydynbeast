# `phydynbeast`

This package allows you to specify arbitrary compartmental models for fitting with genetic data using the BEAST2 PhyDyn package. The model is specified as a list of equations and a list of parameters that specify how different populations should be interpreted within a population-genetics framework. As the model is fitted by Bayesian MCMC within BEAST, all parameters which are estimated must have corresponding priors and operators for sampling. 

Documentation for the BEAST2 PhyDyn package can be found at [https://github.com/mrc-ide/PhyDyn](https://github.com/mrc-ide/PhyDyn).

Note that all PhyDyn analyses require an XML template which can be made in [beauti](https://www.beast2.org/beauti/). An example template for this analysis is included in the package. 
It is easy to create such a template by first configuring an XML to use a Yule tree prior and then removing XML elements that specify the Yule prior and operator. 

# Example

Here is an example which configures a model featuring 

- exponential growth in two demes
- separate birth and death rates
- migration between demes

This model is used to estimate the growth of the BA.2.86 SARS-CoV-2 lineage in the United Kingdom whilst accounting for migration with a global reservoir. 

```
# Configure a model with exponential growth, birth(beta) & death(gamma), in two demes, with migration (mu) between demes
eqns <- list( 
 	confeqn( 'beta*gb', type = 'birth', origin = 'gb', destination = 'gb' )
 	, confeqn( 'beta*oth', type = 'birth', origin = 'oth', destination = 'oth' )
 	, confeqn( 'mu*oth', type = 'migration', origin = 'oth', destination = 'gb' )
 	, confeqn( 'mu*gb', type = 'migration', origin = 'gb', destination = 'oth' )
 	, confeqn( 'gamma*gb', type = 'death', origin = 'gb' )
 	, confeqn( 'gamma*oth', type = 'death', origin = 'oth' )
 )
 parms <- list(
 	      confparm( 'beta'
 	      	       , initial = 1.5*365*(1/5.5)
 	      	       , prior = 'lognormal'
 	      	       , operator = 'realrw' 
 	      	       , lower = .5*365*(1/5.5)
 	      	       , upper = 4*365*(1/5.5)
 	      	       , M = log(  1.5*365*(1/5.5)) , S = 1.0
 	      	       )
 	      , confparm( 'gamma'
 	      		 , initial = 365*(1/5.5)
 	      		 , estimate = FALSE 
 	      		 )
 	      , confparm( 'mu'
 	      		 , initial = 12 
 	      		 , prior = 'exponential'
 	      		 , operator = 'realrw' 
 	      		 , estimate = TRUE
 	      		 , lower = 1/4
 	      		 , upper = 52
 	      		 , mean = 12 
 	      		 )
		# estimate initial deme sizes
 	      , confparm( 'gb'
 	      		 , initial_condition_parameter = TRUE 
 	      		 , initial = 1e-2
 	      		 , estimate = TRUE
 	      		 , prior = 'exponential'
 	      		 , operator = 'realrw'
 	      		 , lower = 0
 	      		 , upper = 1e2 
 	      		 , mean = 10
 	      		 )
 	      , confparm( 'oth'
 	      		 , initial_condition_parameter = TRUE 
 	      		 , initial = 1e-2
 	      		 , estimate = TRUE 
 	      		 , prior = 'exponential'
 	      		 , operator = 'realrw'
 	      		 , lower = 0 
 	      		 , upper = 1e2
 	      		 , mean = 10 
 	      		 )
 	)
 model <- config_phydyn(
  	 system.file( 'extdata', 'ba.2.86_algnWu-Hu-1.1_qc0_beast_template.xml', package = 'phydynbeast' )
  	  , saveto = 'ba.2.86_simodel0.xml'
  	  , t0 = confparm( 't0'
 	      		 , estimate = FALSE
 	      		 , initial = 2023.25
 	      		 )
 	      )
  	  , equations = eqns  
  	  , parameters = parms
  	  , coalescent_approximation = 'PL1'
  	  , integrationSteps = 100 
  	  , minP = 0.001
  	  , penaltyAgtY= 0 
  	  , useStateName = TRUE 
  	  , traj_log_file = 'simodel0-traj.tsv'
  	  , traj_log_frequency = 10000
)
```

# Installation


Install `devtools` and run 

```
devtools::install_github('emvolz/phydynbeast')
```
