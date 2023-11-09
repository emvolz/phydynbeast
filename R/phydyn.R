PRIORS <- c('lognormal','normal', 'exponential', 'beta', '1/x' )
OPERATORS <- c('realrw', 'scale')
EQNTYPES <- c('birth', 'death', 'migration', 'nondeme', 'definition')
COAPPROX <- c('PL1', 'QL', 'PL2')

#' Configure an equation for phydyn analysis 
#' 
#' This function allows for the configuration of a model element, including rates of birth, death, migration between demes, auxliary dynamic state variables (nondeme), and definition of derived quantities. Some model elements (e.g. birth) must specify an *origin* and *destination* since births can happen within or between demes. Other elements (e.g. death) only need to specify *origin*. 
#' 
#' For information about PhyDyn model syntax and supported mathematical syntax and operations, see \url{https://github.com/mrc-ide/PhyDyn/wiki}. 
#' 
#' @param eqn A character string specifying the rate of a given event (birth, death, etc., see *type*)
#' @param type A character string specifying if the given equation gives the rate of birth, death, or migration between demes. Additionally, *eqn* can specify the rate of change of a `nondeme` variable, which is a dynamic variable that is not a deme.
#' @param origin If type is `birth`, `death`, `migration` or `nondeme`, this specifies the name of the deme or variable corresponding to the rate specified in `eqn`
#' @param destination If type is `birth` or `migration` this is a required name of a second deme that is part of the rate specified in `eqn`. For example, if `birth`, the `eqn` specifies the rate that units in deme `origin` generates a copy in `destination`. If `migration`, the `eqn` specifies the rate that units in deme `origin` changes state to deme `destination`. 
#' @return An equation of class `phydynequation`
#' @export 
confeqn <- function(
		  eqn
		  , type = c('birth', 'death', 'migration', 'nondeme', 'definition')
		  , origin = ''
		  , destination = ''
		  )
{
	type = type[1] 
	stopifnot( !is.na( type ))
	stopifnot( is.character( type ))

	origin = origin[1] 
	stopifnot( !is.na( origin ))
	stopifnot( is.character( origin ))

	destination = destination[1] 
	stopifnot( !is.na( destination ))
	stopifnot( is.character( destination ))

	stopifnot( type %in% EQNTYPES )
	if (type %in% c('birth', 'death', 'migration', 'nondeme' ))
		stopifnot( nchar( origin ) > 0 )
	if (type %in% c('birth', 'migration'))
		stopifnot( nchar( destination ) > 0 )

	structure(
		  list( eqn = eqn, type = type, origin = origin, destination = destination
		       , eqnid = glue::glue('{type}.{origin}.{destination}.{rlang::hash(eqn)}')
		       )
		  , class = 'phydynequation'
		  )
}

#' @export 
print.phydynequation <- function( x, ... ) 
{
	stopifnot( inherits( x, 'phydynequation' ))
	
	with(x,{
		if (type %in% c('birth', 'migration') ){
			cat( ifelse(type=='birth', 'Birth:\n', 'Migration:\n') )
			print( data.frame( `Origin deme` = origin
					  , `New deme` = destination
					  , Rate = eqn
		       ), row.names=FALSE, digits=3, right=TRUE)
		} else if ( type == 'death'){
			cat( 'Death or removal:\n' )
			print( data.frame( `Deme` = origin
					  , Rate = eqn
		       ), row.names=FALSE, digits=3, right=TRUE)
		} else if ( type == 'nondeme' ){
			cat('Dynamic state variable, non-deme:\n' )
			print( data.frame( `Variable` = origin
					  , Dynamics = glue::glue('d/dt {origin} = {eqn}')
					  )
			, row.names=FALSE, digits = 3, right = TRUE )
		} else if (type == 'definition'){
			cat( 'Derived variable:\n' )
			print( data.frame( 
					   Definition = eqn
					  )
			, row.names=FALSE, digits = 3, right = TRUE )
		}

	})
	invisible( NULL )
}

#' Configure a numeric parameter for PhyDyn analysis 
#' 
#' This sets the priors, operators (for MCMC sampling), upper and lower bounds, and initial conditions for a PhyDyn model parameter.
#' The names of parameters must correspond to those used in PhyDyn model equations. It is important to carefully specify priors and operators for all parameters since these choices can have a large influence on estimated posteriors and the efficiency of MCMC sampling. 
#' 
#' Supported priors are listed below. Note these may require additional arguments passed to ... 
#' 
#' \itemize{
#' 	\item  normal: This has required arguments `mean` and `sigma`
#' 	\item  lognormal: This has required arguments `M` and `S`
#' 	\item  exponential: This has required argument `mean`
#' 	\item  beta: This has required arguments `alpha` and `beta`. 
#' }
#' 
#' Supported operators which specify MCMC moves are listed below. Note that these may require additional arguments passed to ... 
#' 
#' \itemize{
#' 	\item  realrw: A random walk move. Should specify `windowSize` which is the width of the numeric interval used for the MCMC move. If not set, this will be set to 10% of the initial parameter value.  Optionally pass `useGaussian` set to "true"(default) or "false" to specify if normal or uniform moves should be used. 
#' 	\item  scale: This will inflate or shrink the parameter. Useful for parameters > 0 such as population sizes or rates. Should specify `scaleFactor` (default=0.5) which is the magnitude of the rescaling move. 
#' }
#' 
#' Passing `weight` (default=1) to ... will influence how frequently MCMC moves for this parameter are performed. 
#' 
#' @param name The name of the parameter 
#' @param initial_condition_parameter TRUE if this parameter specifies the initial value of a dynamic variable, such as the size of a deme specified in one of the PhyDyn equations. In that case, the parameter name should be the same as the deme name. 
#' @param initial The starting value (numeric) for this parameter in MCMC sampling.
#' @param estimate If TRUE, this parameter will be estimated by MCMC sampling. If FALSE, the value is fixed to *initial* and there is no need to specify priors or operators. 
#' @param prior A character string specifying the prior distribution. See details. 
#' @param operator A character string specifying the method of MCMC moves for this parameter. See details. 
#' @param lower Numeric lower bound of this parameter. 
#' @param upper Numeric upper bound of this parameter. 
#' @param ... Additional arguments may be required depending on the prior and operator. See details. 
#' @return A `phydynparameter` to be passed to `config_phydyn`
#' @export 
confparm <- function( 
			    name = ''
			    , initial_condition_parameter= FALSE 
			    , initial = NA 
			    , estimate = TRUE 
			    , prior = c('lognormal','normal', 'exponential', 'beta'  )
			    , operator = c('realrw', 'scale')
			    , lower = -Inf
			    , upper = Inf
			    , weight = 1
			    , ...
			    )
{
	#fargs <- c(as.list(environment()), list(...))
	stopifnot( is.character( name ) )
	stopifnot( length( name ) == 1 )
	stopifnot( nchar(name) > 0 )
	initial = initial[1] 
	stopifnot( !is.na( initial ))
	stopifnot( !is.infinite( initial ))
	prior <- prior[1] 
	stopifnot( prior %in% PRIORS )
	operator <- operator[1] 
	stopifnot( operator %in% OPERATORS)

	xargs <- list(...)

	if ( !('useGaussian' %in% names(xargs) ) )
		xargs$useGaussian = 'true' 
	if (!('windowSize' %in% names(xargs)))
		xargs$windowSize <- abs( .1 * initial ) 
	if ( initial == 0  & operator == 'realrw' & xargs$windowSize == 0)
		warning(glue::glue('windowSize for MCMC random walk set to 1.0 for parameter {name}. Set this value to ensure efficient MCMC. ')) 
	if ( xargs$windowSize == 0 ) 
		xargs$windowSize = 1 
	if ( !('scaleFactor'  %in% names(xargs )))
		xargs$scaleFactor <- 0.5 

	structure( list(name = name, prior = prior
			, operator = operator
			, estimate = estimate 
			, lower = lower 
			, upper = upper 
			, initial = initial 
			, weight = weight 
			, initial_condition_parameter =  initial_condition_parameter
			, parmid = glue::glue('phydyn.{name}') 
			, xargs = xargs 
			) |> c( xargs )
		  , class = 'phydynparameter' )
}

#' @export 
print.phydynparameter <- function(x, ... )
{
	stopifnot( inherits( x, 'phydynparameter' ))
	with( x, {
		cat( glue::glue('{name}:'))
		cat('\n')
		print( data.frame(
				  Parameter = name
				  , `Will estimate?` = ifelse(estimate, 'Yes', 'No' )
				  , `Initial guess` = initial 
				  , `Is this an initial conditions parameter?` = ifelse( initial_condition_parameter, 'Yes', 'No' )
				  ) , row.names=FALSE, digits=3, right=TRUE)
		if ( estimate ){
			print( data.frame( 
					  `Prior distribution` = prior
					  , `Lower bound` = lower
					  , `Upper bound` = upper
					  , `Operator type` = operator 
					  , `Operator weight` = weight
					  ), row.names=FALSE, digits=3, right=TRUE)
		}
		cat( 'Additional arguments: \n' )
		print( as.data.frame( xargs ), row.names=FALSE)
	})
	invisible( NULL )
}

.find_add  <- function( xml, searchstr, attrname, attrvalue, addvalue, addattrs )
{
	# find the node matching xpath searchstr and with arrname matching attrvalue 
	nodes0 <- xml2::xml_find_all(xml, searchstr )
	node <- nodes0[ which( xml2::xml_attr( nodes0, attrname ) == attrvalue ) ] 
	# add child with addvalue and addattrs 
	do.call( xml2::xml_add_child
		, list( .x = node, .value = addvalue ) |> c(addattrs) )
}


.find_add_text  <- function( xml, searchstr, attrname, attrvalue, addtext )
{
	# find the node matching xpath searchstr and with arrname matching attrvalue 
	nodes0 <- xml2::xml_find_all(xml, searchstr )
	node <- nodes0[ which( xml2::xml_attr( nodes0, attrname ) == attrvalue ) ] 
	print( node ) 
	xml2::xml_set_text(node, addtext ) 
}

.addeqn <- function( xml, eqn)
{
	stopifnot( inherits(eqn, 'phydynequation') )
	if( eqn$type == 'birth' ){
		.find_add( xml, searchstr = './/popmodel', attrname = 'spec', attrvalue = 'phydyn.model.PopModelODE'
		  	  , addvalue = 'matrixeq' 
		  	  , addattrs = list( id=eqn$eqnid
				, spec = 'phydyn.model.MatrixEquation'
				, type = 'birth'
				, origin = eqn$origin
				, destination = eqn$destination
		  	  )
		)
		.find_add_text( xml, searchstr = './/matrixeq', attrname='id', attrvalue=eqn$eqnid, addtext = eqn$eqn )
	} else if ( eqn$type == 'death' ){
		.find_add( xml, searchstr = './/popmodel', attrname = 'spec', attrvalue = 'phydyn.model.PopModelODE'
		  	  , addvalue = 'matrixeq' 
		  	  , addattrs = list( id=eqn$eqnid
				, spec = 'phydyn.model.MatrixEquation'
				, type = 'death'
				, origin = eqn$origin
		  	  )
		)
		.find_add_text( xml, searchstr = './/matrixeq', attrname='id', attrvalue=eqn$eqnid, addtext = eqn$eqn )
	} else if( eqn$type == 'migration' ){
		.find_add( xml, searchstr = './/popmodel', attrname = 'spec', attrvalue = 'phydyn.model.PopModelODE'
		  	  , addvalue = 'matrixeq' 
		  	  , addattrs = list( id=eqn$eqnid
				, spec = 'phydyn.model.MatrixEquation'
				, type = 'migration'
				, origin = eqn$origin
				, destination = eqn$destination
		  	  )
		)
		.find_add_text( xml, searchstr = './/matrixeq', attrname='id', attrvalue=eqn$eqnid, addtext = eqn$eqn )
	} else if( eqn$type == 'nondeme' ){
		.find_add( xml, searchstr = './/popmodel', attrname = 'spec', attrvalue = 'phydyn.model.PopModelODE'
		  	  , addvalue = 'matrixeq' 
		  	  , addattrs = list( id=eqn$eqnid
				, spec = 'phydyn.model.MatrixEquation'
				, type = 'nondeme'
				, origin = eqn$origin
		  	  )
		)
		.find_add_text( xml, searchstr = './/matrixeq', attrname='id', attrvalue=eqn$eqnid, addtext = eqn$eqn )
	} else if ( eqn$type == 'definition' ){
		.find_add( xml, searchstr = './/popmodel', attrname = 'spec', attrvalue = 'phydyn.model.PopModelODE'
		  	  , addvalue = 'definition' 
		  	  , addattrs = list( id=eqn$eqnid
				, spec = 'phydyn.model.Definition'
				, value = eqn$eqn  
		  	  )
		)
	} else{
		stop('Equation type not supported' )
		stop('Supported equation types: {EQNTYPES}' |> glue::glue() )
	}
}


.addlogger <- function(xml, parm)
{
	stopifnot( inherits(parm, 'phydynparameter') )
	.find_add( xml, searchstr = './/logger', attrname='id', attrvalue='tracelog'
		       , addvalue = 'log'
		       , addattrs = list( idref = parm$parmid )
		       )
	.find_add( xml, searchstr = './/logger', attrname='id', attrvalue='screenlog'
		       , addvalue = 'log'
		       , addattrs = list( idref = parm$parmid )
		       )
}

.addinitialvalue <- function( xml, parm )
{
	stopifnot( inherits(parm, 'phydynparameter') )
	if ( parm$name == 't0' ){ # special case 
		.find_add( xml, searchstr = './/popParams'
			  , attrname = 'id', attrvalue='initValues'
			  , addvalue = 'parameter'
			  , addattrs = list(
				  	    id = parm$parmid
				  	    , spec = 'parameter.RealParameter'
				  	    , name = 't0'
				  	    , estimate = 'false' 
				  	    )
			  )
		.find_add_text( xml, searchstr = './/parameter'
			       , attrname = 'id', attrvalue = parm$parmid
			       , addtext = as.character( parm$initial )
			       )
	} else if ( parm$estimate )
	{
		lbstr <- gsub( x = as.character( parm$lower ), pattern = 'Inf', replacement='Infinity')
		ubstr <- gsub( x = as.character( parm$upper), pattern = 'Inf', replacement='Infinity')
		.find_add( xml, searchstr = './/state', attrname = 'id', attrvalue = 'state'
		  	  , addvalue = 'parameter'
		  	  , addattrs = list(
		  		    	    id=parm$parmid
		  		    	    , spec = 'parameter.RealParameter'
		  		    	    , name = 'stateNode'
		  		    	    , lower = lbstr 
		  		    	    , upper = ubstr 
		  	  )
		  	  )
		.find_add_text(xml, searchstr = './/parameter', attrname='id'
			       , attrvalue=parm$parmid
			       , addtext = as.character( parm$initial ) 
		)
		if( parm$initial_condition_parameter ){ # an initial value for a state variable 
			.find_add( xml, searchstr = './/popParams'
				  , attrname = 'id', attrvalue='initValues'
				  , addvalue = 'initialValue'
				  , addattrs = list(
				  		    id = glue::glue( 'initialValue.{parm$parmid}')
				  		    , spec = 'phydyn.model.ParamValue'
				  		    , pname = parm$name
				  		    , pvalue = glue::glue('@{parm$parmid}')
				  		    )
				  )
		} else{ # a model parameter 
			.find_add( xml, searchstr = './/modelParams'
				  , attrname = 'id', attrvalue='rates'
				  , addvalue = 'param'
				  , addattrs = list(
				  		    id = glue::glue( 'modelParams.{parm$parmid}')
				  		    , spec = 'phydyn.model.ParamValue'
				  		    , pname = parm$name
				  		    , pvalue = glue::glue('@{parm$parmid}')
				  		    )
				  )
		}
	} else{ #not estimating 
		if ( parm$initial_condition_parameter )
		{
			.find_add( xml, searchstr = './/popParams'
			  	  , attrname = 'id', attrvalue='initValues'
			  	  , addvalue = 'initialValue'
			  	  , addattrs = list(
				  	    	    id = glue::glue( 'initialValue.{parm$parmid}')
				  	    	    , spec = 'phydyn.model.ParamValue'
				  	    	    , pname = parm$name
				  	    	    )
			  	  )
			.find_add( xml, searchstr = './/initialValue'
			  	  , attrname = 'id', attrvalue = glue::glue( 'initialValue.{parm$parmid}')
			  	  , addvalue = 'parameter'
			  	  , addattrs = list(
			  		    	    id = parm$parmid
			  		    	    , spec = 'parameter.RealParameter'
			  		    	    , estimate = 'false'
			  		    	    , name = 'pvalue'
			  		    	    , lower = "0.0"
			  		    	    )
			  	  )
			.find_add_text(xml, searchstr = './/parameter', attrname='id'
			       	       , attrvalue=parm$parmid
			       	       , addtext = as.character( parm$initial ) 
			)
		} else{
			.find_add( xml, searchstr = './/modelParams'
			  	  , attrname = 'id', attrvalue='rates'
			  	  , addvalue = 'param'
			  	  , addattrs = list(
				  	    	    id = glue::glue('modelParam.{parm$parmid}')
				  	    	    , spec = 'phydyn.model.ParamValue'
				  	    	    , pname = parm$name
				  	    	    )
			  	  )
			.find_add( xml, searchstr = './/param'
			  	  , attrname = 'id', attrvalue = glue::glue('modelParam.{parm$parmid}')
			  	  , addvalue = 'parameter'
			  	  , addattrs = list(
			  		    	    id = parm$parmid
			  		    	    , spec = 'parameter.RealParameter'
			  		    	    , estimate = 'false'
			  		    	    , name = 'pvalue'
			  		    	    )
			  	  )
			.find_add_text(xml, searchstr = './/parameter', attrname='id'
			       	       , attrvalue=parm$parmid
			       	       , addtext = as.character( parm$initial ) 
			)
		}
				

	}
}

.addoperator <- function( xml, parm )
{
	if (!parm$estimate )
		return( NULL )
	stopifnot( inherits(parm, 'phydynparameter') )
	oid <- glue::glue( '{parm$parmid}.operator' )
	xref <- glue::glue('@{parm$parmid}') 
	if ( parm$operator == 'realrw' ){
		stopifnot( 'windowSize' %in% names(parm) )
		.find_add( xml, searchstr = './/run' , attrname = 'id', attrvalue = 'mcmc'
			  , addvalue = 'operator'
			  , addattrs = list(id = oid, spec = 'RealRandomWalkOperator', parameter=xref, windowSize = as.character( parm$windowSize )
			  		    , weight = as.character(parm$weight)
			  		    , useGaussian=parm$useGaussian)
			  )
	} else if ( parm$operator == 'scale'){
		.find_add( xml, searchstr = './/run' , attrname = 'id', attrvalue = 'mcmc'
			  , addvalue = 'operator'
			  , addattrs = list(id = oid, spec = 'ScaleOperator', parameter=xref
			  		    , weight = as.character(parm$weight)
			  		    , scaleFactor = as.character( parm$scaleFactor )
			  	  )
			  )
	} else{
		stop('Supported operators: {OPERATORS}' |> glue::glue() )
	}
}

.addprior <- function( xml, parm )
{
	if (!parm$estimate )
		return( NULL )
	stopifnot( inherits(parm, 'phydynparameter') )
	xref <- glue::glue('@{parm$parmid}') 
	priorid <- glue::glue('{parm$parmid}.prior')
	.find_add( xml, searchstr = './/distribution', attrname = 'id', attrvalue = 'prior'
		  , addvalue='prior'
		  , addattrs = list( x = xref, name = 'distribution', id = priorid ))
	
	if ( parm$prior == 'lognormal' ){
		lognormid <- glue::glue('{priorid}.lognormal')
		.find_add( xml, searchstr = './/prior', attrname='id', attrvalue = priorid 
			  , addvalue = 'LogNormal' 
			  , addattrs = list(name = 'distr', id = lognormid) 
		)
		stopifnot( 'M' %in% names( parm ))
		stopifnot( 'S' %in% names( parm ))
		parmid1 <- glue::glue('{lognormid}.RealParameter.1')
		parmid2 <- glue::glue('{lognormid}.RealParameter.2')
		.find_add( xml, searchstr = './/LogNormal', attrname='id', attrvalue=lognormid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='M', estimate='false', spec = 'parameter.RealParameter', id = parmid1) #, text = parm$M 
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid1
			       , addtext = as.character( parm$M ))
		.find_add( xml, searchstr = './/LogNormal', attrname='id', attrvalue=lognormid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='S', estimate='false', spec = 'parameter.RealParameter', id = parmid2 ) #text = parm$S 
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid2
			       , addtext = as.character( parm$S ))
	} else if(parm$prior == 'normal' ){
		normid <- glue::glue('{priorid}.normal')
		.find_add( xml, searchstr = './/prior', attrname='id', attrvalue = priorid 
			  , addvalue = 'Normal' 
			  , addattrs = list(name = 'distr', id = normid) 
		)
		stopifnot( 'mean' %in% names( parm ))
		stopifnot( 'sigma' %in% names( parm ))
		parmid1 <- glue::glue('{normid}.RealParameter.1')
		parmid2 <- glue::glue('{normid}.RealParameter.2')
		.find_add( xml, searchstr = './/Normal', attrname='id', attrvalue=normid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='mean', estimate='false', spec = 'parameter.RealParameter', id = parmid1) #, text = parm$M 
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid1
			       , addtext = as.character( parm$mean ))
		.find_add( xml, searchstr = './/Normal', attrname='id', attrvalue=normid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='sigma', estimate='false', spec = 'parameter.RealParameter', id = parmid2 ) #text = parm$S 
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid2
			       , addtext = as.character( parm$sigma ))
	} else if ( parm$prior == 'exponential' ){
		expid <- glue::glue('{priorid}.exponential')
		.find_add( xml, searchstr = './/prior', attrname='id', attrvalue = priorid 
			  , addvalue = 'Exponential' 
			  , addattrs = list(name = 'distr', id = expid) 
		)
		stopifnot( 'mean' %in% names( parm ))
		parmid1 <- glue::glue('{expid}.RealParameter.1')
		.find_add( xml, searchstr = './/Exponential', attrname='id', attrvalue=expid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='mean', estimate='false', spec = 'parameter.RealParameter', id = parmid1)  
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid1
			       , addtext = as.character( parm$mean ))
	} else if( parm$prior == 'beta' ){
		betaid <- glue::glue('{priorid}.Beta')
		.find_add( xml, searchstr = './/prior', attrname='id', attrvalue = priorid 
			  , addvalue = 'Beta' 
			  , addattrs = list(name = 'distr', id = betaid) 
		)
		stopifnot( 'alpha' %in% names( parm ))
		stopifnot( 'beta' %in% names( parm ))
		parmid1 <- glue::glue('{normid}.RealParameter.1')
		parmid2 <- glue::glue('{normid}.RealParameter.2')
		.find_add( xml, searchstr = './/Beta', attrname='id', attrvalue=betaid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='alpha', estimate='false', spec = 'parameter.RealParameter', id = parmid1) #, text = parm$M 
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid1
			       , addtext = as.character( parm$alpha ))
		.find_add( xml, searchstr = './/Beta', attrname='id', attrvalue=betaid
			  , addvalue = 'parameter' 
			  , addattrs = list(name='beta', estimate='false', spec = 'parameter.RealParameter', id = parmid2 ) #text = parm$S 
			  )
		.find_add_text( xml, searchstr='.//parameter', attrname='id', attrvalue=parmid2
			       , addtext = as.character( parm$beta ))
	} else{
		stop('Supported priors: {PRIORS}' |> glue::glue() )
	}
}	

.add_traj_log <- function(xml, logEvery, fileName , popModel, pointFrequency="1", logrates = "all" )
{
	invisible( '
    	    <logger id="popTrajLog" spec="Logger" fileName="seir.algn.traj" logEvery="10000">
        	<log id="trajectoryLogger" spec="phydyn.loggers.TrajectoryLogger" pointFrequency="1" popModel="@pdseirmodel.t:algn" logrates="all"/>
    	    </logger>
	')
	.find_add( xml, searchstr = './/run' , attrname = 'id', attrvalue = 'mcmc'
	  	  , addvalue = 'logger'
	  	  , addattrs = list(id = 'popTrajLog'
	  		    	    , spec = 'Logger'
	  		    	    , fileName=fileName
	  		    	    , logEvery = logEvery
	  	  )
	  )
	.find_add( xml, searchstr = './/logger', attrname = 'id', attrvalue = 'popTrajLog'
		  , addvalue = 'log', addattrs = list( id='trajectoryLogger'
		  				      , spec = 'phydyn.loggers.TrajectoryLogger'
		  				      , pointFrequency = pointFrequency
		  				      , logrates = logrates 
		  				      , popModel = glue::glue( '@{popModel}' ) # @PhyDynPopModel
		  				      )
		  )
}


#' Configure a PhyDyn analysis and generate an XML file for BEAST analysis. 
#' 
#' This function takes a XML template file and includes PhyDyn model elements (equations and parameters). The resulting XML file can then be analysed in BEAST with the PhyDyn package. 
#' 
#' @param xmlfn The filename of the BEAST template XML. 
#' @param saveto Where the file will be saved 
#' @param t0parm A parameter (made using `confparm) for the time (numeric) of intialisation of the demographic or epidemiological simulation. This parameter may be estimated or fixed.
#' @param equations A list of `phydynequation` which can be made using `confeqn`
#' @param parameters A list of `phydynparameter` which can be made using `confparm`
#' @param coalescent_approximation Which approximation to the structured coalescent likelihood should be used? See \url{https://github.com/mrc-ide/PhyDyn/wiki/Structured-Tree-Likelihood}. Note an accuracy/speed tradeoff with different methods. 
#' @param integrationSteps Integer number of time steps to use when solving ODEs. Note an accuracy/speed tradeoff. 
#' @param penaltyAgtY If state variables specify a population size, this parameter can penalise parameter values which predict the number of extant lineages exceeding population size. Usually population sizes are "effective" and no penalty should be applied (default=0). 
#' @param useStateName If TRUE, the deme from which each sequence is sampled is specified in the sequence name using the suffix "_<deme name>". Alternatively, the demes of sampling must be specified in the template XML. 
#' @param traj_log_file The file to which model simulated trajectories will be saved 
#' @param traj_log_frequency Integer number of MCMC steps between logging model trajectories
#' @return  A PhyDyn model. The corresponding XML file is written to disk. 
#' @export 
#' @examples
#' \dontrun{
#' # Configure a model with exponential growth, birth(beta) & death(gamma), in two demes, with migration (mu) between demes
#' eqns <- list( 
#' 	confeqn( 'beta*gb', type = 'birth', origin = 'gb', destination = 'gb' )
#' 	, confeqn( 'beta*oth', type = 'birth', origin = 'oth', destination = 'oth' )
#' 	, confeqn( 'mu*oth', type = 'migration', origin = 'oth', destination = 'gb' )
#' 	, confeqn( 'mu*gb', type = 'migration', origin = 'gb', destination = 'oth' )
#' 	, confeqn( 'gamma*gb', type = 'death', origin = 'gb' )
#' 	, confeqn( 'gamma*oth', type = 'death', origin = 'oth' )
#' )
#' parms <- list(
#' 	      confparm( 'beta'
#' 	      	       , initial = 1.5*365*(1/5.5)
#' 	      	       , prior = 'lognormal'
#' 	      	       , operator = 'realrw' 
#' 	      	       , lower = .5*365*(1/5.5)
#' 	      	       , upper = 4*365*(1/5.5)
#' 	      	       , M = log(  1.5*365*(1/5.5)) , S = 1.0
#' 	      	       )
#' 	      , confparm( 'gamma'
#' 	      		 , initial = 365*(1/5.5)
#' 	      		 , estimate = FALSE 
#' 	      		 )
#' 	      , confparm( 'mu'
#' 	      		 , initial = 12 
#' 	      		 , prior = 'exponential'
#' 	      		 , operator = 'realrw' 
#' 	      		 , estimate = TRUE
#' 	      		 , lower = 1/4
#' 	      		 , upper = 52
#' 	      		 , mean = 12 
#' 	      		 )
#'		# estimate initial deme sizes
#' 	      , confparm( 'gb'
#' 	      		 , initial_condition_parameter = TRUE 
#' 	      		 , initial = 1e-2
#' 	      		 , estimate = TRUE
#' 	      		 , prior = 'exponential'
#' 	      		 , operator = 'realrw'
#' 	      		 , lower = 0
#' 	      		 , upper = 1e2 
#' 	      		 , mean = 10
#' 	      		 )
#' 	      , confparm( 'oth'
#' 	      		 , initial_condition_parameter = TRUE 
#' 	      		 , initial = 1e-2
#' 	      		 , estimate = TRUE 
#' 	      		 , prior = 'exponential'
#' 	      		 , operator = 'realrw'
#' 	      		 , lower = 0 
#' 	      		 , upper = 1e2
#' 	      		 , mean = 10 
#' 	      		 )
#' 	)
#' model <- config_phydyn(
#'  	 system.file( 'extdata', 'ba.2.86_algnWu-Hu-1.1_qc0_beast_template.xml', package = 'phydynbeast' )
#'  	  , saveto = 'ba.2.86_simodel0.xml'
#'  	  , t0 = confparm( 't0'
#' 	      		 , estimate = FALSE
#' 	      		 , initial = 2023.25
#' 	      		 )
#' 	      )
#'  	  , equations = eqns  
#'  	  , parameters = parms
#'  	  , coalescent_approximation = 'PL1'
#'  	  , integrationSteps = 100 
#'  	  , minP = 0.001
#'  	  , penaltyAgtY= 0 
#'  	  , useStateName = TRUE 
#'  	  , traj_log_file = 'simodel0-traj.tsv'
#'  	  , traj_log_frequency = 10000
#' )
#' }
config_phydyn <- function(
			  xmlfn
			  , saveto
			  , t0parm = confparm( name = 't0')
			  , equations = list()  
			  , parameters = list()
			  , coalescent_approximation = c('PL1', 'QL', 'PL2')
			  , integrationSteps = 200 
			  , minP = 0.001
			  , penaltyAgtY= 0 
			  , useStateName = TRUE 
			  , traj_log_file = 'traj.tsv'
			  , traj_log_frequency = 10000
			  )
{
	xml2::read_xml( xmlfn ) -> xml

	coalescent_approximation <- coalescent_approximation[1] 
	stopifnot( !is.na( coalescent_approximation ))
	stopifnot( is.character( coalescent_approximation ))
	stopifnot( coalescent_approximation %in% COAPPROX )
	stopifnot( is.numeric( minP ) & length( minP ) == 1)
	stopifnot( is.numeric( penaltyAgtY ) & length( penaltyAgtY )==1)
	stopifnot( is.numeric( integrationSteps ) & length( integrationSteps )==1)
	stopifnot( is.logical( useStateName ) & length( useStateName )==1)
	stopifnot( class( parameters ) == 'list' )
	stopifnot( class( equations ) == 'list' )

	parameters$t0 <- t0parm 
	
	bdmeqns <- equations[ sapply( equations, '[[', 'type') %in% c('birth', 'death', 'migration' )  ] 
	demes <- sapply( bdmeqns, '[[', 'origin' ) |> c( sapply( bdmeqns, '[[', 'destination' )) |> unique() |> setdiff( '' )
	ndeqns <- equations[ sapply( equations, '[[', 'type') %in% c('nondeme' )  ] 
	nondemes <-sapply( ndeqns, '[[', 'origin' ) |>  unique() |> setdiff( '' ) 
	parmnames <- sapply( parameters, '[[', 'name' ) 
	stopifnot( length( parmnames ) == length(unique(parmnames)))
	if (!('t0' %in% parmnames))
		stop('Must include parameter t0 which is the initial time of the simulation. This can be estimated or fixed. ')
	noic <- setdiff( demes, parmnames )
	if ( length( noic ) > 0 )
		stop('Must include initial conditions for all state variables (demes & nondemes) within *parameters*. The parameter name should match the deme name. These can be estimated or fixed. Missing initial conditions: {noic}' |> glue::glue() )
	
	
	# make STreeLikelihoodODE
	.find_add( xml, searchstr = './/distribution', attrname = 'id', attrvalue = 'prior'
		  , addvalue = 'distribution'
		  , addattrs = list(id='PhyDynModelLikelihood'
			, spec = 'phydyn.distribution.STreeLikelihoodODE'
			, equations = coalescent_approximation
			, minP = as.character( minP )
			, penaltyAgtY = as.character( penaltyAgtY  )
			, useStateName = ifelse( useStateName, 'true', 'false' )
		  )
	)

	# make popmodel 
	.find_add( xml, searchstr = './/distribution', attrname='id', attrvalue='PhyDynModelLikelihood'
		  , addvalue = 'popmodel'
		  , addattrs = list(id='PhyDynPopModel', spec='phydyn.model.PopModelODE' )
		  )


	## popParams
	.find_add( xml, searchstr = './/popmodel', attrname='id', attrvalue = 'PhyDynPopModel'
		  , addvalue = 'popParams'
		  , addattrs = list(id='initValues', spec = 'phydyn.model.TrajectoryParameters', integrationSteps = as.character(integrationSteps))
		  )
	## modelParams
	.find_add( xml, searchstr = './/popmodel', attrname='id', attrvalue = 'PhyDynPopModel'
		  , addvalue = 'modelParams'
		  , addattrs = list(id='rates', spec = 'phydyn.model.ModelParameters')
		  )
	
	#treeIntervals
	nodes0 <- xml2::xml_find_all(xml, './/tree' )
	treenode <- nodes0[ which( xml2::xml_attr( nodes0, 'spec' ) == 'beast.base.evolution.tree.Tree' ) ] 
	treeid <- xml2::xml_attr( treenode, 'id' )[1] 
	.find_add( xml, searchstr = './/distribution', attrname='id', attrvalue='PhyDynModelLikelihood'
		  , addvalue = 'treeIntervals'
		  , addattrs = list(id='STreeIntervals', spec='phydyn.evolution.tree.coalescent.STreeIntervals', tree=glue::glue('@{treeid}' ))
		  )


	# priors 
	for ( parm in parameters )
		.addprior ( xml, parm )

	# operators 
	for ( parm in parameters )
		.addoperator( xml, parm )

	# init cond 
	for( parm in parameters )
		.addinitialvalue( xml, parm )

	# loggers 
	for( parm in parameters )
		.addlogger( xml, parm )

	# equations 
	for ( eqn in equations )
		.addeqn( xml, eqn )

	# # trajectory log 
	.add_traj_log(xml, logEvery=as.character(round(traj_log_frequency))
	 	      , fileName=traj_log_file , popModel='PhyDynPopModel', pointFrequency="1", logrates = "all" )
	 

	# save 
	xml2::write_xml( xml, file = saveto )

	# return structure
	structure( list(xml = xml, equations = equations, parameters = parameters, demes = demes, nondemes = nondemes, parmnames = parmnames, path = saveto   )
		  , class = 'PhyDyn' )
}

#' @export 
print.PhyDyn <- function(x, ... )
{
	stopifnot( inherits( x, 'PhyDyn' ))
	with( x, {
		     cat(glue::glue( '# PhyDyn model with {length(demes)} demes and {length(parmnames)} parameters.' ))
		     cat('\n')
		     cat('\n')
		     cat(glue::glue( "{sapply(parameters,'[[','estimate') |> sum() } parameters to be estimated.") )
		     cat( '\n' )
		     cat( glue::glue( 'XML saved to {path} ' ) )
		     cat( '\n' )
		     cat( '\n' )
		     print( data.frame( `Demes` = demes ), row.names=TRUE)
		     cat('\n')
		     if ( length( nondemes ) > 0 ){
		     	     print( data.frame( `Non-deme state variable` = nondemes), row.names=TRUE)
		     	     cat('\n')
		     }
		     if( length( parmnames ) > 0 ){
		     	     print( data.frame( `Parameters` = parmnames ), row.names=TRUE )
		     	     cat('\n')
		     }
		     cat('\n')
		     cat( '# Equations\n' )
		     for( eqn in equations ){ print (eqn ); cat( '\n' ) }
		     cat('\n' )
		     cat('\n' )
		     cat( '# Parameters\n' )
		     for (parm in parameters) { print( parm ); cat( '\n' ) }
		     cat( '\n' )

		  })
}

# sample_priors <- function( phdynconf )
# {
# }


