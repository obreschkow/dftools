#' Analytical distribution functions
#'
#' This function allows the user to call several pre-defined distribution functions (DFs). In the current implementation, they are all galaxy mass functions (i.e. a type of one-dimensional DFs), such as the Schechter function.
#'
#' @param x Vector of log-masses, typically \code{x = log10(M/Msun)}.
#' @param p Parameters of the analytical function. See argument \code{type}.
#' @param output Specifies what the function is doing. \code{'density'} evaluates the MF at \code{(x,p)}, \code{'npara'} returns the number of parameters of the analytical function, \code{'initial'} returns a vector of typical parameters \code{p} that are used as default initial values when fitting the MF, \code{'equation'} returns a string with the equation of the MF, \code{'function'} returns a function of \code{(x,p)}, \code{'names'} returns a list of strings specifying the parameter names in the expression format.
#' @param type Kind of MF: \code{'Schechter'} for Schechter function; \code{'PL'} for a simple power law; \code{'MRP'} for a Murray-Robotham-Power function, a 4-parameter extension of the Schechter function.
#'
#' @examples
#' # Evaluate and plot a Schechter function
#' x = seq(7,11,length=100)
#' mass = 10^x
#' parameters = c(-2,10,-1.5)
#' phi = dfmodel(x, parameters, type = 'Schechter')
#' plot(mass, phi, type='l', log='xy')
#'
#' @seealso \code{\link{dffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

dfmodel <- function(x = NULL, p = NULL, output = 'density', type = 'Schechter') {

  if (output=='npara') {
    # return number of parameters
    return(length(dfmodel(output = 'initial', type = type)))
  } else if (output=='function') {
    # return function(x,p)
    return(function(x,p) dfmodel(x, p, type = type))
  }

  # call model function
  if (type == 'Schechter') {
    f = .dfmodel.Schechter(x = x, p = p, output = output)
  } else if (type == 'MRP') {
    f = .dfmodel.MRP(x = x, p = p, output = output)
  } else if (type == 'PL') {
    f = .dfmodel.PL(x = x, p = p, output = output)
  } else {
    stop('dfmodel: unknown DF type.')
  }
  
  # finalize
  if (is.null(f)) {
    stop('dfmodel: unknown output type.')
  } else {
    return(f)
  }
  
}

.dfmodel.Schechter = function(x = NULL, p = NULL, output = 'density') {
  # returns the (unlogged) space density of galaxies of mass 10^x
  if (output == 'density') {
    mu = 10^(x-p[2])
    return(log(10)*10^p[1]*mu^(p[3]+1)*exp(-mu))
  } else if (output == 'initial') {
    return(c(-2,11,-1.3))
  } else if (output == 'equation') {
    return('dN/(dVdx) = log(10)*10^p[1]*mu^(p[3]+1)*exp(-mu), where mu=10^(x-p[2])')
  } else if (output == 'names') {
    names = c(expression('log'[10]*'('*phi['*']*')'),
              expression('log'[10]*'(M'['*']*')'),
              expression(alpha))
    return(names)
  }
  return(NULL)
}

.dfmodel.MRP = function(x = NULL, p = NULL, output = 'density') {
  # returns the (unlogged) space density of galaxies of mass 10^x
  # p[1]=normalization, p[2]=Ms, p[3]=alpha, p[4]=beta
  if (output == 'density') {
    mu = 10^(x-p[2])
    return(log(10)*10^p[1]*p[4]*mu^(p[3]+1)*exp(-mu^abs(p[4])))
  } else if (output == 'initial') {
    return(c(-2,11,-1,1))
  } else if (output == 'equation') {
    return('dN/(dVdx) = log(10)*p[4]*10^p[1]*mu^(p[3]+1)*exp(-mu^p[4]), where mu=10^(x-p[2])')
  } else if (output == 'names') {
    names = c(expression('log'[10]*'('*phi['*']*')'),
              expression('log'[10]*'(M'['*']*')'),
              expression(alpha),
              expression(beta))
    return(names)
  }
  return(NULL)
}

.dfmodel.PL = function(x = NULL, p = NULL, sigma = NULL, output = 'density') {
  if (output == 'density') {
    return(10^p[1]*(10^x)^p[2])
  } else if (output == 'initial') {
    return(c(2,-1))
  } else if (output == 'equation') {
    return('dN/(dVdx) = 10^p[1]*(10^x)^p[2]')
  } else if (output == 'names') {
    names = c(expression('log'[10]*'(A)'),
              expression(alpha))
    return(names)
  }
  return(NULL)
}
