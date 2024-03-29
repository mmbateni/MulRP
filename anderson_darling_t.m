## [q,Asq,info] = anderson_darling_t(x,uniform|normal|exponential)		
##		
## Test the hypothesis that x is selected from the given distribution		
## using the Anderson-Darling test.  If the returned q is small, reject		
## the hypothesis at the q*100% level.		
##		
## The Anderson-Darling A^2 statistic is calculated as follows:		
##		
##    A^2_n = -n - \sum_{i=1}^n (2i-1)/n log(z_i (1-z_{n-i+1}))		
##		
## where z_i is the ordered position of the x's in the CDF of the 		
## distribution.  Unlike the Kolmogorov-Smirnov statistic, the 		
## Anderson-Darling statistic is sensitive to the tails of the 		
## distribution.		
##		
## For 'normal' and 'exponential' distributions, estimate the 		
## distribution parameters from the data, convert the values 		
## to CDF values, and compare the result to tabluated critical 		
## values.  This includes an correction for small n which 		
## works well enough for n >= 8, but less so from smaller n.  The		
## returned info.Asq_corrected contains the adjusted statistic.		
##		
## For 'uniform', assume the values are uniformly distributed		
## in (0,1), compute A^2 and return the corresponding p value from		
## 1-anderson_darling_cdf(A^2,n).		
## 		
## If you are selecting from a known distribution, convert your 		
## values into CDF values for the distribution and use 'uniform'.		
## Do not use 'uniform' if the distribution parameters are estimated 		
## from the data itself, as this sharply biases the A^2 statistic		
## toward smaller values.		
##		
## [1] Stephens, MA; (1986), "Tests based on EDF statistics", in		
## D'Agostino, RB; Stephens, MA; (eds.) Goodness-of-fit Techinques.		
## New York: Dekker.		
##		
## Author: Paul Kienzle		
function [q,Asq,info] = anderson_darling_t(x,dist)

  if size(x,1) == 1, x=x(:); end
  x = sort(x);
  n = size(x,1);
  use_cdf = 0;
  # Compute adjustment and critical values to use for stats.
  switch dist
    case 'normal',
       # This expression for adj is used in R.
       # Note that the values from NIST dataplot don't work nearly as well.
       adj = 1 + (.75 + 2.25/n)/n;
       qvals = [ 0.1, 0.05, 0.025, 0.01 ];
       Acrit = [ 0.631, 0.752, 0.873, 1.035];
       x = stdnormal_cdf(zscore(x));
 
     case 'uniform',
       ## Put invalid data at the limits of the distribution
       ## This will drive the statistic to infinity.
       x(x<0) = 0;
       x(x>1) = 1;
       adj = 1.;
       qvals = [ 0.1, 0.05, 0.025, 0.01 ];
       Acrit = [ 1.933, 2.492, 3.070, 3.857 ];
       use_cdf = 1;
 
     case 'XXXweibull',
       adj = 1 + 0.2/sqrt(n);
       qvals = [ 0.1, 0.05, 0.025, 0.01 ];
       Acrit = [ 0.637, 0.757, 0.877, 1.038];
       ## XXX FIXME XXX how to fit alpha and sigma?
       x = weibull_cdf(x,ones(n,1)*alpha,ones(n,1)*sigma);
 
     case 'exponential',
       adj = 1 + 0.6/n;
       qvals = [ 0.1, 0.05, 0.025, 0.01 ];
       # Critical values depend on n.  Choose the appropriate critical set.
       # These values come from NIST dataplot/src/dp8.f.
       Acritn = [0, 1.022, 1.265, 1.515, 1.888
 	    11, 1.045, 1.300, 1.556, 1.927;
 	    21, 1.062, 1.323, 1.582, 1.945;
 	    51, 1.070, 1.330, 1.595, 1.951;
 	    101, 1.078, 1.341, 1.606, 1.957;
 	    ];
      # FIXME: consider interpolating in the critical value table.
      Acrit = Acritn(lookup(Acritn(:,1),n),2:5);
      lambda = 1./mean(x);  # exponential parameter estimation
      x = exponential_cdf(x,ones(n,1)*lambda);
    otherwise
      # FIXME consider implementing more of distributions; a number
      # of them are defined in NIST dataplot/src/dp8.f.
      error("Anderson-Darling test for %s not implemented", dist);
  end

  if any(x<0 | x>1)
     error('Anderson-Darling test requires data in CDF form');
  endif

  i = [1:n]'*ones(1,size(x,2));
  Asq = -n - sum( (2*i-1) .* (log(x) + log(1-x(n:-1:1,:))) )/n;
  # Lookup adjusted critical value in the cdf (if uniform) or in the
  # the critical table.
  if use_cdf
    q = 1-anderson_darling_cdf(Asq*adj, n);
  else
    idx = lookup([-Inf,Acrit],Asq*adj);
    q = [1,qvals](idx); 
  endif
  if nargout > 2,
    info.Asq = Asq;
    info.Asq_corrected = Asq*adj;
    info.Asq_critical = [100*(1-qvals); Acrit]';
    info.p = 1-q;
    info.p_is_precise = use_cdf;
  endif