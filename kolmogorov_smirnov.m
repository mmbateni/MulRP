## Copyright (C) 1995-2012 Kurt Hornik
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{pval}, @var{ks}] =} kolmogorov_smirnov_test (@var{x}, @var{dist}, @var{params}, @var{alt})
## Perform a Kolmogorov-Smirnov test of the null hypothesis that the
## sample @var{x} comes from the (continuous) distribution dist.  I.e.,
## if F and G are the CDFs corresponding to the sample and dist,
## respectively, then the null is that F == G.
##
## The optional argument @var{params} contains a list of parameters of
## @var{dist}.  For example, to test whether a sample @var{x} comes from
## a uniform distribution on [2,4], use
##
## @example
## kolmogorov_smirnov_test (x, "unif", 2, 4)
## @end example
##
## @noindent
## @var{dist} can be any string for which a function @var{dist_cdf}
## that calculates the CDF of distribution @var{dist} exists.
##
## With the optional argument string @var{alt}, the alternative of
## interest can be selected.  If @var{alt} is @code{"!="} or
## @code{"<>"}, the null is tested against the two-sided alternative F
## != G@.  In this case, the test statistic @var{ks} follows a two-sided
## Kolmogorov-Smirnov distribution.  If @var{alt} is @code{">"}, the
## one-sided alternative F > G is considered.  Similarly for @code{"<"},
## the one-sided alternative F > G is considered.  In this case, the
## test statistic @var{ks} has a one-sided Kolmogorov-Smirnov
## distribution.  The default is the two-sided case.
##
## The p-value of the test is returned in @var{pval}.
##
## If no output argument is given, the p-value is displayed.
## @end deftypefn

## Author: KH <Kurt.Hornik@wu-wien.ac.at>
## Description: One-sample Kolmogorov-Smirnov test

function [pval, ks] = kolmogorov_smirnov (x, distname, Params,alt)
names_dist={'exponential','gamma','extreme value' ,'lognormal','normal','weibull'};%'generalized extreme value' 'generalized pareto' 

  if (! isvector (x))
    error ("kolmogorov_smirnov_test: X must be a vector");
  endif

  n = length (x);
  s = sort (x);
##  alt  = "!=";
##  args{1} = s;
##  nvargs = numel (Params);
##  if (nvargs > 0)
##    if (ischar (Params{end}))
##      alt = Params{end};
##      args(2:nvargs) = Params(1:end-1);
##    else
##      args(2:nvargs+1) = Params;
##    endif
##  endif
%%z = reshape (feval (...cdf, args{:}), 1, n);
  
if strcmp(distname,names_dist{1});
        z = expcdf(s,Params);
end
if strcmp(distname,names_dist{2});
        z = gamcdf(s,Params(1),Params(2));
end
if strcmp(distname,names_dist{3});
        z = evcdf(s,Params(1),Params(2));
end
if strcmp(distname,names_dist{4});
        z = logncdf(s,Params(1),Params(2));   
end
if strcmp(distname,names_dist{5});
        z = normcdf(s,Params(1),Params(2));   
end
if strcmp(distname,names_dist{6});
        z = wblcdf(s,Params(1),Params(2));   
end
  if (strcmp (alt, "!=") || strcmp (alt, "<>"))
    ks  = sqrt (n) * max (max ([abs(z - (0:(n-1))/n); abs(z - (1:n)/n)]));
    pval = 1 - kolmogorov_smirnov_cdf (ks);
  elseif (strcmp (alt, ">"))
    ks   = sqrt (n) * max (max ([z - (0:(n-1))/n; z - (1:n)/n]));
    pval = exp (- 2 * ks^2);
  elseif (strcmp (alt, "<"))
    ks   = - sqrt (n) * min (min ([z - (0:(n-1))/n; z - (1:n)/n]));
    pval = exp (- 2 * ks^2);
  else
    error ("kolmogorov_smirnov_test: alternative %s not recognized", alt);
  endif
endfunction


## test for recognition of unifcdf function
%!assert (kolmogorov_smirnov_test (0:100, "unif", 0, 100), 1.0, eps)
## test for recognition of logistic_cdf function
%!assert (kolmogorov_smirnov_test (0:100, "logistic"), 0)
## test for  F < G
%!assert (kolmogorov_smirnov_test (50:100, "unif", 0, 50, "<"))

%!error kolmogorov_smirnov_test (1)
%!error <X must be a vector> kolmogorov_smirnov_test ({}, "unif", 2, 4)
%!error <no not_a_distcdf or not_a_dist_cdf function found>
%!  kolmogorov_smirnov_test (1, "not_a_dist");
%!error <alternative foo not recognized>
%!  kolmogorov_smirnov_test (1, "unif", 2, 4, "foo");

