![](https://travis-ci.org/houtanb/aimad.svg?branch=master)

AIM-AD
=====

Code created for:

Bastani, Houtan and Luca Guerrieri. "On the Application of Automatic Differentiation to the Likelihood Function for Dynamic General Equilibrium Models". **Advances in Automatic Differentiation**, Eds. Bischof, C.H., et. al. Berlin: Springer-Verlag, 2008. 303-314

Available as [Federal Reserve Board International Finance Discussion Paper #920](http://www.federalreserve.gov/pubs/ifdp/2008/920/)

**Abstract**: A key application of automatic differentiation (AD) is to
facilitate numerical optimization problems. Such problems are at the core of
many estimation techniques, including maximum likelihood. As one of the first
applications of AD in the field of economics, we used Tapenade to construct
derivatives for the likelihood function of any linear or linearized general
equilibrium model solved under the assumption of rational expectations.  We
view our main contribution as providing an important check on finite-difference
(FD) numerical derivatives. We also construct Monte Carlo experiments to
compare maximum-likelihood estimates obtained with and without the aid of
automatic derivatives. We find that the convergence rate of our optimization
algorithm can increase substantially when we use AD derivatives.

See also:
- [TAPENADE](http://www-sop.inria.fr/tropics/tapenade.html)
- [Anderson-Moore Algorithm](https://www.federalreserve.gov/econres/ama-index.htm)
