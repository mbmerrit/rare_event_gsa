function rnd = randexp(m, n)
%  Random samples from the standard exponential distribution.
%
%  Here we use a simpler inversion formula as discussed in 
%
%	   L. Devroye, "Non-Uniform Random Variate Generation", 
%	   Springer-Verlag, 1986.

U = rand(m, n);
rnd = -log(U);
