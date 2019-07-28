function [ y N err ] = FD_int_num( eta, j, tol, Nmax )

% Numerical integration of Fermi-Dirac integrals for order j > -1.
% Author: Raseong Kim (Purdue University)
% Date: September 29, 208
% Extended (composite) trapezoidal quadrature rule with variable
% transformation, x = exp( t - exp( t ) )
% Valid for eta ~< 15 with precision ~eps with 60~500 evaluations.
%
% Inputs
% eta: eta_F
% j: FD integral order
% tol: tolerance
% Nmax: number of iterations limit
%
% Note: When "eta" is an array, this function should be executed
% repeatedly for each component.
%
% Outputs
% y: value of FD integral (the "script F" defined by Blakemore (1982))
% N: number of iterations
% err: error
%
% For more information in Fermi-Dirac integrals, see:
% "Notes on Fermi-Dirac Integrals (3rd Edition)" by Raseong Kim and Mark
% Lundstrom at http://nanohub.org/resources/5475
%
% Reference
% [1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
% Numerical recipies: The art of scientific computing, 3rd Ed., Cambridge
% University Press, 2007.

for N = 1 : Nmax
    a = -4.5;                       % limits for t
    b = 5.0;
    t = linspace( a, b, N + 1 );    % generate intervals
    x = exp( t - exp( -t ) );
    f = x .* ( 1 + exp( -t ) ) .* x .^ j ./ ( 1+ exp( x - eta ) );
    y = trapz( t, f );
    
    if N > 1                        % test for convergence
        err = abs( y - y_old );
        if err < tol
            break;
        end
    end
    
    y_old = y;
end

if N == Nmax
    error( 'Increase the maximum number of iterations.')
end

y = y ./ gamma( j + 1 );


