function [y, e, wo, wt] = noiseCancellation( u, d, m, mu, ws )
% NOISECANCELLATION - Implement noise cancellation for exercise 1
%   
% DESCRIPTION
%
%   [Y,E,WO,WT] = NOISECANCELLATION( U,D,M,MU,WS ) performs steepest
%   descent to locate the optimal Wiener coefficients.
%
% INPUT
%
%   U   Input signal u(n)               [n-by-1]
%   D   Desired signal d(n)             [n-by-1]
%   M   Taps (number of Wiener coeff)   [scalar]
%   MU  Steepest descent step           [scalar]
%   WS  Starting value for adaptation   [m-by-1]
%
% OUTPUT
%
%   Y   Output signal e(n)              [n-by-1]
%   E   Error signal d(n) - y(n)        [n-by-1]
%   WO  Optimal Wiener coefficients     [m-by-1]
%   WT  Coefficients at each step       [m-by-n]
%
  
  n = length(u);
  
  % TODO: implement function
  
  y  = zeros( n, 1 );
  e  = zeros( n, 1 );
  wo = zeros( m, 1 );
  wt = zeros( m, n );
  
   %Calculate Autocorrelation u = Ruu
    U = toeplitz( [u; zeros(1,m-1).'], [u(1) zeros(1, m-1)]);
    R_uu = U'*U/n;

    %Calculate Cross Correlation u and d
    d_temp = [d; zeros(m-1,1)];
    p_ud = U'*d_temp / n;


    %Wiener Solution 
    %wo = R_uu \ p_ud;


    wt = ws;
    wo = ws;
    for k=2:n % Iterate 1000 times the adaptation step
        wo = wo +mu *(p_ud-R_uu*wo); % Adaptation Equation ! Quite simple!
        wt =[wt wo]; % Wt records the evolution of vector W
    end 
    
    T = zeros(n,m);
    T(:,1) = u;
    for i = 2:m
        T(i:n,i) = u(1:n - i+1);
    end
    
    y = T*wo;

    e = d - y;
end


% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%
% VERSION
%
%   0.1 - March 15, 2019
%
