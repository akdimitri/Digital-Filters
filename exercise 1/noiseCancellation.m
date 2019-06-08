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
  
  %% 1. Toeplitz Method
  % Calculate Autocorrelation
  %U = toeplitz( [u; zeros(1,m-1).'], [u(1) zeros(1, m-1)]);
  %R_matrix = U'*U/n;


  %Calculate Cross Correlation u and d
  %d_temp = [d; zeros(m-1,1)];
  %p = U'*d_temp / n;
  
  %% 2. Xcorr Method
  % Calculate Autocorrelation using xcorr
  R_vector_xcorr = xcorr(u);
  median = (ceil((2*n-1)/2));
  R_vector_xcorr = R_vector_xcorr(median:( median + m-1));
  R_matrix_xcorr = toeplitz(R_vector_xcorr)/n;
  
  % Calculate cross-correlation using xcorr
  p_xcorr = xcorr(d, u);
  p_xcorr = p_xcorr(median:(median+m-1))/n;
%%  
  R_matrix = R_matrix_xcorr;
  p = p_xcorr;
    
  w = ws;
  wt = ws;
  for( i = 2: n)
      w = w + mu*( p - R_matrix*w);
      wt = [wt w];
  end
  
  wo = w;
  u = [u; zeros(m,1)];
  
  for i= m:n
      y(i) = wo'*u(i:-1:(i-m+1));
      e(i) = d(i) - y(i);
  end
  
  
end
