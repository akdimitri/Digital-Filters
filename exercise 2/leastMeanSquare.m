function [y,w_out] = leastMeanSquare(u, d, m, mu)
% LEASTMEANSQUARE - Online LMS implementation
%
% SYNTAX
%
%   Y = LEASTMEANSQUARE( U, D, M, MU )
%   Y = LEASTMEANSQUARE( U, D, M, MU, D )
%   [Y,W] = LEASTMEANSQUARE( ... )
%
% INPUT
%
%   U   Input signal value              [scalar]
%   D   Desired signal value            [scalar]
%   M   Tap length                      [scalar]
%   MU  Convergence rate                [scalar]
%
% OUTPUT
%
%   Y   Output value of LMS             [scalar]
%   W   Weights (after update)          [m-vector]
%
% DESCRIPTION
%
%   LEASTMEANSQUARE computes one output of LMS algorithm
%


  %% PERSISTENT VaLUES ACROSS CALLS

  persistent w
  persistent u_stream


  %% INITIALIZATIONS

  % initialize w
  if isempty(w)
    w = zeros( m, 1 );
  end

  u_stream(end+1, 1) = u;


  %% COMPUTE NEW VALUE & UPDATE WEIGHTS

  if length( u_stream ) < m             % ----- stream not long enough

    y = u;

  else                                  % ----- calculate y

    % compute output
    y = w' * u_stream(end:-1:end-m+1);

    % get error
    e = d - y;

    if length( u_stream ) > m    % --- stream long enough for
                                 %     weight update

      % update weights
      w = w + mu * e * u_stream( end:-1:end-m+1 );

    end

  end


  %% UPDATE OUTPUT VALUES

  if nargout > 1
    w_out = w;
  end


end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%
% VERSION
%
%   0.1 - March 19, 2018
%
% CHANGELOG
%
%   0.1 (Mar 19, 2018) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------
