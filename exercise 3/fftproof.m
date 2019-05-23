function [yt, yb] = step1( x, w )
% split output to top and bottom
  
  n = length(x);
  
  k = (0:n-1)';
  
  yb = zeros(n/2,1);
  yt = zeros(n/2,1);
  
  for j = 0:n / 2 - 1
    yt(j + 1) = sum(w(n, j * k) .* x(k + 1));
  end

  for j = n / 2:n - 1
    yb(j + 1 - n / 2) = sum(w(n, j * k) .* x(k + 1));
  end
  
end

function [ybe, ybo, yte, yto] = step2( x, w )
% split each output to even and odd
  
  n = length(x);
  
  k = (0:n/2-1)';
  
  ybe = zeros(n/2,1);
  ybo = zeros(n/2,1);
  yte = zeros(n/2,1);
  yto = zeros(n/2,1);
  
  for j = 0:n / 2 - 1
    yte(j + 1) = sum(w(n, j * 2 * k) .* x(2 * k + 1));
    yto(j + 1) = sum(w(n, j * (2 * k + 1)) .* x(2 * k + 1 + 1));
  end
  
  for j = n / 2:n - 1
    ybe(j + 1 - n / 2) = sum(w(n, j * 2 * k) .* x(2 * k + 1));
    ybo(j + 1 - n / 2) = sum(w(n, j * (2 * k + 1)) .* x(2 * k + 1 + 1));
  end

end


function [xk1, xk2] = step3( x, w )
% apply w identity and finish the proof with
% y = [ xk1+xk2; xk1-xk2 ];
  
  n = length(x);
  
  xk1 = zeros(n/2,1);
  xk2 = zeros(n/2,1);
  
  j = (0:n / 2 - 1)';

  for k = 0:n / 2 - 1
    xk1(k + 1) = sum(x(2 * j + 1) .* w(n / 2, j * k));
  end

  for k = 0:n / 2 - 1
    xk2(k + 1) = sum(x(2 * j + 1 + 1) .* w(n / 2, j * k));
  end

  omega = diag(w(n, j));
  xk2 = omega * xk2;
end



