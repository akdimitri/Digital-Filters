function [R,p,w_o] = wiener( a, b, c, sigma_v1, sigma_v2 )

  R = zeros(4);
  p = zeros(4,1);
  w_o = zeros(4,1);
  
  A = zeros(4);
  A(1, :) = [1, a(1), a(2), a(3)];
  A(2, :) = [a(1), 1+a(2), a(3), 0];
  A(3, :) = [a(2), a(1)+a(3), 1, 0];
  A(4, :) = [a(3), a(2), a(1), 1];
  
  S = [sigma_v1, 0, 0, 0]';
  
  R_u = A \ S;
  R = toeplitz(R_u);
  
  C = zeros(4);
  C(1, :) = [1, c(1), c(2), c(3)];
  C(2, :) = [c(1), 1+c(2), c(3), 0];
  C(3, :) = [c(2), c(1)+c(3), 1, 0];
  C(4, :) = [c(3), c(2), c(1), 1];
  
  S = zeros(4, 1);
  R_x = C \ S;
  
  B = zeros(4);
  B(1, :) = [-b(1), -b(2), -b(3), -b(4)];
  B(2, :) = [-b(2), -b(1)-b(3), -b(4), 0];
  B(3, :) = [-b(3), -b(2)-b(4), -b(1), 0];
  B(4, :) = [-b(4), -b(3), -b(2), -b(1)];
  
  R_s = B * R_u;
  p = R_x + R_s;
  w_o = R \ p;
  
end


function [y,w_out] = lms(u, d, m, mu)
% LMS - Least Mean Square adaptation
%   
  
  persistent w
  persistent u_stream
  
  y = 0;
  w_out = zeros( m, 1 );
  
  if isempty(w)
      w = zeros(m, 1);
      u_stream = zeros(m-1, 1);
  end
  
  u_stream = [u_stream; u];
  U = flipud(u_stream(end-m+1:end)); 
  
  y = w' * U;
  e = d - y;
  
  w = w + mu * e * U; 
  
  w_out = w;
end


function [y,w_out] = nlms(u, d, m, mu, alpha)
% LMS - Least Mean Square adaptation
%   
  
  persistent w
  persistent u_stream
  
  y = 0;
  w_out = zeros( m, 1 );
  
  if isempty(w)
      w = zeros(m, 1);
      u_stream = zeros(m-1, 1);
  end
  
  u_stream = [u_stream; u];
  U = flipud(u_stream(end-m+1:end)); 
  
  y = w' * U;
  e = d - y;
  
  m_adapt = mu / (alpha + norm(U)^2);
  w = w + m_adapt * e * U; 
  
  w_out = w;
  
end



function [y,w_out] = rls(u, d, m, lambda, delta)
% LMS - RLS adaptation
%   
  
  persistent w
  persistent u_stream
  persistent P
  
  y = 0;
  w_out = zeros( m, 1 );
  
  if isempty(w)
      w = zeros(m, 1);
      u_stream = zeros(m-1, 1);
      P = 1/delta * eye(m);
  end
  
  u_stream = [u_stream; u];
  U = flipud(u_stream(end-m+1:end));
  
  y = w' * U;
  
  k = (1/lambda) * P * U / (1 + 1/lambda * U' * P * U);
  a = d - U' * w;
  w = w + k * a;
  P = 1 / lambda * (P - k * U' * P);
  
  w_out = w;
  
end