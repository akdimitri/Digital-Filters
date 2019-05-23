function y = fastFIR( u, w )

  n = length(u);
  m = length(w);

  v = [u(m:n); 0; u(1:m-1)];
  y = ifft( fft( v ) .* fft( [w; zeros(n-m+1,1)] ) );
  y = real( y( 1:n-m+1 ) );

end
