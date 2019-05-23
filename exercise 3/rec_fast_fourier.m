function [y, flop] = rec_fast_fourier(x)
    % A recursive implementation of the Fast Fourier Transform
    % 
    % input: 
    %   x: The input signal
    % output:
    %   y: The transformed signal
    %   flop: The flops of the computation
    %

    n = length(x);
    flop = 0;
   
    if n == 1
        y = x;
    else
        [y_top, flop_top] = rec_fast_fourier(x(1:2:n - 1));  % T(n/2) cost.
        [y_bottom, flop_bottom] = rec_fast_fourier(x(2:2:n));  % T(n/2) cost.
        d = exp(-2 * pi * 1i / n) .^ (0:n / 2 - 1);
        z = (d.') .* y_bottom;  % 6n/2 cost.
        y = [y_top + z; y_top - z];  % n complex additions => 2n cost.
        flop = 6*n/2 + 2*n + flop_top + flop_bottom;
    end
end