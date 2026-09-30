function analytic_signal = mhkit_hilbert(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the analytic signal of a real signal using the Hilbert transform
%
% Native MATLAB implementation that needs neither the Signal Processing
% Toolbox nor Python. It returns the analytic signal z = x + i*H[x], where
% H is the Hilbert transform, which is the same quantity returned by the
% Signal Processing Toolbox hilbert function and by scipy.signal.hilbert.
% The algorithm matches scipy.signal.hilbert in SciPy v1.16.2:
% https://github.com/scipy/scipy/blob/b1296b9b4393e251511fe8fdd3e58c22a1124899/scipy/signal/_signaltools.py#L2476-L2594
%
% Parameters
% ------------
% x : double
%   Real, non-empty signal. A vector is transformed along its length. A
%   matrix is transformed column by column, one time series per column.
%
% Returns
% ---------
% analytic_signal : complex double
%   Analytic signal, same size as x. real(analytic_signal) equals x and
%   imag(analytic_signal) is the Hilbert transform of x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        x {mustBeNumeric, mustBeReal}
    end

    % Transform a row vector along its length, like hilbert() and scipy
    [x, was_row] = mhkit_standardize_user_input_to_column_vectors(double(x), ...
        'function_name', 'mhkit_hilbert');

    n = size(x, 1);
    spectrum = fft(x, [], 1);

    h = zeros(n, 1);
    if mod(n, 2) == 0
        h([1, n / 2 + 1]) = 1;
        h(2:n / 2) = 2;
    else
        h(1) = 1;
        h(2:(n + 1) / 2) = 2;
    end

    analytic_signal = ifft(spectrum .* h, [], 1);

    % Always return complex, as scipy does, even when the imaginary part is
    % exactly zero (for example n = 1 or n = 2)
    if isreal(analytic_signal)
        analytic_signal = complex(analytic_signal);
    end

    analytic_signal = mhkit_restore_column_vectors_to_user_input(analytic_signal, was_row);

end
