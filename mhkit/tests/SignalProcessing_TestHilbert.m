classdef SignalProcessing_TestHilbert < matlab.unittest.TestCase
    % Tests for mhkit_hilbert using exact analytic-signal identities.
    %
    % For a signal with an integer number of periods in the record:
    %   cos(w t) -> exp(i w t)          (H[cos] = sin)
    %   sin(w t) -> -i exp(i w t)       (H[sin] = -cos)
    % and the Hilbert transform of a constant (DC) is zero.

    properties (TestParameter)
        n_samples = {63, 64, 1000, 1001};
    end

    methods (Test)
        function test_cosine_and_sine(testCase, n_samples)
            k = 5; % whole periods in the record
            t = (0:n_samples-1)' / n_samples;
            phase = 2 * pi * k * t;

            testCase.verifyEqual(mhkit_hilbert(cos(phase)), exp(1i * phase), 'AbsTol', 1e-12);
            testCase.verifyEqual(mhkit_hilbert(sin(phase)), -1i * exp(1i * phase), 'AbsTol', 1e-12);
        end

        function test_dc_has_zero_hilbert_transform(testCase, n_samples)
            x = 3.5 * ones(n_samples, 1);
            z = mhkit_hilbert(x);
            testCase.verifyEqual(real(z), x, 'AbsTol', 1e-12);
            testCase.verifyEqual(imag(z), zeros(n_samples, 1), 'AbsTol', 1e-12);
        end

        function test_real_part_is_input(testCase, n_samples)
            x = randn(RandStream('twister', 'Seed', n_samples), n_samples, 1);
            testCase.verifyEqual(real(mhkit_hilbert(x)), x, 'AbsTol', 1e-12);
        end

        function test_matrix_is_transformed_column_wise(testCase)
            n = 256;
            t = (0:n-1)' / n;
            x = [cos(2*pi*3*t), 2 * cos(2*pi*7*t), sin(2*pi*11*t)];
            z = mhkit_hilbert(x);
            testCase.verifySize(z, [n, 3]);
            for column = 1:3
                testCase.verifyEqual(z(:, column), mhkit_hilbert(x(:, column)), 'AbsTol', 1e-12);
            end
        end

        function test_row_vector_is_transformed_along_its_length(testCase)
            n = 128;
            phase = 2 * pi * 4 * (0:n-1) / n;
            z = mhkit_hilbert(cos(phase));
            testCase.verifySize(z, [1, n]);
            testCase.verifyEqual(z, exp(1i * phase), 'AbsTol', 1e-12);
        end

        function test_integer_input(testCase)
            x = int16([3; -7; 2; 9; -1; 0; 5]);
            testCase.verifyEqual(mhkit_hilbert(x), mhkit_hilbert(double(x)));
        end

        function test_short_inputs(testCase)
            % n = 1 and n = 2 have no negative frequencies to remove
            testCase.verifyEqual(mhkit_hilbert(4), complex(4, 0));
            testCase.verifyEqual(mhkit_hilbert([1; -2]), complex([1; -2], 0), 'AbsTol', 1e-15);
        end

        function test_invalid_input_errors(testCase)
            testCase.verifyError(@() mhkit_hilbert([1; 1i]), 'MATLAB:validators:mustBeReal');
            testCase.verifyError(@() mhkit_hilbert([]), 'MHKiT:mhkit_hilbert:InvalidInput');
        end
    end
end
