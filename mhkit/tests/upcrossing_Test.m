classdef upcrossing_Test < matlab.unittest.TestCase
    properties
        t
        signal
        zeroCrossApprox
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function setupTestClass(testCase)
            % Define time vector
            testCase.t = linspace(0, 4, 1000);

            % Define signal
            testCase.signal = testCase.exampleWaveform_(testCase.t);

            % Approximate zero crossings
            testCase.zeroCrossApprox = [0, 2.1, 3, 3.8];
        end

    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    % methods (Test)
    %     % Test methods
    %
    %     function unimplementedTest(testCase)
    %         testCase.verifyFail("Unimplemented test");
    %     end
    % end

    methods
        function signal = exampleWaveform_(~, t)
            % Generate a simple waveform form to analyse
            % This has been created to perform
            % a simple independent calcuation that
            % the mhkit functions can be tested against.
            A = [0.5, 0.6, 0.3];
            T = [3, 2, 1];
            w = 2 * pi ./ T;

            signal = zeros(size(t));
            for i = 1:length(A)
                signal = signal + A(i) * sin(w(i) * t);
            end
        end

        function [crests, troughs, heights, periods] = exampleAnalysis_(testCase, signal)
            % NB: This only works due to the construction
            % of our test signal. It is not suitable as
            % a general approach.

            % Gradient-based turning point analysis
            grad = diff(signal);

            % +1 to get the index at turning point
            turningPoints = find(grad(1:end-1) .* grad(2:end) < 0) + 1;

            crestInds = turningPoints(signal(turningPoints) > 0);
            troughInds = turningPoints(signal(turningPoints) < 0);

            crests = signal(crestInds);
            troughs = signal(troughInds);
            heights = crests - troughs;

            % Numerical zero-crossing solution
            zeroCross = zeros(size(testCase.zeroCrossApprox));
            for i = 1:length(testCase.zeroCrossApprox)
                zeroCross(i) = fzero(@(x) testCase.exampleWaveform_(x), ...
                    testCase.zeroCrossApprox(i));
            end

            periods = diff(zeroCross);
        end
    end

    methods (Test)
        %% Test functions without indices (inds)
        function test_peaks(testCase)
            [want, ~, ~, ~] = testCase.exampleAnalysis_(testCase.signal);
            got = uc_peaks(testCase.t, testCase.signal);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        function test_troughs(testCase)
            [~, want, ~, ~] = testCase.exampleAnalysis_(testCase.signal);
            got = uc_troughs(testCase.t, testCase.signal);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        function test_heights(testCase)
            [~, ~, want, ~] = testCase.exampleAnalysis_(testCase.signal);

            got = uc_heights(testCase.t, testCase.signal);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        function test_periods(testCase)
            [~, ~, ~, want] = testCase.exampleAnalysis_(testCase.signal);

            got = uc_periods(testCase.t, testCase.signal);

            testCase.verifyEqual(got, want, 'AbsTol', 2e-3);
        end

        function test_custom(testCase)
            [want, ~, ~, ~] = testCase.exampleAnalysis_(testCase.signal);

            % create a similar function to finding the peaks
            f = @(ind1, ind2) max(testCase.signal(ind1:ind2));

            got = uc_custom(testCase.t, testCase.signal, f);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        %% Test functions with indcies
        function test_peaks_with_inds(testCase)
            [want, ~, ~, ~] = testCase.exampleAnalysis_(testCase.signal);

            inds = upcrossing(testCase.t, testCase.signal);

            got = uc_peaks(testCase.t, testCase.signal, inds);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        function test_trough_with_inds(testCase)
            [~, want, ~, ~] = testCase.exampleAnalysis_(testCase.signal);

            inds = upcrossing(testCase.t, testCase.signal);

            got = uc_troughs(testCase.t, testCase.signal, inds);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        function test_heights_with_inds(testCase)
            [~, ~, want, ~] = testCase.exampleAnalysis_(testCase.signal);

            inds = upcrossing(testCase.t, testCase.signal);

            got = uc_heights(testCase.t, testCase.signal, inds);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

        function test_periods_with_inds(testCase)
            [~, ~, ~, want] = testCase.exampleAnalysis_(testCase.signal);
            inds = upcrossing(testCase.t, testCase.signal);

            got = uc_periods(testCase.t, testCase.signal,inds);

            testCase.verifyEqual(got, want, 'AbsTol', 2e-3);
        end

        function test_custom_with_inds(testCase)
            [want, ~, ~, ~] = testCase.exampleAnalysis_(testCase.signal);
            inds = upcrossing(testCase.t, testCase.signal);

            % create a similar function to finding the peaks
            f = @(ind1, ind2) max(testCase.signal(ind1:ind2));

            got = uc_custom(testCase.t, testCase.signal, f, inds);

            testCase.verifyEqual(got, want, 'AbsTol', 1e-3);
        end

    end
end
