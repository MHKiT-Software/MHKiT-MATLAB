var sourceData33 = {"FileName":"/Users/asimms/Desktop/Programming/mhkit_matlab_simms_dev/MHKiT-MATLAB-2/mhkit/dolfyn/tools/fill_time_gaps.m","RawFileContents":["function out = fill_time_gaps(epoch, sample_rate_hz)","% Fill gaps (NaN values) in the timeseries by simple linear","% interpolation.  The ends are extrapolated by stepping","% forward/backward by 1/sample_rate_hz.","","dt = 1 / sample_rate_hz;","% using fillgaps to interplolate values which uses nan so we convert","% negative values to nan","missing_data = epoch < 0;","epoch(missing_data) = nan;","out = fillgaps(epoch, false);","if isnan(epoch(1))","    i0 = find(~isnan(epoch),1) - 1;","    delta = (-i0:-1)*dt;","    epoch(1:i0) = epoch(i0) + delta;","end","if isnan(epoch(end))","    ie = find(~isnan(epoch),1,'last');","    delta = (1:numel(epoch)-ie) * dt;","    epoch(ie+1:end) = epoch(ie) + delta;","end","","end","",""],"CoverageDisplayDataPerLine":{"Function":{"LineNumber":1,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":52,"ContinuedLine":false},"Statement":[{"LineNumber":6,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":24,"ContinuedLine":false},{"LineNumber":9,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":25,"ContinuedLine":false},{"LineNumber":10,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":26,"ContinuedLine":false},{"LineNumber":11,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":29,"ContinuedLine":false},{"LineNumber":12,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":18,"ContinuedLine":false},{"LineNumber":13,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":35,"ContinuedLine":false},{"LineNumber":14,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":24,"ContinuedLine":false},{"LineNumber":15,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":36,"ContinuedLine":false},{"LineNumber":17,"Hits":2,"StartColumnNumbers":0,"EndColumnNumbers":20,"ContinuedLine":false},{"LineNumber":18,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":38,"ContinuedLine":false},{"LineNumber":19,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":37,"ContinuedLine":false},{"LineNumber":20,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":40,"ContinuedLine":false}]}}