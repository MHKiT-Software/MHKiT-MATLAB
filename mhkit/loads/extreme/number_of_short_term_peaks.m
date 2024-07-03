function n_st = number_of_short_term_peaks(n, t, t_st)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Estimate the number of peaks in a specified period.
%
%     Parameters
%     ----------
%         n : int
%             Number of peaks in analyzed timeseries.
%         t : double
%             Length of time of analyzed timeseries.
%         t_st : double
%             Short-term period for which to estimate the number of peaks.
%
%      Returns
%      -------
%         n_st : double
%             Number of peaks in short term period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_st = n * t_st / t;