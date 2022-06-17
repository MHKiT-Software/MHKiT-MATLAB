function ste = ste_block_maxima_gev(block_maxima)
% Approximate the short-term extreme distribution using the block
% maxima method and the Generalized Extreme Value distribution.
% 
%     Parameters
%     ----------
%     block_maxima: array
%         Block maxima (i.e. largest peak in each block).
% 
%     Returns
%     -------
%     ste: TODO
%         Short-term extreme distribution.

ste_params = gevfit(block_maxima);


ste = 1;
end

