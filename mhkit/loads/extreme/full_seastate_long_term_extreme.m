function lte = full_seastate_long_term_extreme(ste, weights)
% Return the long-term extreme distribution of a response of
% interest using the full sea state approach.
% 
% Parameters
% ----------
% ste: struct (list of ste distribtuions)
%     Short-term extreme distribution of the quantity of interest for
%     each sample sea state.
% weights: list[floats]
%     The weights from the full sea state sampling
% 
% Returns
% -------
% lte: object
%     Long-term extreme distribution.


assert(isstruct(ste), "ste must be a list like structure")
ste_names = fieldnames(ste);
assert(length(ste_names) > 1, "ste must be a list like structure")
assert(isa(ste.(ste_names{1}),"ste_peaks") || ...
    isa(ste.(ste_names{1}),"ste_block_maxima"), "Each item in ste must" + ...
    " be a ste object of type ste_block_maxima or ste_peaks")
assert(isvector(weights),'weights must be array')

lte = long_term_extreme(ste,weights);
end

