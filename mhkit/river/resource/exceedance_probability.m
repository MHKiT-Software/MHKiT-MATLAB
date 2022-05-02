function F=exceedance_probability(Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the exceedance probability
%    
% Parameters
% ----------
%     Q : Array
%       Discharge data [m3/s]
%         
% Returns   
% -------
%     F : Array
%       Exceedance probability [unitless] 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we make a sorted array of discharge values
sort_q = sort(Q,'descend');
% Next we make a ranking array
rank = [1:numel(sort_q)];
% Now each value in Q is assigned a rank where rank is the same as pandas
% ranking algorithm with the 'max' option. For example in the data set 
% [1, 2, 2, 3] the ranks would be [ 1, 2, 3, 4] and the rank of 2 would be
% 3 because 3 is the maximum rank that a 2 has in the total set. 
assigned_rank = arrayfun(@(x) rank(find(sort_q==x,1,'last')), Q);

F = 100 * (assigned_rank / (numel(Q)+1));






