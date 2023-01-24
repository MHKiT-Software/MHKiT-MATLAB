function R = tensorproduct( A , B , subscripts) %ind_R , A , ind_A , B , ind_B)
%TENSORPRODUCT Implementation of Einstein summation convention for
% multidimensional matlab arrays, where repeated indices sum over.
%
% Call:
%     R = tensorproduct(ind_R,A,ind_A,B,ind_B);
%
% Inputs:
%     ind_R : string with indices of the output
%         A : (non-empty) (multidimensional) matlab array
%     ind_A : string with indices of array A
%         B : (non-empty) (multidimensional) matlab array
%     ind_B : string with indices of array B
%
% Supports multiple:
% - Outer products
% - Inner products
% - Singleton dimensions
% - Pages
%
% Does not support (yet):
% - Self-contraction
%
% Example:
%     A = rand(5, 1,4,8);
%     B = rand(4,10,5  );
%     R = tensorproduct('jzgi',A,'gxki',B,'kjg'); %     Outer: 'i','j'
%                                                 %     Inner: 'k'
%                                                 %      Page: 'g'
%                                                 % Singleton: 'x','z'
%                                                 %  size(R) = [10,1,5,8]
%
% Note:
%   Provided indices in (ind_R), (ind_A) and (ind_B) are case-sensitive,    
%   i.e. 'a' and 'A' are different indices.
%
% Version compatibility:
%   This implementation makes use of Matlab built-in function pagemtimes,
%   introduced in Matlab version R2020b. To make use of this implementation
%   in previous Matlab releases, comment out the statement:
%
%     R = pagemtimes(A,B);
%
%   and uncomment the following lines:
%
%     R = zeros(size(A,1),size(B,2),size(A,3));
%     for i = 1:size(A,3)
%         R(:,:,i) = A(:,:,i) * B(:,:,i);
%     end
%   
%   While it will still work, performance (wall-clock time) might be 
%   degraded.
%
% Benchmarking:
%   The performance of tensorproduct can be assessed by calling the
%   function tensorproduct_benchmark.
 
% Copyright 2021 - by David Codony, PhD (dcodony@cimne.upc.edu)
% This software is distributed without any warranty.
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright
% notice is retained, and note is made of any changes that have
% been made.
%% Initialization
[ind_A , ind_B, ind_R] = parse_input_string(subscripts);
[A_permute,B_permute,R_permute,A_reshape,B_reshape,R_reshape] = ...
    tensorproduct_initialize(ind_A,ind_B,ind_R,size(A),size(B));
%% Computation
R = tensorproduct_core(A,B,...
       A_permute,B_permute,R_permute,A_reshape,B_reshape,R_reshape);
%% Parse input string
    function [ind_A , ind_B, ind_R] = parse_input_string(subscripts)
        s=split(subscripts,'->');
        %split input indices
        in=s{1};
        ind_R=s{2};
        in=split(in,',');        
        ind_A=in{1};
        ind_B=in{2};       
    end
end
