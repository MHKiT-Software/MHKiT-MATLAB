function R = tensorproduct_core(A,B,A_permute,B_permute,R_permute,A_reshape,B_reshape,R_reshape)
%% Rearrangement of input
A = reshape( permute( A , A_permute ) , A_reshape );
B = reshape( permute( B , B_permute ) , B_reshape );
%% Multiplication
if ~verLessThan('matlab','9.9')
    % ------ Matlab 2020b and newer releases: ------
    R = pagemtimes(A,B);
else
    % ------ Matlab 2020a and previous releases: ------
    R = zeros(size(A,1),size(B,2),size(A,3));
    for i = 1:size(A,3)
        R(:,:,i) = A(:,:,i) * B(:,:,i);
    end
end
%% Rearrangement of output
R = permute( reshape( R , R_reshape ) , R_permute );
end

