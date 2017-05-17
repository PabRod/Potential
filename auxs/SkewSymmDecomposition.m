function [skew, symm] = SkewSymmDecomposition(x)
%SKEWSYMMDECOMPOSITION Returns the symmetric and skew matrix whose direct 
% sum is equal to x
%
%   Any square matrix can be expressed uniquely as the sum fo a symmetric 
% and a skew matrix

%% Input control
nRows = size(x, 1);
nCols = size(x, 2);
if(nRows ~= nCols)
    msgId = 'SkewSymmDecomposition:WrongArgin';
    errMsg = 'x should be a square matrix';
    error(msgId, errMsg);
end

%% Algorithm
symm = (x + x')/2;
skew = (x - x')/2;
% Note that x = symm + skew

end