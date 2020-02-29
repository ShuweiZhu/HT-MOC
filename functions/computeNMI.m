%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
% This is a demo for the LWEA and LWGP algorithms. If you find this %
% code useful for your research, please cite the paper below.       %
%                                                                   %
% Dong Huang, Chang-Dong Wang, and Jian-Huang Lai.                  %
% "Locally weighted ensemble clustering."                           %
% IEEE Transactions on Cybernetics, to appear, 2017.                %
%                                                                   %
% The code has been tested in Matlab R2014a and Matlab R2015a on a  %
% workstation with Windows Server 2008 R2 64-bit.                   %
%                                                                   %
% https://www.researchgate.net/publication/316681928                %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function allScores = computeNMI(cluster_labels,true_labels)

allScores = zeros(size(cluster_labels,2),1);
for i = 1:size(cluster_labels,2)
    if min(cluster_labels(:,i))>0
        allScores(i) = nmi(cluster_labels(:,i), true_labels);
    end
end

function score = nmi(true_labels, cluster_labels)
%NMI Compute normalized mutual information (NMI) using the true and cluster
%   labels and return the value in 'score'.
%
%   Input    : true_labels    : N-by-1 vector containing true labels
%              cluster_labels : N-by-1 vector containing cluster labels
%
%   Output   : score          : NMI value
%
%   Author   : Shi Zhong, 2003.
%              http://www.cse.fau.edu/~zhong/software/textclust.zip
%
%   Modified : Wen-Yen Chen (wychen@alumni.cs.ucsb.edu)
%			   Chih-Jen Lin (cjlin@csie.ntu.edu.tw)

% Compute the confusion matrix 'cmat', where
%   col index is for true label (CAT),
%   row index is for cluster label (CLS).
n = length(true_labels);
cat = spconvert([(1:n)' true_labels ones(n,1)]);
cls = spconvert([(1:n)' cluster_labels ones(n,1)]);
cls = cls';
cmat = full(cls * cat);

n_i = sum(cmat, 1); % Total number of data for each true label (CAT), n_i
n_j = sum(cmat, 2); % Total number of data for each cluster label (CLS), n_j

% Calculate n*n_ij / n_i*n_j
[row, col] = size(cmat);
product = repmat(n_i, [row, 1]) .* repmat(n_j, [1, col]);
index = find(product > 0);
n = sum(cmat(:));
product(index) = (n*cmat(index)) ./ product(index);
% Sum up n_ij*log()
index = find(product > 0);
product(index) = log(product(index));
product = cmat .* product;
score = sum(product(:));
% Divide by sqrt( sum(n_i*log(n_i/n)) * sum(n_j*log(n_j/n)) )
index = find(n_i > 0);
n_i(index) = n_i(index) .* log(n_i(index)/n);
index = find(n_j > 0);
n_j(index) = n_j(index) .* log(n_j(index)/n);
denominator = sqrt(sum(n_i) * sum(n_j));

% Check if the denominator is zero
if denominator == 0
    score = 0;
else
    score = score / denominator;
end
