function [bcs, baseClsSegs] = getAllSegs(baseCls)

[N,M] = size(baseCls);
% n:    the number of data points.
% M:    the number of base clusterings.
% nCls:     the number of clusters (in all base clusterings).

bcs = baseCls;
nClsOrig = max(bcs,[],1);
C = cumsum(nClsOrig); 
bcs = bsxfun(@plus, bcs,[0 C(1:end-1)]);
nCls = nClsOrig(end)+C(end-1);
baseClsSegs=sparse(bcs(:),repmat([1:N]',1,M), 1,nCls,N); 
