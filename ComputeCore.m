function Core=ComputeCore(X,NumSeed,Subset,UpLimit)

KK2=max(10,ceil(sqrt(UpLimit)));
% KK3=ceil(sqrt(UpLimit));
Core=cell(NumSeed,3);

for ii=1:NumSeed
    % The k-nearest neighborhood (KNN)
    data_ii=find(Subset(ii,:)==1);
    data=X(data_ii,:);
    n=size(data,1);
    if n>KK2
        DD=pdist2(data,data);
        Nk=zeros(n,KK2);
        for i=1:n
            [~,Id]=sort(DD(i,:));
            Nk(i,1:KK2)=Id(2:KK2+1);
        end
        % The reverse nearest neighborhood (RNN)
        Rk_n=zeros(n,1);
        for i=1:n
            Rk_n(i,1)=length(find(Nk(:)==i));
        end
        % To determine core, boundary, and noise points
        Core_points=Rk_n>=KK2;
        Core_ii=data_ii(Core_points);
        Core{ii,1}=Core_ii;   % Core points
        core_all=union(find(Rk_n>=KK2),Nk(Core_points,:))';
        core_all=data_ii(core_all);
        Core{ii,2}=core_all;  % Core points and their reachable points
        Core{ii,3}=setdiff(data_ii,core_all); % Border points
    else
        Core{ii,1}=[];  Core{ii,2}=[]; 
        Core{ii,3}=data_ii;
    end
end
