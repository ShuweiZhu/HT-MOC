function Precomputation=PreCompute(Object,LeafNodes,Subset,L)
% Obtain the object
data=Object.data;
T=Object.T;
locus=Object.locus;
Neighbor=Object.Neighbor;

% Compute the L adjecent neighborhood
Neigh_L=Neighbor(:,2:(L+1));

[n,d]=size(data);
%% Remove interesting edges
t=1; Ct=T;
while t<Ct
    locus_t=locus(t,1);
    locus(t,1)=locus(t,2);  locus(t,3)=0; 
    if locus_t==1 % locus(t,1)==1 (there exist two links between point 1 and another)
        locus(t+1,1)=locus(t+1,2);  locus(t+1,3)=0; 
        t=t+1;Ct=Ct+1;
    end
    t=t+1;
end

% Obtain clustering labeling
Lab_seed=Assignment_seed(locus);
Lab=zeros(n,1);
for ii=1:length(Lab_seed)
    Ind=Subset(ii,:)==1;
    Lab(Ind)=Lab_seed(ii);
end
%% Compute two objectives
% Compute cluster centers and the Variance
K=max(Lab);
tc=zeros(K,d);
Pre_variance=zeros(K,1);
Pre_Nk=zeros(K,1);
Pre_Conn=cell(n,1);
CnnList=zeros(n,1);
for k=1:K
    Ind=find(Lab==k);
    n_k=length(Ind);
    Pre_Nk(k,1)=n_k;
    tc(k,:)=mean(data(Ind,:));
    dist_k=sum((data(Ind,:)-ones(n_k,1)*tc(k,:)).^2,2);  %% remove sqrt
    Pre_variance(k,1)=sum(dist_k);
end
Pre_centers=tc;
Variance=sum(Pre_variance)/n;

% Compute connectedness (Cnn)
Penality=zeros(n,1);
for i=1:n
    Ki=Lab(i,1);
    lab_L=Lab(Neigh_L(i,:),1);
    penality=find(lab_L~=Ki);
    if ~isempty(penality)
        CnnList(i,1)=1;
        Pre_Conn{i,1}=penality;
        Penality(i,1)=sum(1./penality);
    end
end
Cnn=sum(Penality);

% Objectives
Obj=[Variance,Cnn];

% Constructure of precomputation
Precomputation.variance=Pre_variance;
Precomputation.Conn=Pre_Conn;
Precomputation.CnnList=CnnList;
Precomputation.Nk=Pre_Nk;
Precomputation.centers=Pre_centers;
Precomputation.Lab=Lab;

end