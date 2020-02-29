function [pop,fit,fcc,fccSeed,Ks]=initialize_pop_original(Object,local_core,Kmax,L)
%% Obtain the object and parameters
% The Object
data=Object.data;
Locus=Object.locus;
Neighbor=Object.Neighbor;
NeighborSeed=Object.NeighborSeed;

% Global parameters
N=Object.N;
M=2;
L2=Object.L_seed;  % Neighbor size of seed points
Ns=Object.Ns; % Number of seed points

%% Initialization
[n,d]=size(data);
pop=zeros(N,2*Ns);
fit=zeros(N,M);
fcc=zeros(N,n);
fccSeed=zeros(N,Ns);
Ks=zeros(N,1);
if Kmax<60
    K_all=[2:Kmax];
else
    K_all=[2:60];
end
Kn=length(K_all);
K_all_temp=K_all;

%% loop for getting the population
p=1;
while p<=N
    % Set K
    if ~isempty(K_all_temp)
        Kn_temp=length(K_all_temp);
        Kset=K_all_temp(randi(Kn_temp));
        K_all_temp(K_all_temp==Kset)=[];
    else
        K_all_temp=K_all;
        Kset=K_all_temp(randi(Kn));
        K_all_temp(K_all_temp==Kset)=[];
    end
    % Remove interesting edges
    locus=Locus;
    t=1; Ct=Kset;
    while t<Ct
        locus_t=locus(t,1);
        j=locus(t,1);  i=locus(t,2);
        neighbor_i=NeighborSeed(i,1:(L2+1));
        neighbor_i(neighbor_i==j)=[];
        locus(t,1)=neighbor_i(randi(L2));
        if locus_t==1 % there exist two links between point 1 and another
            j=locus(t+1,1);
            locus(t+1,2)=locus(t,1); % the link keeps the same as the previous one
            t=t+1;Ct=Ct+1;
        end
        t=t+1;
    end
    
    % Obtain clustering labeling
    Lab_seed=Assignment_seed(locus);
    Lab=zeros(n,1);
    for ii=1:length(Lab_seed)
        Ind=local_core==ii;
        Lab(Ind)=Lab_seed(ii);
    end
    
    %% Construct the population
    % Compute two objectives
    [Variance,Conn,K]=ObjCompute(data,Lab,Neighbor,L);
    
    if K==1
        disp('Only one class');       
    else
        % obtain population
        pop(p,1:2*Ns)=reshape(locus(:,1:2),1,2*Ns); %% 这里总感觉不是很好
        % Obtain fitness
        fit(p,1:M)=[Variance Conn];
        fcc(p,1:n)=Lab';
        fccSeed(p,1:Ns)=Lab_seed';
        % cluster number
        Ks(p,1)=K;
        p=p+1;
    end
end

end

%% Compute two objectives
function [Variance,Conn,K]=ObjCompute(data,lab,Neighbor,L)

[n,d]=size(data);

% Compute cluster centers and the Variance
K=max(lab);
tc=zeros(K,d);
variance=zeros(K,1);
Nk=zeros(K,1);
for k=1:K
    Ind=find(lab==k);
    n_k=length(Ind);
    Nk(k,1)=n_k;
    tc(k,:)=mean(data(Ind,:));
    dist_k=sum((data(Ind,:)-ones(n_k,1)*tc(k,:)).^2,2);
    variance(k,1)=sum(dist_k);
end
Variance=sum(variance)/n;

% Compute connectedness (Conn)
Penality=zeros(n,1);
for i=1:n
    Ki=lab(i,1);
    neigh_i=Neighbor(i,2:(L+1));
    lab_L=lab(neigh_i,1);
    penality=find(lab_L~=Ki);
    if ~isempty(penality)
        Penality(i,1)=sum(1./penality);
    end
end
Conn=sum(Penality);

end
