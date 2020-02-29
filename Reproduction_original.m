function Pop=Reproduction_original(Object,local_core,pop,InterDist,Kmax,L)
% Obtain the object
data=Object.data;
T=Object.T;
Locus=Object.locus;
Neighbor=Object.Neighbor;
NeighborSeed=Object.NeighborSeed;

% Global parameters
N=Object.N;
L2=Object.L_seed;  % Neighbor size of seed points
Ns=Object.Ns; % Number of seed points

[n,d]=size(data);
m=1;
Kset=zeros(N,1);
Nrand=randperm(N);
N2=30; % population size for ensemble
pop2=pop(Nrand(1:N2)); % population for ensemble
Kset2=zeros(N2,1);

%% %%%%%%%%%%%%%% 1.1.Genetic operation %%%%%%%%%%%%%%%
Nrand2=randperm(N);
N1=N-N2; % population size for genetic clustering
pop1=pop(Nrand2(1:N1));

for p=1:2:N1-1
    %% a) Crossover and Mutation
    P1_temp=pop1(p).pos;
    P1=reshape(P1_temp,Ns,2);
    P2_temp=pop1(p+1).pos;
    P2=reshape(P2_temp,Ns,2);
    Pop1=P1(1:T,:);
    Pop2=P2(1:T,:);
    
    %----------Uniformly crossover----------
    Ind=rand(T,1)>0.5; % Pc=1;
    temp=Pop1(Ind);
    Pop1(Ind)=Pop2(Ind);
    Pop2(Ind)=temp;
    
    Pops=[Pop1 Pop2];
    %----------Neighborhood based mutation----------
    for pp=1:2
        pop_temp=Pops(:,2*pp-1:2*pp);
        NN_ij=zeros(1,T);
        for t=1:T
            Id=pop_temp(t,2);
            NN_ij(1,t)=find(NeighborSeed(Id,:)==pop_temp(t,1));
        end
        Pm=1/T*ones(1,T)+(NN_ij/T).^m; % Pm=1/T;
        Ind2=find(rand(1,T)<Pm);
        % replace the selected links
        if ~isempty(Ind2)
            for i=Ind2
                j=pop_temp(i,1);
                neighbor_i=NeighborSeed(pop_temp(i,2),1:(L2+1));
                if find(neighbor_i==j)
                    neighbor_i(neighbor_i==j)=[];
                end
                pop_temp(i,1)=neighbor_i(randi(L2));
            end
        else % (randomly select one link)
            i=randi(T);
            j=pop_temp(i,1);
            neighbor_i=NeighborSeed(pop_temp(i,2),1:(L2+1));
            if find(neighbor_i==j)
                neighbor_i(neighbor_i==j)=[];
            end
            pop_temp(i,1)=neighbor_i(randi(L2));
        end
        Pops(:,2*pp-1)=pop_temp(:,1);
    end
    
    %----------Obtain the offspring solutions----------
    Pop1=Pops(:,[1,2]); Pop2=Pops(:,[3,4]);
    P1(1:T,:)=Pop1; P2(1:T,:)=Pop2;
    pop1(p).pos=reshape(P1,2*Ns,1)';
    pop1(p+1).pos=reshape(P2,2*Ns,1)';
    
    %----------Assign clustering indices----------
    Lab_seed1=Assignment_seed(P1);
    Lab_seed2=Assignment_seed(P2);
    Labs_seed=[Lab_seed1 Lab_seed2];
    
    Lab1=zeros(n,1);
    Lab2=zeros(n,1);
    for ii=1:length(Lab_seed1)
        Ind=local_core==ii;
        Lab1(Ind)=Lab_seed1(ii);
    end
    for ii=1:length(Lab_seed2)
        Ind=local_core==ii;
        Lab2(Ind)=Lab_seed2(ii);
    end    
    Labs=[Lab1 Lab2];
    
    %% b) Update populations
    for pp=1:2
        lab=Labs(:,pp);
        if length(unique(Labs_seed(:,pp)))==1
            disp('Only one class');
            Kset(p+pp-1,1)=1;
        else
            % Update objectives
            [Variance,Conn,K]=ObjCompute(data,lab,Neighbor,L);
            pop1(p+pp-1).cost=[Variance,Conn];
            % Update labels
            pop1(p+pp-1).fcc=lab';
            pop1(p+pp-1).fccSeed=Labs_seed(:,pp)';
            pop1(p+pp-1).K=K;
            Kset(p+pp-1,1)=K;
        end
    end
end

%%% Delete solution with number of cluster as 1
pop1(Kset==1)=[];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ！！！！！！！！！ Ensemble solutions ！！！！！！！！！ %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clust=zeros(N2,Ns);
for p2=1:N2
    lab_seed=pop2(p2).fccSeed;
    Clust(p2,:)=lab_seed;
end
%%% The numbers of clusters.
% clsNums = [2:40];
Kmax=min(Kmax,60);
K_all=2:Kmax;
for p2=1:N2
    if ~isempty(K_all)
        Kn_temp=length(K_all);
        Kset=K_all(randi(Kn_temp));
        K_all(K_all==Kset)=[];
        clsNums(1,p2)=Kset;
    else
        K_all=2:Kmax;
        Kset=K_all(randi(Kmax-1));
        K_all(K_all==Kset)=[];
        clsNums(1,p2)=Kset;
    end
end
%%% Get all clusters in the ensemble
[bcs, baseClsSegs] = getAllSegs(Clust');
%%% Compute ECI
para_theta = 0.4; % Parameter
ECI = computeECI(bcs, baseClsSegs, para_theta);
%%% Perform LWGP
resultsLWGP = runLWGP(bcs, baseClsSegs, ECI, clsNums);

[~,bb]=sort(Locus(:,2));
Unfixed=Locus(1:T,2);
Fixed=Locus(T+1:end,2);
Unfixed=sort(Unfixed);
Fixed=sort(Fixed);
Locus2=Locus(bb,1:2);
for p2=1:N2
    Clust=resultsLWGP(:,p2);
    K=max(Clust);
    Link=[];
    for k=1:K
        ID_k=find(Clust==k);
        if length(ID_k)==1
            Link_k=Locus2(ID_k,:);
            Link=[Link;Link_k];
            continue
        elseif length(ID_k)==2
            Link_k=[ID_k(2),ID_k(1);ID_k(1),ID_k(2)];
            Link=[Link;Link_k];
            continue
        end
        DD_k=InterDist(ID_k,ID_k);
        % Compute MST links with Prim
        link=PrimMST(DD_k,size(DD_k,1));
        % add the link from 1 to another
        Ind_k=find(link(:,1)==1);
        Id=Ind_k(1);
        if Id<length(ID_k)-1
            link2=[link(1:Id,:); [link(Id,2),1]; link(Id+1:end,:)];
        else
            link2=[link(1:Id,:); [link(Id,2),1]];
        end
        Link_k=[ID_k(link2(:,1)),ID_k(link2(:,2))];
        Link=[Link;Link_k];
    end
    Pop_ens=Link;
    [~,bb]=sort(Pop_ens(:,2));
    Pop_ens=Pop_ens(bb,1:2);
    Pop_ens(Fixed,:)=Locus2(Fixed,:);
    Pop_ens=Pop_ens(Locus(:,2),:);
    Lab_ens=Assignment_seed(Pop_ens);
    Lab=zeros(n,1);
    for ii=1:length(Lab_ens)
        Ind=local_core==ii;
        Lab(Ind)=Lab_ens(ii);
    end  
    if length(unique(Lab_ens))==1
        disp('Only one class');
        Kset2(p2,1)=1;
    else
        % Update population
        pop2(p2).pos=reshape(Pop_ens,2*Ns,1)';
        % Update objectives
        [Variance,Conn,K]=ObjCompute(data,Lab,Neighbor,L);
        pop2(p2).cost=[Variance,Conn];
        % Update labels
        pop2(p2).fcc=Lab';
        pop2(p2).fccSeed=Lab_ens';
        pop2(p2).K=K;
        Kset2(p2,1)=K;
    end
end

%%% Delete solution with number of cluster as 1
pop2(Kset2==1)=[];

Pop=[pop1;pop2];

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

%% Construct the minimum spanning tree (MST) using Prim
function link=PrimMST(DD,n)
% Prim
A=DD;
A(A==0)=Inf;
P=zeros(1,n);
P(1,1)=1;
V=1:n;
V_P=V-P;
link=zeros(n-1,2);
k=1;
while k<n
    p=P(P~=0);
    v=V_P(V_P~=0);
    pv =min(min(A(p,v)));
    [x,y]=find(A==pv);
    for i=1:length(x)
        if  any(P==x(i)) && any(V_P==y(i))
            P(1,y(i))=y(i);
            V_P=V-P;
            link(k,:)=[x(i),y(i)]; %% y direct to x
            k = k+1;
            break;
        end
    end
end

end