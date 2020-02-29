function [Lab,Lab_all,C,Seed,TrainTime,iter]=HL_MOCK(data,Kmax,Neighbor,pop_size,no_gen,eta)

%% Preprocessing
% initialize parameters
X=data;
[n,d]=size(data);
Pm.Xlth=n; Pm.Xwd=d;
Pm.LearnRate1=0.1; % learn rate for the following steps
Pm.ErrorBound=0.01; % learn rate for the following layers
Pm.B=4;
Pm.UpLimit=floor(eta*sqrt(n));

% Train the seed nodes
Timer=tic;
[Seed,DD,NeighborSeed,~,~,~,~]=GrowingTrainNode(Pm,X);
TrainTime = toc(Timer);

% Leaf nodes and the coresponding subset
LeafNodes=Seed.LeafNodes;
local_core=Seed.local_core;
% Subset=Seed.Subset;
% Core=Seed.Core;

% Obtain the links and draw the MST
[link,InterDist]=DrawMST(X,LeafNodes,DD);

%% Construct the object
L_seed=5;
N=pop_size;
Object.data=data;
Object.Ns=size(LeafNodes,1); %Number of seed points
Object.N=N;
Object.T=min(2*Kmax,size(LeafNodes,1));
Object.L_seed=L_seed;
Object.locus=link;
Object.Neighbor=Neighbor;
Object.NeighborSeed=NeighborSeed;

%% Initialize the population
if n<10000
    L=10;
elseif n<20000
    L=15;
else
    L=20;
end

if Kmax>Object.T
    Kmax=Object.T-1;
end
% Precomputation=PreCompute(Object,LeafNodes,Subset,L);
[init_pop,fit_pop,fcc_pop,fccSeed,Ks]=initialize_pop_original(Object,local_core,Kmax,L);

%% Construct the population
empty.pos=[];       % position
empty.cost=[];      % cost value
empty.fcc=[];       % clustering label
empty.fccSeed=[];  % Seed label
empty.K=[];         % cluster number
empty.rank=[];      % dominate ranking
empty.dcount=[];    % dominate count
empty.dset=[];      % dominated set
empty.rank=[];      % dominate ranking
empty.cdis=[];      % crowding distance

pop=repmat(empty,N,1);
for i=1:N
    pop(i).pos=init_pop(i,:);
    pop(i).cost=fit_pop(i,:);
    pop(i).fcc=fcc_pop(i,:);
    pop(i).fccSeed=fccSeed(i,:);
    pop(i).K=Ks(i,:);
end

[pop,F]=non_dominated_sorting(pop);
pop=calculated_crowding_distance(pop,F);
pop=sorting(pop);
[pop,F]=non_dominated_sorting(pop);

pf_t0=cat(1,pop(F{1}).cost);
count=0;
%% Main iteration %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter=1:no_gen
    %% Mating selection (generate mating pool)
    fprintf('The %d-th iteration\n',iter);
    pop2=pop;
    % index=bi_tournament_select(N);
    
    %% Reproduction (genetic operation and clustering)
    pop_new=Reproduction_original(Object,local_core,pop2,InterDist,Kmax,L);
    
    %% Evironmental selection
    [pop,F]=non_dominated_sorting([pop;pop_new]);
    pop=calculated_crowding_distance(pop,F);
    pop=sorting(pop);
    [pop,F]=non_dominated_sorting(pop);
    
    pop=pop(1:N);
    [pop,F]=non_dominated_sorting(pop);
    pop=calculated_crowding_distance(pop,F);
    pop=sorting(pop);
    [pop,F]=non_dominated_sorting(pop);
    
    %% Termination condition
    pf_t1=cat(1,pop(F{1}).cost);
    [~,ie12,ie21] = binary_epsilon_indicator(pf_t1, pf_t0);
    pf_t0=pf_t1;
    if ie12==1&&ie21==1
        count=count+1;
    else
        count=0;
    end
    if count>3
        break
    end
end

C=cat(1,pop(F{1}).cost);
Lab=cat(1,pop(F{1}).fcc);
Lab_all=cat(1,pop.fcc);
