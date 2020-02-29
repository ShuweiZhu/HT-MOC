function result=EnsembleCluster(Clust)

members = Clust';
[N,poolSize]=size(members);

% For each run, M base clusterings will be randomly drawn from the pool.
% Each row in bcIdx corresponds to an ensemble of M base clusterings.
pool=1:poolSize;
Selectpool=randperm(poolSize,10);

%% Construct the ensemble of M base clusterings
% baseCls is an N x M matrix, each row being a base clustering.
baseCls = members(:,Selectpool);

%% Get all clusters in the ensemble
[bcs, baseClsSegs] = getAllSegs(baseCls);

%% Compute ECI
disp('Compute ECI ... ');
ECI = computeECI(bcs, baseClsSegs, para_theta);

%% Compute LWCA
LWCA= computeLWCA(baseClsSegs, ECI, M);

%% Perform LWGP
disp('Run the LWGP algorithm ... ');
resultsLWGP = runLWGP(bcs, baseClsSegs, ECI, clsNums);
disp('--------------------------------------------------------------');

%% Perform LWEA
disp('Run the LWEA algorithm ... ');
resultsLWEA = runLWEA(LWCA, clsNums);

end