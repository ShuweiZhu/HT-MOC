clc
clear

cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
eta=0.8;

for DD=6
    %% 1.1.Load dataset
    DataName=InputData(DD);
    [data,trueclus,Ktrue] = LoadDataSet(DataName);
    [n,d]=size(data);
    
    %% 1.2.Analyze the result of each time
    Evaluation_name=['Evaluation7.7\' DataName '_' num2str(eta) '.txt'];
    delete(Evaluation_name);
    tt=10;
    Evaluation=zeros(tt,3);
    for time=1:tt
        disp([DataName '_iter1']);
        Clust_name=['Result6.28\Clust_' DataName '_' num2str(eta) '_iter' num2str(time) '.mat'];
        load(Clust_name,'Fit','Lab','Lab_all','overallTime');
        %% ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª Ensemble solutions ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª %%
        %%% Get all clusters in the ensemble
        N=size(Lab,1);
        Clust=Lab(randperm(N,fix(0.5*N)),:);
        %         Clust=Lab_all;
        [bcs, baseClsSegs] = getAllSegs(Clust');
        %%% Compute ECI
        disp('Compute ECI ... ');
        para_theta = 0.4; % Parameter
        ECI = computeECI(bcs, baseClsSegs, para_theta);
        %%% Perform LWGP
        disp('Run the LWGP algorithm ... ');
        resultsLWGP = runLWGP(bcs, baseClsSegs, ECI, Ktrue);
        %%% Evaluation the clustering (RI and ARI)
        fc=resultsLWGP;
        % RandIndx and AdjRandIndx
        [A,B]=calculateAB(trueclus,fc);
        [RandIndx,AdjRandIndx] = RandIndices(A,B);
        NMI=nmi(trueclus,fc);
        ACC=accuracy(trueclus,fc)/100;
        Evaluation(time,:)=[AdjRandIndx,NMI,ACC];
        dlmwrite(Evaluation_name,[AdjRandIndx,NMI,ACC],'-append');
    end
    Mean_and_Std=[mean(Evaluation);std(Evaluation)]
end
rmpath(genpath(cd));
