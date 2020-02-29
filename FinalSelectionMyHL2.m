%% 1.Result of NSGA2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
Eta=[0.6 0.8 1.0 1.2 1.4];

for DD=4
    for TT=1:5
        eta=Eta(TT);
        %% 1.1.Load dataset
        DataName=InputData(DD);
        [data,trueclus0,Ktrue] = LoadDataSet(DataName);
        [n,d]=size(data);
        rate=0.1;
        %% 1.2.Analyze the result of each time
        if DD==10||DD==11||DD==12||DD==28
            Evaluation_name=['Evaluation7.7\' DataName '_' num2str(eta) '_' num2str(rate) '.txt'];
        else
            Evaluation_name=['Evaluation7.7\' DataName '_' num2str(eta) '.txt'];
        end
        
%         delete(Evaluation_name);
        tt=5;
        Evaluation=zeros(tt,3);
        for time=1:tt
            disp([DataName '_iter1']);
            if DD==10||DD==11||DD==12||DD==28
                Clust_name=['Results\Clust_' DataName '_' num2str(eta)...
                    '_' num2str(rate)  '_iter' num2str(time) '.mat'];
                load(Clust_name,'Fit','Lab','Lab_all','overallTime','Index');
                trueclus=trueclus0(Index);
            else
                Clust_name=['Results\Clust_' DataName '_' num2str(eta) '_iter' num2str(time) '.mat'];
                load(Clust_name,'Fit','Lab','Lab_all','overallTime');
                trueclus=trueclus0;
            end            
            %% ！！！！！！！！！ Ensemble solutions ！！！！！！！！！ %%
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
            % [ ~,~,~,nmiscore,accurateScore,~,~,~,Fmeasure,adjrand] = evaluate(trueclus,Lab);
            Evaluation(time,:)=[AdjRandIndx,NMI,ACC];
            dlmwrite(Evaluation_name,[AdjRandIndx,NMI,ACC],'-append');
        end
        Mean_and_Std=[mean(Evaluation);std(Evaluation)]
    end
end
rmpath(genpath(cd));
