%% 1.Result of NSGA2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
Eta=[0.6 0.8 1.0 1.2 1.4];

for DD=20
    for tt=1
        eta=Eta(tt);
        %% 1.1.Load dataset
        DataName=InputData(DD);
        [data,trueclus0,Ktrue] = LoadDataSet(DataName);
        [n,d]=size(data);
        
        for time=1:5
            disp([DataName '_iter1']);
            if DD==10||DD==11||DD==12||DD==28
                rate=0.4;eta=1;
                Clust_name=['Result7.7\Clust_' DataName '_' num2str(eta)...
                    '_' num2str(rate)  '_iter' num2str(time) '.mat'];
                load(Clust_name,'Fit','Lab','Lab_all','overallTime','Seed','Index');
                trueclus=trueclus0(Index);
            else
                Clust_name=['Result7.7\Clust_' DataName '_' num2str(eta) '_iter' num2str(time) '.mat'];
                load(Clust_name,'Fit','Lab','Lab_all','overallTime','Seed');
                trueclus=trueclus0;
            end
            [N,n]=size(Lab_all);
            Clust=Lab_all;
            %% ！！！！！！！！！ Ensemble solutions ！！！！！！！！！ %%
            for s=[0.2,0.4,0.6,0.8,1.0]
                Npool=fix(s*N);
                %%% Compute graph-based short path distance
                LeafNodes=Seed.LeafNodes;
                NeighborSeed=Seed.NeighborSeed;
                supk=5;
                short_path=GraphDist(LeafNodes,NeighborSeed,supk);
                
                %%% Get all clusters in the ensemble
                Selectpool=[];
                for p=1:N
                    fc=Clust(p,:);
                    K=length(unique(fc));
                    local_core=zeros(n,1);
                    Subset=Seed.Subset;
                    for i = 1:n
                        local_core(i,1)=find(Subset(:,i)==1);
                    end
                    LCCV(p)=computeMySWC(n,fc,K,LeafNodes,short_path,local_core);
                end
                [aa,bb]=sort(LCCV,'descend');
                Qual=LCCV'./max(LCCV);
                Id=bb(1);
                Selectpool=Id;
                pool=1:N;
                pool(Id)=[];Qual(Id)=[];
                
                % Sampling for a subset
                if n<1000
                    members=Clust;
                else
                    members=Clust(:,randperm(n,1000));
                end
                % Select the remaining subset
                for p=2:Npool
                    Sim=zeros(N-p+1,p-1);
                    for t=1:(N-p+1)
                        Sim(t,:)=computeNMI(members(:,Selectpool),members(:,pool(t)));
                    end
                    Diversity=sum(1-Sim,2)/(p-1);
                    [~,Ind2]=max(Qual./max(Qual)+Diversity./max(Diversity));
                    Selectpool=[Selectpool,pool(Ind2)];
                    pool(Ind2)=[];Qual(Ind2)=[];
                end
                Clust_ens=Clust(Selectpool,:);
                [bcs, baseClsSegs] = getAllSegs(Clust_ens');
                %%% Compute ECI
                disp('Compute ECI ... ');
                para_theta = 0.4; % Parameter
                ECI = computeECI(bcs, baseClsSegs, para_theta);
                %%% Perform LWGP
                disp('Run the LWGP algorithm ... ');
                resultsLWGP = runLWGP(bcs, baseClsSegs, ECI, Ktrue);
                fc=resultsLWGP;
                
                %% Evaluation the clustering (RI and ARI)
                % RandIndx and AdjRandIndx
                [A,B]=calculateAB(trueclus,fc);
                [RandIndx,AdjRandIndx] = RandIndices(A,B);
                NMI=nmi(trueclus,fc);
                ACC=accuracy(trueclus,fc)/100;
                % [ ~,~,~,nmiscore,accurateScore,~,~,~,Fmeasure,adjrand] = evaluate(trueclus,Lab);
                %% 1.2.Analyze the result of each time
                if DD==10||DD==11||DD==12||DD==28
                    Evaluation_name=['EvaluationCVI\' DataName '_' num2str(eta) '_' num2str(Npool) '_' num2str(rate) '.txt'];
                else
                    Evaluation_name=['EvaluationCVI\' DataName '_' num2str(eta) '_' num2str(Npool) '.txt'];
                end
                % delete(Evaluation_name);
                % Evaluation(time,:)=[AdjRandIndx,NMI,ACC];
                dlmwrite(Evaluation_name,[AdjRandIndx,NMI,ACC],'-append');
            end
            % Mean_and_Std=[mean(Evaluation);std(Evaluation)]
        end
    end
end

rmpath(genpath(cd));
