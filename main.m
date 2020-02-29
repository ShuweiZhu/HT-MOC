clc
clear
close all

cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));

%%%%%% Main loop %%%%%%
pop_size=50;
no_gen=50;
NormType=1;
Eta=[0.6 0.8 1.0 1.2 1.4];

for DD=3
    for t=2
        eta=Eta(t);        
       %% 1.Load dataset
        DataName=InputData(DD);
        [data,trueclus,Ktrue] = LoadDataSet(DataName);
        [n,d]=size(data);
        X=data;
        Kmax=fix(sqrt(size(data,1)));
        Kmax=min(Kmax,60);
        
        %normalization
        if NormType==1 % max min normalization
            for i=1:d
                MinV=min(X(:,i));
                MaxV=max(X(:,i));
                X(:,i)=X(:,i)-MinV;
                X(:,i) =X(:,i)/(MaxV-MinV);
            end
        end
        if NormType==2 % standard mean normalization
            for i=1:d
                X(:,i)=(X(:,i)-mean(X(:,i)))./std(X(:,i));
            end
        end
        %% Cancer
        if DD==16
            X=X+0.01*rand(n,d);
        end
        % Find the neighbors of the dataset with KD-tree
        tree = KDTreeSearcher(X);
        [Neighbor,~]=knnsearch( tree, X, 'k', 25);
        
        %% 2.Perform the algorithm for 10 times
        for time=1:10
            disp([DataName '_iter' num2str(time)]);
            startTimer=tic;
            [Lab,Lab_all,Fit,Seed,TrainTime,iter]=HL_MOCK(X,Kmax,Neighbor,pop_size,no_gen,eta);
            % runtime
            overallTime = toc(startTimer);
            % Save the cluster result
            Clust_name=['Clust_' DataName '_' num2str(eta) '_iter' num2str(time) '.mat'];
            save(Clust_name,'Lab','Lab_all','Fit','overallTime','TrainTime','Seed','iter');
            close all
        end
        
    end
end

rmpath(genpath(cd));
