clc
clear

cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));

for DD=11
    %% 1.Load dataset
    DataName=InputData(DD);
    [data,trueclus0,Ktrue] = LoadDataSet(DataName);
    [n,d]=size(data);
    eta=1;
    for rate=[0.1 0.2]
    %% 2.Analyze the result of each time
    for time=1:3
        disp([DataName]);
        if DD==10||DD==11||DD==12||DD==28
            
            Clust_name=['Result7.7\Clust_' DataName '_' num2str(eta)...
                '_' num2str(rate)  '_iter' num2str(time) '.mat'];
            load(Clust_name,'Fit','Lab','Lab_all','overallTime','Index');
            trueclus=trueclus0(Index);
        else
            Clust_name=['Result7.7\Clust_' DataName '_' num2str(eta) '_iter' num2str(time) '.mat'];
            load(Clust_name,'Fit','Lab','Lab_all','overallTime');
            trueclus=trueclus0;
        end
        N=size(Lab,1);
        
        if DD==10||DD==11||DD==12||DD==28
            Evaluation_name=['Evaluation\' DataName '_' num2str(eta) '_' num2str(rate) ...
                             '_iter' num2str(time) '.txt'];
        else
            Evaluation_name=['Evaluation\' DataName '_' num2str(eta) '.txt'];
        end      
        delete(Evaluation_name);
        tic
        for t=1:N
            fc=Lab(t,:);
            K=length(unique(fc));
            % RandIndx and AdjRandIndx
            [A,B]=calculateAB(trueclus',fc);
            [RandIndx,AdjRandIndx] = RandIndices(A,B);
            NMI=nmi(trueclus,fc');
            ACC=accuracy(trueclus,fc')/100;
            dlmwrite(Evaluation_name,[AdjRandIndx,NMI,ACC,K],'-append');
        end
        toc
    end
    end
end

rmpath(genpath(cd));
