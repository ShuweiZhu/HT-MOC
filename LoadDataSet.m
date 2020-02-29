function [data,trueclus,K] = LoadDataSet(dataName)

%% Load the systhetic dataset
dataNameFull = ['data_',dataName,'.mat'];
if exist(dataNameFull)
    gt = [];  fea = [];
    load(['data_',dataName,'.mat'],'fea','gt');
    data=fea; trueclus=gt;
    K = numel(unique(gt)); % The number of clusters
else
    %% dataName is the input argument that determines the desire DataSet
    switch lower(dataName)
        %% Large data sets %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'pendigits'
            gt = []; fea = [];
            load(['data_',dataName,'.mat'],'fea','gt');
            data=fea; trueclus=gt;
            K = numel(unique(gt)); % The number of clusters
            fprintf('Pendigits is loaded.\n');
        case 'usps'
            gt = []; fea = [];
            load(['data_',dataName,'.mat'],'fea','gt');
            data=fea; trueclus=gt;
            K = numel(unique(gt)); % The number of clusters
            fprintf('USPS is loaded.\n');
        case 'letters'
            gt = []; fea = [];
            load(['data_',dataName,'.mat'],'fea','gt');
            data=fea; trueclus=gt;
            K = numel(unique(gt)); % The number of clusters
            fprintf('Letters is loaded.\n');
        case 'mnist'
            gt = []; fea = [];
            load(['data_',dataName,'.mat'],'fea','gt');
            data=fea; trueclus=gt;
            K = numel(unique(gt)); % The number of clusters
            fprintf('MNIST is loaded.\n');
        case 'covertype'
            gt = []; fea = [];
            load(['data_',dataName,'.mat'],'fea','gt');
            data=fea; trueclus=gt;
            K = numel(unique(gt)); % The number of clusters
            fprintf('Covertype is loaded.\n');
        case 'isolet'
            gt = []; fea = [];
            load(['data_',dataName,'.mat'],'fea','gt');
            data=fea; trueclus=gt;
            K = numel(unique(gt)); % The number of clusters
            fprintf('ISOLET is loaded.\n');
       %% Real world data set %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% UCI datasets
        case 'wine'
            S = load('./Datasets/Wine.txt', '-ascii');
            data = S(:,2:end);
            trueclus = S(:,1);
            K = 3;
            fprintf('Wine is loaded.\n');
        case 'cancer'
            S = load('./Datasets/Cancer.txt', '-ascii');
            S = S(:,2:end);
            x = S(:,end)==2;
            S(x,end) =1 ;
            x = S(:,end)==4;
            S(x,end) =2 ;
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 2;
            fprintf('Cancer DataSet is loaded.\n');
        case 'shuttle2'
            S = load('./Datasets/Shuttle.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 7;
            fprintf('Shuttle is loaded.\n');
        case 'shuttle'
            S = load('./Datasets/Shuttle2.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 7;
            fprintf('Shuttle2 is loaded.\n');
        case 'magic'
            S=load('./Datasets/Magic.txt','-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 2;
            fprintf('Magic is loaded.\n');        
        case 'c_cube'
            S = load('./Datasets/c_cube.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 4;
            fprintf('c_cube is loaded.\n');
        case 'opticaldigit'
            S = load('./Datasets/Opticaldigit.txt', '-ascii');
            data = S(:,1:end-1);
            data(:,[1,40])=[];
            trueclus = S(:,end);
            K = 10;
            fprintf('Opticaldigit is loaded.\n');
        case 'pendigit'
            S = load('./Datasets/Pendigit.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 10;
            fprintf('Pendigit is loaded.\n');
            
       %% Sythentic data %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'd31'
            S = load('./Datasets/Artificial_data/D31.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 31;
            fprintf('D31 is loaded.\n');
        case 'r15'
            S = load('./Datasets/Artificial_data/R15.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 15;
            fprintf('R15 is loaded.\n');
        case 'aggregation'
            S = load('./Datasets/Artificial_data/Aggregation.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 7;
            fprintf('Aggregation is loaded.\n');
        case 'flame'
            S = load('./Datasets/Artificial_data/Flame.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            K = 2;
            fprintf('Flame is loaded.\n');
        case 't4.8k'
            S = load('./Datasets/Artificial_data/T4.8k.txt', '-ascii');
            data = S(:,1:end-1);
            trueclus = S(:,end);
            ID=find(trueclus>0);
            data =data(ID,:);
            trueclus =trueclus(ID,:);
            K = 6;
            fprintf('T4.8k is loaded.\n');                         
           
    end
end
