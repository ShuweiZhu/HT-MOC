function [Seed,DD,NeighborSeed,Layer,Node,Child,CorSet]=GrowingTrainNode(Pm,X)

[n,d]=size(X);

% pre-assign some values
Node{1}=mean(X);
Layer(1,1)=1; % map the layer No. and the sequential No. of each nodes
Layer(2,1)=2;
Child{1}=1; % record if have child node for each nodes. 1: has child, 0: no child.
CorSet{1}=true(1,n); % at the beginning, corset of the first node is the whole data set
StopGrow=0; % determiner to decide if to grow the topology or not
PreNode=1; % sequential number of the just processed hypernode
ThisCoarse=1; % sequential number of the coarse node in the just processed layer
ThisLayer=2; % sequential number of the layer to be created
ThisNode=2; % sequential number of the node to be created
Counter=2; % No. of existing hypernode (including several sun nodes)
Reader=(1:n); % reader to read out the sequential number of data objects

%% growing training
% recursively grow and train the topology
LoopTime=0;
TT=0;
while ThisLayer~=1
    if StopGrow==0 % 0: grow deeper layer, 1: stop grow
        SubDataLct=Reader(CorSet{PreNode}(ThisCoarse,:)); % extract the "SubData" for deeper layer growth
        SubData=X(SubDataLct,:);
        % Begain train nodes
        [Node{ThisNode},CorSet{ThisNode},NumSubNode] = TrainNodes(SubData,SubDataLct,Pm);
        Child{ThisNode}=zeros(1,NumSubNode);
    end
    Bad1=find(sum(CorSet{ThisNode},2)>Pm.UpLimit); % find the node with high membership
    Bad2=find(Child{ThisNode}==0); % find the node with no heir
    CoarsePool=intersect(Bad1,Bad2);
    if isempty(CoarsePool)==1 % back to check the upper layer
        ThisLayer=ThisLayer-1;
        ThisNode=PreNode;
        PreNode=find(Layer==ThisLayer-1, 1, 'last' );
        StopGrow=1;
    else % if nodes have no heir and too many data points
        ThisCoarse=CoarsePool(1);
        Counter=Counter+1;
        Child{ThisNode}(ThisCoarse)=1;
        PreNode=ThisNode;
        ThisNode=Counter;
        ThisLayer=ThisLayer+1;
        Layer(ThisNode,1)=ThisLayer;
        StopGrow=0;
    end
end
Pm.TotalNumNode=Counter;

%% Find all leaf nodes (seed) and the coresponding subset %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LeafNodes,LeafLayer,Subset,~]=FindLeafNodes(Pm,Layer,Node,Child,CorSet);
DD=pdist2(LeafNodes,LeafNodes);
[~,NeighborSeed]=sort(DD,2);  % The neighbors of all seed points

figure(3) % Draw the result of nodes in different layers
[~,~]=DrawNodes(X,Layer,Node,CorSet,LeafNodes,LeafLayer,Subset);

%%%%% Re-tained the existing leaf nodes
AdjustTimes=0;
Converge=0;
while Converge==0
    AdjustTimes=AdjustTimes+1;
    Ls=3; % The neighbour size of each seed point
    count=0;  %The number of points to be adjusted for new nodes
    %%% Adjust the border points to the related seed points
    t=1;
    NumNode=size(LeafNodes,1);
    Core=ComputeCore(X,NumNode,Subset,Pm.UpLimit); % Detect core points and border points
    while t<=NumNode
        Point_t=find(Subset(t,:)==1);
        NNt=NeighborSeed(t,2:Ls+1); % (Note: the neighbors do not contain itself)
        if length(Point_t)<=fix(sqrt(sqrt(n))) % Subset size<=n^(1/4), delete this node
            X_t=X(Point_t,:);
            Node_t=LeafNodes(NNt,:);
            [~,lab]=min(pdist2(X_t,Node_t),[],2);
            count=count+length(Point_t);
            Lab2=unique(lab);
            for j=1:length(Lab2)
                Ind=Point_t(lab==Lab2(j));
                Subset(NNt(Lab2(j)),Ind)=1;
            end
            LeafNodes(t,:)=[]; % Remove the node t
            LeafLayer(t,:)=[];
            NeighborSeed(t,:)=[];
            Subset(t,:)=[];
            Core(t,:)=[];
            NumNode=NumNode-1;
            %%%%%% Adjust the neighbor index (important!) %%%%%%%
            NeighborSeed2=NeighborSeed';
            NeighborSeed2(NeighborSeed2==t)=[];
            Ind=NeighborSeed2>t;
            NeighborSeed2(Ind)=NeighborSeed2(Ind)-1;
            NeighborSeed=reshape(NeighborSeed2,NumNode,NumNode)';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % Re-assign the points the correct nodes
            % Border_t=Core{t,3};
            X_t=X(Point_t,:);
            NNt2=[t,NNt]; % (Note: the neighbors contain itself)
            Node_t=LeafNodes(NNt2,:);
            [~,lab]=min(pdist2(X_t,Node_t),[],2);
            Lab2=setdiff(lab,1); % (Note: the diff labels only appear as a single value)
            if ~isempty(Lab2) % Re-assginment of boder points
                for j=1:length(Lab2)
                    count=count+length(find(lab==Lab2(j)));
                    Ind=Point_t(lab==Lab2(j));
                    Subset(t,Ind)=0;
                    Subset(NNt2(Lab2(j)),Ind)=1;
                end
            end
            t=t+1;
        end
    end
    
    %%% Detect the subet with a larger size (Size>1.1*U_L)
    for t=1:NumNode
        SubDataLct=find(Subset(t,:)==1);
        if length(SubDataLct)>1.1*Pm.UpLimit
            SubData=X(SubDataLct,:);
            [SubNode,CorSubSet,NumSubNode]=TrainNodes(SubData,SubDataLct,Pm); % Train nodes
            if min(sum(CorSubSet,2))>fix(sqrt(sqrt(n)))
                LeafLayer=[LeafLayer(1:t-1);
                    ones(NumSubNode,1)*LeafLayer(t);
                    LeafLayer(t+1:NumNode)];
                LeafNodes=[LeafNodes(1:t-1,:); SubNode;
                    LeafNodes(t+1:NumNode,:)];
                Subset=[Subset(1:t-1,:); CorSubSet;
                    Subset(t+1:NumNode,:)];
                NumNode=NumNode+1;
            end
        end
    end
    
    %%% Re-train the seed points
    for i=1:NumNode
        Lab=find(Subset(i,:)==1);
        DataTrain=X(Lab,:);
        for j=1:length(Lab)
            LeafNodes(i,:)=LeafNodes(i,:)+Pm.LearnRate1*(DataTrain(j,:)-LeafNodes(i,:));
        end
    end
    DD=pdist2(LeafNodes,LeafNodes);
    [~,NeighborSeed]=sort(DD,2);
    
    %%% testing the convergence using count
    Count(AdjustTimes,1)=count;
    if count<=0.01*n || AdjustTimes==10
        Converge=1;
    end
end

figure(4) % Draw the result of nodes in different layers
% n=size(X,1);
% n2=randperm(n,fix(0.1*n));
% plot(X(n2,1),X(n2,2),'.');
plot(X(:,1),X(:,2),'.');
hold on
plot(LeafNodes(:,1),LeafNodes(:,2),'ro','MarkerFaceColor','r','MarkerSize',6);
xlim([-0.04,1.04]);
ylim([-0.04,1.04]);
set(gca,'XTick',[],'YTick',[]);
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'LooseInset',[0 0 0 0])
%%% Stage 3: Detect the subet with a larger size
time2=0;
for t=1:NumNode
    aa=sum(Subset(t,:));
    if aa>1*Pm.UpLimit
        time2=time2+1;
%         aa
    end
end

% Core=ComputeCore(X,NumNode,Subset,Pm.UpLimit); % Detect core points and border points
local_core=zeros(n,1);
for i = 1:n
    local_core(i,1)=find(Subset(:,i)==1);
end
Seed.LeafNodes=LeafNodes;
% Seed.NeighborSeed=NeighborSeed;
Seed.LeafLayer=LeafLayer;
Seed.local_core=local_core;
% Seed.dist=DD;
% Seed.Subset=Subset;
% Seed.Core=Core;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all leaf nodes (seed) and the coresponding subset %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LeafNodes,LeafLayer,Subset,Label]=FindLeafNodes(Pm,Layer,Node,Child,CorSet)

for i=1:length(CorSet) % convert CorSet from logic type to double type
    CorSet{i}=double(CorSet{i});
end

%%% Find the leaf nodes
LeafNodes=[];
LeafLayer=[];
Subset=[];
Label=zeros(1,Pm.Xlth);
k=0;
for j=2:length(Node)
    All_Child=Child{j}(:);
    Leaf=find(All_Child==0);
    if isempty(Leaf)==0 % leaf nodes exist, do the following processing
        for h=1:length(Leaf)
            leafnode=Node{j}(Leaf(h),:);
            LeafNodes=[LeafNodes;leafnode];
            LeafLayer=[LeafLayer;Layer(j)];
            sub_set=CorSet{j}(Leaf(h),:);
            Subset=[Subset;sub_set];
            k=k+1;
            Label(sub_set==1)=k;
        end
    end
end

%%% Find the number of leaf nodes in different layers
% UniLeafLayer=unique(LeafLayer);
% nLayer=length(UniLeafLayer);
% NumLeafLayer=zeros(nLayer,1);
% for i=1:nLayer
%       NumLeafLayer(i,1)=length(find(LeafLayer==UniLeafLayer(i)));
% end

end
