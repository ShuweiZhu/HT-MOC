% when the hierarchy structure needs to grow, train the new deeper layer
% with Pm.k new neurons.

function [Node,CorSet,NumNode] = TrainNodes(SubData,DataLct,Pm)
%% preprocessing
[SubXlth,~]=size(SubData);
if SubXlth/(Pm.B-1)>Pm.UpLimit
    NumNode=Pm.B;
else
    NumNode=ceil(SubXlth/Pm.UpLimit);
end
% Node=zeros(NumNode,Pm.Xwd);

% RandomSize=Pm.UpLimit/2;
% for i=1:NumNode
%     NodeOri=SubData(randperm(SubXlth,RandomSize),:);
%     Node(i,:)=mean(NodeOri);
% end
%[~,Node]=kmeans(SubData,NumNode);
% Node=SubData(randperm(SubXlth,NumNode),:);
% for i=1:10
%     MidNode=SubData(randperm(SubXlth,NumNode),:);
%     for j=1:NumNode
%         [~,Winner]=min(dist(MidNode(j,:),Node'));
%         Node(Winner,:)=(Node(Winner,:)*(i-1)+MidNode(j,:))/i;
%     end
% end

% NumSib=length(NodeParentSib(:,1));
% if NumSib ==4 && isempty(NodeGrandParent)==0
%     Node(1,:)=2*NodeGrandParent-NodeParent;
%     for i=1:NumNode
%         Node(i,:)=(2*NodeGrandParent-NodeParentSib(i,:)+Node(1,:))/2;
%     end
% else
%     Node=SubData(randperm(SubXlth,NumNode),:);
% end

%Node=SubData(1:NumNode,:);
Node=SubData(randperm(SubXlth,NumNode),:);
Converge=0;
OriError=0;
ClassLableOri=zeros(SubXlth,1);
ClassLable=zeros(SubXlth,1);
Error=zeros(NumNode,1);
WinTime=ones(NumNode,1);
TotalTime=sum(WinTime);

%% training
TrainTimes=0;
while Converge==0 
    TrainTimes=TrainTimes+1;
    %TrainTimes
    %DataTrain=SubData(randperm(SubXlth,SubXlth),:);
    DataTrain=SubData;
    % training
    for i=1:SubXlth
        [~,Winner]=min(dist(Node,DataTrain(i,:)'));
        WinTime(Winner)=WinTime(Winner)+1;
        TotalTime=TotalTime+1;
        LearnRate=Pm.LearnRate1*(1-WinTime(Winner)/TotalTime);
        Node(Winner,:)=Node(Winner,:)+LearnRate*(DataTrain(i,:)-Node(Winner,:));
    end
    for i=1:SubXlth
        [~,ClassLable(i)]=min(dist(Node,DataTrain(i,:)')); % Euclidian distance
    %or [~,ClassLable(i)]=min(sqrt(sum((DataTrain(i,:)-Node).^2,2)));    
    end
    % testing the convergence using error
    for i=1:NumNode
        Error(i,1)=sum(dist(DataTrain(ClassLable==i,:),Node(i,:)')); % Euclidian distance   
    end
    if abs(mean(Error)-OriError)<Pm.ErrorBound*OriError || TrainTimes==10
        Converge=1;
    else
        OriError=mean(Error);
    end
end

% update the membership of this subset
CorSet=false(NumNode,Pm.Xlth);
for i=1:SubXlth
    [~,Host]=min(dist(Node,SubData(i,:)'));
    CorSet(Host,DataLct(i))=1;
end

end

