function [All_nodes,All_CorSet]=DrawNodes(X,Layer,Node,CorSet,LeafNodes,LeafLayer,Subset)
%% Draw the result of nodes in different layers
%%% Find nodes in different layers
All_nodes=cell(max(Layer),1);
All_CorSet=cell(max(Layer),1);
All_nodes{1}=Node{1};
All_nodes_leaf=cell(length(unique(LeafLayer)),1);
All_nodes_Set =cell(length(unique(LeafLayer)),1);
for t=2:max(Layer)
    Ind=Layer==t;
    All_nodes{t}=cat(1,Node{Ind}); 
    All_CorSet{t}=cat(1,CorSet{Ind}); 
    if t>=min(LeafLayer)
        All_nodes_leaf{t}=LeafNodes(LeafLayer==t,:);
        All_nodes_Set{t}=Subset(LeafLayer==t,:);
    end
end

Y2=[];
for t=1:max(Layer)
%      n=size(X,1);
%      n2=randperm(n,fix(0.1*n));
%      plot(X(n2,1),X(n2,2),'.');
     plot(X(:,1),X(:,2),'.');
     hold on
     Y1=All_nodes{t};
     plot(Y1(:,1),Y1(:,2),'ro','MarkerSize',6);
     if t>=min(LeafLayer)
         Y2=[Y2;All_nodes_leaf{t}];
         plot(Y2(:,1),Y2(:,2),'ro','MarkerFaceColor','r','MarkerSize',6);
     end
     xlim([-0.04,1.04]);
     ylim([-0.04,1.04]);
     set(gca,'XTick',[],'YTick',[]);
     set(gca,'LooseInset',get(gca,'TightInset'))
     hold off     
end

