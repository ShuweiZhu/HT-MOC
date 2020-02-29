function short_path=GraphDist(LeafNodes,NeighborSeed,L)
%Construct the graph
supk=L;
Nnode=size(LeafNodes,1);
weight=zeros(Nnode,Nnode);
dist=pdist2(LeafNodes,LeafNodes);
for i=1:Nnode
    %     if isnoise(i)==0
    for j=2:supk+1
        x=NeighborSeed(i,j);
        weight(i,x)=dist(i,x);
    end
    %     end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the shortest path between cores
fprintf('Start computing the graph-based distance between local cores...\n');
short_path=zeros(Nnode,Nnode);%The shortest path between local cores
% maxd=(exp(alpha*maxd)-1)^(1/alpha);
weight2=sparse(weight);
fprintf('Start clustering...\n');
for i=1:Nnode
    short_path(i,i)=0;
    %      [D,Z]=dijkstra2(weight,cores(i));
    [D,~,Z] = graphshortestpath(weight2,i,'METHOD','Dijkstra');
    for j=i+1:Nnode
        short_path(i,j)=D(j);
        if short_path(i,j)==inf
            short_path(i,j)=0;
        end
        short_path(j,i)=short_path(i,j);
    end
end
maxd=max(max(short_path));
for i=1:Nnode
    for j=1:Nnode
        if short_path(i,j)==0
            short_path(i,j)=maxd;
        end
    end
end
end
