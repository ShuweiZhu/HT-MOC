function [link2,DD3]=DrawMST(X,LeafNodes,DD)
%%% Description: Draw the MST result
% Input: X - data set; 
%           LeafNodes - Seed points
%           DD - DD=pdist2(LeafNodes,LeafNodes)
% Output: MST-graph links

Y=LeafNodes;
N=size(DD,1);
DD2=zeros(N,N);  d_max=zeros(N,1);
for ii=1:N
    dd=DD(ii,:);
    [aa,bb]=sort(dd);
    d_max(ii,1)=aa(end);
    for jj=1:N
        DD2(ii,jj)=find(bb==jj);
    end
end
DD3=zeros(N,N);
MaxDD=max(max(d_max));
for ii=1:N
    for jj=ii+1:N
        DD3(ii,jj)=min(DD2(ii,jj),DD2(jj,ii));
    end
end
DD3=eye(N)+DD3+DD3';
DD3=DD3+DD./MaxDD;

%% Perform PrimMST with different distance matrix
% DD: the Euclidian distance
% DD3: the interesting links combined with DD
link=PrimMST(DD3,N);
[~,Ind_link]=sort(link(:,3),'descend');
link2=link(Ind_link,:);

% add the link from 1 to another
Ind=find(link2(:,1)==1); Id=Ind(1);
if Id<N-1
link2=[link2(1:Id,:);
          [link2(Id,2),1,link2(Id,3)];
           link2(Id+1:end,:)];
else
link2=[link2(1:Id,:);
          [link2(Id,2),1,link2(Id,3)]]; 
end
%% Draw the result
figure(1)
% n=size(X,1);
% n2=randperm(n,fix(0.1*n));
% plot(X(n2,1),X(n2,2),'.');
plot(X(:,1),X(:,2),'.');
hold on
plot(Y(:,1),Y(:,2),'ro','MarkerFaceColor','r','MarkerSize',6);
hold on;
for i=1:N
    plot([Y(link2(i,1),1),Y(link2(i,2),1)],[Y(link2(i,1),2),Y(link2(i,2),2)],'LineWidth',1);
    hold on;
end
xlim([-0.04,1.04]);
ylim([-0.04,1.04]);
set(gca,'XTick',[],'YTick',[]);
set(gca,'LooseInset',get(gca,'TightInset'))
hold off

end

%% Construct the minimum spanning tree (MST) using Prim's algorithm
function link=PrimMST(DD,n)
% Prim
A=DD;
A(A==0)=Inf;
P=zeros(1,n);
P(1,1)=1;
V=1:n;
V_P=V-P;
link=zeros(n-1,3);
k=1;
while k<n
    p=P(P~=0);
    v=V_P(V_P~=0);
    pv =min(min(A(p,v)));
    [x,y]=find(A==pv);
    for i=1:length(x)
        if  any(P==x(i)) && any(V_P==y(i)) %% Key step
            P(1,y(i))=y(i);
            V_P=V-P;
            link(k,:)=[x(i),y(i),DD(y(i),x(i))]; %% y direct to x
            k = k+1;
            break;
        end
    end
end

end
