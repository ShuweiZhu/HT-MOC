function [rho,d_rho,DD,G_ep,dt_set]=Compute_data(data,KK)

[n,d]=size(data);

%% Compute the density rho

DD = pdist2(data,data,'euclidean');
DD2=triu(DD);
vD=DD2(DD2>0);
max_D1=max(vD)+1;

% Defining the cutoff distance dc (see ref. Rodriguez & Laio (2014))
percent=2;
position=(((n*(n-1))/2)*percent)/100;
vD=sort(vD);
dc=vD(fix(position)); % cutoff distance: used to compute the density

% Defining the density rho of each point
rho=zeros(n,1);

% Gaussian Kernel
for i=1:n
    for j=i+1:n
        rho(i,1)=rho(i,1)+exp(-(DD(i,j))^2/(2*dc.^2));
        rho(j,1)=rho(j,1)+exp(-(DD(i,j))^2/(2*dc.^2));
    end
end
% Scaling rho (between 0 and 1) and computing mean and std of rho
rho_max=max(rho);
rho=rho./rho_max;

% Compute mean and std of rho
m_rho=mean(rho);
s_rho=std(rho);
d_rho=m_rho-s_rho;

if d_rho<=0.0
    d_rho=m_rho/2.0;
end

%% Defining the Epistasis Graph (G_ep)
G_ep=zeros(n,KK);
M_ep=zeros(n,n);
for i=1:n
    % Step 1: Connecting each object to the closest object with higher density
    Ind= find(rho>rho(i,1));
    if isempty(Ind)
        % Object with highest density
        dist_i=DD(i,:);
        [~,id2]=min(dist_i(dist_i>0));
        ind_min=id2;
        delta=DD(i,ind_min);
    else
        [~,id]=min(DD(i,Ind));
        ind_min=Ind(id);
        delta=DD(i,ind_min);
    end
    
    % Step 2: Adding edges to each vertex (according to the distance)
    %         until the indegree is equal to K
    G_ep(i,1)=ind_min;
    [~,Ind2]=sort(DD(i,:));
    Ind2(Ind2==ind_min)=[];
    G_ep(i,2:KK)=Ind2(1,2:KK);
    
    M_ep([ind_min,Ind2(2:KK)],i)=1;
end

% Compute parameters dt1, dt2, dt3 according to threeSigmaRule
dt_set=threeSigmaRule(DD,G_ep,n,KK);

end

%% Compute parameters dt1, dt2, dt3
% rho:68.27%, 95.45% and 99.73% of the values lie within one, two
%     and three standard deviations of the mean, respectively
function dt_set=threeSigmaRule(DD,G_ep,n,KK)

xg=zeros(n,KK);

for i=1:n
    xg(i,:)=DD(i,G_ep(i,:));
end
Xg=xg(:);

m_adj=mean(Xg);
s_adj=std(Xg);

dt1=m_adj;
dt2=m_adj+s_adj;
dt3=m_adj+2*s_adj;

dt_set=[dt1,dt2,dt3];

end
