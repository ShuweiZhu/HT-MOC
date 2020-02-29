function ClusterAssignment=Assignment_seed(Locus)

%% Assignment of clustering labeling
K_value=1;
T=size(Locus,1);
previous=zeros(T,1);
ClusterAssignment=zeros(T,1);
[aa,bb]=sort(Locus(:,2));
Locus=Locus(bb,:);

for i=1:T
    ctr=1;
    
    if ClusterAssignment(i,1)==0
        ClusterAssignment(i,1)=K_value;
        previous(ctr)=Locus(i,2);
        ctr=ctr+1;
        next=Locus(i,1);
        
        while ClusterAssignment(next,1)==0
            ClusterAssignment(next,1)=K_value;
            previous(ctr)=next;
            ctr=ctr+1;
            next=Locus(next,1);
        end
        
        if ClusterAssignment(next,1)~=K_value
            ctr=ctr-1;
            while ctr>0
                ClusterAssignment(previous(ctr),1)=ClusterAssignment(next,1);
                ctr=ctr-1;
            end
        else
            K_value=K_value+1;
        end
        
    end
    
end
