function  pop=calculated_crowding_distance(pop,F)

C=cat(1,pop.cost);

nobj=length(pop(1).cost);

NF=length(F);

for i=1:NF
    
    NFM=length(F{i});   % the length of the ranking set
    C0=C(F{i},:);   %the fitness values of the current ranking set
    
    D=zeros(NFM,nobj);
    
    for j=1:nobj
    
        Cj=C0(:,j);  % the fitness values of the j'th object
        
        [value,index]=sort(Cj);  % in ascending order
        
        minc=value(1);
        maxc=value(end);
        
        % for each rank, the D value at the boundery is 10 
        D(index(1),j)=10;   %the min boundery 
        D(index(end),j)=10;   %the max boundery
         
        for k=2:NFM-1            
           D(index(k),j)=abs(value(k+1)-value(k-1))/(maxc-minc); 
        end
        
    end
    
    %Note that the order in D is the same as that in F{i}
    for z=1:NFM
       pop(F{i}(z)).cdis=sum(D(z,:)); 
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                          www.matlabnet.ir                         %
%                   Free Download  matlab code and movie            %
%                          Shahab Poursafary                        %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
