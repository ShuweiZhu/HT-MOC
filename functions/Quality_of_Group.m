function [Qual_K,L_K]=Quality_of_Group(Ks,Quality)

% Compute quality of each group (in terms of K)
Kset=unique(Ks);
SizeK=length(Kset);
Qual_K=zeros(SizeK,1);
L_K=zeros(SizeK,1);
for kk=1:SizeK
    K=Kset(kk);
    Ind=find(Ks==K);
    L_K(kk,1)=length(Ind);
    Qual_K(kk,1)=mean(Quality(Ind));
end

%% Adjust quality of each group based on neighbourhood
% Qual_K2=zeros(SizeK,1);
% for kk=1:SizeK
%     if kk==1
%        Qual_K2(kk,1)=sum(Qual_K(kk:kk+2))/3;
%     elseif kk==SizeK
%        Qual_K2(kk,1)=sum(Qual_K(kk-2:kk))/3;
%     elseif kk==2
%        Qual_K2(kk,1)=sum(Qual_K(kk-1:kk+2))/4;
%     elseif kk==SizeK-1
%        Qual_K2(kk,1)=sum(Qual_K(kk-2:SizeK))/4;   
%     else
%        Qual_K2(kk,1)=sum(Qual_K(kk-2:kk+2))/5;    
%     end
% end
