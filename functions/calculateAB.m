function [A,B]=calculateAB(clus,fc)
class=unique(clus);
kk=length(class);
A=cell(1,kk);
B=cell(1,kk);
for i=1:kk
    A{i}=find(clus==i);
    B{i}=find(fc==i);  
end
