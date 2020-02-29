% function [RandIndx, AdjRandIndx] = RandIndices(A,B)
% This code calculates the rand index between two cluster partitions
% A = input cluster number 1 as structure, example mat.file is given in the folder
% B = Input cluster number 2 as structure, example mat.file is given in the folder
% Date 24 March 2013(Hemanta Medhi)

function [RandIndx, AdjRandIndx] = RandIndices(A,B)

cls = ContingencyTable(A,B);
[m n] = size(cls);
no_objects = sum(sum(cls));
i_dot = sum(cls,2);
j_dot = sum(cls,1);
for i = 1:m
     for j = 1:n
        if cls(i,j)<2
        tr (i,j) = 0;
        else
        tr (i,j) = nchoosek(cls(i,j),2);
        end
     end
end
term1 = sum(sum(tr));
for in = 1:length(i_dot)
    if i_dot(in)<2
       c_i(in) = 0;
    else
       c_i(in) = nchoosek(i_dot(in),2);  
    end
end

for jn = 1:length(j_dot)
    if j_dot(jn)<2
       c_j(jn) = 0;
    else
       c_j(jn) = nchoosek(j_dot(jn),2);  
    end
end

A = sum(sum(tr));
B = sum(c_i) - sum(sum(tr));
C = sum(c_j) - sum(sum(tr));
D = nchoosek(no_objects,2)-A-B-C;
RandIndx =(A+D)/(A+B+C+D);

ARI_t = term1 - (sum(c_i)*sum(c_j))/nchoosek(no_objects,2);
ARI_d = 0.5*(sum(c_i) + sum(c_j)) - (sum(c_i)*sum(c_j))/nchoosek(no_objects,2);
AdjRandIndx = ARI_t/ARI_d;

end