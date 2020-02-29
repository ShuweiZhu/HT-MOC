% function[C_table] = ContingencyTable(A,B)
% This code calculates the Contingency Table for calculating the Rand
% Indices
% A = input cluster number 1
% B = Input cluster number 2
% Date 24 March 2013(Hemanta Medhi)

function[C_table] = ContingencyTable(A,B)

C_table = zeros(length(A),length(B));
for i =1:length(A)
    for j = 1:length(B)
        C_table (i,j) = length(intersect(A{i},B{j}));        
    end
end

end