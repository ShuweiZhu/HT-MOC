function  index=bi_tournament_select(N)

index=zeros(N,1);
for i=1:N
    k=tournoment_binary(N);
    index(i,:)=k;
end

function  j=tournoment_binary(N)
%%% binary tournoment selection

j1=randi([1 N]);
j2=randi([1 N]);
if j1<j2
    j=j1;
elseif j1<j2
    j=j2;
else
    if rand<0.5
        j=j1;
    else
        j=j2;
    end
end