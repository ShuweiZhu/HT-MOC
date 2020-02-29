function  pop=sorting(pop)

% sort the pop according to the crowding distance
[~,index]=sort([pop.cdis],'descend');
pop=pop(index);

% sort the pop according to the ranking value
[~,index]=sort([pop.rank]);
pop=pop(index);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                          www.matlabnet.ir                         %
%                   Free Download  matlab code and movie            %
%                          Shahab Poursafary                        %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%