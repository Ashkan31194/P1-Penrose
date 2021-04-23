function [ List1 ] = Traverse( S,QL5,tQ,n,L_d,L_u,List)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
List1=[List;S];
size(List1,1)
Neigh=Neighbour( S',L_d,L_u );
for i=1:size(Neigh,1)
    S=Neigh(i,:);
    if nnz(ismember(List1(:,1:n),S,'rows'))==0
             if Check_Cut( S',QL5,tQ,n )==1
                 List1=Traverse( S,QL5,tQ,n,L_d,L_u,List1);
             end
    end
end

end

