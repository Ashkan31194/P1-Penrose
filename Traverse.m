function [ List1 ] = Traverse( S,T,n,L_d,L_u,xv,yv,List)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
List1=[List;S'];
size(List1,1)
Neigh=Neighbour( S,L_d,L_u );

for i=1:size(Neigh,1)
    S=Neigh(i,:)';
    if nnz(ismember(List1(:,1:n),S','rows'))==0
            QS=T'*S;
            in = inpolygon(QS(1),QS(2),xv,yv);
             if in==1
                 List1=Traverse( S,T,n,L_d,L_u,xv,yv,List1);
             end
    end
end
end
