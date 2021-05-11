function [ List1 ] = Traverse( S,T,D,n,Lx_d,Lx_u,Ly_d,Ly_u,xv,yv,List)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
List1=[List;S'];
size(List1,1)
Neigh=Neighbour(S);

for i=1:size(Neigh,1)
    S=Neigh(i,:)';
    QS=T'*S;
    QSS=D'*S;
    if QSS(1)>Lx_d && QSS(1)<Lx_u && QSS(2)>Ly_d && QSS(2)<Ly_u
        if nnz(ismember(List1(:,1:n),S','rows'))==0
                in = inpolygon(QS(1),QS(2),xv,yv);
                    if in==1
                         List1=Traverse( S,T,D,n,Lx_d,Lx_u,Ly_d,Ly_u,xv,yv,List1);
                    end
        end
    end
end
