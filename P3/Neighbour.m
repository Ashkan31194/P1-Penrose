function [ Neigh ] = Neighbour( t,L_d,L_u )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=size(t,1);
Neigh=repmat(t',2*n,1);
xx=eye(n);
yy=-eye(n);
dNeigh=[xx;yy];
Neigh=Neigh+dNeigh;
M=max(Neigh,[],2);
m=min(Neigh,[],2);
L1=find(M>L_u);
L2=find(m<L_d);
LL=union(L1,L2);
Neigh(LL,:)=[];

end

