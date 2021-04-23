function [ Neigh ] = Neighbour( t,L_d,L_u )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=size(t,1);
Neigh=repmat(t',n*(n-1),1);
v=[1 -1 zeros(1,n-2)];
P=unique(perms(v),'rows');
dNeigh=P;
Neigh=Neigh+dNeigh;
M=max(Neigh,[],2);
m=min(Neigh,[],2);
L1=find(M>L_u);
L2=find(m<L_d);
LL=union(L1,L2);
Neigh(LL,:)=[];
end
