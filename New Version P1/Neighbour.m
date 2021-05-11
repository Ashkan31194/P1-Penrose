function [ Neigh ] = Neighbour( t )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=size(t,1);
Neigh=repmat(t',n*(n-1),1);
v=[1 -1 zeros(1,n-2)];
P=unique(perms(v),'rows');
dNeigh=P;
Neigh=Neigh+dNeigh;
end

