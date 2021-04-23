function [ check ] = Check_Cut( S,QL5,tQ,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A=zeros(2*n);
A(1:n,1:n)=QL5;
A(n+1:2*n,1:n)=eye(n);
A(n+1:2*n,n+1:2*n)=eye(n);
g=QL5*S;
b=g-tQ;


f=[zeros(1,2*n) ones(1,2*n)];
Aeq=zeros(2*n,4*n);
Aeq(:,1:2*n)=A;
Aeq(:,2*n+1:4*n)=eye(2*n);
beq=ones(1,2*n);
beq(1:n)=b;
lb=zeros(1,4*n);
ub=[];
options=optimoptions('linprog','Display','none');
[~,fval]  = linprog(f',[],[],Aeq,beq,lb,ub,[],options);
if fval<1e-5
    check=1;
else
    check=0;
end


end

