clc;
clear;

n=5; 
m=2; 
u1=sqrt(2/5)*[cos(0);cos(2*pi*(1/5));cos(2*pi*(2/5));cos(2*pi*(3/5));cos(2*pi*(4/5))];
u2=sqrt(2/5)*[sin(0);sin(2*pi*(1/5));sin(2*pi*(2/5));sin(2*pi*(3/5));sin(2*pi*(4/5))];
t=[1;1;1;1;1];

%u1=round(u1,1);
%u2=round(u2,1);

u1=u1/sqrt(dot(u1,u1));
u2=u2/sqrt(dot(u2,u2)); 



PL5=u1*u1'+u2*u2'; 
QL5=eye(n)-PL5; 
tQ=QL5*t;


L_u=5;
L_d=-5;


List=[];
c=1;

value=0;

while value==0
    A1=rand(n,1)<0.5;
    if Check_Cut( A1,QL5,tQ,n )==1
        Value=1;
        break
    end
end
    
S=A1';

tic
List=Traverse( S,QL5,tQ,n,L_d,L_u,List);
toc

PList=PL5*List';
PList=PList';
uu1=repmat(u1',size(PList,1),1);
uu2=repmat(u2',size(PList,1),1);
Fi1=dot(PList,uu1,2);
Fi2=dot(PList,uu2,2);

Final_Points=[Fi1 Fi2];

figure(1)
scatter(Final_Points(:,1),Final_Points(:,2),'filled','sizedata',8)
axis([L_d L_u L_d L_u])

DD = pdist2(List,List,'euclidean');

DD=(DD==1);

figure(3)
gplot(DD,Final_Points,'-*')

figure(4)
scatter3(List(:,1),List(:,2),List(:,3),'filled')

figure(5)
scatter(Final_Points(:,1),Final_Points(:,2),'filled')



save('AperiodicP3TypePenroseTiling.cord')
