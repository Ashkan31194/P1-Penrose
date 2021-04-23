clc;
clear;
close all

n=5; %enter the dimension of the space
m=2; %enter the dimesion of the projection space
% u1 and u2 are the unit vectors of the projection vector - the total
% number of uis mus be equal to m
u1=sqrt(2/5)*[cos(0);cos(2*pi*(1/5));cos(2*pi*(2/5));cos(2*pi*(3/5));cos(2*pi*(4/5))];
u2=sqrt(2/5)*[sin(0);sin(2*pi*(1/5));sin(2*pi*(2/5));sin(2*pi*(3/5));sin(2*pi*(4/5))];
u3=sqrt(2/5)*[cos(0);cos(4*pi*(1/5));cos(4*pi*(2/5));cos(4*pi*(3/5));cos(4*pi*(4/5))];
u4=sqrt(2/5)*[sin(0);sin(4*pi*(1/5));sin(4*pi*(2/5));sin(4*pi*(3/5));sin(4*pi*(4/5))];
u5=sqrt(1/5)*[1;1;1;1;1];
u1=u1/sqrt(dot(u1,u1)); %normalize u1
u2=u2/sqrt(dot(u2,u2)); %normalize u2
u3=u3/sqrt(dot(u3,u3)); %normalize u3
u4=u4/sqrt(dot(u4,u4)); %normalize u4
u5=u5/sqrt(dot(u5,u5)); %normalize u5



D=[[cos(0);cos(2*pi*(1/5));cos(2*pi*(2/5));cos(2*pi*(3/5));cos(2*pi*(4/5))] [sin(0);sin(2*pi*(1/5));sin(2*pi*(2/5));sin(2*pi*(3/5));sin(2*pi*(4/5))]];
T=[[cos(0);cos(4*pi*(1/5));cos(4*pi*(2/5));cos(4*pi*(3/5));cos(4*pi*(4/5))] [sin(0);sin(4*pi*(1/5));sin(4*pi*(2/5));sin(4*pi*(3/5));sin(4*pi*(4/5))] [1;1;1;1;1]];

D=round(D,2);%comment these if aperiodic tilig is required
T=round(T,2);%comment these if aperiodic tilig is required
L_u=30;%determine the max size the search algorithm goes through
L_d=-30;%determine the min size the search algorithm goes through

L = linspace(0,2.*pi,11);
xv = sin(L)';
yv = cos(L)';
taw=(1+sqrt(5))/2;
xv=xv*((taw+2)/2)*1.05;
yv=yv*((taw+2)/2)*1.05;


Layer=[xv yv];


S=[1;0;0;0;0];
QS=T'*S;
[in,on]  = inpolygon(QS(1),QS(2),xv,yv);

List=[];
tic
List = Traverse( S,T,n,L_d,L_u,xv,yv,List);
toc

PList=D'*List';
PList=PList';

figure(2)
scatter(PList(:,1),PList(:,2),'filled','sizedata',8)
axis([L_d L_u L_d L_u])

DD = pdist2(List,List,'euclidean');

AF=find(abs(DD-sqrt(2))<1e-5);
AFF=find(abs(DD-sqrt(2))>1e-5);

DD(AF)=1;
DD(AFF)=0;
tic
for i=1:size(DD,1)
    for j=1:size(DD,1)
        if DD(i,j)==1
            k=List(i,:)-List(j,:);
            if nnz(k+circshift(k',1)')==2
                DD(i,j)=1;
            else
                DD(i,j)=0;
            end
        end
    end
end
counter=1;
for i=1:size(DD,1)
    for j=i+1:size(DD,1)
       if DD(i,j)==1
           Edges(counter,:)=[i j];
           counter=counter+1;
       end
    end
end
R=rotx(11);
RR=R(2:3,2:3);
QQ=RR*PList';
QQ=QQ';

%this for can be removed from the code! the speed becomes much faster. But
%a few inappropriate edges produce.
for i=1:size(Edges,1)
    for j=i+1:size(Edges,1);
        A=QQ(Edges(i,1),:);
        B=QQ(Edges(i,2),:);
        C=QQ(Edges(j,1),:);
        D=QQ(Edges(j,2),:);
        m1=(B(2)-A(2))/(B(1)-A(1));
        m2=(D(2)-C(2))/(D(1)-C(1));
        x_int=(m1*A(1)-m2*C(1)-A(2)+C(2))/(m1-m2);
        y_int=m1*(x_int-A(1))+A(2);
        Int=[x_int y_int];
        check1=dot((Int-A),(Int-B));
        check2=dot((Int-C),(Int-D));
        if check1<-1e-5 && check2<-1e-5
            DD(Edges(i,1),Edges(i,2))=0;
            DD(Edges(i,2),Edges(i,1))=0;
            DD(Edges(j,1),Edges(j,2))=0;
            DD(Edges(j,2),Edges(j,1))=0;  
        end
    end
end
Edges=[];
counter=1;
for i=1:size(DD,1)
    for j=i+1:size(DD,1)
       if DD(i,j)==1
           Edges(counter,:)=[i j];
           counter=counter+1;
       end
    end
end
      
toc
 figure(3)
 gplot(DD,PList,'k-')
 
[X,Y]=gplot(DD,PList);

figure(4)
plot(X,Y,'-mo',...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',4);
axis([L_d L_u L_d L_u])

figure(5)
scatter(PList(:,1),PList(:,2),'filled','sizedata',10)

save('P11.cord')
