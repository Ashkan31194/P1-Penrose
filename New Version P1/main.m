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
%t=[0.3;0.3;0.3;0.3;0.3]; %enter the t vector which moves the unit hypercube to point t


u1=u1/sqrt(dot(u1,u1)); %normalize u1
u2=u2/sqrt(dot(u2,u2)); %normalize u2
u3=u3/sqrt(dot(u3,u3)); %normalize u3
u4=u4/sqrt(dot(u4,u4)); %normalize u4
u5=u5/sqrt(dot(u5,u5)); %normalize u5

PL=u1*u1'+u2*u2'; %Projection Matrix - total number of sums=m
QL=eye(n)-PL; 


D=[[cos(0);cos(2*pi*(1/5));cos(2*pi*(2/5));cos(2*pi*(3/5));cos(2*pi*(4/5))] [sin(0);sin(2*pi*(1/5));sin(2*pi*(2/5));sin(2*pi*(3/5));sin(2*pi*(4/5))]];
T=[[cos(0);cos(4*pi*(1/5));cos(4*pi*(2/5));cos(4*pi*(3/5));cos(4*pi*(4/5))] [sin(0);sin(4*pi*(1/5));sin(4*pi*(2/5));sin(4*pi*(3/5));sin(4*pi*(4/5))] [1;1;1;1;1]];
n1=norm(T(:,1));
n2=norm(T(:,2));
n3=norm(T(:,3));
%D=[u1 u2];
%T=[u3 u4 u5];
%D=round(D,1);
order1=6; %input variable of fibonacci estimation 1
order2=8; %input variable of fibonacci estimation 2
p1=fibonacci(order1+1);
q1=fibonacci(order1);
p2=fibonacci(order2+1);
q2=fibonacci(order2);
T=[[2*q1;-p1;p1-q1;p1-q1;-p1] [0;q2;-p2;p2;-q2] [1;1;1;1;1]] ;
T(:,1)=T(:,1)*(n1/norm(T(:,1)));
T(:,2)=T(:,2)*(n2/norm(T(:,2)));
T(:,3)=T(:,3)*(n3/norm(T(:,3)));
L = linspace(0,2.*pi,11);
xv = sin(L)';
yv = cos(L)';
taw=(1+sqrt(5))/2;
xv=xv*(((p1/q1)+2)/2)*1.05;
yv=yv*(((p2/q2)+2)/2)*1.05;

Period1=[1 0]*5*(taw^order1); %x length of box
Period2=[0 2*sin(pi/5)]*(taw^(order2+1)); %y length of box
AA=[1 0];
AA(1)=AA(1)-Period1(1)/2;
AA(2)=AA(2)-Period2(2)/2;

figure(100)
rectangle('Position',[AA Period1(1) Period2(2)])

figure(14)
plot(xv,yv) % polygon
axis([-2.5 2.5 -2.5 2.5])


L_u=30;
L_d=-30;

S=[1;0;0;0;0];
QS=T'*S;
in  = inpolygon(QS(1), QS(2),xv,yv);



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

for i=1:size(Edges,1)
        A=QQ(Edges(i,1),:);
        B=QQ(Edges(i,2),:);
        A=repmat(A,size(Edges,1),1);
        B=repmat(B,size(Edges,1),1);
        C=QQ(Edges(:,1),:);
        D=QQ(Edges(:,2),:);
        m1=(B(:,2)-A(:,2))./(B(:,1)-A(:,1));
        m2=(D(:,2)-C(:,2))./(D(:,1)-C(:,1));
        x_int=(m1.*A(:,1)-m2.*C(:,1)-A(:,2)+C(:,2))./(m1-m2);
        y_int=m1.*(x_int-A(:,1))+A(:,2);
        Int=[x_int y_int];
        check1=dot((Int-A),(Int-B),2);
        check2=dot((Int-C),(Int-D),2);
        j=find((check1<-1e-5) & (check2<-1e-5));
        if size(j,1)>0
            DD(Edges(i,1),Edges(i,2))=0;
            DD(Edges(i,2),Edges(i,1))=0;
            DD(Edges(j,1),Edges(j,2))=0;
            DD(Edges(j,2),Edges(j,1))=0;  
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

Gr=graph(DD);
cycles2 = allcycles(Gr,'MinCycleLength',10,'MaxCycleLength',10);
reqcyc=[];
counter=1;
for i=1:size(cycles2,1)
    tenP=cell2mat(cycles2(i));
    Cords=PList(tenP,:);
    check=checkConvex(Cords(:,1)',Cords(:,2)');
    if check==1
        reqcyc(counter,:)=i;
        counter=counter+1;
    end
end

tenP=cell2mat(cycles2(reqcyc(1)));
Cords=PList(tenP,:);
SS=0.99*(Cords-mean(Cords))+mean(Cords);
pgn=polyshape(SS(:,1)',SS(:,2)');
    
for i=2:size(reqcyc,1)
    tenP=cell2mat(cycles2(reqcyc(i)));
    Cords=PList(tenP,:);
    SS=0.99*(Cords-mean(Cords))+mean(Cords);
    pgn1=polyshape(SS(:,1)',SS(:,2)');
    pgn=union(pgn,pgn1);
end
%[INPOLY, ONPOLY]=isinterior(pgn,PList(:,1),PList(:,2));
%TFin=INPOLY-ONPOLY;
TFin=isinterior(pgn,PList(:,1),PList(:,2));
nds=find(TFin==1);
PList(nds,:)=[];
DD(nds,:)=[];
DD(:,nds)=[];
      
toc
 figure(3)
 gplot(DD,PList,'k-')
 hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
[X,Y]=gplot(DD,PList);

figure(4)
plot(X,Y,'-mo',...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',4);
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)

figure(5)
scatter(PList(:,1),PList(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)

xvv=[AA(1);AA(1)+Period1(1);AA(1)+Period1(1);AA(1);AA(1)];
yvv=[AA(2);AA(2);AA(2)+Period2(2);AA(2)+Period2(2);AA(2)];
figure(16)
plot(xvv,yvv) % polygon
inn  = inpolygon(PList(:,1), PList(:,2),xvv,yvv);

PBList=PList(inn,:);
ExtraList=PList(~inn,:);
figure(6)
scatter(PBList(:,1),PBList(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
%scatter(ExtraList(:,1),ExtraList(:,2),'filled','r')

figure(7)
gplot(DD(inn,inn),PList(inn,:),'k-')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
[XX,YY]=gplot(DD(inn,inn),PList(inn,:));
figure(8)
plot(XX,YY,'-mo',...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',4);
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)



%save('P1.mat')
