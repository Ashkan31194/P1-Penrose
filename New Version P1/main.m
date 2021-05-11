clc;
clear;
close all

n=5; 
D=[[cos(0);cos(2*pi*(1/5));cos(2*pi*(2/5));cos(2*pi*(3/5));cos(2*pi*(4/5))] [sin(0);sin(2*pi*(1/5));sin(2*pi*(2/5));sin(2*pi*(3/5));sin(2*pi*(4/5))]];
T=[[cos(0);cos(4*pi*(1/5));cos(4*pi*(2/5));cos(4*pi*(3/5));cos(4*pi*(4/5))] [sin(0);sin(4*pi*(1/5));sin(4*pi*(2/5));sin(4*pi*(3/5));sin(4*pi*(4/5))] [1;1;1;1;1]];
Offset=5;%the width between periodic Box and search space (5 is good)
order1=5;%enter the value of approximation for golen ratio in x direction Lx_Box=5*(taw^order1)
order2=8;%enter the value of approximation for golen ratio in y direction Ly_Box=(2*sin(pi/5))*(taw^(order2+1))

n1=norm(T(:,1));
n2=norm(T(:,2));
n3=norm(T(:,3));
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
xv=xv*(((p1/q1)+2)/2)*1.05;% be careful about the coefficient 1.05 (it might needs some modification case by case_1.05 works for most cases
yv=yv*(((p2/q2)+2)/2)*1.05;% be careful about the coefficient 1.05 (it might needs some modification case by case_1.05 works for most cases

Period1=[1 0]*5*(taw^order1);
Period2=[0 2*sin(pi/5)]*(taw^(order2+1));
AA=[1 0];
AA(1)=AA(1)-Period1(1)/2;
AA(2)=AA(2)-Period2(2)/2;

Dx_d=AA(1);
Ux_u=AA(1)+Period1(1);
Dy_d=AA(2);
Uy_u=AA(2)+Period2(2);

Lx_d=Dx_d-Offset;
Lx_u=Ux_u+Offset;
Ly_d=Dy_d-Offset;
Ly_u=Uy_u+Offset;

S=[1;0;0;0;0];
QS=T'*S;
in  = inpolygon(QS(1), QS(2),xv,yv);


tic
List=[];
List = Traverse( S,T,D,n,Lx_d,Lx_u,Ly_d,Ly_u,xv,yv,List);
toc

PList=D'*List';
PList=PList';

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
TFin=isinterior(pgn,PList(:,1),PList(:,2));
nds=find(TFin==1);
PList(nds,:)=[];
DD(nds,:)=[];
DD(:,nds)=[];

Grr=graph(DD);
cycles1 = allcycles(Grr,'MinCycleLength',4,'MaxCycleLength',4);
deldel=[];
counter=1;
for i=1:size(cycles1,1)
    fourP=cell2mat(cycles1(i));
    AAA=sum(DD(fourP,:),2);
    chch=find(AAA==2);
    cccc=find(AAA==5);
    if size(chch,1)==2 && size(cccc,1)==2
         deldel(2*counter-1:2*counter)=fourP(chch);
        counter=counter+1;
    end
end
PList(deldel,:)=[];
DD(deldel,:)=[];
DD(:,deldel)=[];

Grrr=graph(DD);
[bins,binsizes]=conncomp(Grrr);
[~,indx]=max(binsizes);
delCC=find(bins~=indx);
PList(delCC,:)=[];
DD(delCC,:)=[];
DD(:,delCC)=[]; 

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

inn  = inpolygon(PList(:,1), PList(:,2),xvv,yvv);

PBList=PList(inn,:);
ExtraList=PList(~inn,:);
figure(6)
scatter(PBList(:,1),PBList(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(ExtraList(:,1),ExtraList(:,2),'filled','r')

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



save('P1.mat')
