clc;
clear;
close all

n=5; 
D=[[cos(0);cos(2*pi*(1/5));cos(2*pi*(2/5));cos(2*pi*(3/5));cos(2*pi*(4/5))] [sin(0);sin(2*pi*(1/5));sin(2*pi*(2/5));sin(2*pi*(3/5));sin(2*pi*(4/5))]];
T=[[cos(0);cos(4*pi*(1/5));cos(4*pi*(2/5));cos(4*pi*(3/5));cos(4*pi*(4/5))] [sin(0);sin(4*pi*(1/5));sin(4*pi*(2/5));sin(4*pi*(3/5));sin(4*pi*(4/5))] [1;1;1;1;1]];
Offset=6;%the width between periodic Box and search space (5 is good)
order1=4;%enter the value of approximation for golen ratio in x direction Lx_Box=5*(taw^order1)
order2=6;%enter the value of approximation for golen ratio in y direction Ly_Box=(2*sin(pi/5))*(taw^(order2+1))
%if you want a box with same lengths in x and y directions the following relation must holds: order2=order1+2
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
if size(reqcyc,1)~=0
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
end
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

FG=graph(DD);
cycles5 = allcycles(FG,'MinCycleLength',5,'MaxCycleLength',5);
reqcyc5=[];
NewList=[];
counter=1;
for i=1:size(cycles5,1)
    fiveP=cell2mat(cycles5(i));
    Cords=PList(fiveP,:);
    check=checkConvex(Cords(:,1)',Cords(:,2)');
    if check==1
        NewList(counter,:)=mean(PList(fiveP,:));
        reqcyc5(counter,:)=i;
        counter=counter+1;
    end
end
DDnew=zeros(size(NewList,1));
DDSS=cell2mat(cycles5(reqcyc5));
%scatter(NewList(:,1),NewList(:,2),'filled')
for i=1:size(NewList,1)
    for j=i+1:size(NewList,1)
        if size(intersect(DDSS(i,:),DDSS(j,:)),2)==2
            DDnew(i,j)=1;
            DDnew(j,i)=1;
        end
    end
end


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
figure(9)
gplot(DDnew,NewList,'k-')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
[XN,YN]=gplot(DDnew,NewList);
figure(10)
plot(XN,YN,'-mo',...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',4);
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
iinn  = inpolygon(NewList(:,1), NewList(:,2),xvv,yvv);

PBnewList=NewList(iinn,:);
ExtranewList=NewList(~iinn,:);
figure(11)
gplot(DDnew(iinn,iinn),NewList(iinn,:),'k-')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
[XXN,YYN]=gplot(DDnew(iinn,iinn),NewList(iinn,:));
figure(12)
plot(XXN,YYN,'-mo',...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',4);
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
figure(13)
scatter(PBnewList(:,1),PBnewList(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(ExtranewList(:,1),ExtranewList(:,2),'filled','r')
%save('P1.mat')





xvv=[AA(1);AA(1)+Period1(1);AA(1)+Period1(1);AA(1);AA(1)];
yvv=[AA(2);AA(2);AA(2)+Period2(2);AA(2)+Period2(2);AA(2)];
epsln=0.01;
PXout=[AA(1)-epsln;AA(1)+Period1(1)+epsln;AA(1)+Period1(1)+epsln;AA(1)-epsln;AA(1)-epsln];
PYout=[AA(2)-epsln;AA(2)-epsln;AA(2)+Period2(2)+epsln;AA(2)+Period2(2)+epsln;AA(2)-epsln];
PXin=[AA(1)+epsln;AA(1)+Period1(1)-epsln;AA(1)+Period1(1)-epsln;AA(1)+epsln;AA(1)+epsln];
PYin=[AA(2)+epsln;AA(2)+epsln;AA(2)+Period2(2)-epsln;AA(2)+Period2(2)-epsln;AA(2)+epsln];

inPin=inpolygon(PList(:,1), PList(:,2),PXin,PYin);
inPout=inpolygon(PList(:,1), PList(:,2),PXout,PYout);
inPboundary=logical(inPout-inPin);


PListin=PList(inPin,:);
PListboundary=PList(inPboundary,:);
figure(14)
scatter(PListin(:,1),PListin(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(PListboundary(:,1),PListboundary(:,2),'filled','r')


FALeft=find(PList(:,1)>AA(1)-epsln & PList(:,1)<AA(1)+epsln);
FARight=find(PList(:,1)>AA(1)+Period1(1)-epsln & PList(:,1)<AA(1)+Period1(1)+epsln);
if size(FARight,1)<size(FALeft,1)
    change=PList(FARight,:)-Period1;
    UR=PList(FALeft,:);
    delD=[];
    counter=1;
    for i=1:size(UR,1)
        TFT=find(abs(UR(i,2)-change(:,2))<1e-5);
        if size(TFT,1)<1
            delD(counter,:)=i;
            counter=counter+1;
        end
        inPboundary(FARight,:)=0;
        inPboundary(FALeft(delD),:)=0;   
    end
    
elseif size(FARight,1)>size(FALeft,1)
    inPboundary(FARight,:)=0;
else
    inPboundary(FARight,:)=0;
end

FADown=find(PList(:,2)>AA(2)-epsln & PList(:,2)<AA(2)+epsln);
FAUp=find(PList(:,2)>AA(2)+Period2(2)-epsln & PList(:,2)<AA(2)+Period2(2)+epsln);
if size(FADown,1)<size(FAUp,1)
    inPboundary(FAUp,:)=0;
elseif size(FADown,1)>size(FAUp,1)
    change1=PList(FAUp,:)-Period2;
    UR1=PList(FADown,:);
    delD1=[];
    counter=1;
    for i=1:size(UR1,1)
        TFT1=find(abs(UR1(i,2)-change1(:,2))<1e-5);
        if size(TFT1,1)<1
            delD1(counter,:)=i;
            counter=counter+1;
        end
        inPboundary(FAUp,:)=0;   
        inPboundary(FAdown(delD1),:)=0;   
    end
else
    inPboundary(FAUp,:)=0;
end

inPBOX=logical(inPin+inPboundary);
PListinBOX=PList(inPBOX,:);
PListoutBOX=PList(~inPBOX,:);
figure(15)
scatter(PListinBOX(:,1),PListinBOX(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(PListoutBOX(:,1),PListoutBOX(:,2),'filled','r')

effEdge=find(PList(:,1)>AA(1)-epsln & PList(:,2)>AA(2)-epsln);
ineffEdge=zeros(size(PList,1),1);
ineffEdge(effEdge)=1;
ineffEdge=logical(ineffEdge);
infEdge=and(ineffEdge,~inPBOX);
figure(16)
scatter(PList(inPBOX,1),PList(inPBOX,2),'filled','b')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(PList(infEdge,1),PList(infEdge,2),'filled','r')

PEdges=[];
counter=1;
for i=1:size(DD,1)
    for j=i+1:size(DD,1)
       if DD(i,j)==1 && inPBOX(i)==1
           iPer=find(abs(PListinBOX(:,1)-PList(i,1))<1e-5 & abs(PListinBOX(:,2)-PList(i,2))<1e-5);
           if inPBOX(j)==1
               jPer=find(abs(PListinBOX(:,1)-PList(j,1))<1e-5 & abs(PListinBOX(:,2)-PList(j,2))<1e-5);
               PEdges(counter,:)=[iPer jPer 0 0];
               counter=counter+1;
           elseif inPBOX(j)==0 && infEdge(j)==1
               C1=floor((PList(j,1)-AA(1))/Period1(1));
               C2=floor((PList(j,2)-AA(2))/Period2(2));
               if abs((PList(j,2)-AA(2))/Period2(2)-(C2+1))<1e-5 
                   C2=C2+1;
               end
               if abs((PList(j,2)-AA(2))/Period2(2)-(C2-1))<1e-5
                   C2=C2-1;
               end
               if abs((PList(j,1)-AA(1))/Period1(1)-(C1+1))<1e-5
                   C1=C1+1;
               end
               if abs((PList(j,1)-AA(1))/Period1(1)-(C1-1))<1e-5
                   C1=C1-1;
               end 

               W1=PList(j,1)-C1*Period1(1);
               W2=PList(j,2)-C2*Period2(2);
               WFW=find(abs(PListinBOX(:,1)-W1)<1e-5 & abs(PListinBOX(:,2)-W2)<1e-5);
               PEdges(counter,:)=[iPer WFW C1 C2];
               counter=counter+1;
           end
       end
    end
end







inPNin=inpolygon(NewList(:,1), NewList(:,2),PXin,PYin);
inPNout=inpolygon(NewList(:,1), NewList(:,2),PXout,PYout);
inPNboundary=logical(inPNout-inPNin);


NewListin=NewList(inPNin,:);
NewListboundary=NewList(inPNboundary,:);
figure(17)
scatter(NewListin(:,1),NewListin(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(NewListboundary(:,1),NewListboundary(:,2),'filled','r')


FANLeft=find(NewList(:,1)>AA(1)-epsln & NewList(:,1)<AA(1)+epsln);
FANRight=find(NewList(:,1)>AA(1)+Period1(1)-epsln & NewList(:,1)<AA(1)+Period1(1)+epsln);
if size(FANRight,1)<size(FANLeft,1)
    changeN=NewList(FANRight,:)-Period1;
    URN=NewList(FANLeft,:);
    delDN=[];
    counter=1;
    for i=1:size(URN,1)
        TFT=find(abs(URN(i,2)-changeN(:,2))<1e-5);
        if size(TFT,1)<1
            delDN(counter,:)=i;
            counter=counter+1;
        end
        inPNboundary(FANRight,:)=0;
        inPNboundary(FANLeft(delDN),:)=0;   
    end
    
elseif size(FANRight,1)>size(FANLeft,1)
    inPNboundary(FANRight,:)=0;
else
    inPNboundary(FANRight,:)=0;
end

FANDown=find(NewList(:,2)>AA(2)-epsln & NewList(:,2)<AA(2)+epsln);
FANUp=find(NewList(:,2)>AA(2)+Period2(2)-epsln & NewList(:,2)<AA(2)+Period2(2)+epsln);
if size(FANDown,1)<size(FANUp,1)
    inPNboundary(FANUp,:)=0;
elseif size(FANDown,1)>size(FANUp,1)
    changeN1=NewList(FANUp,:)-Period2;
    URN1=NewList(FANDown,:);
    delDN1=[];
    counter=1;
    for i=1:size(URN1,1)
        TFT1=find(abs(URN1(i,2)-changeN1(:,2))<1e-5);
        if size(TFT1,1)<1
            delDN1(counter,:)=i;
            counter=counter+1;
        end
        inPNboundary(FANUp,:)=0;   
        inPNboundary(FANDown(delDN1),:)=0;   
    end
else
    inPNboundary(FANUp,:)=0;
end

inPNBOX=logical(inPNin+inPNboundary);
NewListinBOX=NewList(inPNBOX,:);
NewListoutBOX=NewList(~inPNBOX,:);
figure(18)
scatter(NewListinBOX(:,1),NewListinBOX(:,2),'filled')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(NewListoutBOX(:,1),NewListoutBOX(:,2),'filled','r')

effNEdge=find(NewList(:,1)>AA(1)-epsln & NewList(:,2)>AA(2)-epsln);
ineffNEdge=zeros(size(NewList,1),1);
ineffNEdge(effNEdge)=1;
ineffNEdge=logical(ineffNEdge);
infNEdge=and(ineffNEdge,~inPNBOX);
figure(19)
scatter(NewList(inPNBOX,1),NewList(inPNBOX,2),'filled','b')
hold on
rectangle('Position',[AA Period1(1) Period2(2)],'linewidth',2)
scatter(NewList(infNEdge,1),NewList(infNEdge,2),'filled','r')

PNEdges=[];
counter=1;
for i=1:size(DDnew,1)
    for j=i+1:size(DDnew,1)
       if DDnew(i,j)==1 && inPNBOX(i)==1
           iPer=find(abs(NewListinBOX(:,1)-NewList(i,1))<1e-5 & abs(NewListinBOX(:,2)-NewList(i,2))<1e-5);
           if inPNBOX(j)==1
               jPer=find(abs(NewListinBOX(:,1)-NewList(j,1))<1e-5 & abs(NewListinBOX(:,2)-NewList(j,2))<1e-5);
               PNEdges(counter,:)=[iPer jPer 0 0];
               counter=counter+1;
           elseif inPNBOX(j)==0 && infNEdge(j)==1
               C1=floor((NewList(j,1)-AA(1))/Period1(1));
               C2=floor((NewList(j,2)-AA(2))/Period2(2));
               if abs((NewList(j,2)-AA(2))/Period2(2)-(C2+1))<1e-5 
                   C2=C2+1;
               end
               if abs((NewList(j,2)-AA(2))/Period2(2)-(C2-1))<1e-5
                   C2=C2-1;
               end
               if abs((NewList(j,1)-AA(1))/Period1(1)-(C1+1))<1e-5
                   C1=C1+1;
               end
               if abs((NewList(j,1)-AA(1))/Period1(1)-(C1-1))<1e-5
                   C1=C1-1;
               end 

               W1=NewList(j,1)-C1*Period1(1);
               W2=NewList(j,2)-C2*Period2(2);
               WFW=find(abs(NewListinBOX(:,1)-W1)<1e-5 & abs(NewListinBOX(:,2)-W2)<1e-5);
               PNEdges(counter,:)=[iPer WFW C1 C2];
               counter=counter+1;
           end
       end
    end
end


TypeP = ones(size(PListinBOX,1),1);
outptFile = fopen('out1.txt','W');
fprintf(outptFile,'%d Vertices\n\n%d Edges\n\n%d Vertex types\n\nBOX\n\n%f %f xlo xhi\n%f %f ylo yhi \n\n',...
                    size(PListinBOX,1),size(PEdges,1),size(unique(TypeP),1),AA(1),AA(1)+Period1(1),AA(2),AA(2)+Period2(2));
fprintf(outptFile,'Vertices\n\n');
fprintf(outptFile,'%d %d %f %f\n',[1:size(PListinBOX,1);TypeP';PListinBOX(:,1)';PListinBOX(:,2)']);
fprintf(outptFile,'\nEdges\n\n');
fprintf(outptFile,'%d %d %d %d %d\n',[1:size(PEdges,1);PEdges(:,1)';PEdges(:,2)';PEdges(:,3)';PEdges(:,4)']);
fclose(outptFile);

TypePN = ones(size(NewListinBOX,1),1);
outptFile1 = fopen('out2.txt','W');
fprintf(outptFile1,'%d Vertices\n\n%d Edges\n\n%d Vertex types\n\nBOX\n\n%f %f xlo xhi\n%f %f ylo yhi \n\n',...
                    size(NewListinBOX,1),size(PNEdges,1),size(unique(TypePN),1),AA(1),AA(1)+Period1(1),AA(2),AA(2)+Period2(2));
fprintf(outptFile1,'Vertices\n\n');
fprintf(outptFile1,'%d %d %f %f\n',[1:size(NewListinBOX,1);TypePN';NewListinBOX(:,1)';NewListinBOX(:,2)']);
fprintf(outptFile1,'\nEdges\n\n');
fprintf(outptFile1,'%d %d %d %d %d\n',[1:size(PNEdges,1);PNEdges(:,1)';PNEdges(:,2)';PNEdges(:,3)';PNEdges(:,4)']);
fclose(outptFile1);

save('P1.mat')
