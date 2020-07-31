%% RC_shandong.m is used for parameter estimation of shandong
clear all;
load baseloadshandong.mat

aircustomer=airsummer_es;
airtrue=aircustomer;
tempcustomer=tempsummer;

%% onoffstate segment
Pmax=max(max(aircustomer));
onoffcustomer=zeros(size(aircustomer,1),size(aircustomer,2));
onoffsumcustomer=zeros(1,size(aircustomer,2));
onoffprobability=zeros(1,size(aircustomer,2));

for i=1:size(aircustomer,1)
    for j=1:size(aircustomer,2)
        
            if aircustomer(i,j)>=Pmax/100
               onoffcustomer(i,j)=1;
            end
       
         
    end
end

%% onoffprobability

for i=1:1
    for j=1:size(aircustomer,2)
    onoffprobability(i,j)=sum(onoffcustomer(:,j))./size(aircustomer,1);
onoffsumcustomer(i,j)=onoffprobability(i,j)*size(aircustomer,3);
    end
end

figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,6])
boxplot(onoffsumcustomer)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Number of customers')

%% static-dynamic state segment
sdcustomer=zeros(size(aircustomer,1),size(aircustomer,2));

 for k=1:size(airtrue,1)
     
air_1st(k,:)=diff(airtrue(k,:));
temp_1st(k,:)=diff(tempcustomer(k,:));
airtemp(k,:)=airtrue(k,:)./tempcustomer(k,:);
airtemp_1st(k,:)=diff(airtemp(k,:));

if find(abs(airtemp_1st(k,:))<=4)
    
    %staticpoint{k}=intersect(intersect(find(abs(airtemp_1st(k,:))<=0.5),find(airtemp(k,:)>0)),find(onoffcustomer(k,:)==1));
    staticpoint{k}=intersect(intersect(find(abs(airtemp_1st(k,:))<=4),find(airtemp(k,:)>0)),find(onoffcustomer(k,:)==1));
  end

if ~isempty(intersect(intersect(find(abs(airtemp_1st(k,:))<=4),find(airtemp(k,:)>0)),find(onoffcustomer(k,:)==1)))
    
a1=[staticpoint{k}' staticpoint{k}'+1];
b1=[a1(1:end-1,2)+1 a1(2:end,1)-1 ];


if find((a1(:,1)-a1(:,2))>=0)
    a1(find((a1(:,1)-a1(:,2))>=0),:)=[];
end

if find((b1(:,1)-b1(:,2))>=0)
   b1(find((b1(:,1)-b1(:,2))>=0),:)=[];
end

for i=2:size(a1,1)
    if a1(i,1)<=a1(i-1,2)
        a1(i,1)=a1(i-1,2);
    end
end

for i=2:size(b1,1)
    if b1(i,1)<=b1(i-1,2)
        b1(i,1)=b1(i-1,2);
    end
end

if find(b1(:,2)>=24)
   b1(find(b1(:,2)>=24),:)=[];
end

if find(a1(:,2)>=24)
   a1(find(a1(:,2)>=24),:)=[];
end

c1=unique(a1(:,1));

for i=1:size(c1,1)
    sdcustomer(k,c1(i,1):c1(i,1)+1)=ones(1,2);
end

static{k}=a1;
dynamic{k}=b1;

end

 end
 %% onoffsd segment

onoffsdcustomer=zeros(size(airtrue,1),size(airtrue,2));
for i=1:size(airtrue,1)
 for j=1:size(airtrue,2)

    if (onoffcustomer(i,j)==1)&&(sdcustomer(i,j)==1)
        onoffsdcustomer(i,j)=1;
        
    end  
    if (onoffcustomer(i,j)==1)&&(sdcustomer(i,j)==0)
        onoffsdcustomer(i,j)=2;
    end 
    

 end
end

figure(2)
set(gcf,'unit','centimeters','position',[0,0,16,12])
bar(sum(sdcustomer)')
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Time (h)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Number of customers') 
a=onoffsdcustomer;
figure(3)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot(tempcustomer(find(a==1)),airtrue(find(a==1)),'ro','MarkerSize',4)
hold on
plot(tempcustomer(find(a==2)),airtrue(find(a==2)),'bo','MarkerSize',2)
hold on
plot(tempcustomer(find(a==0)),airtrue(find(a==0)),'go','MarkerSize',2)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Ambient Temperature (℃)')
ylabel('\fontsize{10.5}\fontname{Times new roman}ACLs (kW)')

b=onoffcustomer;
figure(4)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot(tempcustomer(find(b==1)),airtrue(find(b==1)),'ro','MarkerSize',4)
hold on
plot(tempcustomer(find(b==2)),airtrue(find(b==2)),'bo','MarkerSize',2)
hold on
plot(tempcustomer(find(b==0)),airtrue(find(b==0)),'go','MarkerSize',2)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Ambient Temperature (℃)')
ylabel('\fontsize{10.5}\fontname{Times new roman}ACLs (kW)')


%% static parameter of linear regression 1

rsindex=zeros(1,1);
RT=zeros(1,1);
for i=1:size(airtrue,3)
    
a=onoffsdcustomer;
temptotal=tempcustomer(find(a==1));
airtotal=airtrue;
airtotal=airtotal(find(a==1));
x=temptotal;
Y=airtotal;
X=[ones(size(x,1),1),x];
[b1,bint1,r1,rint1,s1]=regress(Y,X);

RT(i,1)=1./b1(2,1);
RT(i,2)=-b1(1,1)./b1(2,1);
rsindex(i,1)=s1(1);

end

% figure
% plot(temptotal,airtotal,'.');
airnew=airtotal;

%% static parameter of linear regression 2
rsindex1=zeros(1,1);
RT1=zeros(1,1);
for i=1:size(airtrue,3)
    
a=onoffcustomer;
temptotal1=tempcustomer(find(a==1));
airtotal=airtrue;
airtotal1=airtotal(find(a==1));
x1=temptotal1;
Y1=airtotal1;
X1=[ones(size(x1,1),1),x1];
[b2,bint2,r2,rint2,s2]=regress(Y1,X1);

RT1(i,1)=1./b2(2,1);
RT1(i,2)=-b2(1,1)./b2(2,1);
rsindex1(i,1)=s2(1);

end

%% static parameter of constrainted regression

RT2=zeros(1,1);
Tin2=zeros(size(airtrue,1),size(airtrue,2));
rsindex2=zeros(1,1);
fit2=cell(1,1);

ox=onoffcustomer;
[a,b]=find(ox==1);
tx=tempcustomer(find(ox==1));
ax=airtrue;
ax=ax(find(ox==1));
T=tx;
HVAC=ax;
T_shandong=T;
HVAC_shandong=HVAC;

%% optimization

%% decision variable
timescale=length(HVAC);
b1=sdpvar(1,timescale,'full');
k1=sdpvar(1,1,'full');
%% constraint boundary
Tmin=16*ones(1,timescale);
Tmax=27*ones(1,timescale);
kmin=0.0005;
kmax=0.002;
%% constraints
Constraints = [];

for i = 1:timescale
  Constraints = [Constraints, Tmin(1,i)*k1<= -b1(1,i) <=Tmax(1,i)*k1];
end

for i = 1:timescale
  Constraints = [Constraints, k1*kmin<= 1 <=k1*kmax];
end

for i = 2:timescale
  Constraints = [Constraints, (T(i,1)-T(i-1,1))*(b1(1,i)-b1(1,i-1))>=0 ];
  
end
%% objective

Objective = 0;
for i = 1:timescale

  Objective = Objective +(0.75*HVAC(i,1)-k1*T(i,1)-b1(1,i))*(0.75*HVAC(i,1)-k1*T(i,1)-b1(1,i));

end
ops = sdpsettings('solver','cplex');
ops.cplex= cplexoptimset('cplex');
ops.cplex.mip.tolerances.absmipgap = 0.01;
optimize(Constraints,Objective,ops);

RT2(1,1)=1./value(k1);
Tx=-value(b1)./value(k1); 
Tx=Tx';
Tb=Tin2;

for i=1:timescale
    
Tb(a(i),b(i))=Tx(i);
end

Tin=Tb;
Tinmean=zeros(1,1);
Tinmin=zeros(1,1);
Tinmax=zeros(1,1);
Tintotal=cell(1,1);
Touttotal=cell(1,1); 
tx=Tin;   
tx1=tx(find(tx~=0));
Tintotal{1,1}=tx1;
Touttotal{1,1}=tempcustomer(find(tx~=0));      
Tinmean(1,1)=mean(tx1);
Tinmin(1,1)=min(min(tx1));
Tinmax(1,1)=max(max(tx1));
RT2(1,2)=Tinmean(1,1);
Tx_shandong=Tx;

figure(5)
set(gcf,'unit','centimeters','position',[0,0,8,6])
boxplot(Tx_shandong)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Identification results of air conditioning setpoint temperature')
ylabel('\fontsize{8}\fontname{Times new roman}Temperature (℃)')

figure(6)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(Touttotal{1,1},Tintotal{1,1} ,'bo','MarkerSize',1)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Ambient Temperature (℃)')
ylabel('\fontsize{8}\fontname{Times new roman}Setpoint Temperature (℃)')
% k=9.65, yitaR=0.1036,Tinmean=21.09;

%%  C estimation

Cfinalmax=1.293*5*1e3*1000*max(max(airtrue))*1.2/(100*3.6e6);
Cfinalmin=1.293*5*1e3*1000*max(max(airtrue))*1/(100*3.6e6);

% h=5; Q=150;

%RC=6.787
Rfinal=RT2(1,1);
Tfinal=RT2(1,2);
timescale=size(tx,2);

%% PSO
%粒子群算法中的两个参数
c1 = 0.001;
c2 = 0.001;

Vmax=0.01;
Vmin=-0.01;

popmax=Cfinalmax;
popmin=Cfinalmin;

sizepop=100;
%% initialization
pop=zeros(sizepop,1,'single');
V=zeros(sizepop,1,'single');
fitness=zeros(1,sizepop,'single');

for i=1:sizepop
    
      pop(i)=Cfinalmin(1)+(Cfinalmax(1)-Cfinalmin(1))*rand(1,1);
  
    if rand(1,1)<0.5
        
        V(i,:)=0.01*rand(1,1);  %初始化速度
    else
        V(i,:)=-0.01*rand(1,1);
    end
    
% fitness
   Tin_base=zeros(size(Tin,1),size(Tin,2));
   Tout_base=tempcustomer;
   Tin_base(1,1)=Tout_base(1); %initialization
   R_base=Rfinal;
   M_base=exp(0.25./(R_base.*pop(i,1)));
   PACL_base=aircustomer;
  %% Tin simulation
for k=1:size(Tin,1)
    
if k==1
    for j=2:size(Tin,2)
    Tin_base(k,j)=Tin_base(k,j-1)./M_base+(1-1./M_base)*(Tout_base(k,j-1)-R_base*0.75*PACL_base(k,j-1));
    end
end 

if k~=1
    for j=1:size(Tin,2)
 if j==1       
    Tin_base(k,j)=Tin_base(k-1,size(Tin,2))./M_base+(1-1./M_base)*(Tout_base(k-1,size(Tin,2))-R_base*0.75*PACL_base(k-1,size(Tin,2)));
 end
 
 if j~=1       
    Tin_base(k,j)=Tin_base(k,j-1)./M_base+(1-1./M_base)*(Tout_base(k,j-1)-R_base*0.75*PACL_base(k,j-1));
 end 
    end
end  

end

Tins_base=Tin;

fitness(i)=sqrt(sum((Tin_base(find(ox==1))-Tins_base(find(ox==1))).^2.)/length(find(ox==1)));
        
end

%% optimization
[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);   %全局最佳
gbest=pop;    %个体最佳
fitnessgbest=fitness;   %个体最佳适应度值
fitnesszbest=bestfitness;   %全局最佳适应度值
maxgen=100;

%% iteration
yy=zeros(1,maxgen,'single');
xx=zeros(1,maxgen,'single');

for i=1:maxgen
    
    for j=1:sizepop
        %速度更新
         V(j,:) = V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        %种群更新
        pop(j,:)=pop(j,:)+V(j,:);
      
        if pop(j,1)>popmax
            pop(j,1)=popmax;
        end
        if pop(j,1)<popmin
            pop(j,1)=popmin;
        end 
        
 % fitness
   Tin_base=zeros(size(Tin,1),size(Tin,2));
   Tout_base=tempcustomer;
   Tin_base(1,1)=Tout_base(1); %initialization
   R_base=Rfinal;
   M_base=exp(0.25./(R_base.*pop(j,1)));
   PACL_base=aircustomer;
  %% Tin simulation
for k=1:size(Tin,1)
    
if k==1
    for s=2:size(Tin,2)
    Tin_base(k,s)=Tin_base(k,s-1)./M_base+(1-1./M_base)*(Tout_base(k,s-1)-R_base*0.75*PACL_base(k,s-1));
    end
end 

if k~=1
    for s=1:size(Tin,2)
 if s==1       
    Tin_base(k,s)=Tin_base(k-1,size(Tin,2))./M_base+(1-1./M_base)*(Tout_base(k-1,size(Tin,2))-R_base*0.75*PACL_base(k-1,size(Tin,2)));
 end
 
 if s~=1       
    Tin_base(k,s)=Tin_base(k,s-1)./M_base+(1-1./M_base)*(Tout_base(k,s-1)-R_base*0.75*PACL_base(k,s-1));
 end 
    end
end  

end

Tins_base=Tin;

fitness(j)=sqrt(sum((Tin_base(find(ox==1))-Tins_base(find(ox==1))).^2.)/length(find(ox==1)));        
     
    end
    
    for j=1:sizepop
        
        %个体最优更新
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %群体最优更新
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);     %%群体最优粒子
            fitnesszbest = fitness(j);
        end
    end 
    
    
    yy(i)=fitnesszbest;     %%最优值对应的适应度值
    xx(i)=zbest;
    
end
% C/yita=68.5991;

figure(7)
set(gcf,'unit','centimeters','position',[0,0,8,6])

plot(Tin_base(find(ox==1)),'r-.','linewidth',1)
hold on
plot(Tins_base(find(ox==1)),'b-','linewidth',1)

set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Point')
ylabel('\fontsize{8}\fontname{Times new roman}Temperature (℃)')
ylim([10 35])

b=onoffcustomer;
figure(8)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(tempcustomer(find(b==1)),airtrue(find(b==1)),'bo','MarkerSize',2)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Ambient Temperature (℃)')
ylabel('\fontsize{8}\fontname{Times new roman}ACLs (kW)')

R_shandong=Rfinal;
Tin_shandong=Tfinal;
C_shandong=zbest;
Cmax_shandong=Cfinalmax;
Cmin_shandong=Cfinalmin;

result_shandong=[R_shandong Tin_shandong RT1(1,1) RT1(1,2) C_shandong Cmax_shandong Cmin_shandong];



save RC_shandong.mat

tr=rand(size(temptotal,1),1);
tr=tr/4;
temptotal=temptotal+tr;


figure(9)
set(gcf,'unit','centimeters','position',[0,0,16,8])
subplot(1,2,1)
a=plot(temptotal,airnew,'.');
set(gca,'FontName','Times New Roman','FontSize',8)

hold on
x1=20.6318:0.01:40;
tempa=1/RT2(1,1);
y1=tempa*(x1-RT2(1,2));
b=plot(x1,y1);
hold on
%c=plot('k=833.3 KW/℃','Tset=20.63℃');
%text(30,1,'k=833.3 KW/℃','FontName','Times New Roman');

legend([a(1),b(1)],'Static Load data','Fitted Line','fontname','Times new roman','fontsize',8);
xlabel('\fontsize{8}\fontname{Times new roman}Temperature(℃)')
ylabel('\fontsize{8}\fontname{Times new roman}Load(KW)')
xlim([10 40])
%legend()

subplot(1,2,2)

boxplot(Tx_shandong)
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel({'\fontsize{8}\fontname{Times new roman}Identification results of ';'\fontsize{8}\fontname{Times new roman}air conditioning setpoint temperature'})
ylabel('\fontsize{8}\fontname{Times new roman}Temperature (℃)')
 