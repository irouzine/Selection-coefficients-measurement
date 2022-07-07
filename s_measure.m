function sdis = s_measure(genome1,genome2,genome3,C,appr,apprR)
global tsec1 tsec2 tsec3 r s0 M L N tf f0 run run2 mu ac bc
 
%Uncomment if you want to import data from name2.dat
%Gdata = {'genome1','genome2','genome3'};
%load('name2.dat',Gdata{:});

scon = {};
[~,l] = size(genome1{1,1});
%l - quantity of the alleles 
[R,~] = size(genome1);
%R - quantity of runs
f1 = zeros(R,l); % array for fi at 1 time (for R runs)
f2 = zeros(R,l); % array for fi at 2 time (for R runs)
f3 = zeros(R,l); % array for fi at 3 time (for R runs)
eps = 10^-3; % difference between zero and value of the y coordinate of the triangle center, when the loop stops
yccord = 1; % basic value of the y coordinate of the triangle center
step =10^-3; % value of the step
findc = true;


% calculation of fi

for rs = 1:R
        f1(rs,:)=mean(genome1{rs,1});
    f2(rs,:)=mean(genome2{rs,1});
    f3(rs,:)=mean(genome3{rs,1});
end

    %{
for i = 1:l
    pr =0; 
    pr2 =0; 
    pr3 =0; 
   for k = 1:ns
       pr = pr +  genome1{rs,1}(k,i);
       pr2 = pr2 + genome2{rs,1}(k,i);
       pr3 = pr3 + genome3{rs,1}(k,i);
   end
   f1(rs,i) = pr/ns;
   f2(rs,i) = pr2/ns;
   f3(rs,i) = pr3/ns;
end
    %}

%finding average of each fi
%{
for j = 1:l
    fpr1 = 0;
    fpr2 = 0;
    fpr3 = 0;
for r2 = 1:R
    fpr1 = fpr1 + f1(r2,j);
    fpr2 = fpr2 + f2(r2,j);
    fpr3 = fpr3 + f3(r2,j);
end
f1f(j)= fpr1/R;
f2f(j)= fpr2/R;
f3f(j)= fpr3/R;
end
%}
f1f = mean(f1);
f2f = mean(f2);
f3f = mean(f3);
%save(strcat('fi1','.mat'),'f1f');


%fnorm =5;
while abs(yccord) > eps
%calculation of si
bsi1 = -1*log(f1f)-C;
bsi2 = -1*log(f2f)-C;
bsi3 = -1*log(f3f)-C;
%sorting and calculating
[B1,I1] = sort(bsi1,'descend');
[B2,I2] = sort(bsi2,'descend');
[B3,I3] = sort(bsi3,'descend');



xdots = {};
ydots = {};
xdotss = {};
ydotss = {};
switch appr
    case 'poly'
%polynomial approximation
p1 = polyfit([1:1:length(I1)],B1,apprR);
p2 = polyfit([1:1:length(I1)],B2,apprR);
p3 = polyfit([1:1:length(I1)],B3,apprR);
inter1 = p1 - p2;
inter2 = p1 - p3;
inter3 = p2 - p3;
xdots{1,1} = roots(inter1);
xdots{1,2} = roots(inter2);
xdots{1,3} = roots(inter3);
ydots{1,1} = polyval(p1,xdots{1,1});
ydots{1,2} = polyval(p1,xdots{1,2});
ydots{1,3} = polyval(p2,xdots{1,3});
ybord{1,1} = polyval(p1,1);
ybord{1,2} = polyval(p2,1);
ybord{1,3} = polyval(p3,1);
ybord{2,1} = polyval(p1,length(I1));
ybord{2,2} = polyval(p2,length(I1));
ybord{2,3} = polyval(p3,length(I1));
    case 'spline'
%spline approximation
sp1 = spline([1:1:length(I1)],B1);
sp2 = spline([1:1:length(I1)],B2);
sp3 = spline([1:1:length(I1)],B3);

sinter1 = @(x) ppval(sp1,x)-ppval(sp2,x);
sinter2 = @(x) ppval(sp1,x)-ppval(sp3,x);
sinter3 = @(x) ppval(sp2,x)-ppval(sp3,x);
xdots{1,1} = fzero(sinter1,mean([1:1:length(I1)]));
xdots{1,2} = fzero(sinter2,mean([1:1:length(I1)]));
xdots{1,3} = fzero(sinter3,mean([1:1:length(I1)]));
ydots{1,1} = ppval(sp1,xdots{1,1});
ydots{1,2} = ppval(sp1,xdots{1,2});
ydots{1,3} = ppval(sp2,xdots{1,3});
ybord{1,1} = ppval(sp1,1);
ybord{1,2} = ppval(sp2,1);
ybord{1,3} = ppval(sp3,1);
ybord{2,1} = ppval(sp1,length(I1));
ybord{2,2} = ppval(sp2,length(I1));
ybord{2,3} = ppval(sp3,length(I1));
    case 'pchip'
%pchip approximation   
pp1 = pchip([1:1:length(I1)],B1);
pp2 = pchip([1:1:length(I1)],B2);
pp3 = pchip([1:1:length(I1)],B3);
pinter1 = @(x) ppval(pp1,x)-ppval(pp2,x);
pinter2 = @(x) ppval(pp1,x)-ppval(pp3,x);
pinter3 = @(x) ppval(pp2,x)-ppval(pp3,x);
xdots{1,1} = fzero(pinter1,mean([1:1:length(I1)]));
xdots{1,2} = fzero(pinter2,mean([1:1:length(I1)]));
xdots{1,3} = fzero(pinter3,mean([1:1:length(I1)]));       
ydots{1,1} = ppval(pp1,xdots{1,1});
ydots{1,2} = ppval(pp1,xdots{1,2});
ydots{1,3} = ppval(pp2,xdots{1,3});
ybord{1,1} = ppval(pp1,1);
ybord{1,2} = ppval(pp2,1);
ybord{1,3} = ppval(pp3,1);
ybord{2,1} = ppval(pp1,length(I1));
ybord{2,2} = ppval(pp2,length(I1));
ybord{2,3} = ppval(pp3,length(I1));
    case "test"
findc=false;
yccord = 0;
C = 0;
end

if findc
%coordinates checking
for times = 1: length(xdots)
    ncount1 = 0;
for param1= 1: length(xdots{1,times})
    if xdots{1,times}(param1)<=length(I1) && xdots{1,times}(param1)>=1 && imag(xdots{1,times}(param1))==0 && ydots{1,times}(param1)<=ybord{1,times} && ydots{1,times}(param1)>=ybord{2,times}
       ncount1 = ncount1 +1; 
        xdotss{times,ncount1} = xdots{1,times}(param1);
        ydotss{times,ncount1} = ydots{1,times}(param1);
    end
end
%&& polyval(p1,xdots1(param1))<=max(polyval(p1,[1:1:length(I1)]))&& polyval(p1,xdots1(param1))>=max(-1*polyval(p1,[1:1:length(I1)]))
end
%coordinates of triangle tops
xdcord = [];
ydcord = [];
[shr,dl] = size(xdotss);
for h =1:shr
    for c = 1:dl
        xdcord = [xdcord xdotss{h,c}];
        ydcord = [ydcord ydotss{h,c}];
    end
end


%xdcord = [xdots1sorted xdots2sorted xdots3sorted];
%ydcord = [ydots1 ydots2 ydots3];
%center finding
polyin = polyshape({xdcord},{ydcord});
[xccord,yccord] = centroid(polyin);

if yccord>0
    C = C+ step;
end
if yccord<0
    C = C - step;
end

disp(yccord);
end
end

%plotting the Srank curves

plot([1:1:length(I1)],B1,'m-');
title('Srank curves')
xlabel('m_i')
ylabel('Bs_i')
hold on
plot([1:1:length(I2)],B2,'b-');
hold on
plot([1:1:length(I3)],B3,'r-');
txt1 = sprintf(' t=%g',tsec1);
text(1,B1(1),strcat('\leftarrow',txt1));
txt2 = sprintf(' t=%g',tsec2);
text(1,B2(1),strcat('\leftarrow',txt2));
txt3 = sprintf(' t=%g',tsec3);
text(1,B3(1),strcat('\leftarrow',txt3));
txt4 = sprintf('r=%g,s0=%g,M=%g,L=%g,N=%g,tf=%g,f0=%g,runs=%g,attempt=%g,MuL=%g,ac=%g,bc=%g',r,s0,M,L,N,tf,f0,run,run2,mu*L,ac,bc);
%text(90,B3(90)+3,txt4)
title(strcat('Srank curves',txt4));
hold on
%plot([1:1:length(I1)],polyval(p1,[1:1:length(I1)]),'g-',[1:1:length(I1)],polyval(p2,[1:1:length(I1)]),'g-',[1:1:length(I1)],polyval(p3,[1:1:length(I1)]),'g-')
hold on
%plot(polyin);
%plot(xdotss{1,1},ydotss{1,1},'ko',xdotss{1,2},ydotss{1,2},'ko',xdotss{1,3},ydotss{1,3},'ko');
hold on
if findc
plot(xccord,yccord,'ko');
hold on
plot([1:1:length(I1)],zeros(1,l),'k--');

end
scon{1,4} = 0;
if findc
scon{1,4} = xccord;
end
scon{1,1} = bsi1;
scon{1,2} = bsi2;
scon{1,3} = bsi3;
%plotting the bars (Figure 4)
figure
bar(bsi1);
txt5 = strcat(txt4,txt1);
title(strcat('Bsi(i)',txt5));
figure
bar(bsi2);
txt5 = strcat(txt4,txt2);
title(strcat('Bsi(i)',txt5));
figure
bar(bsi3);
txt5 = strcat(txt4,txt3);
title(strcat('Bsi(i)',txt5));
sdis=scon;