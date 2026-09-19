function [Gansu_SF,BETA,new_cases,dt,n,time,SS,NN] = function_beta(Label,AA,BB,i,Lambda_new,Total_P,Total_P_0_14)
% Find the row indices where the 16th column contains 'Gansu'
row_indices = find(cellfun(@(x) contains(x, Label{i}, 'IgnoreCase', true), BB(:,16)));
Gansu_SF=AA(row_indices-1,17)+1e-3;
n=length(Gansu_SF);
time=[1:1:n]; % number of months

dt=0.001;
[YY, dy] = B_spline_der(time,log(Gansu_SF),dt);
yy_der=exp(YY).*dy;
yy=exp(YY);

%% Compute E(t)
Lambda=Lambda_new(i);
sigma=30/3;
gamma=30/6.25;
g=1/(15*12);
d=1/(76*12);
E=yy/sigma;
E_der=yy_der/sigma;

%% Solve the discrete inverse equation
y0 =[Gansu_SF(1),0,Total_P(i)-Total_P_0_14(i),Total_P(i)]; % initial values
y= Inverse_model([1:dt:n],y0,E,dt,Lambda);

I=y(1,:); 
R=y(2,:);
A=y(3,:);
N=y(4,:);
S=N-I-E-R-A;

NN=N(1:1/dt:end);
SS=S(1:1/dt:end);


%% Compute beta(t)
% for i=1:length([1:dt:n])
%     beta(i)=max(0,-(S_der(i)-(Lambda-(g+d)*S(i)))*N(i)/(S(i)*I(i)));
% end

for j=1:length([1:dt:n])
    beta(j)=max(0,(E_der(j)+sigma*E(j)+d*E(j)+g*E(j))*N(j)/(S(j)*I(j)));
end

BETA=beta(1:1/dt:end);

%% Validate by solving the differential equation
y00 =[Total_P_0_14(i),E(1),Gansu_SF(1),0,Total_P(i)-Total_P_0_14(i)]; % initial values
out = Model([1:dt:n],y00,beta,Lambda,dt);
new_cases=sigma*out(2,:);

end