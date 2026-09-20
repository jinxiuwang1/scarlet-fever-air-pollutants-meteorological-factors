function out = Model(t,y,beta,Lambda,dt)


%% Known parameters
sigma=30/3;
gamma=30/6.25;
g=1/(15*12);
d=1/(76*12);


%% Model variable assignment
S(1)=y(1);
E(1)=y(2);
I(1)=y(3);
R(1)=y(4);
A(1)=y(5);
N(1)=S(1)+E(1)+I(1)+R(1)+A(1);

for i=1:length(t)-1
    S(i+1)=S(i)+(Lambda-beta(i)*S(i)*I(i)/N(i)-(d+g)*S(i))*dt;
    E(i+1)=E(i)+(beta(i)*S(i)*I(i)/N(i)-sigma*E(i)-(d+g)*E(i))*dt;
    I(i+1)=I(i)+(sigma*E(i)-gamma*I(i)-g*I(i)-d*I(i))*dt;
    R(i+1)=R(i)+(gamma*I(i)-g*R(i)-d*R(i))*dt;
    A(i+1)=A(i)+(g*(S(i)+E(i)+I(i)+R(i))-d*A(i))*dt;
    N(i+1)=S(i+1)+E(i+1)+I(i+1)+R(i+1)+A(i+1);
end

out=[S;E;I;R;A];