function out = Inverse_model(t,y,E,dt,Lambda)


%% Known parameters
sigma=30/3;
gamma=30/6.25;
g=1/(15*12);
d=1/(76*12);


%% Model variable assignment
I(1)=y(1);
R(1)=y(2);
A(1)=y(3);
N(1)=y(4);

for i=1:length(t)-1
    I(i+1)=I(i)+(sigma*E(i)-gamma*I(i)-g*I(i)-d*I(i))*dt;
    R(i+1)=R(i)+(gamma*I(i)-g*R(i)-d*R(i))*dt;
    A(i+1)=A(i)+(g*(N(i)-A(i))-d*A(i))*dt;
    N(i+1)=N(i)+(Lambda-d*N(i))*dt;
end

out=[I;R;A;N];