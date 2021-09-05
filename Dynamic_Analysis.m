clear
close all
clc

%%Parameter
%Mass
m= 1592;            %kg

%Spring and damper constant
k = 15000;             %N/m pegas roda depan
             %N/m pegas roda belakang
c = 0;             %Ns/m damper roda depan
             %Ns/m damper roda belakang

%Inertia of Body
I1xx = 477.60;
I1yy = 1326.67;
I1zz = 944.59;
Ihat1 = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];  %kgm2

h=0.01;     %step
t=0:h:5;    %time
iter=0; 

%Jarak titik pegas dan spring
ua1 = [-1.3;-0.3;0.9];               %Jarak ke lokal coordinate body 1
ub1 = [1.3;-0.3;0.9];
uc1 = [1.3;-0.3;-0.9];
ud1 = [-1.3;-0.3;-0.9];


g = 9.81;           %m/s^2
L0 = 0.25;           %m
Lawalf = L0+(m*g/(4*k));
Lawalr = L0+(m*g/(4*k));

%Initial Condition
q           =zeros(7,1);
qdt         =zeros(7,1);
qddt_lamda  =zeros(7,1);
q_dep       =0;

q(1)= 0;
q(2)= 0.55;
q(3)= 0;
q(4)= 0.000001;
q(5)= 0;
q(6)= 0;
q(7)= 0;
qdt(3)= 0.1;



%Data Size
DataSize    =size(t,2);
q_alltime   =zeros(7,DataSize);
v_alltime   =zeros(7,DataSize);
a_alltime   =zeros(7,DataSize);

for j=1:DataSize
           
    epsilon          =1e-12;
    delta_q_dep_norm =1;
    i_newton         =1;
    
    while abs(delta_q_dep_norm)>epsilon
        
        Config;
        
        delta_q_dep         =inv(Cq_dep)*(-C);
        delta_q_dep_norm    =norm(delta_q_dep);
        
        q_dep =q_dep+delta_q_dep;
        q     =[q(1);q(2);q(3);q_dep;q(5);q(6);q(7)];
        
        i_newton = i_newton+1;
        
        if i_newton > 100
           break
        end
    end
    
        qdt_dep =-(Cq_dep)\(Cq_ind*[qdt(1);qdt(2);qdt(3);qdt(5);qdt(6);qdt(7)]);
        qdt     =[qdt(1);qdt(2);qdt(3);qdt_dep;qdt(5);qdt(6);qdt(7)];
        
    Lamda;
        
    %Data tiap waktu
    q_alltime(:,j) = q;
    v_alltime(:,j) = qdt;
    a_alltime(:,j) = qddt_lamda(1:7);
        
        
    Runge_Kutta;

    iter = iter+1;   

end

figure(1)
plot(t,q_alltime(1,:))
title('Posisi X')
figure(2)
plot(t,q_alltime(2,:))
title('Posisi Y')
figure(3)
plot(t,q_alltime(3,:))
title('Posisi Z')
%figure(4)
%plot(t,q_alltime(5,:))
%figure(5)
%plot(t,q_alltime(6,:))
%figure(6)
%plot(t,q_alltime(7,:))
