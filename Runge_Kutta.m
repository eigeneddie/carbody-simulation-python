%Runge kutta integration
q_ind=[q(1);q(2);q(3);q(5);q(6);q(7)];
qdt_ind=[qdt(1);qdt(2);qdt(3);qdt(5);qdt(6);qdt(7)];
gi=[qddt_lamda(1);qddt_lamda(2);qddt_lamda(3);qddt_lamda(5);qddt_lamda(6);qddt_lamda(7)];
X=[q_ind; qdt_ind];

f=[X(7:12);gi];
k1=h*f;
t1=t+0.5*h;
X1=X+0.5*k1;
q(1)=X1(1);
q(2)=X1(2);
q(3)=X1(3);
q(5)=X1(4);
q(6)=X1(5);
q(7)=X1(6);
qdt(1)=X1(7);
qdt(2)=X1(8);
qdt(3)=X1(9);
qdt(5)=X1(10);
qdt(6)=X1(11);
qdt(7)=X1(12);


Config;
Lamda;

q_ind=[q(1);q(2);q(3);q(5);q(6);q(7)];
qdt_ind=[qdt(1);qdt(2);qdt(3);qdt(5);qdt(6);qdt(7)];
gi=[qddt_lamda(1);qddt_lamda(2);qddt_lamda(3);qddt_lamda(5);qddt_lamda(6);qddt_lamda(7)];

f1=[X1(7:12);gi];
k2=h*f1;
t2=t+0.5*h;
X2=X+0.5*k2;
q(1)=X2(1);
q(2)=X2(2);
q(3)=X2(3);
q(5)=X2(4);
q(6)=X2(5);
q(7)=X2(6);
qdt(1)=X2(7);
qdt(2)=X2(8);
qdt(3)=X2(9);
qdt(5)=X2(10);
qdt(6)=X2(11);
qdt(7)=X2(12);


Config;
Lamda;

q_ind=[q(1);q(2);q(3);q(5);q(6);q(7)];
qdt_ind=[qdt(1);qdt(2);qdt(3);qdt(5);qdt(6);qdt(7)];
gi=[qddt_lamda(1);qddt_lamda(2);qddt_lamda(3);qddt_lamda(5);qddt_lamda(6);qddt_lamda(7)];

f2=[X2(7:12);gi];
k3=h*f2;
t3=t+0.5*h;
X3=X+0.5*k3;
q(1)=X3(1);
q(2)=X3(2);
q(3)=X3(3);
q(5)=X3(4);
q(6)=X3(5);
q(7)=X3(6);
qdt(1)=X3(7);
qdt(2)=X3(8);
qdt(3)=X3(9);
qdt(5)=X3(10);
qdt(6)=X3(11);
qdt(7)=X3(12);


Config;
Lamda;

q_ind=[q(1);q(2);q(3);q(5);q(6);q(7)];
qdt_ind=[qdt(1);qdt(2);qdt(3);qdt(5);qdt(6);qdt(7)];
gi=[qddt_lamda(1);qddt_lamda(2);qddt_lamda(3);qddt_lamda(5);qddt_lamda(6);qddt_lamda(7)];

f3=[X3(7:12);gi];
k4=h*f3;
X=X+1/6*(k1+2*k2+2*k3+k4);
q(1)=X(1);
q(2)=X(2);
q(3)=X(3);
q(5)=X(4);
q(6)=X(5);
q(7)=X(6);
qdt(1)=X(7);
qdt(2)=X(8);
qdt(3)=X(9);
qdt(5)=X(10);
qdt(6)=X(11);
qdt(7)=X(12);