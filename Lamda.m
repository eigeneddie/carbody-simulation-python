%Kalkulasi untuk acceleration lamda
MCq=[M Cq'; Cq 0];                                         %Matrix M dan Cq

%%
%%Ketidakrataan Jalan

%if j>max(size(t))/2
        
    
   ktdf = 0;
   ktdr = 0;
   
   ua0 = [-1.3;ktdr;0.9];                  %Jarak ke lokal coordinate ground
   ub0 = [1.3;ktdf;0.9];
   uc0 = [1.3;ktdf;-0.9];
   ud0 = [-1.3;ktdr;-0.9];
   %radt = -2*pi*0.01*sin(2*pi*1*t(j))-qdt(1:3)-A1*(cross(wlokal,ua1));
   %rbdt = 2*pi*0.01*cos(2*pi*1*t(j))-qdt(1:3)-A1*(cross(wlokal,ub1));
   %rcdt = 2*pi*0.01*cos(2*pi*1*t(j))-qdt(1:3)-A1*(cross(wlokal,uc1));
   %rddt = -2*pi*0.01*sin(2*pi*1*t(j))-qdt(1:3)-A1*(cross(wlokal,ud1));
    
   %else
    %ua0 = [-1.3;0;0.9];                  %Jarak ke lokal coordinate ground
    %ub0 = [1.3;0;0.9];
    %uc0 = [1.3;0;-0.9];
    %ud0 = [-1.3;0;-0.9];
    ra0dt = [0; 0; 0]; % Tanpa Simpangan
    rb0dt = [0; 0; 0];
    rc0dt = [0; 0; 0];
    rd0dt = [0; 0; 0];
    
    ra1dt = qdt(1:3)+A1*wlokal_skew*ua1;
    rb1dt = qdt(1:3)+A1*wlokal_skew*ub1;
    rc1dt = qdt(1:3)+A1*wlokal_skew*uc1;
    rd1dt = qdt(1:3)+A1*wlokal_skew*ud1;
    
    radt = ra0dt-ra1dt;
    rbdt = rb0dt-rb1dt;
    rcdt = rc0dt-rc1dt;
    rddt = rd0dt-rd1dt;

 %end

%% Spring
ra0 = ua0;
rb0 = ub0;
rc0 = uc0;
rd0 = ud0;

ra = ra0-q(1:3)-A1*ua1;
rb = rb0-q(1:3)-A1*ub1;
rc = rc0-q(1:3)-A1*uc1;
rd = rd0-q(1:3)-A1*ud1;

La = norm(ra);
Lb = norm(rb);
Lc = norm(rc);
Ld = norm(rd);

rau = ra/La;
rbu = rb/Lb;
rcu = rc/Lc;
rdu = rd/Ld;
          
%Fa = kr*(La-Lawalr)+cr*Ladt;
%Fb = kf*(Lb-Lawalf)+cf*Lbdt;
%Fc = kf*(Lc-Lawalf)+cf*Lcdt;
%Fd = kr*(Ld-Lawalr)+cr*Lddt;

Fak = k*(La-Lawalr)*rau;
Fbk = k*(Lb-Lawalf)*rbu;
Fck = k*(Lc-Lawalf)*rcu;
Fdk = k*(Ld-Lawalr)*rdu;

F1k = Fak+Fbk+Fck+Fdk;

%%
%%Damper
Ladt = norm(radt);
Lbdt = norm(rbdt);
Lcdt = norm(rcdt);
Lddt = norm(rddt);

Fac = c*Ladt*rau;
Fbc = c*Lbdt*rbu;
Fcc = c*Lcdt*rcu;
Fdc = c*Lddt*rdu;

F1c = Fac+Fbc+Fcc+Fdc;

%%
%%Qteta

uglobala = A1*ua1;
uglobalb = A1*ub1;
uglobalc = A1*uc1;
uglobald = A1*ud1;

uags = [0 -uglobala(3) uglobala(2); 
        uglobala(3) 0 -uglobala(1); 
        -uglobala(2) uglobala(1) 0];
ubgs = [0 -uglobalb(3) uglobalb(2); 
        uglobalb(3) 0 -uglobalb(1); 
        -uglobalb(2) uglobalb(1) 0];
ucgs = [0 -uglobalc(3) uglobalc(2); 
        uglobalc(3) 0 -uglobalc(1); 
        -uglobalc(2) uglobalc(1) 0];
udgs = [0 -uglobald(3) uglobald(2); 
        uglobald(3) 0 -uglobald(1); 
        -uglobald(2) uglobald(1) 0];

Fatk = -k*(La-Lawalr)*G1'*uags'*rau;
Fbtk = -k*(Lb-Lawalr)*G1'*ubgs'*rbu;
Fctk = -k*(Lc-Lawalr)*G1'*ucgs'*rcu;
Fdtk = -k*(Ld-Lawalr)*G1'*udgs'*rdu;

Fatc = -c*Ladt*G1'*uags'*rau;
Fbtc = -c*Lbdt*G1'*ubgs'*rbu;
Fctc = -c*Lcdt*G1'*ucgs'*rcu;
Fdtc = -c*Lddt*G1'*udgs'*rdu;

%Qteta = zeros(4,1);
Qtetak = Fatk+Fbtk+Fctk+Fdtk;
Qtetac = Fatc+Fbtc+Fctc+Fdtc;
Qteta = Qtetak+Qtetac;

%%
%%Qe
Qev1  = [0;-m*g;0]+F1k+F1c;


%%Qv
Qv = -Ghat1'*(cross(wlokal,(Ihat1*wlokal))+Ihat1*Ghat1dt*[qdt(4);qdt(5);qdt(6);qdt(7)]);

%Q
Qe_R=Qev1;
Qe_t=Qteta+Qv;
Qe=[Qe_R;Qe_t];        %Matrix gaya external
Qd=-2*(qdt(4)^2 + qdt(5)^2 + qdt(6)^2 + qdt(7)^2);                                                %Matrix Qd
Q=[Qe;Qd];                                                          %Matrix Q total
qddt_lamda=MCq\Q;                                                   %Matrix percepatan lamda
