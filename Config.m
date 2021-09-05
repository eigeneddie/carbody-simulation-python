%%Transformation Matrix
%Car body

G1 = 2*[-q(5) q(4) -q(7) q(6); -q(6) q(7) q(4) -q(5); -q(7) -q(6) q(5) q(4)];
Ghat1 = 2*[-q(5) q(4) q(7) -q(6); -q(6) -q(7) q(4) q(5); -q(7) q(6) -q(5) q(4)];
Ghat1dt = 2*[-qdt(5) qdt(4) qdt(7) -qdt(6); -qdt(6) -qdt(7) qdt(4) qdt(5); -qdt(7) qdt(6) -qdt(5) qdt(4)];


A1 = [1-2*q(6)^2-2*q(7)^2 2*(q(5)*q(6)-q(4)*q(7)) 2*(q(5)*q(7)+q(4)*q(6));
      2*(q(5)*q(6)+q(4)*q(7)) 1-2*q(5)^2-2*q(7)^2 2*(q(6)*q(7)-q(4)*q(5));
      2*(q(5)*q(7)-q(4)*q(6)) 2*(q(6)*q(7)+q(4)*q(5)) 1-2*q(5)^2-2*q(6)^2;];


wlokal = Ghat1*[qdt(4);qdt(5);qdt(6);qdt(7)];
wlokal_skew = [0 -wlokal(3) wlokal(2); wlokal(3) 0 -wlokal(1); -wlokal(2) wlokal(1) 0]; 

%%Mass Matrix
mrr = [m zeros(1,2);0 m 0; zeros(1,2) m];

mtt = Ghat1'*Ihat1*Ghat1;

M = [mrr zeros(3,4); zeros(4,3) mtt];


%%Constraint Equation

C1     = q(4)^2+q(5)^2+q(6)^2+q(7)^2-1; 
 
C      = C1;                          %Matrix Constraint
                      

%Constraint Jacobian Matrix
Cq_eu1  =[zeros(1,3) 2*q(4) 2*q(5) 2*q(6) 2*q(7)];

Cq=Cq_eu1;                          %Constraint Jacobian Matrix

Cq_dep=Cq(:,4);                                 %Constraint Jacobian Matrix (Dependent Variable)
Cq_ind=Cq(:,[1:3 5:6 7]);                         %Constraint Jacobian Matrix (Independent Variable)
       
%if Cq_ind(1,6)<1e-12
    %Cq_ind(1,6) = 0;
%else
    %Cq_ind(1,6) = Cq_ind(1,6);
%end
       