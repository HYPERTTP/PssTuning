function [F ,X]= osimb (X)


T1=X(1);
T2=X(2);
Kpss=X(3);


% ================
Re=0;
Xe=0.4;
Vt=cos(15*pi/180)+sqrt(-1)*sin(15*pi/180);
Vinf=1.05;
H=3.25;
Tpdo=9.6;
KA=400;
TA=0.2;
Rs=0;
Xq=2.1;
Xd=2.5;
Xpd=0.39;
D=20;
ws=377;

% ========================>>> Find Initial Value <<<========================
i=sqrt(-1);
IG=(Vt-Vinf)/(Re+sqrt(-1)*Xe);
Vt0=Vt;
gama=angle(IG);
teta=angle(Vt0);
delta0=angle(Vt0+(Rs+sqrt(-1)*Xq)*IG);
Idq= abs(IG)*(cos(pi/2+gama-delta0)+sqrt(-1)*sin(pi/2+gama-delta0));
Vdq=abs(Vt0)*(cos(pi/2+teta-delta0)+sqrt(-1)*sin(pi/2+teta-delta0));
Id0=real(Idq);
Iq0=imag(Idq);
Vd0=real(Vdq);
Vq0=imag(Vdq);
Epd0=Vd0+Rs*Id0-Xpd*Iq0;
Epq0=Vq0+Rs*Iq0+Xpd*Id0;
Efd0=Epq0+(Xd-Xpd)*Id0;
Vref=abs(Vt0)+Efd0/KA;
Omega=ws;
TM=Epq0*Iq0+(Xq-Xpd)*Id0*Iq0;
% =======================================
DELTA=Re^2+(Xe+Xq)*(Xe+Xpd);
K1=(-1/DELTA)*(Iq0*Vinf*(Xpd-Xq)*((Xq+Xe)*sin(delta0)-Re*cos(delta0))+...
Vinf*((Xpd-Xq)*Id0-Epq0)*((Xpd+Xe)*cos(delta0)+Re*sin(delta0)));
K2=(1/DELTA)*(Iq0*DELTA-Iq0*(Xpd-Xq)*(Xq+Xe)-Re*(Xpd-Xq)*Id0+Re*Epq0);
K3=1/(1+((Xd-Xpd)*(Xq+Xe)/DELTA));
K4=(Vinf*(Xd-Xpd)/DELTA)*((Xq+Xe)*sin(delta0)-Re*cos(delta0));
K5=(1/DELTA)*((Vd0/abs(Vt0))*Xq*(Re*Vinf*sin(delta0)+Vinf*cos(delta0)*...
(Xpd+Xe))+(Vq0/abs(Vt0))*(Xpd*(Re*Vinf*cos(delta0)-Vinf*(Xq+Xe)*...
sin(delta0))));
K6=(1/DELTA)*((Vd0/abs(Vt0))*Xq*Re-(Vq0/abs(Vt0))*Xpd*(Xq+Xe))+...
(Vq0/abs(Vt0));

%==================================

x1=[-1/(K3*Tpdo)  -K4/Tpdo 0   1/Tpdo 0];
x2=[ 0              0      ws    0    0];
x3=[-K2/(2*H)  -K1/(2*H)   -D*ws/(2*H)  0 0];
x4=[-KA*K6/TA  -KA*K5/TA  0 -1/TA  KA/TA];
x5=[-K2*T1*Kpss/(2*H*T2)  -K1*T1*Kpss/(2*H*T2)  Kpss/T2  0 -1/T2];
Apss = [ x1 ;x2 ; x3 ; x4; x5];
eig0 = -1;
eta0 = 0.2;
Time = 20;
dVref = 0.0;
dTm= 0.20;
tsim= 1.0005;
[~,a2,a3 ] = damp(Apss);
eigi= real (a3);

dd1=length(a3);

for no = 1:dd1
    j1 (no)= (eig0-eigi(no))^2;
    j2 (no)= (eta0-a2(no))^2;
    
end

 F = 1*sum(j1)+ 1000*sum(j2);
% T=3;
% PSS= 1;
% simopt = simset('solver','ode45','SrcWorkspace','Current','DstWorkspace','Current');  % Initialize sim options
% [tout xout yout]=sim('SMIB_sim',[0 Time],simopt);  %#ok
% 
% F=1e2*trapz(tout.*(abs(teta).^2));
end