clear all
clc
format short g
re=0;
xe=0.5;
vt=(cos((15*pi)/180)+sqrt(-1)*sin((15*pi)/180));
vinf=1.05;
h=3.2;
tdop=9.6;
ka=400;
ta=0.2;
rs=0.0;
xq=2.1;
xd=2.5;
xdp=0.39;
d=0;
ws=377;
ig = (vt-vinf)/(sqrt(-1)*xe);
igangle=angle(ig);
igangledegree=rad2deg(igangle);
ed= (vt+ (sqrt(-1)*xq)*(ig));
edangle=angle(ed);
edangledegree= rad2deg(edangle);  
igd=abs(ig)*(cos(-igangle-(pi/2)+edangle)-sin(-igangle-(pi/2)+edangle)*sqrt(-1));
id=real(igd);
iq=imag(igd);
 (cos(-(pi/12)-(pi/2)+edangle)-sin(-(pi/12)-(pi/2)+edangle)*sqrt(-1));
vd=real(cos(-(pi/12)-(pi/2)+edangle)-sin(-(pi/12)-(pi/2)+edangle)*sqrt(-1));
vq=imag((cos(-(pi/12)-(pi/2)+edangle)-sin(-(pi/12)-(pi/2)+edangle)*sqrt(-1)));
eqp= vq+ xdp*id;
efd= eqp+ (xd-xdp)*id;
vref= abs(vt)+ (efd/ka);
tm= eqp*iq + (xq-xdp)*id*iq;
del= re*re +(xe+xq)*(xe+xdp);
 1 + ((xd-xdp)*(xq+xe)/del);
 k3=1/(1 + ((xd-xdp)*(xq+xe)/del));
 k4= (vinf*(xd-xdp)/del)*((xq+xe)*sin(edangle) - re*cos(edangle) );
 k2=(1/del)*(iq*del -iq*(xdp-xq)*( xq+ xe)- re*(xdp-xq)*id +re*eqp );
 k6=abs(1/del*((vd/vt)*xq*re- (vq/vt)*xdp*(xq+xe)) + vq/vt);
k5=abs((1/del)*((vd/vt)*xq*( re*vinf*sin(edangle) + vinf*cos(edangle)*(xdp+xe)) +(vq/vt)*(xdp*(re*vinf*cos(edangle)-vinf*(xe+xq)*sin(edangle)))));
k1= (-1/del)*(iq*vinf*(xdp-xq)*((xq+xe)*sin(edangle)-re*cos(edangle)) + vinf*((xdp-xq)*id-eqp)*((xdp+xe)*cos(edangle) + re*sin(edangle) ));

%%without exciter 
ka=400;
kl = [-1/(k3*tdop) -k4/(tdop) 0 1/tdop; 0 0 ws 0 ; -k2/(2*h) -k1/(2*h) (-d*ws)/2*h 0; (-ka*k6)/ta  -(ka*k5)/ta 0 -1/ta ];
eval=eig(kl);

%%with exciter
kpss=0.5;
t1=0.5;
t2=0.1;
ml = [-1/(k3*tdop) -k4/(tdop) 0 1/tdop 0;
       0 0 ws 0 0 ; 
       -k2/(2*h) -k1/(2*h) 0 0 0; 
       (-ka*k6)/ta  -(ka*k5)/ta  0  -1/ta ka/ta; 
       (-k2*t1*kpss)/(t2*2*h) (-k1*t1*kpss)/(t2*2*h) kpss/t2 0 -1/t2];
eval2=eig(ml)
kpss=0.5;
tw=10;
t3=1;
t4=1;
t5=1;
t6=1;
%% without second washout block
gl = [-1/(k3*tdop) -k4/(tdop) 0 1/tdop 0 0;
       0 0 ws 0 0 0; 
       -k2/(2*h) -k1/(2*h) 0 0 0 0; 
       (-ka*k6)/ta  -(ka*k5)/ta  0  -1/ta ka/ta 0; 
       (-k2*t1*kpss)/(t2*2*h) (-k1*t1*kpss)/(t2*2*h) kpss/t2 0 -1/t2 0;
        (-kpss*t1*k2)/(2*h*t2)  (-kpss*k1*t1)/(2*h*t2) (-ws*kpss*t1)/2*h*t2 0 ((tw+t1)/(t2*tw)) -1/t2   ];   
    eig(gl);
    d=1;
    t2=1;
    
%%with second washout block
    pss = [      -1/(k3*tdop) -k4/(tdop) 0 1/tdop 0 0 0;
                    0 0 ws 0 0 0 0; 
                  -k2/(2*h) -k1/(2*h) 0 0 0 0 0; 
                  (-ka*k6)/ta  -(ka*k5)/ta  0  -1/ta ka/ta 0 0; 
                  
                    (-k2*kpss)/(2*h) (-k1*kpss)/(2*h) -(kpss*d*ws)/2*h 0 -1/tw 0 0;
            
                    (-kpss*t1*k2)/(2*h*t2)  (-kpss*k1*t1)/(2*h*t2) (-ws*d*kpss*t1)/(2*h*t2) 0 (1/t2)-t1/(tw*t2)   -1/t2  0 ;
         -(t1*t3*k2*kpss)/(2*h*t2*t4)  -(t1*t3*k1*kpss)/(2*h*t2*t4) (-d*t3*kpss*t1*ws)/(2*h*t2*t4)  0  t3/(t2*t4)-(t3*t1)/(t2*t4*tw)  (-1/(t4))+(t3)/t4*t2  -1/t4 ]  
    eig(pss)
    
%% new pss filter%%
    
kp1=1;
kp2 =1;
kp3=1;
t3=1;
t4=1;
t5=1;
t6=1;
fl = [-1/(k3*tdop) -k4/(tdop) 0 1/tdop 0 0 0 0 0 0 0;
       0 0 ws 0 0 0 0 0 0 0 0; 
       -k2/(2*h) -k1/(2*h) 0 0 0 0 0 0 0 0 0; 
       (-ka*k6)/ta  -(ka*k5)/ta  0  -1/ta ka/ta 0 0 0 0 0 0;
        0 0 0 0   (-1/t1)     0     0   0     0   0    0            ;
        0 0 0 0    -kp1/t1  -1/t2   0          0       0     0     0        ;
        0 0 0 0    0         0    (-1/t3)      0       0     0      0      ;
        0 0 0 0   0          0     -kp2/t3    -1/t4     0     0     0      ;
        0 0 0 0    0         0       0          0     -1/t5   0     0      ;
        0 0 0 0    0         0       0          0    -kp3/t5 -1/t6  0      ;
        0 0 0 0 -(1+kp1)/t1 -1/t2  -(1+kp2)/t3 -1/t4   -(1+kp3)/t5   -1/t6 0 ];
      eig(fl);
 
