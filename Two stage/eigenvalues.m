clear all
clc
 i=sqrt(-1); 
R=[   -2.5231 +     10.984i
      -2.5231 -     10.984i
     -0.11418 +     5.1072i
     -0.11418 -     5.1072i




 ]

a=real(R);
b=imag(R);
ii=size(R);
for ff=1:ii(1)
damping(ff)=(-a(ff))/(sqrt(a(ff)^2+b(ff)^2));
freq(ff)=abs(b(ff))/(2*pi);
end
damping'
freq'

