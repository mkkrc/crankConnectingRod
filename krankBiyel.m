clc
clear

F=1000;
m2=4;
is2=0.1;
m3=2;
is3=0.35;
m4=1;

l=0.8;
r=0.25;
n=l^2/r^2;
qd=120;
qdd=0;

i=1;
for q=0:0.05*pi:2*pi;
    
%analitik çözümleme s3 için
s3=asin(r*sin(q)/l);
z31=-r*cos(q);
z32=l*cos(s3);
g0t=(z31/z32);
s3d=g0t*qd; %s3 hýzý
z31t=r*sin(q);
z32t=-l*g0t*sin(q);
g0tt=(z31t*z32-z31*z32t)/z32^2;
s3dd=g0tt*qd^2+g0t*qdd; %s3 ivmesi

%analitik çözümleme s4 için
A=r*cos(q);
B=r^2-l^2;
s4=A+(A^2-B)^(1/2);
z41=r*-s4*sin(q);
z42=s4-r*cos(q);
gxt=(z41/z42);
s4d=gxt*qd; %s4 hýzý
z41t=-r*(gxt*sin(q)+s4*cos(q));
z42t=gxt+r*sin(q);
gxtt=(z41t*z42-z41*z42t)/z42^2;
s4dd=gxtt*q^2+gxt*qdd; %s4 ivmesi

xS3=(r.*cos(q)+sqrt(r.^2.*cos(q)^2-r^2+l^2)/2);
yS3=r*sin(q)/2;
xS3d=-r.*sin(q).*qd+(1/4).*(-2.*r.^2.*sin(q).*cos(q).*qd).*(r.^2.*cos(q).^2-r.^2+l.^2).^(-1/2);
yS3d=r.*cos(q).*qd./2;
x3dd=-r*(cos(q)+(8*cos(2*q)*sqrt(n^2-sin(q)^2)+4*sin(2*q)/sqrt(n^2-sin(q)^2))/(16*(n^2-sin(q)^2)))*qd^2-r*(sin(q)+sin(2*q)/(4*sqrt(n^2-sin(q)^2)))*qdd;
y3dd=-r*sin(q)*qd^2/2+r*cos(q)*qdd;

rx32=-l*cos(s3)/2;
ry32=-l*sin(s3)/2;
rx34=l*cos(s3)/2;
ry34=l*sin(s3)/2;
rx21=0;
ry21=0;
rx23=-r*cos(q)/2;
ry23=-r*sin(q)/2;

%dinamik analiz ve isletme momenti
A=[ 1 0 -1 0 0 0 0 0 0;
    0 1 0 -1 0 0 0 0 0;
    -ry32 rx32 ry34 -rx34 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    -1 0 0 0 0 0 1 0 0;
    0 -1 0 0 0 0 0 1 0;
    ry23 -rx23 0 0 0 0 -ry21 rx21 1];

%B=[Rx23; Ry23; Rx34; Ry34; R14; M14; Rx12; Ry12; Mc]
C=[m3*x3dd; m3*y3dd; m3*is3^2*s3dd; m4*s4dd+F; 0; 0; 0; 0; 0];
B=A^-1*C;

    Rx23(i)=B(1);
    Ry23(i)=B(2);
    Rx34(i)=B(3);
    Ry34(i)=B(4);
    R14(i)=B(5);
    M14(i)=B(6);
    Rx12(i)=B(7);
    Ry12(i)=B(8);
    Mc(i)=B(9);
i=i+1;

end

q=0:0.05*pi:2*pi;
plot(q,Mc)

clear g0t g0tt gxt gxtt z31 z31t z32 z32t z41 z41t z42 z42t r l F g0t g0tt...
    gxt gxtt is2 is3 l m2 m3 m4 r rx21 rx23 rx32 rx34 ry21 ry23 ry34 z31...
    z31t z32 z32t z41 z41t z42 z42t ry32 R14 Rx12 Rx23 Rx34 Ry12 Ry23 Ry34 M14