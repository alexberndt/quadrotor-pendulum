% PD controlled quadcopter
% coordinates q=[x y z phi theta psi]

clear
% close all;
clc

% quadcopter parameters
m=1/2;
g=9.81;

% setpoints ( phi_r=0, theta_r=0 )
x_r=2;
y_r=4;
z_r=5;
psi_r=pi/4;

% sim parameters
tspan=[0:0.05:10];
Xr=[x_r;y_r;z_r;0;0;psi_r;zeros(6,1)];

X0=Xr+0.2*randn(12,1);

% visualization parameters
w=1.5;
Ax=[-w+x_r w+x_r -w+y_r w+y_r -w+z_r w+z_r];

% control parameters
kz=3;
dz=2;

kpsi=3;
dpsi=3;

k1=3;
d1=4;

k2=25;
d2=10;

% saturation
sl=1.5;
sat=@(x)[min(max(x,-sl),sl)];

% vectorfield
f=@(t,x)[x(7:12);
         (cos(x(4))*cos(x(6))*sin(x(5))+sin(x(4))*sin(x(6)))*(g-kz*(x(3)-z_r)-dz*x(9))/(cos(x(4))*cos(x(5)));
         (cos(x(4))*sin(x(6))*sin(x(5))-sin(x(4))*cos(x(6)))*(g-kz*(x(3)-z_r)-dz*x(9))/(cos(x(4))*cos(x(5)));
         cos(x(4))*cos(x(5))*(g-kz*(x(3)-z_r)-dz*x(9))/(cos(x(4))*cos(x(5)))-g;
         (sin(x(6))*(-k1*(x(1)-x_r)-d1*x(7))-cos(x(6))*(-k1*(x(2)-y_r)-d1*x(8)))-k2*x(4)-d2*x(10);
         (cos(x(6))*(-k1*(x(1)-x_r)-d1*x(7))+sin(x(6))*(-k1*(x(2)-y_r)-d1*x(8)))-k2*x(5)-d2*x(11);
         -kpsi*(x(6)-psi_r)-dpsi*x(12)];
  
% ODE simulation     
[T,X]=ode45(@(t,x)f(t,x),tspan,X0);
figure(1);
clf;
subplot 211
plot(T,X(:,1:6));
h=legend('$x$','$y$','$z$','$\phi$','$\theta$','$\psi$');set(h,'interpreter','latex')
subplot 212
plot(T,X(:,7:12));
h=legend('$\dot x$','$\dot y$','$\dot z$','$\dot \phi$','$\dot \theta$','$\dot\psi$');set(h,'interpreter','latex')
xlabel('$t$','interpreter','latex')

%% visualisation
l=0.3;
rc=0.1;
rx=rc*cos(linspace(0,2*pi,1e1));
rx=[rx rx(1)];
ry=rc*sin(linspace(0,2*pi,1e1));
ry=[ry ry(1)];


%%
figure(2)
clf;
grid();

X = X(:,1:6);

[N,~] = size(X);

for j=1:N
   Rt=R(X(j,4:6)); 
    R1=Rt*([l+rc+rx;ry;zeros(size(rx))])+X(j,1:3)';
    R2=Rt*([rx;l+rc+ry;zeros(size(rx))])+X(j,1:3)';
    R3=Rt*([-l-rc+rx;ry;zeros(size(rx))])+X(j,1:3)';
    R4=Rt*([rx;-l-rc+ry;zeros(size(rx))])+X(j,1:3)';

    A1=Rt*([-l l;0 0;0 0])+X(j,1:3)';
    A2=Rt*([0 0; -l l;0 0])+X(j,1:3)';


plot3(x_r,y_r,z_r,'r*',A1(1,:),A1(2,:),A1(3,:),'k',A2(1,:),A2(2,:),A2(3,:),'k',R1(1,:),R1(2,:),R1(3,:),'r',R2(1,:),R2(2,:),R2(3,:),'b',R3(1,:),R3(2,:),R3(3,:),'b',R4(1,:),R4(2,:),R4(3,:),'b');

axis(Ax);
set(gca,'box','on')
drawnow
end


function y=R(Xrot)
phi=Xrot(1);
theta=Xrot(2);
psi=Xrot(3);

Rpsi=[cos(psi)  -sin(psi)   0;
      sin(psi)  cos(psi)    0;
      0         0           1];
  
% rotation around y with theta
Rtheta=[cos(theta)    0       sin(theta);
      0             1       0
      -sin(theta)    0       cos(theta)];
  
% rotation around x with phi 
Rphi=[1       0           0;
        0       cos(phi)    -sin(phi);
        0       sin(phi)    cos(phi)];

y=Rpsi*Rtheta*Rphi;
end
