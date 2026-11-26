function UAV001()
%无人机参数输入
m=1.0;Ix=0.0025;Iy=0.0025;Iz=0.0045;d=0.17;g=9.8;
cR=1250;omegab=0;Tm=0.05;cT=5*1e-6;cM=2.5*1e-7;
x0=zeros(13,1);%13个状态量
%x,y,z,q0,q1,q2,q3,x_d,y_d,z_d,p,q,r
x0(4)=1;%q0=1;
omega0=zeros(4,1);
X0=[x0;omega0];
scenario=2;

switch scenario
    case 1
    fprintf('场景1：先上升再悬停\n');
    Sigma=@hover;
    
    case 2
    fprintf('场景2：恒定滚转力矩\n');
    Sigma=@(t)[0.5604755;0.5604755;0.5600000;0.5600000];

    case 3
    fprintf('场景3：恒定俯仰力矩\n');
    Sigma=@(t)[0.56023; 0.55977; 0.55977; 0.56023];

    case 4
    fprintf('场景4：恒定偏航力矩\n');
    Sigma=@(t)[0.56057; 0.55943; 0.56057; 0.55943];


end

Tspan=0:0.01:5;
opts=odeset('RelTol',1e-8,'AbsTol',1e-10);%给ode求解器的误差参数

[t,X]=ode45(@UAV_ode,Tspan,X0,opts);



% 图1：三维轨迹
figure('Name','三维轨迹','Color','w');
x=X(:,1:13);
plot3(x(:,1),x(:,2),x(:,3),'LineWidth',1.5);
grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('三维轨迹');
axis equal; % 让x、y、z三个坐标轴的单位长度相同
set(gca,'ZDir','reverse');
%图2：位置
figure('Name','位置','Color','w');
plot(t, x(:,1:3),'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('位置 [m]');
legend('x','y','z','Location','best');
title('位置');

%图3：角度
% 先计算角度
q0=x(:,4);q1=x(:,5);q2=x(:,6);q3=x(:,7);
phi=atan2(2*(q0.*q1+q2.*q3),1-2*(q1.^2+q2.^2));
theta=asin(2*(q0.*q2-q3.*q1));
psi=atan2(2*(q0.*q3+q1.*q2),1-2*(q2.^2+q3.^2));
phi=unwrap(phi);
theta=unwrap(theta);
psi=unwrap(psi);

figure('Name','角度','Color','w');
plot(t, [phi,theta,psi],'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('角度 [rad]');
legend('phi','theta','psi','Location','best');
title('角度');

%图4：速度
figure('Name','速度','Color','w');
plot(t, x(:,8:10),'LineWidth',1.2); grid on;
xlabel('t [s]'); ylabel('速度 [m/s]');
legend('x\_d','y\_d','z\_d','Location','best'); 
title('速度');
fprintf('       t           x         y         z        phi      theta      psi\n');


disp([t,x(:,1:3),phi*180/pi,theta*180/pi,psi*180/pi])


    function dxdt=UAV_ode(t,X)

q0=X(4);q1=X(5);q2=X(6);q3=X(7);
x_d=X(8);y_d=X(9);z_d=X(10);
p=X(11);q=X(12);r=X(13);
omega1=X(14);omega2=X(15);omega3=X(16);omega4=X(17);
omega=[omega1;omega2;omega3;omega4];

    B=[cT cT cT cT 
        -sqrt(2)*d*cT/2 -sqrt(2)*d*cT/2 sqrt(2)*d*cT/2 sqrt(2)*d*cT/2
        sqrt(2)*d*cT/2 -sqrt(2)*d*cT/2 -sqrt(2)*d*cT/2 sqrt(2)*d*cT/2
        cM -cM cM -cM];
    omegass=cR*Sigma(t)+omegab;
    squareomega=omega.^2;
    U=B*squareomega;
    omega_d=(-omega+omegass)/Tm;%一阶惯性环节
omega1_d=omega_d(1);omega2_d=omega_d(2);omega3_d=omega_d(3);omega4_d=omega_d(4);
Q=[q0;q1;q2;q3];
mode=sqrt(Q'*Q);
q0=q0/mode;
q1=q1/mode;
q2=q2/mode;
q3=q3/mode;

U1=U(1);U2=U(2);U3=U(3);U4=U(4);
Q=[q0;q1;q2;q3];
Q_d=[0 -p -q -r
    p 0 r -q
    q -r 0 p
    r q -p 0]*Q/2;

q0_d=Q_d(1);
q1_d=Q_d(2);
q2_d=Q_d(3);
q3_d=Q_d(4);


x_dd=-U1*2*(q1*q3+q0*q2)/m;

y_dd=-U1*2*(q2*q3-q0*q1)/m;


p_d=(U2+(Iy-Iz)*q*r)/Ix;
q_d=(U3+(Iz-Ix)*p*r)/Iy;
r_d=(U4+(Ix-Iy)*p*q)/Iz;

if X(3) >= 0 && (U1*(1-2*q1^2-2*q2^2) <= m*g)

z_dd=0;
else
z_dd=g-U1*(1-2*q1^2-2*q2^2)/m;
end
    dxdt=[x_d;y_d;z_d;
    q0_d;q1_d;q2_d;q3_d;
    x_dd;y_dd;z_dd;
    p_d;q_d;r_d;omega1_d;omega2_d;omega3_d;omega4_d];
end

    function sigmahover=hover(t)
        if t<2.85
            s=0.58;
        elseif t<3.0
            s=0.58-2/15*(t-2.85);
        else
            s=0.56;
        end
        sigmahover=[s;s;s;s];
    end
end
