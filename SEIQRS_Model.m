% SEIQRS Model analysis and optimal control with time delays diffusion effect
% % % case 1:Disease-free equilibrium point
% clc
% clear
% hold on
% b=33;
% p=0.025;
% d=0.01;
% theta=0.015;
% beta=0.0002;
% sigma=0.0014;
% omega=0.0155;
% xi=0.01;
% eplison=0.014;
% gamma=0.03;
% eta=0.01;
% delta=0.003;
% a1=0.01;
% a2=0.018;
% t1=36.9;
% t2=36.9;
% t3=35;
% 
% Z1=theta+d+omega+xi;
% Z2=d+a1+eplison+gamma;
% Z3=d+a2+eta;
% Z4=d+delta;
% ddex1dez=@(t,y,z)[(1-p)*b+theta*y(2)-d*y(1)+delta*z(5,3)-beta*z(1,1)*z(3,2)/(1+sigma*z(1,1));
%     beta*z(1,1)*z(3,2)/(1+sigma*z(1,1))-Z1*y(2);
%     omega*y(2)-Z2*y(3);
%     eplison*y(3)-Z3*y(4);
%     p*b+xi*y(2)+gamma*y(3)+eta*y(4)-d*y(5)-delta*z(5,3)];
% lags=[t1 t2 t3];
% S0=b*(d-d*p+delta)/(d*Z4);%S0 in disease-free equilibrium
% R0=(p*b)/Z4;%R0 in disease-free equilibrium
% R00=(beta*omega*S0)/(Z1*Z2*(1+sigma*S0));%Basic regeneration number
% %Positive equilibrium point values
% S1=(Z1*Z2)/(beta*omega-Z1*Z2*sigma);
% E1=(Z2*Z3*(Z1*Z2*d^2+Z1*Z2*d*delta-b*beta*d*omega-b*beta*delta*omega+Z1*Z2*b*d*sigma+Z1*Z2*b*delta*sigma+b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% I1=(Z3*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% Q1=(eplison*omega*E1)/(Z2*Z3);
% R1=(Z1*Z2^2*Z3*d*xi - Z3*b*beta*gamma*omega^2 - b*beta*eplison*eta*omega^2 + Z1*Z2^2*Z3*b*sigma*xi + Z3*b*beta*gamma*omega^2*p + b*beta*eplison*eta*omega^2*p + Z1^2*Z2^2*Z3*b*p*sigma + Z1*Z2*Z3*d*gamma*omega + Z1*Z2*d*eplison*eta*omega - Z2*Z3*b*beta*omega*xi - Z1*Z2*Z3*b*beta*omega*p + Z1*Z2*Z3*b*gamma*omega*sigma + Z1*Z2*b*eplison*eta*omega*sigma + Z2*Z3*b*beta*omega*p*theta + Z2*Z3*b*beta*omega*p*xi - Z1*Z2^2*Z3*b*p*sigma*theta - Z1*Z2^2*Z3*b*p*sigma*xi - Z1*Z2*Z3*b*gamma*omega*p*sigma - Z1*Z2*b*eplison*eta*omega*p*sigma)/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% 
% sol=dde23(ddex1dez,lags,[3000 1250 600 600 1000],[0,2000]);
% plot(sol.x,sol.y,'LineWidth',1.6);
% set(gca,'FontSize',16)
% xlabel('t');
% ylabel('Number of hosts/Unit');
% legend('S','E','I','Q','R')
% tint=linspace(0,5,10);
% sint=deval(sol,tint);
% hold on
% grid on
% % % for i=1:50
% % %    while 1
% % %         ii=60+8*rand()*rand*i;
% % %         ie=40+10*rand*i;
% % %         iq=30+10*rand()*i;
% % %         ir=60+60*rand*i;
% % %         is=2000-ii-ir-ie-iq;
% % %         if(is > 0)
% % %             break;
% % %         end
% % %    end
% % %     sol=dde23(ddex1dez,lags,[is ie ii iq  ir],[0,600]);
% % % %     plot3(sol.y(3,:),sol.y(1,:),sol.y(2,:),'-.');
% % %    plot3(sol.y(1,:),sol.y(5,:),sol.y(4,:),'-.');
% % % end
% % % xlabel('\itI(\itt)');
% % % ylabel('\itR(\itt)');
% % % zlabel('\itQ(\itt)');
% % % set(gca,'FontSize',10,'LineWidth',1);
% % % grid on;
% % % hold off;
% % % % % Hopf Bifurcation
% % plot3(sol.y(5,:),sol.y(4,:),sol.y(1,:),'-.');
% % 
% % xlabel('R(t)');
% % ylabel('Q(t)');
% % zlabel('S(t)');
% % set(gca,'FontSize',10,'LineWidth',1);
% % grid off;
% % hold off;
% % % % 
% % plot(sol.y(5,:),sol.y(1,:),'-.');
% % xlabel('R(t)');
% % ylabel('S(t)');
% % set(gca,'FontSize',10,'LineWidth',1);
% % grid off;
% % hold off;



% % % %The SEIQRS model in different situations
% % %Case 2:t1>0,t2=0,The t1 is 82.84184
% % clc
% % clear
% % b=33;
% % p=0.025;
% % d=0.01;
% % theta=0.0015;
% % beta=0.0005;
% % sigma=0.0014;
% % omega=0.0155;
% % xi=0.01;
% % eplison=0.014;
% % gamma=0.03;
% % eta=0.01;
% % delta=0.03;
% % a1=0.01;
% % a2=0.018;
% % % % t1=88.24;
% % % % t2=88.24;
% % t1=90;t2=90;
% % % % t1=70;t2=70;
% % Z1=theta+d+omega+xi;
% % Z2=d+a1+eplison+gamma;
% % Z3=d+a2+eta;
% % Z4=d+delta;
% % ddex1dez=@(t,y,z)[(1-p)*b+theta*y(2)-d*y(1)+delta*y(5)-beta*z(1,1)*z(3,2)/(1+sigma*z(1,1));
% %     beta*z(1,1)*z(3,2)/(1+sigma*z(1,1))-Z1*y(2);
% %     omega*y(2)-Z2*y(3);
% %     eplison*y(3)-Z3*y(4);
% %     p*b+xi*y(2)+gamma*y(3)+eta*y(4)-d*y(5)-delta*y(5)];
% % lags=[t1 t2];
% % S0=(b*d + b*delta - b*d*p)/(d^2 + delta*d);%S0 in disease-free equilibrium
% % R0=(b*p)/(d+delta);%R0 in disease-free equilibrium
% % R00=(beta*omega*S0)/(Z1*Z2*(1+sigma*S0));%Basic regeneration number
% % %  Positive equilibrium point values
% % S1=(Z1*Z2)/(beta*omega - Z1*Z2*sigma);
% % E1=(Z2*Z3*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % I1=(Z3*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % Q1=(eplison*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % R1=(Z1*Z2^2*Z3*d*xi - Z3*b*beta*gamma*omega^2 - b*beta*eplison*eta*omega^2 + Z1*Z2^2*Z3*b*sigma*xi + Z3*b*beta*gamma*omega^2*p + b*beta*eplison*eta*omega^2*p + Z1^2*Z2^2*Z3*b*p*sigma + Z1*Z2*Z3*d*gamma*omega + Z1*Z2*d*eplison*eta*omega - Z2*Z3*b*beta*omega*xi - Z1*Z2*Z3*b*beta*omega*p + Z1*Z2*Z3*b*gamma*omega*sigma + Z1*Z2*b*eplison*eta*omega*sigma + Z2*Z3*b*beta*omega*p*theta + Z2*Z3*b*beta*omega*p*xi - Z1*Z2^2*Z3*b*p*sigma*theta - Z1*Z2^2*Z3*b*p*sigma*xi - Z1*Z2*Z3*b*gamma*omega*p*sigma - Z1*Z2*b*eplison*eta*omega*p*sigma)/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % 
% % sol=dde23(ddex1dez,lags,[500 1250 300 150 600],[0,6000]);
% % plot(sol.x,sol.y,'LineWidth',1.6);
% % set(gca,'FontSize',16)
% % xlabel('t');
% % ylabel('Number of hosts/Unit');
% % legend('S','E','I','Q','R')
% % tint=linspace(0,5,10);
% % sint=deval(sol,tint);
% % hold on
% % grid on
% % % SER phase diagram
% % plot3(sol.y(5,:),sol.y(4,:),sol.y(1,:),'-.');
% % xlabel('R(t)');
% % ylabel('Q(t)');
% % zlabel('S(t)');
% % set(gca,'FontSize',10,'LineWidth',1);
% % grid off;
% % hold off;
% % % % SER phase diagram
% % plot3(sol.y(5,:),sol.y(2,:),sol.y(1,:),'-.');
% % hold on
% % xlabel('R(t)');
% % ylabel('E(t)');
% % zlabel('S(t)');
% % set(gca,'FontSize',10,'LineWidth',1);
% % grid off;
% % hold off;



% % %case3:t1=0,t2>0, t1 is 70
% clc
% clear
% b=33;
% p=0.025;
% d=0.01;
% theta=0.0015;
% beta=0.0005;
% sigma=0.0014;
% omega=0.0155;
% xi=0.01;
% eplison=0.014;
% gamma=0.03;
% eta=0.01;
% delta=0.03;
% a1=0.01;
% a2=0.018;
% % t=65;
% t=73;
% Z1=theta+d+omega+xi;
% Z2=d+a1+eplison+gamma;
% Z3=d+a2+eta;
% Z4=d+delta;
% ddex1dez=@(t,y,z)[(1-p)*b+theta*y(2)-d*y(1)+delta*z(5,1)-beta*y(1)*y(3)/(1+sigma*y(1));
%     beta*y(1)*y(3)/(1+sigma*y(1))-Z1*y(2);
%     omega*y(2)-Z2*y(3);
%     eplison*y(3)-Z3*y(4);
%     p*b+xi*y(2)+gamma*y(3)+eta*y(4)-d*y(5)-delta*z(5,1)];
% 
% lags=[t t];
% S0=(b*d + b*delta - b*d*p)/(d^2 + delta*d);%S0 in disease-free equilibrium
% R0=(b*p)/(d+delta);%R0 in disease-free equilibrium
% R00=(beta*omega*S0)/(Z1*Z2*(1+sigma*S0));%Basic regeneration number
% % % Positive equilibrium point values
% S1=(Z1*Z2)/(beta*omega - Z1*Z2*sigma);
% E1=(Z2*Z3*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% I1=(Z3*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% Q1=(eplison*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% R1=(Z1*Z2^2*Z3*d*xi - Z3*b*beta*gamma*omega^2 - b*beta*eplison*eta*omega^2 + Z1*Z2^2*Z3*b*sigma*xi + Z3*b*beta*gamma*omega^2*p + b*beta*eplison*eta*omega^2*p + Z1^2*Z2^2*Z3*b*p*sigma + Z1*Z2*Z3*d*gamma*omega + Z1*Z2*d*eplison*eta*omega - Z2*Z3*b*beta*omega*xi - Z1*Z2*Z3*b*beta*omega*p + Z1*Z2*Z3*b*gamma*omega*sigma + Z1*Z2*b*eplison*eta*omega*sigma + Z2*Z3*b*beta*omega*p*theta + Z2*Z3*b*beta*omega*p*xi - Z1*Z2^2*Z3*b*p*sigma*theta - Z1*Z2^2*Z3*b*p*sigma*xi - Z1*Z2*Z3*b*gamma*omega*p*sigma - Z1*Z2*b*eplison*eta*omega*p*sigma)/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% 
% 
% sol=dde23(ddex1dez,lags,[500 1250 300 150 600],[0,6000]);
% plot(sol.x,sol.y,'LineWidth',1.6);
% set(gca,'FontSize',16)
% xlabel('t');
% ylabel('Number of hosts/Unit');
% legend('S','E','I','Q','R')
% tint=linspace(0,5,10);
% sint=deval(sol,tint);
% hold on
% grid on
% % % % % 
% % % %Hopf bifurcation(SER phase diagram)
% plot3(sol.y(5,:),sol.y(2,:),sol.y(1,:),'-.');
% hold on
% xlabel('R(t)');
% ylabel('E(t)');
% zlabel('S(t)');
% set(gca,'FontSize',10,'LineWidth',1);
% grid off;
% hold off;



% % % % %case4:t1=t2=t3
% clc
% clear
% b=33;
% p=0.025;
% d=0.01;
% theta=0.0015;
% beta=0.0005;
% sigma=0.0014;
% omega=0.0155;
% xi=0.01;
% eplison=0.014;
% gamma=0.03;
% eta=0.01;
% delta=0.03;
% a1=0.01;
% a2=0.018;
% % t1=36;
% % t2=36;
% % t3=36;
% t1=30;
% t2=30;
% t3=30;
% Z1=theta+d+omega+xi;
% Z2=d+a1+eplison+gamma;
% Z3=d+a2+eta;
% Z4=d+delta;
% ddex1dez=@(t,y,z)[(1-p)*b+theta*y(2)-d*y(1)+delta*z(5,3)-beta*z(1,1)*z(3,2)/(1+sigma*z(1,1));
%     beta*z(1,1)*z(3,2)/(1+sigma*z(1,1))-Z1*y(2);
%     omega*y(2)-Z2*y(3);
%     eplison*y(3)-Z3*y(4);
%     p*b+xi*y(2)+gamma*y(3)+eta*y(4)-d*y(5)-delta*z(5,3)];
% 
% lags=[t1 t2 t3];
% S0=(b*d + b*delta - b*d*p)/(d^2 + delta*d);%S0 in disease-free equilibrium
% R0=(b*p)/(d+delta);%R0 in disease-free equilibrium
% R00=(beta*omega*S0)/(Z1*Z2*(1+sigma*S0));%Basic regeneration number
% % Positive equilibrium point values
% S1=(Z1*Z2)/(beta*omega - Z1*Z2*sigma);
% E1=(Z2*Z3*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% I1=(Z3*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% Q1=(eplison*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% R1=(Z1*Z2^2*Z3*d*xi - Z3*b*beta*gamma*omega^2 - b*beta*eplison*eta*omega^2 + Z1*Z2^2*Z3*b*sigma*xi + Z3*b*beta*gamma*omega^2*p + b*beta*eplison*eta*omega^2*p + Z1^2*Z2^2*Z3*b*p*sigma + Z1*Z2*Z3*d*gamma*omega + Z1*Z2*d*eplison*eta*omega - Z2*Z3*b*beta*omega*xi - Z1*Z2*Z3*b*beta*omega*p + Z1*Z2*Z3*b*gamma*omega*sigma + Z1*Z2*b*eplison*eta*omega*sigma + Z2*Z3*b*beta*omega*p*theta + Z2*Z3*b*beta*omega*p*xi - Z1*Z2^2*Z3*b*p*sigma*theta - Z1*Z2^2*Z3*b*p*sigma*xi - Z1*Z2*Z3*b*gamma*omega*p*sigma - Z1*Z2*b*eplison*eta*omega*p*sigma)/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% sol=dde23(ddex1dez,lags,[500 1250 300 150 600],[0,6000]);
% plot(sol.x,sol.y,'LineWidth',1.6);
% set(gca,'FontSize',16)
% xlabel('t');
% ylabel('Number of hosts/Unit');
% legend('S','E','I','Q','R')
% tint=linspace(0,5,10);
% sint=deval(sol,tint);
% hold on
% grid on

% % % % % %Hopf bifurcation
% plot3(sol.y(5,:),sol.y(3,:),sol.y(1,:),'-.');
% hold on
% xlabel('R(t)');
% ylabel('I(t)');
% zlabel('S(t)');
% set(gca,'FontSize',10,'LineWidth',1);
% grid off;
% hold off
% xlabel('\itS(\itt)');
% ylabel('\itR(\itt)');
% zlabel('\itI(\itt)');
% set(gca,'FontSize',10,'LineWidth',1);
% grid off;
% hold off;



% % %case5 and case6:t1>0,t2>0
% % clc
% % clear
% % b=33;
% % p=0.025;
% % d=0.01;
% % theta=0.0015;
% % beta=0.0005;
% % sigma=0.0014;
% % omega=0.0155;
% % xi=0.01;
% % eplison=0.014;
% % gamma=0.03;
% % eta=0.01;
% % delta=0.03;
% % a1=0.01;
% % a2=0.018;
% % t1=36.9;
% % t2=36.9;
% % t3=35;
% % Z1=theta+d+omega+xi;
% % Z2=d+a1+eplison+gamma;
% % Z3=d+a2+eta;
% % Z4=d+delta;
% % ddex1dez=@(t,y,z)[(1-p)*b+theta*y(2)-d*y(1)+delta*z(5,3)-beta*z(1,1)*z(3,2)/(1+sigma*z(1,1));
% %     beta*z(1,1)*z(3,2)/(1+sigma*z(1,1))-Z1*y(2);
% %     omega*y(2)-Z2*y(3);
% %     eplison*y(3)-Z3*y(4);
% %     p*b+xi*y(2)+gamma*y(3)+eta*y(4)-d*y(5)-delta*z(5,3)];
% % 
% % lags=[t1 t2 t3];
% % S0=(b*d + b*delta - b*d*p)/(d^2 + delta*d);%S0 in disease-free equilibrium
% % R0=(b*p)/(d+delta);%R0 in disease-free equilibrium
% % R00=(beta*omega*S0)/(Z1*Z2*(1+sigma*S0));%Basic regeneration number
% % % Positive equilibrium point values
% % S1=(Z1*Z2)/(beta*omega - Z1*Z2*sigma);
% % E1=(Z2*Z3*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % I1=(Z3*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % Q1=(eplison*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % R1=(Z1*Z2^2*Z3*d*xi - Z3*b*beta*gamma*omega^2 - b*beta*eplison*eta*omega^2 + Z1*Z2^2*Z3*b*sigma*xi + Z3*b*beta*gamma*omega^2*p + b*beta*eplison*eta*omega^2*p + Z1^2*Z2^2*Z3*b*p*sigma + Z1*Z2*Z3*d*gamma*omega + Z1*Z2*d*eplison*eta*omega - Z2*Z3*b*beta*omega*xi - Z1*Z2*Z3*b*beta*omega*p + Z1*Z2*Z3*b*gamma*omega*sigma + Z1*Z2*b*eplison*eta*omega*sigma + Z2*Z3*b*beta*omega*p*theta + Z2*Z3*b*beta*omega*p*xi - Z1*Z2^2*Z3*b*p*sigma*theta - Z1*Z2^2*Z3*b*p*sigma*xi - Z1*Z2*Z3*b*gamma*omega*p*sigma - Z1*Z2*b*eplison*eta*omega*p*sigma)/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
% % 
% % sol=dde23(ddex1dez,lags,[500 1250 300 150 600],[0,6000]);
% % plot(sol.x,sol.y,'LineWidth',1.6);
% % set(gca,'FontSize',16)
% % xlabel('t');
% % ylabel('Number of hosts/unit');
% % legend('S','E','I','Q','R')
% % tint=linspace(0,5,10);
% % sint=deval(sol,tint);
% % hold on
% % grid on
% % 
% % % % % %Hopf bifurcation
% % % plot3(sol.y(5,:),sol.y(2,:),sol.y(1,:),'-.');
% % % hold on
% % % xlabel('R(t)');
% % % ylabel('E(t)');
% % % zlabel('S(t)');
% % % set(gca,'FontSize',10,'LineWidth',1);
% % % grid off;
% % % hold off



% %Comparative tests
clc
clear
hold on
b=3;
p=0.025;
d=0.01;
theta=0.015;
beta=0.00022;
sigma=0.0014;
omega=0.0155;
xi=0.01;
eplison=0.014;
gamma=0.03;
eta=0.01;
delta=0.003;
a1=0.01;
a2=0.018;
t1=5;
t2=5;
t3=15;
Z1=theta+d+omega+xi;
Z2=d+a1+eplison+gamma;
Z3=d+a2+eta;
Z4=d+delta;
% %SIRS
ddex1dez=@(t,y,z)[(1-p)*b-d*y(1)+delta*y(3)-beta*z(1,1)*z(2,2)/(1+sigma*z(1,1));
    beta*z(1,1)*z(2,2)/(1+sigma*z(1,1))-(d+a1+gamma)*y(2);
    p*b+gamma*y(2)-(d+delta)*y(3)];
lags=[t1 t2]
sol=dde23(ddex1dez,lags,[700 100 700 ],[0,1200]);
plot(sol.x,sol.y(2,:),'LineWidth',1.6);
set(gca,'FontSize',16)
xlabel('t');
ylabel('Number of hosts/Unit');
% legend('I1-SIRS')
tint=linspace(0,5,10);
sint=deval(sol,tint);
hold on
grid on

%SIQRS
ddex1dez=@(t,y,z)[(1-p)*b-d*y(1)+delta*y(4)-beta*z(1,1)*z(2,2)/(1+sigma*z(1,1));
    beta*z(1,1)*z(2,2)/(1+sigma*z(1,1))-(d+a1+gamma+eplison)*y(2);
    eplison*y(2)-(d+a2+eta)*y(3);
    p*b+gamma*y(2)-d*y(4)-delta*z(4,3)];
lags=[t1 t2 t3]
sol=dde23(ddex1dez,lags,[700 100 300 900],[0,1200]);
plot(sol.x,sol.y(2,:),'LineWidth',1.6);
set(gca,'FontSize',16)
xlabel('t');
ylabel('Number of hosts/Unit');
% legend('I2-SIQRS')
tint=linspace(0,5,10);
sint=deval(sol,tint);
hold on
grid on

% SEIQRS
ddex1dez=@(t,y,z)[(1-p)*b+theta*y(2)-d*y(1)+delta*z(5,3)-beta*z(1,1)*z(3,2)/(1+sigma*z(1,1));
    beta*z(1,1)*z(3,2)/(1+sigma*z(1,1))-Z1*y(2);
    omega*y(2)-Z2*y(3);
    eplison*y(3)-Z3*y(4);
    p*b+xi*y(2)+gamma*y(3)+eta*y(4)-d*y(5)-delta*z(5,3)];
lags=[t1 t2 t3];

S0=(b*d + b*delta - b*d*p)/(d^2 + delta*d);%S0 in disease-free equilibrium
R0=(b*p)/(d+delta);%R0 in disease-free equilibrium
R00=(beta*omega*S0)/(Z1*Z2*(1+sigma*S0));%Basic regeneration number
% Positive equilibrium point values
S1=(Z1*Z2)/(beta*omega - Z1*Z2*sigma);
E1=(Z2*Z3*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
I1=(Z3*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
Q1=(eplison*omega*(Z1*Z2*d^2 + Z1*Z2*d*delta - b*beta*d*omega - b*beta*delta*omega + Z1*Z2*b*d*sigma + Z1*Z2*b*delta*sigma + b*beta*d*omega*p - Z1*Z2*b*d*p*sigma))/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);
R1=(Z1*Z2^2*Z3*d*xi - Z3*b*beta*gamma*omega^2 - b*beta*eplison*eta*omega^2 + Z1*Z2^2*Z3*b*sigma*xi + Z3*b*beta*gamma*omega^2*p + b*beta*eplison*eta*omega^2*p + Z1^2*Z2^2*Z3*b*p*sigma + Z1*Z2*Z3*d*gamma*omega + Z1*Z2*d*eplison*eta*omega - Z2*Z3*b*beta*omega*xi - Z1*Z2*Z3*b*beta*omega*p + Z1*Z2*Z3*b*gamma*omega*sigma + Z1*Z2*b*eplison*eta*omega*sigma + Z2*Z3*b*beta*omega*p*theta + Z2*Z3*b*beta*omega*p*xi - Z1*Z2^2*Z3*b*p*sigma*theta - Z1*Z2^2*Z3*b*p*sigma*xi - Z1*Z2*Z3*b*gamma*omega*p*sigma - Z1*Z2*b*eplison*eta*omega*p*sigma)/(Z1^2*Z2^2*Z3*d*sigma + Z1^2*Z2^2*Z3*delta*sigma + Z3*beta*delta*gamma*omega^2 + beta*delta*eplison*eta*omega^2 - Z1*Z2^2*Z3*d*sigma*theta - Z1*Z2^2*Z3*delta*sigma*theta - Z1*Z2^2*Z3*delta*sigma*xi - Z1*Z2*Z3*beta*d*omega - Z1*Z2*Z3*beta*delta*omega + Z2*Z3*beta*d*omega*theta + Z2*Z3*beta*delta*omega*theta + Z2*Z3*beta*delta*omega*xi - Z1*Z2*Z3*delta*gamma*omega*sigma - Z1*Z2*delta*eplison*eta*omega*sigma);

sol=dde23(ddex1dez,lags,[700 700 100 600 900],[0,1200]);
plot(sol.x,sol.y(3,:),'LineWidth',1.6);
set(gca,'FontSize',16)
xlabel('t');
ylabel('Number of hosts/Unit');
legend('I1-SIRS','I2-SIQRS','I3-SEIQRS')
tint=linspace(0,5,10);
sint=deval(sol,tint);
hold on
grid on