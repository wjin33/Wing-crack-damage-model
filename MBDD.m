function [r_p,r_q,r_epsT,r_epsel,r_epsed,r_epsE,r_epsid,r_Omega,r_a_i,r_f,r_Yd,r_energy] = MBDD()
%%     Damage Consitutive Model (one element test)
clc
close all
clear all

global nu0 E0 N_V0 a_0  Ki Qsym
global FTOL c0 c1 Iter

%% Orientation integral points and weights
BETA=33.2699078510*pi/180;
T1=cos(pi/4)*cos(pi/2-BETA);
T2=cos(BETA);
n_cosines=[1 0 0;
           0 1 0;
           0 0 1;
           sqrt(2)/2   sqrt(2)/2  0;
           sqrt(2)/2  -sqrt(2)/2  0;
           sqrt(2)/2   0          sqrt(2)/2;
           sqrt(2)/2   0         -sqrt(2)/2;
           0           sqrt(2)/2  sqrt(2)/2;
           0           sqrt(2)/2 -sqrt(2)/2;
           T1          T1         T2;
           T1          T1        -T2;
           T1         -T1         T2;
           T1         -T1        -T2;
           T1          T2         T1;
           T1          T2        -T1;
           T1         -T2         T1;
           T1         -T2        -T1;
           T2          T1         T1;
           T2          T1        -T1;
           T2         -T1         T1;
           T2         -T1        -T1;];
n_cosines(22:42,1:3)=-n_cosines(1:21,1:3);   
n_weight(1:3,1)=0.0265214244093;
n_weight(4:9,1)=0.0199301476312;
n_weight(10:21,1)=0.0250712367487;
n_weight(22:42,1)=n_weight(1:21,1);

%% Input parameters for the model




%% Parameters for  computation
FTOL = 1e-6;
Iter=50;




%% Store the parameters
% Pars=zeros(1,6);
% Pars(1)=E0;
% Pars(2)=nu0;
% Pars(3)=N_V0;
% Pars(4)=a_0;
% Pars(5)=Ki;
% Pars(6)=Qsym;

load('calied_MBDD_2.mat');
a= Result';
Pars = a.*[(10^10) 1 1 1 10^7 10^6];
E0   = Pars(1);
nu0  = Pars(2);
N_V0 = Pars(3);
a_0  = Pars(4);
Ki   = Pars(5);
Qsym = Pars(6);
% E0 = 3.6552E10;  % unit=Pa
% nu0 = 0.32;
% N_V0= 690;
% a_0 = 4.86E-2;
% Ki=93.416E7;
% Qsym=10000E6;
c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;

%% Load path illustration:
%  eta = %%  0 = iso, 1 = uniaxial
%% Load path example
%  %tensile-compressive
% eta = [2 2 2];
% steps = 3;
% chargs = [-45e6 45e6 60e6 -200e6];     % s11 vertical-tension
% ninc = [60 1 100 1];


% cyclic-compressive
eta = [2 2 2 2];
steps = 4;
chargs = [200e6 -200e6 250e6 -250e6];     % s11 vertical-tension
ninc = [100 1 100 1];



%%         INITIALIZATIONS:
sigmaT = zeros(3,3); % non-load
% sigmaT = [10E6 0 0;0 0 0;0 0 0;];
sigmaT_v = mat2_mat1(sigmaT);
% A = ai*ones(1,9);
Omega0= N_V0*a_0^3*eye(3);
Omega = Omega0; % undamaged
%========================================================================
%         Strain decomposition:
%         epsT = epsel + epsed + epsid
%========================================================================
epsT = zeros(3,3);
epsel = zeros(3,3);
epsed = zeros(3,3);
epsE = zeros(3,3);
epsid = zeros(3,3);
Yd1 = zeros(3,3);
Yd = zeros(42,1);
a_i=zeros(42,1);
energy = zeros(5,1);

% transfer matrix to vector
epsT_v = mat2_mat1(epsT);
epsel_v = mat2_mat1(epsel); 
epsed_v = mat2_mat1(epsed);
epsE_v = mat2_mat1(epsE);
epsid_v = mat2_mat1(epsid);
Omega_v = mat2_mat1(Omega-Omega0);
%% Storage
r_p(1,:) = sum(sigmaT_v)/3;
r_q(1,:) = sigmaT(1,1)-sigmaT(2,2);
r_sigmaT(1,:) = sigmaT_v;
r_epsel(1,:) = epsel_v;
r_epsed(1,:) = epsed_v;
r_epsid(1,:) = epsid_v;
r_epsE(1,:) = epsE_v;
r_epsT(1,:) = epsT_v;
r_Omega(1,:) = Omega_v;
r_Yd(1,:) = Yd';
r_a_i(1,:) = a_i';
r_energy(1,:) = energy';
r_f(1,:)=-ones(42,1);

%% SIMULATION (for a single stress path component):

% tinc = 2;
for icharg = 1:steps 
    disp(['============= load step #',num2str(icharg),' ============='])
    
    for inc = 1:ninc(icharg) % load increments
        disp(['              increments #',num2str(inc),'              '])
        
        if eta(icharg) == 1 %iso
            dsig = chargs(icharg)/ninc(icharg)*[1 0 0;0 1 0;0 0 1;];
        elseif eta(icharg) == 2 %uniaxial 
            dsig = chargs(icharg)/ninc(icharg)*[1 0 0;0 0 0;0 0 0;];    
        end
        
        sigmaTT = sigmaT + dsig ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        sigma_nn=zeros(42,1);
        tau_nm=zeros(42,1);
        m_cosines=zeros(42,3);
        sigma_mm=zeros(42,1);
        theta=zeros(42,1);
        sigma_max=zeros(42,1);
        sigma_min=zeros(42,1);
        n_theta=zeros(42,3);
        a_thetai=zeros(42,1);
        fd=-ones(42,1);
        
        for i_point=1:42
            sigma_nn(i_point,1)=n_cosines(i_point,:)*sigmaTT*n_cosines(i_point,:)';
            Temp=sigmaTT*n_cosines(i_point,:)'-(n_cosines(i_point,:)*sigmaTT*n_cosines(i_point,:)')*n_cosines(i_point,:)';
            tau_nm(i_point,1)=sqrt(sum(Temp.*Temp));
            if tau_nm(i_point,1)==0
                [V,D]= eig(sigmaTT);
                flag=0;
                for k_D=1:3
                    if D(k_D,k_D)~=sigma_nn(i_point,1)
                        if flag==0
                            D_max=D(k_D,k_D);
                            flag=1;
                            m_cosines(i_point,:)=V(:,k_D)';
                        end
                        if flag==1 && D(k_D,k_D)>D_max
                           m_cosines(i_point,:)=V(:,k_D)';
                        end
                    end
                end
            else
                m_cosines(i_point,:)=Temp'/tau_nm(i_point,1);
            end
            sigma_mm(i_point,1)=m_cosines(i_point,:)*sigmaTT*m_cosines(i_point,:)';
            stress_local=[sigma_mm(i_point,1)  tau_nm(i_point,1); tau_nm(i_point,1)   sigma_nn(i_point,1);];
            [V,D]= eig(stress_local);
            min_local=D(1,1);
            max_local=D(2,2);
            k_D=1;
            if D(2,2)<min_local
                min_local=D(2,2);
                max_local=D(1,1);
                k_D=1;
            end
            sigma_min(i_point,1)=min_local;
            sigma_max(i_point,1)=max_local;
            if V(2,k_D)>=0
                theta(i_point,1)=acos(V(1,k_D));
                if theta(i_point,1)>=pi/2
                    theta(i_point,1)=theta(i_point,1)-pi;
                end
            else
                theta(i_point,1)=-acos(V(1,k_D));
                if theta(i_point,1)<=-pi/2
                    theta(i_point,1)=theta(i_point,1)+pi;
                end       
            end
            a_i(i_point,1) = (n_cosines(i_point,:)*Omega*n_cosines(i_point,:)'/N_V0)^(1/3);
            n_theta(i_point,:)=m_cosines(i_point,:).*cos(theta(i_point,1))+n_cosines(i_point,:).*sin(theta(i_point,1));
            a_thetai(i_point,1)=(n_theta(i_point,:)*Omega*n_theta(i_point,:)'/N_V0)^(1/3);

            if sigma_nn(i_point,1)<=0 %%tensile crack
                C_0 = sigma_nn(i_point,1)^2;
                C_1 = Aij_Bij((sigmaTT*sigmaTT),Ai_Aj(n_cosines(i_point,:)));
                C_2 = Aij_Bij(sigmaTT, Aijkl_Bkl( Ai_Aj_Ak_Al(n_cosines(i_point,:)),sigmaTT));
                Yd(i_point,1)=n_weight(i_point,1)*(c0*C_0 + c1*(C_1-C_2)); 
                fd(i_point,1)=-sigma_nn(i_point,1)-a_i(i_point,1)/(1/Ki+a_i(i_point,1)/Qsym);
            else
                C_1 = Aij_Bij((sigmaTT*sigmaTT),Ai_Aj(n_cosines(i_point,:)));
                C_2 = Aij_Bij(sigmaTT, Aijkl_Bkl( Ai_Aj_Ak_Al(n_cosines(i_point,:)),sigmaTT));
                Yd(i_point,1)=n_weight(i_point,1)*(c1*(C_1-C_2));
                fd(i_point,1)=(cos(theta(i_point,1))*tau_nm(i_point,1)*(a_i(i_point,1)/a_thetai(i_point,1))^2-sigma_min(i_point,1))...
                                -a_thetai(i_point,1)/(1/Ki+a_thetai(i_point,1)/Qsym);
            end
        end
        
%         end_r = length(r_p);
%         if end_r==41
%             fd=fd
%         end
        
        if fd <= FTOL % elastic increment
            % Elastic trial:
            [ depsel, epsed ] = Elas_strain(Omega,Pars,dsig,sigmaT,n_cosines,n_weight);
            depsid = zeros(3,3);
            end_r = length(r_p);
            depsed = epsed-mat1_mat2(r_epsed(end_r,:));
            depsE = depsel + depsed;
            depsT= depsE + depsid;
            dOmega = zeros(3,3);

            sigmaT = sigmaT + dsig;
            sigEner = sigmaT-0.5*dsig;
            Omega = Omega + dOmega;
            epsE = epsE + depsE;
            epsT = epsT + depsT ;
            epsid = epsid + depsid;
            epsel = epsel + depsel; 

            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            % transfer matrix to vector
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsE_v = mat2_mat1(epsE);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsid_v = mat2_mat1(epsid);
            Omega_v = mat2_mat1(Omega-Omega0);

            energy(1,1)=Aij_Bij(sigEner,depsT);
            energy(2,1)=Aij_Bij(sigEner,depsE);
            energy(3,1)=Aij_Bij(sigEner,depsel);
            energy(4,1)=Aij_Bij(sigEner,depsid);
            for i_point=1:21
                Yd1 = Yd1 + Yd(i_point,1)*Ai_Aj(n_cosines(i_point,:));
            end 
            energy(5,1)=Aij_Bij(Yd1,dOmega);
            
            %Start to store components
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            r_sigmaT(end_r,:)=sigmaT_v;
            r_epsT(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsid(end_r,:) = epsid_v;
            r_Omega(end_r,:) = Omega_v;
            r_f(end_r,:) = fd';
            r_Yd(end_r,:) = Yd';
            r_a_i(end_r,:) = a_i';
            r_energy(end_r,:) =r_energy(end_r-1,1:5)+energy';
        else
%             incinc=1;
%             while sum(fd>0) > 0 %&& (incinc<Iter)
                [ depsid, dOmega,d_a_i,fd] = InElas_strain(Pars,dsig,sigmaT,n_cosines,n_weight,sigma_nn,...
                                               tau_nm,theta,sigma_min,a_i,a_thetai,fd,n_theta);
%                  if  fd > FTOL
%                      a_i = a_i + d_a_i;
%                      Omega = Omega + dOmega;
%                      epsid = epsid + depsid;
%                      incinc=incinc+1;
%                  end  
%             end
            
%             dOmega
%             Omega
            
            Omega = Omega + dOmega;
            
            [depsel, epsed] = Elas_strain(Omega,Pars,dsig,sigmaT,n_cosines,n_weight);
            end_r = length(r_p);
            depsed = epsed-mat1_mat2(r_epsed(end_r,:));
            depsE = depsel + depsed;
            depsT= depsE + depsid;
             
            sigmaT = sigmaT + dsig;
            sigEner = sigmaT-0.5*dsig;
            epsE = epsE + depsE;
            epsT = epsT + depsT ;
            epsid = epsid + depsid;
            epsel = epsel + depsel;
            a_i = a_i + d_a_i;

            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            % transfer matrix to vector
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsE_v = mat2_mat1(epsE);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsid_v = mat2_mat1(epsid);
            Omega_v = mat2_mat1(Omega-Omega0);

            energy(1,1)=Aij_Bij(sigEner,depsT);
            energy(2,1)=Aij_Bij(sigEner,depsE);
            energy(3,1)=Aij_Bij(sigEner,depsel);
            energy(4,1)=Aij_Bij(sigEner,depsid);
            
            for i_point=1:21
                Yd1 = Yd1 + Yd(i_point,1)*Ai_Aj(n_cosines(i_point,:));
            end 
            energy(5,1)=Aij_Bij(Yd1,dOmega);
            
            %Start to store components
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            r_sigmaT(end_r,:)=sigmaT_v;
            r_epsT(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsid(end_r,:) = epsid_v;
            r_Omega(end_r,:) = Omega_v;
            r_f(end_r,:) = fd';
            r_Yd(end_r,:) = Yd';
            r_a_i(end_r,:) = a_i';
            r_energy(end_r,:) =r_energy(end_r-1,1:5)+energy';                                
                                           

        end
%         tinc = tinc + 1;% just counts total increments (if several loadings)       
    end % increments

end % loading parts


% return



figure('Name','q(eps1^T,eps3^T)','NumberTitle','off');
plot(r_epsT(:,1),r_q/1000000,'-b',r_epsT(:,3),r_q/1000000,'-r','Linewidth',2)
xlabel('Strain, \epsilon^T','FontSize',20)
ylabel('Deviatoric stress, \sigma_1-\sigma_3 (MPa)','FontSize',20)
legend('\epsilon_1^T', '\epsilon_3^T','Location','Best')
grid
set(gca,'FontSize',20)
% set(gca,'XDir','reverse')
% set(gca,'YDir','reverse')
print -depsc 'MBDD_q_eps1-3-cyc.eps'


% 
% figure('Name','q(eps1^{ed},eps3^{ed})','NumberTitle','off');
% plot(r_epsed(:,1),r_q/1000000,'-b',r_epsed(:,3),r_q/1000000,'-r','Linewidth',2)
% xlabel('Axial strain, \epsilon^{ed}','FontSize',20)
% ylabel('Deviatoric stress, q (MPa)','FontSize',20)
% legend('\epsilon_1^{ed}', '\epsilon_3^{ed}','Location','Best')
% grid
% set(gca,'FontSize',20)
% print -depsc 'MBDD_q_eps_ed_1-3.eps'
% 
% 
% figure('Name','q(eps1^{id},eps3^{id})','NumberTitle','off');
% plot(r_epsid(:,1),r_q/1000000,'-b',r_epsid(:,3),r_q/1000000,'-r','Linewidth',2)
% xlabel('Axial strain, \epsilon^{id}','FontSize',20)
% ylabel('Deviatoric stress, q (MPa)','FontSize',20)
% legend('\epsilon_1^{id}', '\epsilon_3^{id}','Location','Best')
% grid
% set(gca,'FontSize',20)
% print -depsc 'MBDD_q_eps_id_1-3.eps'
% 



% figure('Name','eps33_11(q)','NumberTitle','off')
% plot(r_epsT(:,3)-r_epsT(:,1),(r_sigmaT(:,3)-r_sigmaT(:,1))./10^6,'-bo','Linewidth',1)
% xlabel('\epsilon_{33}-\epsilon_{11}','FontSize',20)
% ylabel('\sigma_{33}-\sigma_{11}','FontSize',20)
% grid
% set(gca,'FontSize',20)
% % print -depsc 'stress_strain.eps'
% 
% 
% 
% 
figure('Name','Omega1(eps1)','NumberTitle','off')
plot(r_q/1000000,r_Omega(:,1),'b-','Linewidth',3)
hold on
plot(r_q/1000000,r_Omega(:,2),'r--','Linewidth',4)
plot(r_q/1000000,r_Omega(:,3),'k-','Linewidth',1.5)
legend('\Omega_1','\Omega_2','\Omega_3','Location','Northwest')
xlabel('Axial stess, \sigma_1 (MPa)','FontSize',20)
ylabel('Damage variable, \Omega','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'MBDD_q_omegas-cyc.eps'
%


figure1 = figure('NumberTitle','off','Name','Energy-tension');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_sigmaT(:,1)/1000000,r_energy(:,1),'Parent',axes1,'LineWidth',2);
plot2 = plot(r_sigmaT(:,1)/1000000,r_energy(:,2),'Parent',axes1,'Linewidth',2);
plot3 = plot(r_sigmaT(:,1)/1000000,r_energy(:,4),'Parent',axes1,'Linewidth',2);
% plot4 = plot(r_sigmaT(:,1)/1000000,r_energy(:,5),'Parent',axes1,'Linewidth',3);

set(plot1,'LineStyle','-','Color',[0 0 1],'DisplayName','External work');
set(plot2,'LineStyle','--','Color',[1 0 0],'DisplayName','Elastic strain energy');
set(plot3,'LineStyle','-.','Color',[0 0 0.5],'DisplayName','Inelastic strain energy');
% set(plot4,'LineStyle','-.','Color',[0 1 1],...
%      'DisplayName','Crack debonding');
% Create xlabel
xlabel('Axial stress, \sigma_{1} (MPa)','FontSize',20);
% Create ylabel
ylabel('Energy (J)','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','Best','FontSize',15);
print -depsc 'MBDD_energy-cyc.eps'


end




%%                  Functions
function [ depsel, epsed ] = Elas_strain(Omega,Pars,dsig,sigmaT,n_cosines,n_weight)
    epsed = zeros(3,3);
    E0  = Pars(1);  % unit=Pa
    nu0 = Pars(2);
    N_V0= Pars(3);
    a_0 = Pars(4);
    global  c0 c1
    b1 = (1+nu0)/E0;
    b2 = nu0/E0;

    depsel = b1*dsig-b2*trace(dsig)*eye(3);
    sigmaTT = sigmaT + dsig ;
    
    for i_point=1:42
        rho = n_cosines(i_point,:)*(Omega-N_V0*a_0^3*eye(3))*n_cosines(i_point,:)';
        C_0 = n_cosines(i_point,:)*sigmaTT*n_cosines(i_point,:)';
        C_1 = sigmaTT*Ai_Aj(n_cosines(i_point,:));
        C_2 = C_0*Ai_Aj(n_cosines(i_point,:));
        if C_0 < 0
            epsed = epsed + n_weight(i_point,1)*rho*(2*c0*C_0*Ai_Aj(n_cosines(i_point,:)) + 2*c1*(C_1-C_2));
        else
            epsed = epsed + n_weight(i_point,1)*rho*2*c1*(C_1-C_2);
        end
    end   

end


function [ depsid, dOmega, d_a_i, fd] = InElas_strain(Pars,dsig,sigmaT,n_cosines,n_weight,sigma_nn,...
                                                tau_nm,theta,sigma_min,a_i,a_thetai,fd,n_theta)

    global  c0 c1 Ki Qsym FTOL Iter
    
    d_a_i=zeros(42,1);
    d_a_thetai=zeros(42,1);
    dOmega=zeros(3,3);
    depsid=zeros(3,3);   
    N_V0=Pars(3);
    
    sigmaTT = sigmaT + dsig ;
    
    for i_point=1:42
        if fd(i_point,1)>=0
            if sigma_nn(i_point,1)>0 %%compressive shearing crack
                a_thetainit=a_thetai(i_point,1);
                incinc=0;
                while (abs(fd(i_point,1)) > FTOL) && (incinc<Iter)                      
                    partial_fd_a = cos(theta(i_point,1))*tau_nm(i_point,1)*(-2)*(a_i(i_point,1)^2/a_thetai(i_point,1)^3)-1/Ki/(1/Ki+a_i(i_point,1)/Qsym)^2;
                    a_thetai(i_point,1)= a_thetai(i_point,1) -fd(i_point,1)/partial_fd_a;
                    fd(i_point,1)=(cos(theta(i_point,1))*tau_nm(i_point,1)*(a_i(i_point,1)/a_thetai(i_point,1))^2-sigma_min(i_point,1))...
                                -a_thetai(i_point,1)/(1/Ki+a_thetai(i_point,1)/Qsym);
                    incinc=incinc+1;           
                end
                d_a_thetai(i_point,1)=a_thetai(i_point,1)-a_thetainit;
            else  %%tensile shearing crack
                a_init=a_i(i_point,1);
                a_i(i_point,1)=-sigma_nn(i_point,1)/Ki/(1+sigma_nn(i_point,1)/Qsym);
                fd(i_point,1)=-sigma_nn(i_point,1)-a_i(i_point,1)/(1/Ki+a_i(i_point,1)/Qsym);
                d_a_i(i_point,1)= a_i(i_point,1)-a_init;
            end
        end
    end
    

    for i_point=1:42
        dOmega = dOmega + N_V0*d_a_i(i_point,1)^3*Ai_Aj(n_cosines(i_point,:)) + N_V0*d_a_thetai(i_point,1)^3*Ai_Aj(n_theta(i_point,:));
    end
    
    for i_point=1:42
        drho=n_cosines(i_point,:)*dOmega*n_cosines(i_point,:)';
        C_0 = n_cosines(i_point,:)*sigmaTT*n_cosines(i_point,:)';
        C_1 = sigmaTT*Ai_Aj(n_cosines(i_point,:));
        C_2 = C_0*Ai_Aj(n_cosines(i_point,:));
        if C_0 < 0
            depsid = depsid + n_weight(i_point,1)*drho*(2*c0*C_0*Ai_Aj(n_cosines(i_point,:)) + 2*c1*(C_1-C_2));
        else
            depsid = depsid + n_weight(i_point,1)*drho*2*c1*(C_1-C_2);
        end
        d_a_i(i_point,1)= (n_cosines(i_point,:)*dOmega*n_cosines(i_point,:)'/N_V0)^(1/3);
    end   

end


function [ C ] = Ai_Aj( A )

C(1:3,1:3)=0;
for i= 1:3
    for j=1:3 
        C(i,j) = C(i,j) + A(i)*A(j); 
    end
end

end

function [ C ] = Ai_Aj_Ak_Al( A )

C(1:3,1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C(i,j,k,l) = C(i,j,k,l) + A(i)*A(j)*A(k)*A(l);
            end
        end
    end
end

end

function [ C ] = Ai_Aj_Ak( A )

C(1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            C(i,j,k) = C(i,j,k) + A(i)*A(j)*A(k);
        end
    end
end

end


function [ C ] = Ai_Bjk( A, B )

C(1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            C(i,j,k) = C(i,j,k) + A(i)*B(j,k);
        end
    end
end

end


function [ C ] = Aijk_Bk( A,B )

C(1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            C(i,j) = C(i,j) + A(i,j,k)*B(k);
        end
    end
end

end

function [ C ] = Ai_Bijk( A,B )

C(1:3,1:3) = 0;
for j = 1:3
    for k = 1:3
        for i = 1:3
            C(j,k) = C(j,k) + A(i)*B(i,j,k);
        end
    end
end

end




function [ C ] = Ai_Bj_Ck( A, B, C)

C(1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            C(i,j,k) = C(i,j,k) + A(i)*B(j)*C(k);
        end
    end
end

end




function [ scalar ] = Aij_Bij( A,B )
%   Double contraction: scalar = A_ij*B_ij

scalar = 0;
for i= 1:3
    for j=1:3 
        scalar = scalar + A(i,j)*B(i,j); 
    end
end
   
end

function [ C ] = Aij_Bkl( A,B )

C(1:3,1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C(i,j,k,l) = A(i,j)*B(k,l);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bij( A,B )

C = zeros(3,3);
for k = 1:3
    for l = 1:3
        for i = 1:3
            for j = 1:3
              C(k,l) = C(k,l)+A(i,j,k,l)*B(i,j);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bkl( A,B )

C = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
              C(i,j) = C(i,j)+A(i,j,k,l)*B(k,l);
            end
        end
    end
end

end


function [ f, Yd ] = fdDP( sigmaT,dsig,A,Pars,omega)

    f=zeros(1,9);
    Yd=zeros(1,9);
    
    E0  = Pars(1);  % unit=Pa
    nu0 = Pars(2);
    miu   = Pars(6);
    K_Ic  = Pars(7);
    K_IIc = Pars(8);
    alpha = Pars(9);
    beta  = Pars(10);

    sigmaTT = sigmaT + dsig ;

    tau12=heaviside((sigmaTT(1)-sigmaTT(2))- miu*(sigmaTT(1)+sigmaTT(2)))*((sigmaTT(1)-sigmaTT(2))/2- miu*(sigmaTT(1)+sigmaTT(2))/2);
    tau13=heaviside((sigmaTT(1)-sigmaTT(3))- miu*(sigmaTT(1)+sigmaTT(3)))*((sigmaTT(1)-sigmaTT(3))/2- miu*(sigmaTT(1)+sigmaTT(3))/2);
    tau23=heaviside((sigmaTT(2)-sigmaTT(3))- miu*(sigmaTT(2)+sigmaTT(3)))*((sigmaTT(2)-sigmaTT(3))/2- miu*(sigmaTT(2)+sigmaTT(3))/2);
    f(1,1)=((sqrt(2)/2)*tau12*A(4)^2/A(1)^2-sigmaTT(2))*sqrt(pi*A(1))/(K_Ic+alpha*A(1))+(sqrt(2)/2)*tau12*A(4)^2*sqrt(pi*A(1))/A(1)^2/(K_IIc+beta*A(1))-1;
    f(1,2)=((sqrt(2)/2)*tau13*A(5)^2/A(2)^2-sigmaTT(3))*sqrt(pi*A(2))/(K_Ic+alpha*A(2))+(sqrt(2)/2)*tau13*A(5)^2*sqrt(pi*A(2))/A(2)^2/(K_IIc+beta*A(2))-1;
    f(1,3)=((sqrt(2)/2)*tau23*A(6)^2/A(3)^2-sigmaTT(3))*sqrt(pi*A(3))/(K_Ic+alpha*A(3))+(sqrt(2)/2)*tau23*A(6)^2*sqrt(pi*A(3))/A(3)^2/(K_IIc+beta*A(3))-1;
    f(1,4)=tau12*sqrt(pi*A(4))-(K_IIc+beta*A(4));
    f(1,5)=tau13*sqrt(pi*A(5))-(K_IIc+beta*A(5));
    f(1,6)=tau23*sqrt(pi*A(6))-(K_IIc+beta*A(6));
    f(1,7)=f(1,4);
    f(1,8)=f(1,5);
    f(1,9)=f(1,6);
    
    gamma= (1- nu0^2)/((2-nu0)*E0);
    zeta = (1-nu0^2)/E0;
    
    temp1=heaviside(tau12*A(4)^2*cos(pi/4)/A(1)^2-sigmaTT(2));
    temp2=heaviside(tau13*A(5)^2*cos(pi/4)/A(2)^2-sigmaTT(3));
    temp3=heaviside(tau23*A(6)^2*cos(pi/4)/A(3)^2-sigmaTT(3));

    Yd(1,1)= 2*(4*pi+8/3)*((-1/3)*omega(4)^(4/3)*omega(1)^(-4/3)*gamma*tau12^2+temp1*zeta*((-1/6)*omega(4)^(4/3)*omega(1)^(-4/3)*tau12^2-...
             (sqrt(2)/3)*omega(4)^(2/3)*omega(1)^(-2/3)*tau12*sigmaTT(2) + sigmaTT(2)^2));
    Yd(1,2)= 2*(4*pi+8/3)*((-1/3)*omega(5)^(4/3)*omega(2)^(-4/3)*gamma*tau13^2+temp2*zeta*((-1/6)*omega(5)^(4/3)*omega(2)^(-4/3)*tau13^2-...
             (sqrt(2)/3)*omega(5)^(2/3)*omega(2)^(-2/3)*tau13*sigmaTT(3) + sigmaTT(3)^2));
    Yd(1,3)= 2*(4*pi+8/3)*((-1/3)*omega(6)^(4/3)*omega(3)^(-4/3)*gamma*tau23^2+temp3*zeta*((-1/6)*omega(6)^(4/3)*omega(3)^(-4/3)*tau23^2-...
             (sqrt(2)/3)*omega(6)^(2/3)*omega(3)^(-2/3)*tau23*sigmaTT(3) + sigmaTT(3)^2));
    Yd(1,4)= (32/3)*gamma*tau12^2 + (4*pi+8/3)*((4/3)*omega(4)^(1/3)*omega(1)^(-1/3)*gamma*tau12^2+...
             temp1*zeta*((2/3)*omega(4)^(1/3)*omega(1)^(-1/3)*tau12^2-(2*sqrt(2)/3)*omega(4)^(-1/3)*omega(1)^(1/3)*tau12*sigmaTT(2)));
    Yd(1,5)= (32/3)*gamma*tau13^2 + (4*pi+8/3)*((4/3)*omega(5)^(1/3)*omega(2)^(-1/3)*gamma*tau13^2+...
             temp2*zeta*((2/3)*omega(5)^(1/3)*omega(2)^(-1/3)*tau13^2-(2*sqrt(2)/3)*omega(5)^(-1/3)*omega(2)^(1/3)*tau13*sigmaTT(3)));
    Yd(1,6)= (32/3)*gamma*tau23^2 + (4*pi+8/3)*((4/3)*omega(6)^(1/3)*omega(3)^(-1/3)*gamma*tau23^2+...
             temp3*zeta*((2/3)*omega(6)^(1/3)*omega(3)^(-1/3)*tau23^2-(2*sqrt(2)/3)*omega(6)^(-1/3)*omega(3)^(1/3)*tau23*sigmaTT(3)));    
    Yd(1,7)= Yd(1,4);
    Yd(1,8)= Yd(1,5);
    Yd(1,9)= Yd(1,6);
end


function [ dY_dsig ] = dY_dsigf( sigmaT )
%
global E0 nu0 a1 a2 a3 a4
E = eye(3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                dY_dsig(i,j,k,l)=2*a1*trace(sigmaT)*E(i,j)*E(k,l)+1/2*a2*(E(i,k)*sigmaT(l,j)+...
                    E(i,l)*sigmaT(j,k)+E(j,l)*sigmaT(i,k)+E(j,k)*sigmaT(i,l))+a3*(E(k,l)*sigmaT(i,j)+...
                1/2*trace(sigmaT)*(E(i,k)*E(j,l)+E(i,l)*E(j,k)))+2*a4*sigmaT(k,l)*E(i,j);
            end
        end
    end
end

end



function [ invA ] = invmat4( A )
% calculate inverse of 4th-order tensor
E6=eye(6);
A_2=mat4_mat2(A,2);
invA_2=A_2\E6;
invA=mat2_mat4(invA_2,1);
end

function [ vector ] = mat2_mat1( matrix )
%==========================================================================
%
%    MAT2_MAT1 
%    Transfer a 3*3 matrix to 6*1 vextor (default format in ABAQUS)
%
%    sig_11 sig_12 sig_13
%    sig_21 sig_22 sig_23 ==> [sig_11 sig_22 sig_33 sig_12 sig_13 sig_23]^T
%    sig_31 sig_32 sig_33
%
%==========================================================================
vector = zeros(1,6);
for i = 1:3
    vector(i) = matrix(i,i);
    for j = i+1:3
        vector(i+j+1) = matrix(i,j);
    end
end

end


function [ matrix ] = mat1_mat2( vector )
matrix = zeros(3,3);
for i = 1:3
    matrix(i,i) = vector(i);
end
matrix(1,2)=vector(4);
matrix(2,1)=vector(4);
matrix(1,3)=vector(5);
matrix(3,1)=vector(5);
matrix(2,3)=vector(6);
matrix(3,2)=vector(6);
end


function [ tensor ] = mat2_mat4( matrix,coe )
      tensor(1:3,1:3,1:3,1:3) = 0;
      if coe==1
          coe1=1;
          coe2=1;
      elseif coe==2
          coe1=2;
          coe2=4;
      end
      tensor(1,1,1,1) = matrix(1,1);
      tensor(1,1,2,2) = matrix(1,2);
      tensor(1,1,3,3) = matrix(1,3);
      tensor(1,1,1,2) = matrix(1,4);
      tensor(1,1,2,1) = matrix(1,4)/coe1;
      tensor(1,1,2,3) = matrix(1,5)/coe1;
      tensor(1,1,3,2) = matrix(1,5)/coe1;
      tensor(1,1,1,3) = matrix(1,6)/coe1;
      tensor(1,1,3,1) = matrix(1,6)/coe1;

      tensor(2,2,1,1) = matrix(2,1);
      tensor(2,2,2,2) = matrix(2,2);
      tensor(2,2,3,3) = matrix(2,3);
      tensor(2,2,1,2) = matrix(2,4)/coe1;
      tensor(2,2,2,1) = matrix(2,4)/coe1;
      tensor(2,2,2,3) = matrix(2,5)/coe1;
      tensor(2,2,3,2) = matrix(2,5)/coe1;
      tensor(2,2,1,3) = matrix(2,6)/coe1;
      tensor(2,2,3,1) = matrix(2,6)/coe1;

      tensor(3,3,1,1) = matrix(3,1);
      tensor(3,3,2,2) = matrix(3,2);
      tensor(3,3,3,3) = matrix(3,3);
      tensor(3,3,1,2) = matrix(3,4)/coe1;
      tensor(3,3,2,1) = matrix(3,4)/coe1;
      tensor(3,3,2,3) = matrix(3,5)/coe1;
      tensor(3,3,3,2) = matrix(3,5)/coe1;
      tensor(3,3,1,3) = matrix(3,6)/coe1;
      tensor(3,3,3,1) = matrix(3,6)/coe1;

      tensor(1,2,1,1) = matrix(4,1)/coe1;
      tensor(1,2,2,2) = matrix(4,2)/coe1;
      tensor(1,2,3,3) = matrix(4,3)/coe1;
      tensor(1,2,1,2) = matrix(4,4)/coe2;
      tensor(1,2,2,1) = matrix(4,4)/coe2;
      tensor(1,2,2,3) = matrix(4,5)/coe2;
      tensor(1,2,3,2) = matrix(4,5)/coe2;
      tensor(1,2,1,3) = matrix(4,6)/coe2;
      tensor(1,2,3,1) = matrix(4,6)/coe2;

      tensor(2,3,1,1) = matrix(5,1)/coe1;
      tensor(2,3,2,2) = matrix(5,2)/coe1;
      tensor(2,3,3,3) = matrix(5,3)/coe1;
      tensor(2,3,1,2) = matrix(5,4)/coe2;
      tensor(2,3,2,1) = matrix(5,4)/coe2;
      tensor(2,3,2,3) = matrix(5,5)/coe2;
      tensor(2,3,3,2) = matrix(5,5)/coe2;
      tensor(2,3,1,3) = matrix(5,6)/coe2;
      tensor(2,3,3,1) = matrix(5,6)/coe2;

      tensor(1,3,1,1) = matrix(6,1)/coe1;
      tensor(1,3,2,2) = matrix(6,2)/coe1;
      tensor(1,3,3,3) = matrix(6,3)/coe1;
      tensor(1,3,1,2) = matrix(6,4)/coe2;
      tensor(1,3,2,1) = matrix(6,4)/coe2;
      tensor(1,3,2,3) = matrix(6,5)/coe2;
      tensor(1,3,3,2) = matrix(6,5)/coe2;
      tensor(1,3,1,3) = matrix(6,6)/coe2;
      tensor(1,3,3,1) = matrix(6,6)/coe2;
      
      tensor(2,1,1,1) = matrix(4,1)/coe1;
      tensor(2,1,2,2) = matrix(4,2)/coe1;
      tensor(2,1,3,3) = matrix(4,3)/coe1;
      tensor(2,1,1,2) = matrix(4,4)/coe2;
      tensor(2,1,2,1) = matrix(4,4)/coe2;
      tensor(2,1,2,3) = matrix(4,5)/coe2;
      tensor(2,1,3,2) = matrix(4,5)/coe2;
      tensor(2,1,1,3) = matrix(4,6)/coe2;
      tensor(2,1,3,1) = matrix(4,6)/coe2;

      tensor(3,2,1,1) = matrix(5,1)/coe1;
      tensor(3,2,2,2) = matrix(5,2)/coe1;
      tensor(3,2,3,3) = matrix(5,3)/coe1;
      tensor(3,2,1,2) = matrix(5,4)/coe2;
      tensor(3,2,2,1) = matrix(5,4)/coe2;
      tensor(3,2,2,3) = matrix(5,5)/coe2;
      tensor(3,2,3,2) = matrix(5,5)/coe2;
      tensor(3,2,1,3) = matrix(5,6)/coe2;
      tensor(3,2,3,1) = matrix(5,6)/coe2;

      tensor(3,1,1,1) = matrix(6,1)/coe1;
      tensor(3,1,2,2) = matrix(6,2)/coe1;
      tensor(3,1,3,3) = matrix(6,3)/coe1;
      tensor(3,1,1,2) = matrix(6,4)/coe2;
      tensor(3,1,2,1) = matrix(6,4)/coe2;
      tensor(3,1,2,3) = matrix(6,5)/coe2;
      tensor(3,1,3,2) = matrix(6,5)/coe2;
      tensor(3,1,1,3) = matrix(6,6)/coe2;
      tensor(3,1,3,1) = matrix(6,6)/coe2;

end

function [ matrix ] = mat4_mat2( tensor,coe )

      matrix = zeros (6,6);
      if coe==1
          coe1=1;
          coe2=1;
      elseif coe==2
          coe1=2;
          coe2=4;
      end
      matrix(1,1)=tensor(1,1,1,1);
      matrix(1,2)=tensor(1,1,2,2);
      matrix(1,3)=tensor(1,1,3,3);
      matrix(1,4)=tensor(1,1,1,2)*coe1;
      matrix(1,5)=tensor(1,1,2,3)*coe1;
      matrix(1,6)=tensor(1,1,1,3)*coe1;

      matrix(2,1)=tensor(2,2,1,1);
      matrix(2,2)=tensor(2,2,2,2);
      matrix(2,3)=tensor(2,2,3,3);
      matrix(2,4)=tensor(2,2,1,2)*coe1;
      matrix(2,5)=tensor(2,2,2,3)*coe1;
      matrix(2,6)=tensor(2,2,1,3)*coe;

      matrix(3,1)=tensor(3,3,1,1);
      matrix(3,2)=tensor(3,3,2,2);
      matrix(3,3)=tensor(3,3,3,3);
      matrix(3,4)=tensor(3,3,1,2)*coe1;
      matrix(3,5)=tensor(3,3,2,3)*coe1;
      matrix(3,6)=tensor(3,3,1,3)*coe1;

      matrix(4,1)=tensor(1,2,1,1)*coe1;
      matrix(4,2)=tensor(1,2,2,2)*coe1;
      matrix(4,3)=tensor(1,2,3,3)*coe1;
      matrix(4,4)=tensor(1,2,1,2)*coe2;
      matrix(4,5)=tensor(1,2,2,3)*coe2;
      matrix(4,6)=tensor(1,2,1,3)*coe2;

      matrix(5,1)=tensor(2,3,1,1)*coe1;
      matrix(5,2)=tensor(2,3,2,2)*coe1;
      matrix(5,3)=tensor(2,3,3,3)*coe1;
      matrix(5,4)=tensor(2,3,1,2)*coe2;
      matrix(5,5)=tensor(2,3,2,3)*coe2;
      matrix(5,6)=tensor(2,3,1,3)*coe2;

      matrix(6,1)=tensor(1,3,1,1)*coe1;
      matrix(6,2)=tensor(1,3,2,2)*coe1;
      matrix(6,3)=tensor(1,3,3,3)*coe1;
      matrix(6,4)=tensor(1,3,1,2)*coe2;
      matrix(6,5)=tensor(1,3,2,3)*coe2;
      matrix(6,6)=tensor(1,3,1,3)*coe2;
end

function [C] = Aijpq_Bpqkl(A,B)
C = zeros(3,3,3,3);
  for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,p,q)*B(p,q,k,l);
                    end
                end
            end
        end
    end
  end

end

function [ C ] = Aijklpq_Bpq(A,B)
C = zeros(3,3,3,3);
  for i=1:3
      for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,k,l,p,q)*B(p,q);
                    end
                end
            end
        end
    end
  end

end

function [matDz, matS] = matDO1(Omega)
global E0 nu0 a1 a2 a3 a4
b1 = (1+nu0)/E0/2;
b2 = nu0/E0;
E = eye(3);
E6 = eye(6);
trOmega = trace(Omega);
matS(1:3,1:3,1:3,1:3) = 0;
matS_2(1:6,1:6) = 0;
matDz(1:3,1:3,1:3,1:3) = 0;
matDz_2(1:6,1:6) = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                matS(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-...
                    b2*E(i,j)*E(k,l)+2*a1*trOmega*E(i,j)*E(k,l)+...
                    0.5*a2*(E(i,k)*Omega(j,l)+E(i,l)*Omega(j,k)+...
                    Omega(i,k)*E(j,l)+Omega(i,l)*E(j,k))+...
                    a3*(E(i,j)*Omega(k,l)+Omega(i,j)*E(k,l))+...
                    a4*trOmega*(E(i,k)*E(j,l)+E(i,l)*E(j,k));
            end
        end
    end
end   
matS_2 = mat4_mat2(matS,2);
matDz_2 = matS_2\E6;
matDz = mat2_mat4(matDz_2,1);

end

function [ P_1 ] = matP_1(sigmaT)
%   P_1 projection tensor
%   P_1=H(sigmaT)-H(-sigmaT)

P_1(1:3,1:3,1:3,1:3)=0;

n1=[1;0;0];
n2=[0;1;0];
n3=[0;0;1];
[V,D] = eig(sigmaT);
n1=V(:,1);
n2=V(:,2);
n3=V(:,3);
P_1=(heaviside(D(1,1))-heaviside(-D(1,1)))*ni_nj_nk_nl(n1)+...
    (heaviside(D(2,2))-heaviside(-D(2,2)))*ni_nj_nk_nl(n2)+...
    (heaviside(D(3,3))-heaviside(-D(3,3)))*ni_nj_nk_nl(n3);

end

function [ P_2 ] = matP_2(sigmaT)

P_2(1:3,1:3,1:3,1:3)=0;

n1=[1;0;0];
n2=[0;1;0];
n3=[0;0;1];
[V,D] = eig(sigmaT);
n1=V(:,1);
n2=V(:,2);
n3=V(:,3);
s(1:3)=0;
s(1)=D(1,1)-max(D(1,1),max(D(2,2),D(3,3)));
s(2)=D(2,2)-max(D(1,1),max(D(2,2),D(3,3)));
s(3)=D(3,3)-max(D(1,1),max(D(2,2),D(3,3)));
for ii=1:3
    if s(ii)==0
        s(ii)=0;
    else
        s(ii)=heaviside(-s(ii));
    end
end

P_2=s(1)*ni_nj_nk_nl(n1)+...
    s(2)*ni_nj_nk_nl(n2)+...
    s(3)*ni_nj_nk_nl(n3);

end

function [ M ] = ni_nj_nk_nl( n )

M(1:3,1:3,1:3,1:3)=0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                M(i,j,k,l) = n(i)*n(j)*n(k)*n(l);
            end
        end
    end
end 
end

function [ fdm,sigm,Omegam,depsidm ] = fd_lam( lambdam, Omega, sigmaT, epsT, epsid, Pars,deps)
%    Iteration solving for lambdan

a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

E=eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
F2F2=Aij_Bij(F2ij,F2ij);
if F2F2==0
    dg_dY=zeros(3,3);
else    
    dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
end
dY_dsig = dY_dsigf(sigmaT);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
% updating depsid
depsidm = lambdam*Aijkl_Bij(dY_dsig,df_dY);
% depsidn = lambdan*Aijkl_Bij(dY_dsig,df_dY);
% updating Omega
dOmegam = lambdam*dg_dY;
% dOmegan = lambdan*dg_dY;
Omegam = Omega+dOmegam;
% Omegan = Omega+dOmegan;
[matD,S] = matDO1(Omega);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end
temp1 = -invmat4(Aijpq_Bpqkl(S,S));
dmatD_dOmega=zeros(3,3,3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        for p=1:3
                            for q=1:3
    dmatD_dOmega(i,j,k,l,m,n)=dmatD_dOmega(i,j,k,l,m,n)+temp1(i,j,p,q)*dS_dOmega(p,q,k,l,m,n);
                            end
                        end
                    end
                end
            end
        end
    end
end

matDem = matD + lambdam*Aijklpq_Bpq(dmatD_dOmega,dg_dY);
epsm = epsT+deps-epsid-depsidm;
sigm = Aijkl_Bkl(matDem,epsm);

[fdm,Yd1] = fdDP(sigm,zeros(3,3),Omegam,Pars);
end

function [ fdm,epsm,Omegam,depsidm ] = fd_lam2( lambdam, Omega, sigmaT, epsid, Pars,dsig)
%    Iteration solving for lambdan
%

a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

E=eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
F2F2=Aij_Bij(F2ij,F2ij);
if F2F2==0
    dg_dY=zeros(3,3);
else    
    dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
end
dY_dsig = dY_dsigf(sigmaT);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
% updating depsid
depsidm = lambdam*Aijkl_Bij(dY_dsig,df_dY);
% depsidn = lambdan*Aijkl_Bij(dY_dsig,df_dY);
% updating Omega
dOmegam = lambdam*dg_dY;
% dOmegan = lambdan*dg_dY;
Omegam = Omega+dOmegam;
% Omegan = Omega+dOmegan;
[matD,S] = matDO1(Omega);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end

Sm=zeros(3,3,3,3);
Sm=S+Aijklpq_Bpq(dS_dOmega,dOmegam);

sigm = sigmaT+dsig;
epsm = Aijkl_Bkl(Sm,sigm)+epsid+depsidm;
[fdm,Yd1] = fdDP(sigm,zeros(3,3),Omegam,Pars);

end
