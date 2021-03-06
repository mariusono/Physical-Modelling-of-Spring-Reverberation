clc;
clear all;
close all;

pathSave = cd;
pathSave = fullfile(pathSave, '..\Data');

pathSaveFig = cd;
pathSaveFig = fullfile(pathSaveFig, '..\Figures');


global fs;
fs = 44100;
% fs = 44100*3; % Oversampling ! 

global k_val K_val q_val gamma_val alpha_val q_slash_val analyticSol locCheck ;
global beta_num;

k_val = 1/fs;
alpha_val = 0.5;


E = 210000000000; % young's modulus for steel [Pa] [N/m^2]
rho = 7850; % steel density [kg/m^3]
% % % 
% % % r = 0.001/2; % radius of the wire. half the spring thickness [m]
% % % R = 0.010; % coil radius [m]
% % % 
% % % % Derived parameters
% % % epsilon = 1/R;
% % % K_val = sqrt(E/rho)*(r/(2*L^2));
% % % gamma_val = (1/L) * sqrt(E/rho);
% % % q_val = epsilon*L;

K_vec = [0.02,0.04,0.08,0.16];
q_vec = [200, 400, 600 ,800, 1000, 1500, 2000, 3000];
gamma_vec = [1000,1800];

% Single cases 
K_vec = [0.05];
q_vec = [800];
gamma_vec = [2000];
% % 
% % K_vec = [0.0246];
% % q_vec = [1905];
% % gamma_vec = [1190];


fc_vec = zeros(length(K_vec),length(q_vec),length(gamma_vec));
r_vec = zeros(length(K_vec),length(q_vec),length(gamma_vec));
R_vec = zeros(length(K_vec),length(q_vec),length(gamma_vec));
L_vec = zeros(length(K_vec),length(q_vec),length(gamma_vec));
for iG = 1:length(gamma_vec)
    for iK = 1:length(K_vec)
        for iQ = 1:length(q_vec)

            K_val = K_vec(iK);
            q_val = q_vec(iQ);
            gamma_val = gamma_vec(iG);
            alpha_val = 0.5;

            r_vec(iK,iQ,iG) = 2*K_val*sqrt(E/rho)/gamma_val^2;
            R_vec(iK,iQ,iG) = sqrt(E/rho)/(gamma_val*q_val);
            L_vec(iK,iQ,iG) = sqrt(E/rho)/gamma_val;

            fc_vec(iK,iQ,iG) = 3*K_val*q_val^2/(8*pi*sqrt(5));

        end
    end
end

% % % figure(23455432);
% % % plot(q_vec,fc_vec(:,:,1).','-x','linewidth',2)
% % % grid on
% % % hold all
% % % xlims = get(gca,'xlim')
% % % plot(xlims,[fs/2,fs/2],'k-')
% % % xlabel('$q$ [-]','interpreter','latex')
% % % ylabel('$f$ [Hz]','interpreter','latex')
% % % legend('$\kappa = 0.02$','$\kappa = 0.04$','$\kappa = 0.08$','$\kappa = 0.16$','$f_s/2$','interpreter','latex','location','best')
% % % saveas(figure(23455432),'ir_database_overview.png')

for iG = 1:length(gamma_vec)
    for iK = 1:length(K_vec)
        for iQ = 1:length(q_vec)

            if fc_vec(iK,iQ) > 10000
                continue;
            end

            K_val = K_vec(iK);
            q_val = q_vec(iQ);
            gamma_val = gamma_vec(iG);
            alpha_val = 0.5;

            K_val_str = num2str(K_val);
            K_val_str = strrep(K_val_str,'.','p');

            q_val_str = num2str(q_val);
            q_val_str = strrep(q_val_str,'.','p');

            gamma_val_str = num2str(gamma_val);
            gamma_val_str = strrep(gamma_val_str,'.','p');        

            beta_num = linspace(0,40000,300);
            analyticSol = zeros(size(beta_num));
            for iBeta = 1:length(beta_num)
                analyticSol(iBeta) = sqrt(2)*sqrt(K_val.^2*beta_num(iBeta).^4 - 2*K_val.^2*beta_num(iBeta).^2*q_val.^2 + K_val.^2*q_val.^4 + beta_num(iBeta).^2*gamma_val.^2 + gamma_val.^2*q_val.^2 - sqrt(K_val.^4*beta_num(iBeta).^8 - 4*K_val.^4*beta_num(iBeta).^6*q_val.^2 + 6*K_val.^4*beta_num(iBeta).^4*q_val.^4 - 4*K_val.^4*beta_num(iBeta).^2*q_val.^6 + K_val.^4*q_val.^8 - 2*K_val.^2*beta_num(iBeta).^6*gamma_val.^2 + 6*K_val.^2*beta_num(iBeta).^4*gamma_val.^2*q_val.^2 - 6*K_val.^2*beta_num(iBeta).^2*gamma_val.^2*q_val.^4 + 2*K_val.^2*gamma_val.^2*q_val.^6 + beta_num(iBeta).^4*gamma_val.^4 + 2*beta_num(iBeta).^2*gamma_val.^4*q_val.^2 + gamma_val.^4*q_val.^4))/2;
            end

            analyticSol = abs(analyticSol);

            loc = find(analyticSol/(2*pi)<fs/2,1,'last')+1;

            beta_num =linspace(0,beta_num(loc),200);
            analyticSol = zeros(size(beta_num));
            for iBeta = 1:length(beta_num)
                analyticSol(iBeta) = sqrt(2)*sqrt(K_val.^2*beta_num(iBeta).^4 - 2*K_val.^2*beta_num(iBeta).^2*q_val.^2 + K_val.^2*q_val.^4 + beta_num(iBeta).^2*gamma_val.^2 + gamma_val.^2*q_val.^2 - sqrt(K_val.^4*beta_num(iBeta).^8 - 4*K_val.^4*beta_num(iBeta).^6*q_val.^2 + 6*K_val.^4*beta_num(iBeta).^4*q_val.^4 - 4*K_val.^4*beta_num(iBeta).^2*q_val.^6 + K_val.^4*q_val.^8 - 2*K_val.^2*beta_num(iBeta).^6*gamma_val.^2 + 6*K_val.^2*beta_num(iBeta).^4*gamma_val.^2*q_val.^2 - 6*K_val.^2*beta_num(iBeta).^2*gamma_val.^2*q_val.^4 + 2*K_val.^2*gamma_val.^2*q_val.^6 + beta_num(iBeta).^4*gamma_val.^4 + 2*beta_num(iBeta).^2*gamma_val.^4*q_val.^2 + gamma_val.^4*q_val.^4))/2;
            end
            analyticSol = abs(analyticSol);

            locCheck = find(analyticSol/(2*pi)<fs/2+2000,1,'last') ; % add some bumper freq zones to make sure that the optimal solution has full bandwidth

            % x0 = [0.3,0.0003];
            x0 = [0.1,0.0001];
% %             options = optimset('Display','iter','PlotFcns',@optimplotfval); % plots the loss function
            options = optimset('Display','iter');
            fun = @w_numerical_function;
            
            [x,fval,exitflag,output] = fminsearch(fun,x0,options);
% %             [x fval history] = myproblem([0.1,0.0001]); % history will contain all iterations of x ! 
            
            eta_val = x(1);
            theta_val = x(2);
            
            
            free_params_scheme.eta = eta_val;
            free_params_scheme.theta = theta_val;
            save(fullfile(pathSave,['free_params_scheme_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.mat']),'free_params_scheme')

            theta_plus = (theta_val+abs(theta_val))/2;
            eta_plus = (eta_val+abs(eta_val))/2;

            h1 = 2*gamma_val*k_val*sqrt(theta_plus);
            hVec = [0:0.000001:1];
            val = zeros(length(hVec),1);
            for i = 1:length(hVec)
                val(i) = hVec(i) - sqrt(K_val*k_val*(2*eta_plus + sqrt(4*eta_plus^2+(1+abs(cos(q_val*hVec(i))))^2)));
            end
            h2 = hVec(find(val>0,1,'first'));
            h_val = max([h1,h2]);

%             h_val = 0.01;
            
            q_slash_val = (2/h_val)*sin(q_val*h_val/2);    
%             q_slash_val = q_val;
            
            % f_vec = linspace(0,fs/2-2000,600);
            f_vec = linspace(0,fs/2+2000,600);
            w_vec = f_vec*2*pi;

            w_sel = zeros(size(beta_num));
            legendStr = {};
            p1 = [];
            count = 1;
            for iBeta = 1:length(beta_num)
                diff_val = zeros(size(w_vec));    
                for iW = 1:length(w_vec)
                    diff_val(iW) = -gamma_val.^4*q_val.^2*(-0.5*alpha_val + 0.5*(alpha_val - 1)*exp(1i*beta_num(iBeta)*h_val) + (-0.5*alpha_val*exp(1i*k_val*w_vec(iW)) + alpha_val + 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*k_val*w_vec(iW)) + (0.5*alpha_val*exp(1i*k_val*w_vec(iW)) - alpha_val - 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + 0.5)*(-0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val) + 0.5*alpha_val + (-0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + alpha_val*exp(1i*beta_num(iBeta)*h_val) + 0.5*alpha_val*exp(1i*k_val*w_vec(iW)) - alpha_val + 0.5*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) - 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*k_val*w_vec(iW)) + 0.5*exp(1i*beta_num(iBeta)*h_val) - 0.5)*exp(-1i*beta_num(iBeta)*h_val)*exp(-2*1i*k_val*w_vec(iW))./h_val.^2 + (gamma_val.^2*k_val.^2*(0.5*alpha_val + theta_val - 0.5) + gamma_val.^2*k_val.^2*(0.5*alpha_val*exp(1i*k_val*w_vec(iW)) - alpha_val + theta_val*exp(1i*k_val*w_vec(iW)) - 2*theta_val - 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*k_val*w_vec(iW)) + gamma_val.^2*k_val.^2*(0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val) - 1.0*alpha_val + theta_val*exp(1i*beta_num(iBeta)*h_val) - 2*theta_val - 0.5*exp(1i*beta_num(iBeta)*h_val) + 1.0)*exp(1i*beta_num(iBeta)*h_val) + gamma_val.^2*k_val.^2*(0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) - alpha_val*exp(1i*beta_num(iBeta)*h_val) - 1.0*alpha_val*exp(1i*k_val*w_vec(iW)) + 2*alpha_val + theta_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) - 2*theta_val*exp(1i*beta_num(iBeta)*h_val) - 2*theta_val*exp(1i*k_val*w_vec(iW)) + 4*theta_val - 0.5*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + 1.0*exp(1i*k_val*w_vec(iW)))*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + h_val.^2*(exp(1i*k_val*w_vec(iW)) - 2)*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + h_val.^2*exp(1i*beta_num(iBeta)*h_val))*(K_val.^2*q_slash_val.^4 + 2*K_val.^2*q_slash_val.^2*exp(1i*beta_num(iBeta)*h_val)./h_val.^2 - 4*K_val.^2*q_slash_val.^2./h_val.^2 + 2*K_val.^2*q_slash_val.^2*exp(-1i*beta_num(iBeta)*h_val)./h_val.^2 + K_val.^2*exp(2*1i*beta_num(iBeta)*h_val)./h_val.^4 - 4*K_val.^2*exp(1i*beta_num(iBeta)*h_val)./h_val.^4 + 6*K_val.^2./h_val.^4 - 4*K_val.^2*exp(-1i*beta_num(iBeta)*h_val)./h_val.^4 + K_val.^2*exp(-2*1i*beta_num(iBeta)*h_val)./h_val.^4 + K_val*eta_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 2*K_val*eta_val*exp(1i*beta_num(iBeta)*h_val)./(h_val.^2*k_val) + K_val*eta_val*exp(1i*beta_num(iBeta)*h_val)*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 2*K_val*eta_val*exp(1i*k_val*w_vec(iW))./(h_val.^2*k_val) + 4*K_val*eta_val./(h_val.^2*k_val) - 2*K_val*eta_val*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val) + K_val*eta_val*exp(-1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 2*K_val*eta_val*exp(-1i*beta_num(iBeta)*h_val)./(h_val.^2*k_val) + K_val*eta_val*exp(-1i*beta_num(iBeta)*h_val)*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 0.5*alpha_val*gamma_val.^2*q_val.^2*exp(1i*k_val*w_vec(iW)) + alpha_val*gamma_val.^2*q_val.^2 - 0.5*alpha_val*gamma_val.^2*q_val.^2*exp(-1i*k_val*w_vec(iW)) + 0.5*gamma_val.^2*q_val.^2*exp(1i*k_val*w_vec(iW)) + 0.5*gamma_val.^2*q_val.^2*exp(-1i*k_val*w_vec(iW)) + exp(1i*k_val*w_vec(iW))./k_val.^2 - 2./k_val.^2 + exp(-1i*k_val*w_vec(iW))./k_val.^2)*exp(-1i*beta_num(iBeta)*h_val)*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val.^2);
                    %             diff_val[iW] = -gamma_val**4*q_val**2*(-0.5*alpha_val + 0.5*(alpha_val - 1)*np.exp(1j*beta_num[iBeta]*h_val) + (-0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) + alpha_val + 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*k_val*w_vec[iW]) + (0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) - alpha_val - 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + 0.5)*(-0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val) + 0.5*alpha_val + (-0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + alpha_val*np.exp(1j*beta_num[iBeta]*h_val) + 0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) - alpha_val + 0.5*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) - 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*k_val*w_vec[iW]) + 0.5*np.exp(1j*beta_num[iBeta]*h_val) - 0.5)*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(-2*1j*k_val*w_vec[iW])/h_val**2 + (K_val*k_val*theta_val*(np.exp(1j*beta_num[iBeta]*h_val) - 2)*np.exp(1j*beta_num[iBeta]*h_val) + K_val*k_val*theta_val*(np.exp(1j*k_val*w_vec[iW]) - 2)*np.exp(1j*k_val*w_vec[iW]) + K_val*k_val*theta_val*(np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) - 2*np.exp(1j*beta_num[iBeta]*h_val) - 2*np.exp(1j*k_val*w_vec[iW]) + 4)*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + K_val*k_val*theta_val + 0.5*gamma_val**2*k_val**2*(alpha_val - 1) + gamma_val**2*k_val**2*(0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) - alpha_val - 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*k_val*w_vec[iW]) + gamma_val**2*k_val**2*(0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val) - 1.0*alpha_val - 0.5*np.exp(1j*beta_num[iBeta]*h_val) + 1.0)*np.exp(1j*beta_num[iBeta]*h_val) + gamma_val**2*k_val**2*(0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) - alpha_val*np.exp(1j*beta_num[iBeta]*h_val) - 1.0*alpha_val*np.exp(1j*k_val*w_vec[iW]) + 2*alpha_val - 0.5*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + 1.0*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + h_val**2*(np.exp(1j*k_val*w_vec[iW]) - 2)*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + h_val**2*np.exp(1j*beta_num[iBeta]*h_val))*(K_val**2*q_slash_val**4 + 2*K_val**2*q_slash_val**2*np.exp(1j*beta_num[iBeta]*h_val)/h_val**2 - 4*K_val**2*q_slash_val**2/h_val**2 + 2*K_val**2*q_slash_val**2*np.exp(-1j*beta_num[iBeta]*h_val)/h_val**2 + K_val**2*np.exp(2*1j*beta_num[iBeta]*h_val)/h_val**4 - 4*K_val**2*np.exp(1j*beta_num[iBeta]*h_val)/h_val**4 + 6*K_val**2/h_val**4 - 4*K_val**2*np.exp(-1j*beta_num[iBeta]*h_val)/h_val**4 + K_val**2*np.exp(-2*1j*beta_num[iBeta]*h_val)/h_val**4 + K_val*eta_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW])/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(1j*beta_num[iBeta]*h_val)/(h_val**2*k_val) + K_val*eta_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(1j*k_val*w_vec[iW])/(h_val**2*k_val) + 4*K_val*eta_val/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val) + K_val*eta_val*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW])/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(-1j*beta_num[iBeta]*h_val)/(h_val**2*k_val) + K_val*eta_val*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val) - 0.5*alpha_val*gamma_val**2*q_val**2*np.exp(1j*k_val*w_vec[iW]) + alpha_val*gamma_val**2*q_val**2 - 0.5*alpha_val*gamma_val**2*q_val**2*np.exp(-1j*k_val*w_vec[iW]) + 0.5*gamma_val**2*q_val**2*np.exp(1j*k_val*w_vec[iW]) + 0.5*gamma_val**2*q_val**2*np.exp(-1j*k_val*w_vec[iW]) + np.exp(1j*k_val*w_vec[iW])/k_val**2 - 2/k_val**2 + np.exp(-1j*k_val*w_vec[iW])/k_val**2)*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val**2)        
                end
            %     [~,ic] = min(abs(diff_val));
            %     w_sel(iBeta) = w_vec(ic);     

                ana = diff(abs(diff_val));
                ana = ana(1:end-1).*ana(2:end);
                locs = find(ana<0); 
                if ~isempty(locs)
                    w_sel(iBeta) = w_vec(locs(1)+1);     
                else
                    [~,ic] = min(abs(diff_val));
                    w_sel(iBeta) = w_vec(ic);     
                end    

% %                 if ismember(iBeta,[20,120,150])
% %                     figure(11);
% % %                     plot(w_vec/(2*pi),abs(diff_val),'linewidth',2)
% %                     p1(count) = plot(w_vec/(2*pi),diff_val,'linewidth',2)
% %                     grid on
% %                     hold all
% % %                     plot(w_vec(1:end-1)/(2*pi),diff(abs(diff_val)),'linewidth',2)           
% %                     p3 = plot(w_vec(locs(1))/(2*pi),diff_val(locs(1)),'ro','linewidth',2) ;         
% %                     legendStr = cat(1,legendStr,['$\beta = ',num2str(beta_num(iBeta)),'$']);
% %                     count = count +1;
% %                 end
            end
% % %             xlabel('$f$ [Hz]','interpreter','latex')
% % %             ylabel('det($\mathbf{D}$) [-]','interpreter','latex')
% % %             legendStr = cat(1,legendStr,'Solutions')
% % %             xlim([0,20000])
% % %             legend([p1(:);p3],legendStr,'interpreter','latex');

            figure(11112)
            plot(beta_num,analyticSol/(2*pi),'-','linewidth',2)
            hold all
            plot(beta_num,w_sel/(2*pi),'-','linewidth',2)
            legend('Model Dispersion Rel',...
                   'Numerical Dispersion Rel',...
                   'location','northwest')               
            ylim([0,fs/2])
            grid on
            ylabel('$f$ [Hz]','interpreter','latex')
            xlabel('$\beta$ [rad/m]','interpreter','latex')   
            saveas(figure(11112),fullfile(pathSaveFig,['model_vs_num_disp_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.png']));            
            close(figure(11112));

            %% Finite difference part

            k = k_val;        
            h = h_val;        
            q_slash = q_slash_val;
            K = K_val;
            q = q_val;
            gamma = gamma_val;
            alpha = alpha_val;
            eta = eta_val;
            theta = theta_val;
            
            % HARDCODED
            sig_t = 1.65;
            sig_l = 1.65;

            N = floor(1/h);
            % grid step
            h = 1/N;

            % Matrix operators
            Dxmin = (1/h).*toeplitz([1 -1 zeros(1,N-3)],[1 zeros(1,N-2)]);
            Dxpl = (1/h).*toeplitz([-1 zeros(1,N-2)],[-1 1 zeros(1,N-3)]);
            Dxxxx = (1/h^4).*toeplitz([6 -4 1 zeros(1,N-4)]);
            Dxxxx(1,:) = (1/h^4).*[5 -4 1 zeros(1,N-4)];
            Dxxxx(end,:) = (1/h^4).*fliplr([5 -4 1 zeros(1,N-4)]);
            Dxx = (1/h^2).*toeplitz([-2 1 zeros(1,N-3)]);
            Mx = (1/2).*toeplitz([0 1 zeros(1,N-3)],[0 1 zeros(1,N-3)]); % Mx+ 
            
            % Matrices for update equation
            A1 = (eye(N-1) + eye(N-1)*0.5*gamma.^2*k.^2*q.^2*(1 - alpha) + k*(Dxx*K*eta + eye(N-1)*sig_t))/k.^2;
            B1 = 2*Dxx*K.^2*q_slash.^2 - 2*Dxx*K*eta/k + Dxxxx*K.^2 + eye(N-1)*K.^2*q_slash.^4 - 2*eye(N-1)/k.^2 + eye(N-1)*alpha*gamma.^2*q.^2;
            C1 = (eye(N-1) + eye(N-1)*0.5*gamma.^2*k.^2*q.^2*(1 - alpha) + k*(Dxx*K*eta - eye(N-1)*sig_t))/k.^2;
            D1 = 0.5*Dxmin*gamma.^2*q.^2*(alpha - 1);
            E1 = -Dxmin*alpha*gamma.^2*q.^2;
            F1 = 0.5*Dxmin*gamma.^2*q.^2*(alpha - 1);

            A2 = 0.5*Dxpl*gamma.^2*(1 - alpha);
            B2 = Dxpl*alpha*gamma.^2;
            C2 = 0.5*Dxpl*gamma.^2*(1 - alpha);
            D2 = (Dxx*gamma.^2*k.^2*(0.5*alpha + theta - 0.5) + eye(N-1) + eye(N-1)*k*sig_l)/k.^2;
            E2 = -(Dxx*gamma.^2*k.^2*(alpha + 2*theta) + 2*eye(N-1))/k.^2;
            F2 = (Dxx*gamma.^2*k.^2*(0.5*alpha + theta - 0.5) + eye(N-1) - eye(N-1)*k*sig_l)/k.^2;

            A = sparse([A1 D1; A2 D2]);
            B = sparse([B1 E1; B2 E2]);
            C = sparse([C1 F1; C2 F2]);

            % Sweep tone Excitation signal
            excitation = sweeptone(0.5,3,fs,'SweepFrequencyRange',[20 20000])./10;
            Ntime = length(excitation);

% % % %             Dirac Excitation signal 
% % %             dur = 2; % s
% % %             Ntime = floor(dur*fs);
% % %             excitation = zeros(Ntime,1);
% % %             excitation(1) = 0.1;

            t = [0:1/fs:length(excitation).*1/fs - 1/fs];

% %             figure(444);
% %             plot(t,excitation,'linewidth',1)
% %             grid on
% %             xlabel('Time [s]')
% %             ylabel('Amplitude [m]')
    
            % % % Finite Difference Scheme !
            % Domain is from u_{1} to u_{N-1}
            u2 = zeros(N-1,1); 
            u2(1) = excitation(1);

            u1 = u2; 
            u = zeros(N-1,1); 

            zeta2 = zeros(N-1,1);
            zeta1 = zeros(N-1,1);
            zeta = zeros(N-1,1);

            w2 = [u2;zeta2];
            w1 = [u1;zeta1];
            w = [u;zeta];

            out = zeros(Ntime,1);
            out2 = zeros(Ntime,1);
            out3 = zeros(Ntime,1);

            x = linspace(0,1,N+1); % N+1 points for N segments in x scaled.
    %         tic
            count = 1;
            for n=2:Ntime

                w = A\(-B*w1-C*w2);

                u = w(1:length(w)/2);
                zeta = w(1+length(w)/2:end);

                % include excitation part:
                u(1) = excitation(n) + u(1);
                w = [u;zeta];

                uPlot = [0;u;0];
                
                w2 = w1; % update
                w1 = w; % update

                out(n) = u(N-1);
                out2(n) = sign(u(N-1))*sqrt(u(N-1)^2+zeta(N-1)^2);                 
                out3(n) = zeta(N-1);
                
            end
    %         toc


            irEstimate = impzest(excitation,out);
            irEstimate_zeta = impzest(excitation,out3);
            irEstimate_zeta_u_comb = impzest(excitation,out2);
            
            save(fullfile(pathSave,['irEstimate_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.mat']),'irEstimate')
            save(fullfile(pathSave,['irEstimate_zeta_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.mat']),'irEstimate_zeta')
            save(fullfile(pathSave,['irEstimate_zeta_u_comb_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.mat']),'irEstimate_zeta_u_comb')

            t_ir = [0:1/fs:length(irEstimate).*1/fs - 1/fs];

%             soundsc(irEstimate,fs)

            figure(222);
            hold all
            plot(t_ir,irEstimate,'-','linewidth',2)
            grid on
            xlabel('Time [s]')
            ylabel('Amplitude [m]')
            saveas(figure(222),fullfile(pathSaveFig,['irEstimate_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.png']));
            close(figure(222));

            irEstimate_stereo = [irEstimate,irEstimate]./max(abs(irEstimate));
            irEstimate_zeta_stereo = [irEstimate_zeta,irEstimate_zeta]./max(abs(irEstimate_zeta));
            irEstimate_zeta_u_comb_stereo = [irEstimate_zeta_u_comb,irEstimate_zeta_u_comb]./max(abs(irEstimate_zeta_u_comb));

            audiowrite(fullfile(pathSave,['irEstimate_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.wav']),irEstimate_stereo,fs)
            audiowrite(fullfile(pathSave,['irEstimate_zeta_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.wav']),irEstimate_stereo,fs)
            audiowrite(fullfile(pathSave,['irEstimate_zeta_u_comb_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.wav']),irEstimate_stereo,fs)

            M = 256;
            % overlap = M/2; % overlap
            overlap = round(M*3/4); % overlap
            K = M*2;

            % win = ones(M,1);
            win = window(@hann,M);

            [stft,freq,time] = spectrogram(irEstimate,win,overlap,K,fs);
%             [stft,freq,time] = spectrogram(irEstimate_zeta,win,overlap,K,fs);

            figure(8989987);
            surf(time,freq,20*log10(abs(stft)))
            view([0 90])
            shading interp
            xlabel('Time [s]')
            ylabel('Freq [Hz]')
            xlim([0,1])
            ylim([0,20000])
            caxis([-60,20])
            set(gca,'yscale','linear');
            saveas(figure(8989987),fullfile(pathSaveFig,['spectrogram_ir_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.png']));
            hold all

            group_vel_num = diff(w_sel)./diff(beta_num);
            phase_vel_num = w_sel./beta_num;
            
            group_vel = diff(analyticSol)./diff(beta_num);
            phase_vel = analyticSol./beta_num;  
            
            Td = (h.*N)./(group_vel);
            Td_num = (h.*N)./(group_vel_num);

            p3 = plot3(Td,analyticSol(1:end-1)./(2*pi),ones(size(Td)).*100,'k-','linewidth',1);
            plot3(Td.*3,analyticSol(1:end-1)./(2*pi),ones(size(Td)).*100,'k-','linewidth',1)
            saveas(figure(8989987),fullfile(pathSaveFig,['spectrogram_ir_with_echoes_K_',K_val_str,'_q_',q_val_str,'_gamma_',gamma_val_str,'.png']));
            legend(p3,'Dispersive Echoes - Model')
            close(figure(8989987));


% % %             [sound,fs_sound] = audioread('repeated_E.wav');
% % %             sound = sound(:,1);
% % % %             figure;
% % % %             plot(sound)
% % %             soundsc(sound,fs_sound)
% % % 
% % %             sound_reverb = conv(irEstimate,sound);
% % % 
% % %             soundsc(sound_reverb,fs_sound)
            
        end
    end
end




