function f = w_numerical_function(x)
        
    global fs beta_num k_val K_val q_val gamma_val alpha_val h_val q_slash_val analyticSol locCheck;

    eta_val = x(1);
    theta_val = x(2);
    
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
    
    q_slash_val = (2/h_val)*sin(q_val*h_val/2);    
    
%     f_vec = [0:100:44100/2-1000];
    f_vec = linspace(0,fs/2-1000,200);
    w_vec = f_vec*2*pi;
    
    w_sel = zeros(size(beta_num));
    for iBeta = 1:length(beta_num)
        diff_val = zeros(size(w_vec));    
        for iW = 1:length(w_vec)
            diff_val(iW) = -gamma_val.^4*q_val.^2*(-0.5*alpha_val + 0.5*(alpha_val - 1)*exp(1i*beta_num(iBeta)*h_val) + (-0.5*alpha_val*exp(1i*k_val*w_vec(iW)) + alpha_val + 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*k_val*w_vec(iW)) + (0.5*alpha_val*exp(1i*k_val*w_vec(iW)) - alpha_val - 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + 0.5)*(-0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val) + 0.5*alpha_val + (-0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + alpha_val*exp(1i*beta_num(iBeta)*h_val) + 0.5*alpha_val*exp(1i*k_val*w_vec(iW)) - alpha_val + 0.5*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) - 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*k_val*w_vec(iW)) + 0.5*exp(1i*beta_num(iBeta)*h_val) - 0.5)*exp(-1i*beta_num(iBeta)*h_val)*exp(-2*1i*k_val*w_vec(iW))./h_val.^2 + (gamma_val.^2*k_val.^2*(0.5*alpha_val + theta_val - 0.5) + gamma_val.^2*k_val.^2*(0.5*alpha_val*exp(1i*k_val*w_vec(iW)) - alpha_val + theta_val*exp(1i*k_val*w_vec(iW)) - 2*theta_val - 0.5*exp(1i*k_val*w_vec(iW)))*exp(1i*k_val*w_vec(iW)) + gamma_val.^2*k_val.^2*(0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val) - 1.0*alpha_val + theta_val*exp(1i*beta_num(iBeta)*h_val) - 2*theta_val - 0.5*exp(1i*beta_num(iBeta)*h_val) + 1.0)*exp(1i*beta_num(iBeta)*h_val) + gamma_val.^2*k_val.^2*(0.5*alpha_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) - alpha_val*exp(1i*beta_num(iBeta)*h_val) - 1.0*alpha_val*exp(1i*k_val*w_vec(iW)) + 2*alpha_val + theta_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) - 2*theta_val*exp(1i*beta_num(iBeta)*h_val) - 2*theta_val*exp(1i*k_val*w_vec(iW)) + 4*theta_val - 0.5*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + 1.0*exp(1i*k_val*w_vec(iW)))*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + h_val.^2*(exp(1i*k_val*w_vec(iW)) - 2)*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW)) + h_val.^2*exp(1i*beta_num(iBeta)*h_val))*(K_val.^2*q_slash_val.^4 + 2*K_val.^2*q_slash_val.^2*exp(1i*beta_num(iBeta)*h_val)./h_val.^2 - 4*K_val.^2*q_slash_val.^2./h_val.^2 + 2*K_val.^2*q_slash_val.^2*exp(-1i*beta_num(iBeta)*h_val)./h_val.^2 + K_val.^2*exp(2*1i*beta_num(iBeta)*h_val)./h_val.^4 - 4*K_val.^2*exp(1i*beta_num(iBeta)*h_val)./h_val.^4 + 6*K_val.^2./h_val.^4 - 4*K_val.^2*exp(-1i*beta_num(iBeta)*h_val)./h_val.^4 + K_val.^2*exp(-2*1i*beta_num(iBeta)*h_val)./h_val.^4 + K_val*eta_val*exp(1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 2*K_val*eta_val*exp(1i*beta_num(iBeta)*h_val)./(h_val.^2*k_val) + K_val*eta_val*exp(1i*beta_num(iBeta)*h_val)*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 2*K_val*eta_val*exp(1i*k_val*w_vec(iW))./(h_val.^2*k_val) + 4*K_val*eta_val./(h_val.^2*k_val) - 2*K_val*eta_val*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val) + K_val*eta_val*exp(-1i*beta_num(iBeta)*h_val)*exp(1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 2*K_val*eta_val*exp(-1i*beta_num(iBeta)*h_val)./(h_val.^2*k_val) + K_val*eta_val*exp(-1i*beta_num(iBeta)*h_val)*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val) - 0.5*alpha_val*gamma_val.^2*q_val.^2*exp(1i*k_val*w_vec(iW)) + alpha_val*gamma_val.^2*q_val.^2 - 0.5*alpha_val*gamma_val.^2*q_val.^2*exp(-1i*k_val*w_vec(iW)) + 0.5*gamma_val.^2*q_val.^2*exp(1i*k_val*w_vec(iW)) + 0.5*gamma_val.^2*q_val.^2*exp(-1i*k_val*w_vec(iW)) + exp(1i*k_val*w_vec(iW))./k_val.^2 - 2./k_val.^2 + exp(-1i*k_val*w_vec(iW))./k_val.^2)*exp(-1i*beta_num(iBeta)*h_val)*exp(-1i*k_val*w_vec(iW))./(h_val.^2*k_val.^2);
            %             diff_val[iW] = -gamma_val**4*q_val**2*(-0.5*alpha_val + 0.5*(alpha_val - 1)*np.exp(1j*beta_num[iBeta]*h_val) + (-0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) + alpha_val + 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*k_val*w_vec[iW]) + (0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) - alpha_val - 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + 0.5)*(-0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val) + 0.5*alpha_val + (-0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + alpha_val*np.exp(1j*beta_num[iBeta]*h_val) + 0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) - alpha_val + 0.5*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) - 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*k_val*w_vec[iW]) + 0.5*np.exp(1j*beta_num[iBeta]*h_val) - 0.5)*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(-2*1j*k_val*w_vec[iW])/h_val**2 + (K_val*k_val*theta_val*(np.exp(1j*beta_num[iBeta]*h_val) - 2)*np.exp(1j*beta_num[iBeta]*h_val) + K_val*k_val*theta_val*(np.exp(1j*k_val*w_vec[iW]) - 2)*np.exp(1j*k_val*w_vec[iW]) + K_val*k_val*theta_val*(np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) - 2*np.exp(1j*beta_num[iBeta]*h_val) - 2*np.exp(1j*k_val*w_vec[iW]) + 4)*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + K_val*k_val*theta_val + 0.5*gamma_val**2*k_val**2*(alpha_val - 1) + gamma_val**2*k_val**2*(0.5*alpha_val*np.exp(1j*k_val*w_vec[iW]) - alpha_val - 0.5*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*k_val*w_vec[iW]) + gamma_val**2*k_val**2*(0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val) - 1.0*alpha_val - 0.5*np.exp(1j*beta_num[iBeta]*h_val) + 1.0)*np.exp(1j*beta_num[iBeta]*h_val) + gamma_val**2*k_val**2*(0.5*alpha_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) - alpha_val*np.exp(1j*beta_num[iBeta]*h_val) - 1.0*alpha_val*np.exp(1j*k_val*w_vec[iW]) + 2*alpha_val - 0.5*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + 1.0*np.exp(1j*k_val*w_vec[iW]))*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + h_val**2*(np.exp(1j*k_val*w_vec[iW]) - 2)*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW]) + h_val**2*np.exp(1j*beta_num[iBeta]*h_val))*(K_val**2*q_slash_val**4 + 2*K_val**2*q_slash_val**2*np.exp(1j*beta_num[iBeta]*h_val)/h_val**2 - 4*K_val**2*q_slash_val**2/h_val**2 + 2*K_val**2*q_slash_val**2*np.exp(-1j*beta_num[iBeta]*h_val)/h_val**2 + K_val**2*np.exp(2*1j*beta_num[iBeta]*h_val)/h_val**4 - 4*K_val**2*np.exp(1j*beta_num[iBeta]*h_val)/h_val**4 + 6*K_val**2/h_val**4 - 4*K_val**2*np.exp(-1j*beta_num[iBeta]*h_val)/h_val**4 + K_val**2*np.exp(-2*1j*beta_num[iBeta]*h_val)/h_val**4 + K_val*eta_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW])/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(1j*beta_num[iBeta]*h_val)/(h_val**2*k_val) + K_val*eta_val*np.exp(1j*beta_num[iBeta]*h_val)*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(1j*k_val*w_vec[iW])/(h_val**2*k_val) + 4*K_val*eta_val/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val) + K_val*eta_val*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(1j*k_val*w_vec[iW])/(h_val**2*k_val) - 2*K_val*eta_val*np.exp(-1j*beta_num[iBeta]*h_val)/(h_val**2*k_val) + K_val*eta_val*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val) - 0.5*alpha_val*gamma_val**2*q_val**2*np.exp(1j*k_val*w_vec[iW]) + alpha_val*gamma_val**2*q_val**2 - 0.5*alpha_val*gamma_val**2*q_val**2*np.exp(-1j*k_val*w_vec[iW]) + 0.5*gamma_val**2*q_val**2*np.exp(1j*k_val*w_vec[iW]) + 0.5*gamma_val**2*q_val**2*np.exp(-1j*k_val*w_vec[iW]) + np.exp(1j*k_val*w_vec[iW])/k_val**2 - 2/k_val**2 + np.exp(-1j*k_val*w_vec[iW])/k_val**2)*np.exp(-1j*beta_num[iBeta]*h_val)*np.exp(-1j*k_val*w_vec[iW])/(h_val**2*k_val**2)        
        end
%         [~,ic] = min(abs(diff_val));
%         w_sel(iBeta) = w_vec(ic);  
        
        ana = diff(abs(diff_val));
        ana = ana(1:end-1).*ana(2:end);
        locs = find(ana<0);  
        if ~isempty(locs)
            w_sel(iBeta) = w_vec(locs(1)+1);     
        else
            [~,ic] = min(abs(diff_val));
            w_sel(iBeta) = w_vec(ic);     
        end
    end
    f = sum((analyticSol(1:locCheck)-w_sel(1:locCheck)).^2)/length(w_sel(1:locCheck));

%     figure;
%     plot(w_vec,abs(diff_val));
%     grid on
    
end
