clear; clc;
rng(3);
tic;

% Parameters
nMonte = 5000; %5000;
NIs = 8:8:64; %8:8:64;
NGs = [1,2,4,8,0]; %[1,2,4,8,0];
NT = 2; % 2 or 8
PT = 10e-3; % Transmit power [W]
K = 1; % Rician factor (1 or 10)


% Main loop
PR = nan(nMonte,length(NIs),length(NGs),2);
parfor iMonte = 1:nMonte
    if mod(iMonte,100) == 0
        fprintf(['iMonte: ',num2str(iMonte),'\n'])
    end

    PR_iter = nan(length(NIs),length(NGs),2);

    for iNI = 1:length(NIs)
        NI = NIs(iNI);

        for iNG = 1:length(NGs)
            NG = NGs(iNG);

            % Generate channels hRT, hIT and hRI     
            [GRT,GRI,GIT] = func_path_gain();
            hRI = sqrt(GRI) * sqrt(1/2) * (randn(1,NI) + 1i * randn(1,NI)); % Rayleigh
            HIT_LoS = exp(1i * 2 * pi * rand(NI,NT));
            HIT_NLoS = sqrt(1/2) * (randn(NI,NT) + 1i * randn(NI,NT));
            HIT = sqrt(GIT) * (sqrt(K/(1+K)) * HIT_LoS + sqrt(1/(1+K)) * HIT_NLoS); % Rician

            % Initialize w (3diag)
            w = randn(NT,1) + 1i * randn(NT,1);
            w = w / norm(w);

            % Optimize Theta (3diag)
            P_tmp = zeros(1,1e3);
            for iIter = 1:1e3
                hRIeff = hRI;
                hITeff = HIT * w;

                hRIeff_norm = hRIeff / norm(hRIeff);
                hITeff_norm = hITeff / norm(hITeff);

                [~, ~, T_tmp] = func_OPT_3diag(0, hITeff_norm, hRIeff_norm, NG);

                w = (hRI*T_tmp*HIT)' / norm(hRI*T_tmp*HIT);

                P_tmp(iIter+1) = (norm(hRI*T_tmp*HIT))^2; % same as abs(hRI*T_tmp*hIT*w)^2;
                % Stopping condition
                if (P_tmp(iIter+1) - P_tmp(iIter))/P_tmp(iIter) < 1e-4
                    break;
                end
            end
            Theta_3diag = T_tmp;

            % Initialize w (group)
            w = randn(NT,1) + 1i * randn(NT,1);
            w = w / norm(w);

            % Optimize Theta (group)
            P_tmp = zeros(1,1e3);
            for iIter = 1:1e3
                hRIeff = hRI;
                hITeff = HIT * w;

                hRIeff_norm = hRIeff / norm(hRIeff);
                hITeff_norm = hITeff / norm(hITeff);

                T_tmp = func_theta(hRIeff_norm,hITeff_norm,NG);

                w = (hRI*T_tmp*HIT)' / norm(hRI*T_tmp*HIT);

                P_tmp(iIter+1) = (norm(hRI*T_tmp*HIT))^2; % same as abs(hRI*T_tmp*hIT*w)^2;
                % Stopping condition
                if (P_tmp(iIter+1) - P_tmp(iIter))/P_tmp(iIter) < 1e-4
                    break;
                end
            end
            Theta_group = T_tmp;

            % Compute received signal power
            PR_iter(iNI,iNG,1) = PT * norm(hRI*Theta_3diag*HIT) ^ 2;
            PR_iter(iNI,iNG,2) = PT * norm(hRI*Theta_group*HIT) ^ 2;
        end
    end
    PR(iMonte,:,:,:,:) = PR_iter;
end

PR_av = squeeze(mean(PR));
toc;

%% Plot
figure('defaultaxesfontsize',11)
LineW = 2; MarkS = 8;
hold on;
plot(NIs,PR_av(:,5,1)*1e9,'-p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Tree-conn.')
plot(NIs,PR_av(:,4,2)*1e9,'-v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Group-conn., group size 8')
plot(NIs,PR_av(:,4,1)*1e9,'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Forest-conn., group size 8')
plot(NIs,PR_av(:,3,2)*1e9,'->','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Group-conn., group size 4')
plot(NIs,PR_av(:,3,1)*1e9,'--<','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Forest-conn., group size 4')
plot(NIs,PR_av(:,2,2)*1e9,'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Group-conn., group size 2')
plot(NIs,PR_av(:,2,1)*1e9,'--d','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Forest-conn., group size 2')
plot(NIs,PR_av(:,1,2)*1e9,'-*','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Single-conn.')
grid on; box on;
if K == 0
    title('Rayleigh fading')
else
    title(['Rician factor = ',num2str(10*log10(K)),' dB'])
end
xlabel('Number of RIS elements');
ylabel('Average received signal power [nW]');
set(gca,'GridLineStyle',':','GridAlpha',0.8,'LineWidth',1.5);
legend('location','northwest');
ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
if NT == 8
    ax.YTick = 0:0.2:1.6;
    ax.YLim = [0 1.6];
elseif NT == 2
    ax.YTick = 0:0.2:1.2;
    ax.YLim = [0 1.2];
end
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);