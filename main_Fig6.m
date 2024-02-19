clear; clc;
rng(3);
tic;

% Parameters
nMonte = 5000; %5000
NIs = 8:8:64;
NTs = [2,4,6,8];
PT = 10e-3; % Transmit power [W]
K = 1; % Rician factor (1 or 10)

% Main loop
PR = nan(nMonte,length(NIs),length(NTs),2);
for iMonte = 1:nMonte
    if mod(iMonte,100) == 0
        fprintf(['iMonte: ',num2str(iMonte),'\n'])
    end

    for iNI = 1:length(NIs)
        NI = NIs(iNI);
        
        for iNTNR = 1:length(NTs)
            NT = NTs(iNTNR);

            % Generate channels hRT, hIT and hRI
            [GRT,GRI,GIT] = func_path_gain();
            HRT = zeros(1,NT);
            HRI = sqrt(GRI) * sqrt(1/2) * (randn(1,NI) + 1i * randn(1,NI)); % Rayleigh
            HIT_LoS = exp(1i * 2 * pi * rand(NI,NT));
            HIT_NLoS = sqrt(1/2) * (randn(NI,NT) + 1i * randn(NI,NT));
            HIT = sqrt(GIT) * (sqrt(K/(1+K)) * HIT_LoS + sqrt(1/(1+K)) * HIT_NLoS); % Rician
    
            % Normalize channels hIT and hRI
            [URI,~,~] = svd(HRI');
            uRI = URI(:,1);
            [UIT,~,~] = svd(HIT);
            uIT = UIT(:,1);

            % Compute Theta
            [~, ~, Theta] = func_OPT_3diag(0, uIT, uRI', NI);
            
            % Compute received signal power
            PR(iMonte,iNI,iNTNR,1) = PT * norm(HRT + HRI*Theta*HIT) ^ 2;
            PR(iMonte,iNI,iNTNR,2) = PT * norm(HRI) ^ 2 * norm(HIT) ^ 2;
        end
    end
end

PR_av = squeeze(mean(PR));
toc;

%% Plot
figure('defaultaxesfontsize',11)
LineW = 2; MarkS = 8;
hold on;
plot(NIs,PR_av(:,4,2)*1e9,'-h','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,4,1)*1e9,'--p','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,3,2)*1e9,'->','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,3,1)*1e9,'--<','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,2,2)*1e9,'-s','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,2,1)*1e9,'--d','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,1,2)*1e9,'-*','linewidth',LineW,'MarkerSize',MarkS)
plot(NIs,PR_av(:,1,1)*1e9,'--o','linewidth',LineW,'MarkerSize',MarkS)
grid on; box on;
if K == 0
    title('Rayleigh fading')
else
    title(['Rician factor = ',num2str(10*log10(K)),' dB'])
end
xlabel('Number of RIS elements');
ylabel('Average received signal power [nW]')
set(gca,'GridLineStyle',':','GridAlpha',0.8,'LineWidth',1.5);
legend('M = 8 - Fully-conn.','Tree-conn.',...
       'M = 6 - Fully-conn.','Tree-conn.',...
       'M = 4 - Fully-conn.','Tree-conn.',...
       'M = 2 - Fully-conn.','Tree-conn.',...
       'location','northwest','numColumns',2);
plots=get(gca, 'Children');
legend(plots([8,6,4,2,7,5,3,1]));
ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
ax.YTick = 0:0.2:1.6;
ax.YLim = [0 1.6];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);