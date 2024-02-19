clear; clc;

% Parameters
NI = 1:64;

NI8 = NI(rem(NI, 8) == 0);
NI4 = NI(rem(NI, 4) == 0);
NI2 = NI(rem(NI, 2) == 0);

% Fully, group and single connected
FC = NI .* (NI + 1) ./ 2;
GC8 = NI8 .* (8 + 1) ./ 2;
GC4 = NI4 .* (4 + 1) ./ 2;
GC2 = NI2 .* (2 + 1) ./ 2;
SC = NI;

% Tree connected
TC = 2 .* NI - 1;
TGC8 = NI8 .* (2 - 1/8);
TGC4 = NI4 .* (2 - 1/4);

%% Plot
figure('defaultaxesfontsize',11)
hold on;
LineW = 2; MarkS = 8;
plot(NI,FC,'-k','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-conn.')
plot(NI8,GC8,'-^','linewidth',LineW,'MarkerSize',MarkS,'Color',[0.0000, 0.4470, 0.7410],'DisplayName','Group-conn., group size 8')
plot(NI4,GC4,'->','linewidth',LineW,'MarkerSize',MarkS,'Color',[0.8500, 0.3250, 0.0980],'DisplayName','Group-conn., group size 4')
plot(NI2,GC2,'-v','linewidth',LineW,'MarkerSize',MarkS,'Color',[0.4940, 0.1840, 0.5560],'DisplayName','Group-conn., group size 2')
plot(NI,TC,'--k','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Tree-conn.')
plot(NI8,TGC8,'--^','linewidth',LineW,'MarkerSize',MarkS,'Color',[0.0000, 0.4470, 0.7410],'DisplayName','Forest-conn., group size 8')
plot(NI4,TGC4,'-->','linewidth',LineW,'MarkerSize',MarkS,'Color',[0.8500, 0.3250, 0.0980],'DisplayName','Forest-conn., group size 4')
plot(NI,SC,'-','linewidth',LineW,'MarkerSize',MarkS,'Color',[0.4660, 0.6740, 0.1880],'DisplayName','Single-conn.')
grid on; box on;
xlabel('Number of RIS elements');
ylabel('Circuit topology complexity');
set(gca,'GridLineStyle',':','GridAlpha',0.8,'LineWidth',1.5);
legend('location','northwest');
ax = gca;
ax.XTick = 0:8:64;
ax.XLim = [0 64];
ax.YTick = 0:50:300;
ax.YLim = [0 300];
set(gca, 'GridLineStyle', ':','GridAlpha', 0.5);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);