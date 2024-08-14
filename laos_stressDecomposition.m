%carry out method from A geometrical interpretation of large amplitude oscillatory shear response. J. Rheol. 2005

% s(t) = s'(t) + s''(t)

% s'(t) = (s(t)-s(-t))/2
% s"(t) = (s(t)+s(-t))/2

% s'(t) = Gamma'*gamma(t)
% s"(t) = Gamma"*gammdot(t)/omega
load('laosdata.mat')
sp = cell(numel(datas),1);
spp=sp;
Gam_p = sp;
Gam_pp = sp;
for k = 1:numel(datas)
    stress = datas{k}(:,4);
    time = datas{k}(:,1);
    
    nt = numel(stress);
    sp{k} = [];
    spp{k} = [];
    for kk = 1:nt-1
        sp{k}(kk) = (stress(kk)-stress(nt-kk))/2;
        spp{k}(kk) = (stress(kk)+stress(nt-kk))/2;
    end
    Gam_p{k} = sp{k}./datas{k}(2:end,4)';
    Gam_pp{k}= spp{k}.*(2*pi/(datas{1}(end,1)))./(datas{k}(:,3));
end

%% plotting
figNum = 3;
cmap_juarez = (1/256)*[181 1 0; 0 142 198; 253 173 0; 
    224 82 70; 93 122 30; 0 56 107];
cmap_hiroshige = (1/256)*[255 80 77; 252 133 51; 254 168 69;
    255 205 90; 253 231 173; 153 226 217;
    85 192 212; 58 142 176; 28 104 156;
    7 69 116];
%fn = fieldnames(datas);no 
numCurves = size(datas,1);
strainAmps = zeros(numCurves,1);
strainLabs = cell(numCurves,1);
figure(figNum); hold on;

dts = zeros(numel(numCurves:-1:2),2);
dts_max_gampr = zeros(numel(numCurves:-1:2),2);
h=1;
for i = numCurves:-1:1
    
    k = i-2;
    td = datas{i};
    
    dts(h,:) = [max(td(1:end-1,2)),max(sp{i})]';
    
    %plot(max(td(1:end-1,2)),max(sp{i}),'-o','Color',cr,'MarkerFaceColor',cr,'LineWidth',2);
    %plot(td(:,2),td(:,4),'LineWidth',2,'Color',cmap_juarez(4,:));
    %strainAmps(i) = max(td(:,2));
    %strainLabs{i} = [num2str(round(strainAmps(i)*100)),'%'];
    [~,idx_gam0] = min(abs(td(:,2)-max(td(:,2))));
    %dts_max_gampr(h,:) = [td(idx_gam0,2),Gam_p{i}(1,idx_gam0-50)];
    est_slp = (sp{i}(1,idx_gam0-50) -  sp{i}(1,idx_gam0-100)) ./ (td(idx_gam0-50,2) - td(idx_gam0-100,2));
    dts_max_gampr(h,:) = [td(idx_gam0,2),est_slp];
    h=h+1;

end
col = cmap_hiroshige(7,:);

% Paper figure 3d
%first find index at gamma_0
plot(dts_max_gampr(:,1),dts_max_gampr(:,2)./(dts_max_gampr(end,2)),'-o','Color',col,'MarkerFaceColor',col,'LineWidth',0.5);

%plot(dts(:,1),dts(:,2),'-o','Color',col,'MarkerFaceColor',col,'LineWidth',2);
%plot(dts(:,1),dts(:,2)./dts(1,2),'-o','Color',col,'MarkerFaceColor',col);
%% Plot sigma' vs gamma
for i = numCurves:-1:1
    k = i-2;
    td = datas{i};
    %plot(td(1:end-1,2),sp{i}','LineWidth',2,'Color',cmap_hiroshige(10 - k,:));

    plot(td(1:end-1,2),Gam_p{i}','LineWidth',2,'Color',cmap_hiroshige(10 - k,:));
    %plot(td(:,2),td(:,4),'LineWidth',2,'Color',cmap_juarez(4,:));
    %strainAmps(i) = max(td(:,2));
    %strainLabs{i} = [num2str(round(strainAmps(i)*100)),'%'];
end

xlabel('\gamma (-)');
ylabel('\sigma'' (Pa)');
set(gca,'FontSize',13);
set(gca,'linewidth',2)
set(gcf,'Position',[50,50,480,340]);
set(gca,'FontName','Arial');
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'Box','off');
set(gcf, 'Color', 'w')
set(gca,'TickDir','out');

%% Plot sigma'' vs gammadot/omega
% figure(figNum+1); hold on;
% freq = (2*pi/(datas{1}(end,1)));
% for i = numCurves:-1:1
%     k = i-2;
%     td = datas{i};
%     plot(td(1:end-1,3)./freq,spp{i}','LineWidth',2,'Color',cmap_hiroshige(10 - k,:));
% end
% xlabel('$\dot{\gamma}/\omega [-]$','interpreter','latex')
% ylabel('\sigma'''' [Pa]');

%figure(figNum+2); hold on;
%plot max of e
%l = legend(strainLabs);
%legend('boxoff');
% set(l,'Location','northeast');
% set(l,'FontSize',8);
%set(l,'labels',strainLabs);
% Standard figure code:
set(gca,'FontSize',13);
set(gca,'linewidth',2)
set(gcf,'Position',[50,50,580,440]);
set(gca,'FontName','Arial');
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'Box','off');
set(gcf, 'Color', 'w')
set(gca,'TickDir','out');
figure(figNum);
set(gca,'FontSize',13);
set(gca,'linewidth',2)
set(gcf,'Position',[550,50,580,440]);
set(gca,'FontName','Arial');
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'Box','off');
set(gcf, 'Color', 'w')
set(gca,'TickDir','out');