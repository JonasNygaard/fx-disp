%% Main script: "Cross-sectional return dispersion and currency momentum"
% ########################################################################### %
% Script for replicating the main results in Eriksen (2019): "Cross-sectional %
% return dispersion and currency momentum".                                   %
%                                                                             %
% Written by:                                                                 %
% Jonas N. Eriksen                                                            %
% Department of Economics and Business Economics                              %
% Aarhus University and CREATES                                               %
%                                                                             %
% Encoding: UTF8                                                              %
% Last modified: June, 2017                                                   %
% ########################################################################### %

clear; clc; tStart = tic; close all; format shortg; c = clock;
disp('# ****************************************************************** #');
disp('# Running the mainRep.m script                                        ');
fprintf('# Code initiated: %.0f:%.0f on %.0f / %0.f - %0.f \n',c([4:5 3 2 1]));
disp('# ****************************************************************** #');

%% Loading and preparing data
% ########################################################################### %
%{
    Data are constructed in the data/ folder using the scripts stored therein. 
    Here we load the processed data from matfiles stored in the data/matfiles/ 
    folder. For more details on the data and the data construction, we refer to 
    the Markdown file located in the data/ folder.
%}
% ########################################################################### %

disp('<*> Loading and preparing data');

% Reading in data from matfiles
load('data/businessCycles.mat');
load('data/currencyFactors.mat')
load('data/currencyPortfolios.mat')

% Create date and datenum vectors for graphs
datenumIndx         = datenum({'01-Nov-1983';'02-Jan-2017'});
datevecIndx         = datevec(datenumIndx(1):1:datenumIndx(2));
datevecIndx         = datevecIndx(datevecIndx(2:end,2) ~= datevecIndx(1:end-1,2),1:3);
datenumIndx         = datenum(datevecIndx);
nberDates(nberDates==0) = NaN;

% Building DOL and CAR factors following Lustig et al. (2011).  
dolFactor           = nanmean(carPortfolios,2); 
dolFactorBA         = nanmean(carPortfoliosBA,2); 
hmlFactor           = carPortfolios(:,end) - carPortfolios(:,1);
hmlFactorBA         = carPortfoliosBA(:,end) - carPortfoliosBA(:,1);
dolFactorMom        = nanmean(momPortfolios,2); 
dolFactorMomBA      = nanmean(momPortfoliosBA,2); 
momFactor           = momPortfolios(:,end) - momPortfolios(:,1);
momFactorBA         = momPortfoliosBA(:,end) - momPortfoliosBA(:,1);
dolFactorGap        = nanmean(gapPortfolios,2); 
dolFactorGapBA      = nanmean(gapPortfoliosBA,2); 
gapFactor           = gapPortfolios(:,end) - gapPortfolios(:,1);
gapFactorBA         = gapPortfoliosBA(:,end) - gapPortfoliosBA(:,1);

statMomentum        = cell(3,1);
statMomentum{1,1}   = portfolioStatistics([momPortfolios dolFactorMom  momFactor]);
statMomentum{2,1}   = portfolioStatistics([momPortfoliosSP nanmean(momPortfoliosSP,2)  momPortfoliosSP(:,end)-momPortfoliosSP(:,1)]);
statMomentum{3,1}   = portfolioStatistics([momPortfoliosFD nanmean(momPortfoliosFD,2)  momPortfoliosFD(:,end)-momPortfoliosFD(:,1)]);

statMomentumBA      = cell(3,1);
statMomentumBA{1,1} = portfolioStatistics([momPortfoliosBA dolFactorMomBA  momFactorBA]);
statMomentumBA{2,1} = portfolioStatistics([momPortfoliosBASP nanmean(momPortfoliosBASP,2)  momPortfoliosBASP(:,end)-momPortfoliosBASP(:,1)]);
statMomentumBA{3,1} = portfolioStatistics([momPortfoliosBAFD nanmean(momPortfoliosBAFD,2)  momPortfoliosBAFD(:,end)-momPortfoliosBAFD(:,1)]);

statCarry           = cell(3,1);
statCarry{1,1}      = portfolioStatistics([carPortfolios dolFactor  hmlFactor]);
statCarry{2,1}      = portfolioStatistics([carPortfoliosSP nanmean(carPortfoliosSP,2)  carPortfoliosSP(:,end)-carPortfoliosSP(:,1)]);
statCarry{3,1}      = portfolioStatistics([carPortfoliosFD nanmean(carPortfoliosFD,2)  carPortfoliosFD(:,end)-carPortfoliosFD(:,1)]);

statCarryBA         = cell(3,1);
statCarryBA{1,1}    = portfolioStatistics([carPortfoliosBA dolFactorBA  hmlFactorBA]);
statCarryBA{2,1}    = portfolioStatistics([carPortfoliosBASP nanmean(carPortfoliosBASP,2)  carPortfoliosBASP(:,end)-carPortfoliosBASP(:,1)]);
statCarryBA{3,1}    = portfolioStatistics([carPortfoliosBAFD nanmean(carPortfoliosBAFD,2)  carPortfoliosBAFD(:,end)-carPortfoliosBAFD(:,1)]);

statGAP         = cell(3,1);
statGAP{1,1}    = portfolioStatistics([gapPortfolios dolFactorGap  gapFactor]);
statGAP{2,1}    = portfolioStatistics([gapPortfoliosSP nanmean(gapPortfoliosSP,2)  gapPortfoliosSP(:,end)-gapPortfoliosSP(:,1)]);
statGAP{3,1}    = portfolioStatistics([gapPortfoliosFD nanmean(gapPortfoliosFD,2)  gapPortfoliosFD(:,end)-gapPortfoliosFD(:,1)]);

statGAPBA       = cell(3,1);
statGAPBA{1,1}  = portfolioStatistics([gapPortfoliosBA dolFactorGapBA  gapFactorBA]);
statGAPBA{2,1}  = portfolioStatistics([gapPortfoliosBASP nanmean(gapPortfoliosBASP,2)  gapPortfoliosBASP(:,end)-gapPortfoliosBASP(:,1)]);
statGAPBA{3,1}  = portfolioStatistics([gapPortfoliosBAFD nanmean(gapPortfoliosBAFD,2)  gapPortfoliosBAFD(:,end)-gapPortfoliosBAFD(:,1)]);


% Writing descriptive statistics to table
headerFile  = {'P1','P2','P3','P4','P5','DOL','HML'};
fid = fopen('output/table1.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','',headerFile{:,:});
fprintf(fid,'%s \\\\\\midrule\n','\multicolumn{8}{c}{Panel A: Momentum portfolios}');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Mean',statMomentumBA{1,1}(1,:));
fprintf(fid,'%s & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] \\\\\n','',statMomentumBA{1,1}(2,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Std',statMomentumBA{1,1}(3,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','SR',statMomentumBA{1,1}(4,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Skew',statMomentumBA{1,1}(5,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Kurt',statMomentumBA{1,1}(6,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s \\\\\n','spot',statMomentumBA{2,1}(1,1:5),'','');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s \\\\\\midrule\n','fd',statMomentumBA{3,1}(1,1:5),'','');
fprintf(fid,'%s \\\\\\midrule\n','\multicolumn{8}{c}{Panel B: Carry portfolios}');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Mean',statCarryBA{1,1}(1,:));
fprintf(fid,'%s & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] \\\\\n','',statCarryBA{1,1}(2,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Std',statCarryBA{1,1}(3,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','SR',statCarryBA{1,1}(4,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Skew',statCarryBA{1,1}(5,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Kurt',statCarryBA{1,1}(6,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s \\\\\n','spot',statCarryBA{2,1}(1,1:5),'','');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s \\\\\\midrule\n','fd',statCarryBA{3,1}(1,1:5),'','');
fprintf(fid,'%s \\\\\\midrule\n','\multicolumn{8}{c}{Panel C: Output gap portfolios}');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Mean',statGAPBA{1,1}(1,:));
fprintf(fid,'%s & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] & [%.2f] \\\\\n','',statGAPBA{1,1}(2,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Std',statGAPBA{1,1}(3,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','SR',statGAPBA{1,1}(4,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Skew',statGAPBA{1,1}(5,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Kurt',statGAPBA{1,1}(6,:));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s \\\\\n','spot',statGAPBA{2,1}(1,1:5),'','');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s \\\\\n','fd',statGAPBA{3,1}(1,1:5),'','');
fprintf(fid,'\n');
fclose(fid);

% Plotting cumulative returns currency strategies
figure;
hold on
b1 = bar(datenumIndx,nberDates.*199);
b2 = bar(datenumIndx,-nberDates.*49);
b1.EdgeColor = [0.8 0.8 0.8];
b1.FaceColor = [0.8 0.8 0.8];
b2.EdgeColor = [0.8 0.8 0.8];
b2.FaceColor = [0.8 0.8 0.8];
b1.ShowBaseLine = 'Off';
b2.ShowBaseLine = 'Off';
p1 = plot(datenumIndx,nancumsum(momFactorBA));
p1.LineStyle = '-'; p1.Color = colorBrewer(1); p1.LineWidth = 1.2;
p2 = plot(datenumIndx,nancumsum(hmlFactorBA));
p2.LineStyle = '-'; p2.Color = colorBrewer(2); p2.LineWidth = 1.2;
p3 = plot(datenumIndx,nancumsum(gapFactorBA));
p3.LineStyle = '-'; p3.Color = colorBrewer(3); p3.LineWidth = 1.2;
p4 = plot(datenumIndx,nancumsum(dolFactorBA));
p4.LineStyle = '-'; p4.Color = colorBrewer(4); p4.LineWidth = 1.2;
hold off
datetick('x','yyyy');
box on
axis([-inf inf -50 200]);
ylabel('Cumulative net excess returns (in % p.a.)');
leg = legend([p1 p2 p3 p4],'MOM','CAR','GAP','DOL');
set(leg,'box','off','Location','NorthWest','FontSize',8);
set(gcf, 'PaperUnits', 'Centimeters','PaperSize',[21 11],'PaperPosition',[0.5 0.5 20 10]);
print(gcf,'-depsc','output/figure1.eps');


%% Computing cross-sectional return dispersion
% ########################################################################### %
%{
    In this part, we compute and check the impact of currency return dispersion
    on currency momentum portfolios. We compute currency return dispersion in
    the following way: RD = (1/5) * sum |RX - DOL|. 
%}
% ########################################################################### %

disp('<*> Computing cross-sectional return dispersion');

dispFX              = nanstd(carPortfolios,1,2);

% Computing global FX volatility innovations
dispARcoef          = [ones(size(dispFX,1)-2,1) dispFX(2:end-1,1) ]\dispFX(3:end,1);
dispFXinno          = dispFX(2:end,1) - [ones(size(dispFX,1)-1,1) dispFX(1:end-1,1)]*dispARcoef;
dispFXinno          = [NaN(1,1);dispFXinno];

figure;
subplot(2,1,1);
hold on
b1 = bar(datenumIndx,nberDates.*3.98);
b1.EdgeColor = [0.8 0.8 0.8];
b1.FaceColor = [0.8 0.8 0.8];
b1.ShowBaseLine = 'Off';
p1 = plot(datenumIndx,dispFX);
p1(1).Color = colorBrewer(1);
p1(1).LineWidth = 1.2;
hold off
datetick('x','yyyy');
axis([-inf inf 0 4]);
box on
title('Cross-sectional dispersion','FontSize',10,'FontWeight','Normal');
ylabel('Level (in %)');
leg = legend(p1,'$CSV_{t}$');
set(leg,'Interpreter','latex');
set(leg,'box','off','Location','NorthEast','FontSize',8);
subplot(2,1,2);
hold on
b1 = bar(datenumIndx,nberDates.*3.98);
b2 = bar(datenumIndx,-nberDates.*1.98);
b1.EdgeColor = [0.8 0.8 0.8];
b1.FaceColor = [0.8 0.8 0.8];
b2.EdgeColor = [0.8 0.8 0.8];
b2.FaceColor = [0.8 0.8 0.8];
b1.ShowBaseLine = 'Off';
b2.ShowBaseLine = 'Off';
p1 = plot(datenumIndx,dispFXinno);
p1(1).Color = colorBrewer(2);
p1(1).LineStyle = '--';
p1(1).LineWidth = 1.2;
hold off
datetick('x','yyyy');
axis([-inf inf -2 4]);
box on
title('Cross-sectional dispersion innovations','FontSize',10,'FontWeight','Normal');
ylabel('Innovations (in %)');
leg = legend(p1,'$\Delta CSV_{t}$');
set(leg,'Interpreter','latex');
set(leg,'box','off','Location','NorthEast','FontSize',8);
set(gcf, 'PaperUnits', 'Centimeters','PaperSize',[21 16],'PaperPosition',[0.5 0.5 20 15]);
print(gcf,'-depsc','output/figure2.eps');

%% Cross-sectional asset pricing: Return dispersion
% ########################################################################### %
%{
    In this part, we conduct simple cross-sectional asset pricing tests. 
    We consider linear GMM estimation following Hansen (1982) and supplement 
    with Fama-MacBeth two-pass cross-sectional regressions. 
%}
% ########################################################################### %

disp('<*> Cross-sectional asset pricing: Return dispersion');
gmmMomDisp      = linearGMM(momPortfoliosBA(3:end,:),[dolFactor(3:end,1) dispFXinno(3:end,1)]);
gmmMomDispU     = linearGMM(momPortfolios(3:end,:),[dolFactor(3:end,1) dispFXinno(3:end,1)]);

% Fama-MacBeth asset pricing tests for all countries
[tsMomDisp,fmMomDisp]   = famaMacBeth(momPortfoliosBA(3:end,:),[dolFactor(3:end,1) dispFXinno(3:end,:)]);
[tsMomDispU,fmMomDispU] = famaMacBeth(momPortfolios(3:end,:),[dolFactor(3:end,1) dispFXinno(3:end,:)]);

coefMom     = linRegNWA(dispFXinno(3:end,1),momPortfoliosBA(3:end,:),1);
coefMomU   = linRegNWA(dispFXinno(3:end,1),momPortfolios(3:end,:),1);

% Scaling coefficients for easier interpretation
coefMom                 = coefMom.b(2:end)./(coefMom.r2./100);
coefMomU                = coefMomU.b(2:end)./(coefMomU.r2./100);

% Building factor-mimicking portfolios
mimickMom               = momPortfoliosBA(3:end,:)*coefMom;
mimickMomU              = momPortfolios(3:end,:)*coefMomU;

% Running GMM estimations
gmmMimickDisp           = linearGMM(momPortfoliosBA(3:end,:),[dolFactor(3:end,1) mimickMom]);
gmmMimickDispU          = linearGMM(momPortfolios(3:end,:),[dolFactor(3:end,1) mimickMomU]);

% Fama-MacBeth asset pricing tests
[tsMimickDisp,fmMimickDisp]     = famaMacBeth(momPortfoliosBA(3:end,:),[dolFactor(3:end,1) mimickMom]);
[tsMimickDispU,fmMimickDispU]   = famaMacBeth(momPortfolios(3:end,:),[dolFactor(3:end,1) mimickMomU]);


headerFile1 = {'b$_{DOL}$','b$_{\Delta CSV}$',...
'$\lambda_{DOL}$','$\lambda_{\Delta CSV}$','R$^{2}$','HJ'};
headerFile2 = {'b$_{DOL}$','b$_{\Delta CSV^{\text{FM}}}$',...
'$\lambda_{DOL}$','$\lambda_{\Delta CSV^{\text{FM}}}$','R$^{2}$','HJ'};
headerFile3 = {'$\beta_{\text{DOL}}$','$\beta_{\Delta CSV}$','R$^{2}\left(\%\right)$'};
headerFile4 = {'$\beta_{\text{DOL}}$','$\beta_{\Delta CSV^{\text{FM}}}$','R$^{2}\left(\%\right)$'};

% Writing results for factor-mimicking portfolios to latex table
fid = fopen('output/table2.tex','w');
fprintf(fid,'%s & %s & %s & %s \\\\\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n','','\multicolumn{2}{c}{Factor loadings}','\multicolumn{2}{c}{Risk prices}','\multicolumn{2}{c}{Fit}');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','',headerFile1{:,:});
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','GMM',gmmMomDisp.b,gmmMomDisp.lambda,gmmMomDisp.r2,gmmMomDisp.hj);
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & %s & [%.2f] \\\\\\midrule\n','',gmmMomDisp.btstat,gmmMomDisp.tlambda,'',gmmMomDisp.hjp);
fprintf(fid,'%s & %s & %s \\\\\\cmidrule(lr){2-6}\\cmidrule(lr){7-7}\n','','\multicolumn{5}{c}{Factor betas}','FMB');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','','P1','P2','P3','P4','P5','$\lambda$');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',headerFile3{1},tsMomDisp.beta(2,:),fmMomDisp.lambda(1));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\\n','',tsMomDisp.tstat(2,:),fmMomDisp.tstat(1));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',headerFile3{2},tsMomDisp.beta(3,:),fmMomDisp.lambda(2));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\\n','',tsMomDisp.tstat(3,:),fmMomDisp.tstat(2));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s \\\\\n',headerFile3{3},tsMomDisp.r2,'');
fprintf(fid,'\n');
fclose(fid);


% Writing results for factor-mimicking portfolios to latex table
fid = fopen('output/table3.tex','w');
fprintf(fid,'%s & %s & %s & %s \\\\\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n','','\multicolumn{2}{c}{Factor loadings}','\multicolumn{2}{c}{Risk prices}','\multicolumn{2}{c}{Fit}');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','',headerFile2{:,:});
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','GMM',gmmMimickDisp.b,gmmMimickDisp.lambda,gmmMimickDisp.r2,gmmMimickDisp.hj);
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & %s & [%.2f] \\\\\\midrule\n','',gmmMimickDisp.btstat,gmmMimickDisp.tlambda,'',gmmMimickDisp.hjp);
fprintf(fid,'%s & %s & %s \\\\\\cmidrule(lr){2-6}\\cmidrule(lr){7-7}\n','','\multicolumn{5}{c}{Factor betas}','FMB');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','','P1','P2','P3','P4','P5','$\lambda$');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',headerFile4{1},tsMimickDisp.beta(2,:),fmMimickDisp.lambda(1));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\\n','',tsMimickDisp.tstat(2,:),fmMimickDisp.tstat(1));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',headerFile4{2},tsMimickDisp.beta(3,:),fmMimickDisp.lambda(2));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\\n','',tsMimickDisp.tstat(3,:),fmMimickDisp.tstat(2));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s \\\\\n',headerFile4{3},tsMimickDisp.r2,'');
fprintf(fid,'\n');
fclose(fid);

figure;
scatter(fmMomDisp.fit.*12,fmMomDisp.mean.*12,[],'white');
axis([-3 6 -3 6]);
ref = refline(1,0);
ref.Color = 'k';
xlabel('Model implied returns (in %)');
ylabel('Average realized returns (in %)');
box on
labels = num2str((1:size(fmMomDisp.fit,1))','P%d');
text(fmMomDisp.fit.*12, fmMomDisp.mean.*12, labels,'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color',colorBrewer(1));
set(gcf, 'PaperUnits', 'Centimeters','PaperSize',[21 13],'PaperPosition',[0.5 0.5 20 12]);
print(gcf,'-depsc','output/figure3.eps');


%% Cross-sectional asset pricing: Carry, VOL, CROWD, and SKEW
% ########################################################################### %
%{
    In this part, we conduct simple cross-sectional asset pricing tests. 
    We consider linear GMM estimation following Hansen (1982) and supplement 
    with Fama-MacBeth two-pass cross-sectional regressions.   
%}
% ########################################################################### %

disp('<*> Cross-sectional asset pricing: Carry, VOL, CROWD, and SKEW');


testAssets      = momPortfoliosBA(3:end,:);
gmmDisp         = linearGMM(testAssets,[dolFactor(3:end,1) dispFXinno(3:end,1)]);
gmmHML          = linearGMM(testAssets,[dolFactor(3:end,1) hmlFactor(3:end,1)]);
gmmVol          = linearGMM(testAssets,[dolFactor(3:end,1) volFXinno(3:end,1)]);
gmmHMLDisp      = linearGMM(testAssets,[dolFactor(3:end,1) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
gmmVolDisp      = linearGMM(testAssets,[dolFactor(3:end,1) volFXinno(3:end,1) dispFXinno(3:end,1)]);

testAssets      = [momPortfoliosBA(3:end,:) carPortfoliosBA(3:end,:)];
gmmDispC        = linearGMM(testAssets,[dolFactor(3:end,1) dispFXinno(3:end,1)]);
gmmHMLC         = linearGMM(testAssets,[dolFactor(3:end,1) hmlFactor(3:end,1)]);
gmmVolC         = linearGMM(testAssets,[dolFactor(3:end,1) volFXinno(3:end,1)]);
gmmHMLDispC     = linearGMM(testAssets,[dolFactor(3:end,1) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
gmmVolDispC     = linearGMM(testAssets,[dolFactor(3:end,1) volFXinno(3:end,1) dispFXinno(3:end,1)]);


% [tsHMLDispC,fmHMLDispC]   = famaMacBeth(testAssets,[dolFactor(3:end,1) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
% [tsVolDispC,fmVolDispC]   = famaMacBeth(testAssets,[dolFactor(3:end,1) volFXinno(3:end,1) dispFXinno(3:end,1)]);


headerFile = {'b$_{DOL}$','b$_{CAR}$','b$_{\Delta VOL}$','b$_{\Delta CSV}$',...
    '$\lambda_{DOL}$','$\lambda_{CAR}$','$\lambda_{\Delta VOL}$',...
    '$\lambda_{\Delta CSV}$','R$^{2}$','HJ'};
fid = fopen('output/table4.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n',headerFile{:,:});
fprintf(fid,'%s \\\\\\cmidrule(lr){1-10}\n','\multicolumn{10}{c}{Panel A: Five momentum portfolios}');
fprintf(fid,'%.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmmDisp.b(1),'','',gmmDisp.b(2),gmmDisp.lambda(1),'','',gmmDisp.lambda(2),gmmDisp.r2,gmmDisp.hj);
fprintf(fid,'(%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmDisp.btstat(1),'','',gmmDisp.btstat(2),gmmDisp.tlambda(1),'','',gmmDisp.tlambda(2),'',gmmDisp.hjp);
fprintf(fid,'%.2f & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmmHML.b,'','',gmmHML.lambda,'','',gmmHML.r2,gmmHML.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmmHML.btstat,'','',gmmHML.tlambda,'','','',gmmHML.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %.2f & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmmVol.b(1),'',gmmVol.b(2),'',gmmVol.lambda(1),'',gmmVol.lambda(2),'',gmmVol.r2,gmmVol.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & (%.2f) & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmmVol.btstat(1),'',gmmVol.btstat(2),'',gmmVol.tlambda(1),'',gmmVol.tlambda(2),'','',gmmVol.hjp);
fprintf(fid,'%.2f & %.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmmHMLDisp.b(1:2),'',gmmHMLDisp.b(3),gmmHMLDisp.lambda(1:2),'',gmmHMLDisp.lambda(3),gmmHMLDisp.r2,gmmHMLDisp.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmHMLDisp.btstat(1:2),'',gmmHMLDisp.btstat(3),gmmHMLDisp.tlambda(1:2),'',gmmHMLDisp.tlambda(3),'',gmmHMLDisp.hjp);
fprintf(fid,'%.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmmVolDisp.b(1),'',gmmVolDisp.b(2:3),gmmVolDisp.lambda(1),'',gmmVolDisp.lambda(2:3),gmmVolDisp.r2,gmmVolDisp.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\\cmidrule(lr){1-10}\n',gmmVolDisp.btstat(1),'',gmmVolDisp.btstat(2:3),gmmVolDisp.tlambda(1),'',gmmVolDisp.tlambda(2:3),'',gmmVolDisp.hjp);
fprintf(fid,'%s \\\\\\cmidrule(lr){1-10}\n','\multicolumn{10}{c}{Panel B: Five momentum and five carry portfolios}');
fprintf(fid,'%.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmmDispC.b(1),'','',gmmDispC.b(2),gmmDispC.lambda(1),'','',gmmDispC.lambda(2),gmmDispC.r2,gmmDispC.hj);
fprintf(fid,'(%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmDispC.btstat(1),'','',gmmDispC.btstat(2),gmmDispC.tlambda(1),'','',gmmDispC.tlambda(2),'',gmmDispC.hjp);
fprintf(fid,'%.2f & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmmHMLC.b,'','',gmmHMLC.lambda,'','',gmmHMLC.r2,gmmHMLC.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmmHMLC.btstat,'','',gmmHMLC.tlambda,'','','',gmmHMLC.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %.2f & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmmVolC.b(1),'',gmmVolC.b(2),'',gmmVolC.lambda(1),'',gmmVolC.lambda(2),'',gmmVolC.r2,gmmVolC.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & (%.2f) & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmmVolC.btstat(1),'',gmmVolC.btstat(2),'',gmmVolC.tlambda(1),'',gmmVolC.tlambda(2),'','',gmmVolC.hjp);
fprintf(fid,'%.2f & %.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmmHMLDispC.b(1:2),'',gmmHMLDispC.b(3),gmmHMLDispC.lambda(1:2),'',gmmHMLDispC.lambda(3),gmmHMLDispC.r2,gmmHMLDispC.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmHMLDispC.btstat(1:2),'',gmmHMLDispC.btstat(3),gmmHMLDispC.tlambda(1:2),'',gmmHMLDispC.tlambda(3),'',gmmHMLDispC.hjp);
fprintf(fid,'%.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmmVolDispC.b(1),'',gmmVolDispC.b(2:3),gmmVolDispC.lambda(1),'',gmmVolDispC.lambda(2:3),gmmVolDispC.r2,gmmVolDispC.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\n',gmmVolDispC.btstat(1),'',gmmVolDispC.btstat(2:3),gmmVolDispC.tlambda(1),'',gmmVolDispC.tlambda(2:3),'',gmmVolDispC.hjp);
fprintf(fid,'\n');
fclose(fid);



%% Cross-sectional asset pricing: Ouptut gap portfolios
% ########################################################################### %
%{
    In this part, we conduct simple cross-sectional asset pricing tests. 
    We consider linear GMM estimation following Hansen (1982) and supplement 
    with Fama-MacBeth two-pass cross-sectional regressions.   
%}
% ########################################################################### %

disp('<*> Cross-sectional asset pricing: Ouptut gap portfolios');


testAssets  = [momPortfoliosBA(3:end,:) gapPortfoliosBA(3:end,:)];
gmm7        = linearGMM(testAssets,[dolFactor(3:end,:) dispFXinno(3:end,1)]);
gmm8        = linearGMM(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1)]);
gmm9        = linearGMM(testAssets,[dolFactor(3:end,:) hmlFactor(3:end,1)]);
gmm10       = linearGMM(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1) dispFXinno(3:end,1)]);
gmm11       = linearGMM(testAssets,[dolFactor(3:end,:) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
gmm12       = linearGMM(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
[tsDisp,fmDisp]     = famaMacBeth(testAssets,[dolFactor(3:end,1) gapFactor(3:end,1) dispFXinno(3:end,1)]);

testAssets  = [momPortfoliosBA(3:end,:) carPortfoliosBA(3:end,:) gapPortfoliosBA(3:end,:)];
gmm13       = linearGMM(testAssets,[dolFactor(3:end,:) dispFXinno(3:end,1)]);
gmm14       = linearGMM(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1)]);
gmm15       = linearGMM(testAssets,[dolFactor(3:end,:) hmlFactor(3:end,1)]);
gmm16       = linearGMM(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1) dispFXinno(3:end,1)]);
gmm17       = linearGMM(testAssets,[dolFactor(3:end,:) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
gmm18       = linearGMM(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1) hmlFactor(3:end,1) dispFXinno(3:end,1)]);
[tsAll,fmAll]       = famaMacBeth(testAssets,[dolFactor(3:end,:) gapFactor(3:end,1) hmlFactor(3:end,1) dispFXinno(3:end,1)]);


headerFile = {'b$_{DOL}$','b$_{GAP}$','b$_{CAR}$','b$_{\Delta CSV}$',...
    '$\lambda_{DOL}$','$\lambda_{GAP}$','$\lambda_{CAR}$',...
    '$\lambda_{\Delta CSV}$','R$^{2}$','HJ'};


fid = fopen('output/table5.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n',headerFile{:,:});

fprintf(fid,'%s \\\\\\cmidrule(lr){1-10}\n','\multicolumn{10}{c}{Panel A: Five momentum and five output gap portfolios}');
fprintf(fid,'%.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmm7.b(1),'','',gmm7.b(2),gmm7.lambda(1),'','',gmm7.lambda(2),gmm7.r2,gmm7.hj);
fprintf(fid,'(%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmm7.btstat(1),'','',gmm7.btstat(2),gmm7.tlambda(1),'','',gmm7.tlambda(2),'',gmm7.hjp);
fprintf(fid,'%.2f & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmm8.b,'','',gmm8.lambda,'','',gmm8.r2,gmm8.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmm8.btstat,'','',gmm8.tlambda,'','','',gmm8.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %.2f & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmm9.b(1),'',gmm9.b(2),'',gmm9.lambda(1),'',gmm9.lambda(2),'',gmm9.r2,gmm9.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & (%.2f) & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmm9.btstat(1),'',gmm9.btstat(2),'',gmm9.tlambda(1),'',gmm9.tlambda(2),'','',gmm9.hjp);
fprintf(fid,'%.2f & %.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmm10.b(1:2),'',gmm10.b(3),gmm10.lambda(1:2),'',gmm10.lambda(3),gmm10.r2,gmm10.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmm10.btstat(1:2),'',gmm10.btstat(3),gmm10.tlambda(1:2),'',gmm10.tlambda(3),'',gmm10.hjp);
fprintf(fid,'%.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmm11.b(1),'',gmm11.b(2:3),gmm11.lambda(1),'',gmm11.lambda(2:3),gmm11.r2,gmm11.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\n',gmm11.btstat(1),'',gmm11.btstat(2:3),gmm11.tlambda(1),'',gmm11.tlambda(2:3),'',gmm11.hjp);
fprintf(fid,'%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',gmm12.b,gmm12.lambda,gmm12.r2,gmm12.hj);
fprintf(fid,'(%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & %s & [%.2f] \\\\\\cmidrule(lr){1-10}\n',gmm12.btstat,gmm12.tlambda,'',gmm12.hjp);
fprintf(fid,'%s \\\\\\cmidrule(lr){1-10}\n','\multicolumn{10}{c}{Panel B: Five momentum, five carry, and five output gap portfolios}');
fprintf(fid,'%.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmm13.b(1),'','',gmm13.b(2),gmm13.lambda(1),'','',gmm13.lambda(2),gmm13.r2,gmm13.hj);
fprintf(fid,'(%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmm13.btstat(1),'','',gmm13.btstat(2),gmm13.tlambda(1),'','',gmm13.tlambda(2),'',gmm13.hjp);
fprintf(fid,'%.2f & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmm14.b,'','',gmm14.lambda,'','',gmm14.r2,gmm14.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmm14.btstat,'','',gmm14.tlambda,'','','',gmm14.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %.2f & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmm15.b(1),'',gmm15.b(2),'',gmm15.lambda(1),'',gmm15.lambda(2),'',gmm15.r2,gmm15.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & (%.2f) & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmm15.btstat(1),'',gmm15.btstat(2),'',gmm15.tlambda(1),'',gmm15.tlambda(2),'','',gmm15.hjp);
fprintf(fid,'%.2f & %.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmm16.b(1:2),'',gmm16.b(3),gmm16.lambda(1:2),'',gmm16.lambda(3),gmm16.r2,gmm16.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmm16.btstat(1:2),'',gmm16.btstat(3),gmm16.tlambda(1:2),'',gmm16.tlambda(3),'',gmm16.hjp);
fprintf(fid,'%.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmm17.b(1),'',gmm17.b(2:3),gmm17.lambda(1),'',gmm17.lambda(2:3),gmm17.r2,gmm17.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\n',gmm17.btstat(1),'',gmm17.btstat(2:3),gmm17.tlambda(1),'',gmm17.tlambda(2:3),'',gmm17.hjp);
fprintf(fid,'%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n',gmm18.b,gmm18.lambda,gmm18.r2,gmm18.hj);
fprintf(fid,'(%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & %s & [%.2f] \\\\\n',gmm18.btstat,gmm18.tlambda,'',gmm18.hjp);
fprintf(fid,'\n');
fclose(fid);






%% Cross-sectional asset pricing: Liquidity and volatility factors
% ########################################################################### %
%{
    In this part, we conduct simple cross-sectional asset pricing tests. 
    We consider linear GMM estimation following Hansen (1982) and supplement 
    with Fama-MacBeth two-pass cross-sectional regressions.   
%}
% ########################################################################### %

disp('<*> Cross-sectional asset pricing: Liquidity and volatility factors');

% Computing TED sprad innovations
tedSpread1      = tedSpread(27:end,:);
tedARcoef       = [ones(size(tedSpread1,1)-1,1) tedSpread1(1:end-1,1)]\tedSpread1(2:end,1);
tedFXinno       = tedSpread1(2:end,1) - [ones(size(tedSpread1,1)-1,1) tedSpread1(1:end-1,1)]*tedARcoef;
tedFXinno       = [NaN(27,1);tedFXinno];

% Computing innovations to NVIX
vix1            = vix(75:end,:);
vixARcoef       = [ones(size(vix1,1)-1,1) vix1(1:end-1,1)]\vix1(2:end,1);
vixFXinno       = vix1(2:end,1) - [ones(size(vix1,1)-1,1) vix1(1:end-1,1)]*vixARcoef;
vixFXinno       = [NaN(75,1); vixFXinno;];

% Computing innovations to NVIX
nvix1           = nvix(1:end-9,:);
nvixARcoef      = [ones(size(nvix1,1)-1,1) nvix1(1:end-1,1)]\nvix1(2:end,1);
nvixFXinno      = nvix1(2:end,1) - [ones(size(nvix1,1)-1,1) nvix1(1:end-1,1)]*nvixARcoef;
nvixFXinno      = [NaN(1,1); nvixFXinno; NaN(9,1)];

% Computing innovations to FX illiquidity
fxIlliquid1     = fxIlliquid(87:end-2,:);
illiqARcoef     = [ones(size(fxIlliquid1,1)-1,1) fxIlliquid1(1:end-1,1)]\fxIlliquid1(2:end,1);
illiqFXinno     = fxIlliquid1(2:end,1) - [ones(size(fxIlliquid1,1)-1,1) fxIlliquid1(1:end-1,1)]*illiqARcoef;
illiqFXinno     = [NaN(87,1); illiqFXinno; NaN(2,1)];

% GMM asset pricing horse races: Dispersion innovations
testAssets      = momPortfoliosBA(3:end,:);
gmmBASU         = linearGMM(testAssets,[dolFactor(3:end,1) basFXinno(3:end,1)]);
gmmBASDispU     = linearGMM(testAssets,[dolFactor(3:end,1) basFXinno(3:end,1) dispFXinno(3:end,1)]);
gmmTEDU         = linearGMM(testAssets(26:end,:),[dolFactor(28:end,1) tedFXinno(28:end,1)]);
gmmTEDDispU     = linearGMM(testAssets(26:end,:),[dolFactor(28:end,1) tedFXinno(28:end,1) dispFXinno(28:end,1)]);
gmmNVIXU        = linearGMM(testAssets(1:end-9,:),[dolFactor(3:end-9,1) nvixFXinno(3:end-9,1)]);
gmmNVIXDispU    = linearGMM(testAssets(1:end-9,:),[dolFactor(3:end-9,1) nvixFXinno(3:end-9,1) dispFXinno(3:end-9,1)]);
gmmILLU         = linearGMM(testAssets(86:end-2,:),[dolFactor(88:end-2,1) illiqFXinno(88:end-2,1)]);
gmmILLDispU     = linearGMM(testAssets(86:end-2,:),[dolFactor(88:end-2,1) illiqFXinno(88:end-2,1) dispFXinno(88:end-2,1)]);
gmmVIXU         = linearGMM(testAssets(74:end,:),[dolFactor(76:end,1) vixFXinno(76:end,1)]);
gmmVIXDispU     = linearGMM(testAssets(74:end,:),[dolFactor(76:end,1) vixFXinno(76:end,1) dispFXinno(76:end,1)]);

% Creating first table with results for cross-sectional horse races
headerFile = {'b$_{DOL}$','b$_{\Delta BAS}$','b$_{\Delta LFI}$','b$_{\Delta TED}$','b$_{\Delta NVIX}$','b$_{\Delta VIX}$','b$_{\Delta CSV}$',...
    '$\lambda_{DOL}$','$\lambda_{\Delta BAS}$','$\lambda_{\Delta LFI}$','$\lambda_{\Delta TED}$','$\lambda_{\Delta NVIX}$','$\lambda_{\Delta VIX}$','$\lambda_{\Delta CSV}$','R$^{2}$','HJ'};

fid = fopen('output/table6.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n',headerFile{:,:});
fprintf(fid,'%.2f & %.2f & %s & %s & %s & %s & %s & %.2f & %.2f & %s & %s & %s & %s & %s & %.2f & %.2f \\\\\n',gmmBASU.b(1:2),'','','','','',gmmBASU.lambda(1:2),'','','','','',gmmBASU.r2,gmmBASU.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & %s & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & %s & %s & %s & [%.2f] \\\\\n',gmmBASU.btstat(1:2),'','','','','',gmmBASU.tlambda(1:2),'','','','','','',gmmBASU.hjp);
fprintf(fid,'%.2f & %.2f & %s & %s & %s & %s & %.2f & %.2f & %.2f & %s & %s & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmmBASDispU.b(1:2),'','','','',gmmBASDispU.b(3),gmmBASDispU.lambda(1:2),'','','','',gmmBASDispU.lambda(3),gmmBASDispU.r2,gmmBASDispU.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & %s & %s & (%.2f) & (%.2f) & (%.2f) & %s & %s & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmBASDispU.btstat(1:2),'','','','',gmmBASDispU.btstat(3),gmmBASDispU.tlambda(1:2),'','','','',gmmBASDispU.tlambda(3),'',gmmBASDispU.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %s & %s & %s & %.2f & %s & %.2f & %s & %s & %s & %s & %.2f & %.2f \\\\\n',gmmILLU.b(1),'',gmmILLU.b(2),'','','','',gmmILLU.lambda(1),'',gmmILLU.lambda(2),'','','','',gmmILLU.r2,gmmILLU.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & %s & %s & %s & (%.2f) & %s & (%.2f) & %s & %s & %s & %s & %s & [%.2f] \\\\\n',gmmILLU.btstat(1),'',gmmILLU.btstat(2),'','','','',gmmILLU.tlambda(1),'',gmmILLU.tlambda(2),'','','','','',gmmILLU.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %s & %s & %.2f & %.2f & %s & %.2f & %s & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmmILLDispU.b(1),'',gmmILLDispU.b(2),'','','',gmmILLDispU.b(3),gmmILLDispU.lambda(1),'',gmmILLDispU.lambda(2),'','','',gmmILLDispU.lambda(3),gmmILLDispU.r2,gmmILLDispU.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & %s & %s & (%.2f) & (%.2f) & %s & (%.2f) & %s & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmILLDispU.btstat(1),'',gmmILLDispU.btstat(2),'','','',gmmILLDispU.btstat(3),gmmILLDispU.tlambda(1),'',gmmILLDispU.tlambda(2),'','','',gmmILLDispU.tlambda(3),'',gmmILLDispU.hjp);
fprintf(fid,'%.2f & %s & %s & %.2f & %s & %s & %s & %.2f & %s & %s & %.2f & %s & %s & %s & %.2f & %.2f \\\\\n',gmmTEDU.b(1),'','',gmmTEDU.b(2),'','','',gmmTEDU.lambda(1),'','',gmmTEDU.lambda(2),'','','',gmmTEDU.r2,gmmTEDU.hj);
fprintf(fid,'(%.2f) & %s & %s & (%.2f) & %s & %s & %s & (%.2f) & %s & %s & (%.2f) & %s & %s & %s & %s & [%.2f] \\\\\n',gmmTEDU.btstat(1),'','',gmmTEDU.btstat(2),'','','',gmmTEDU.tlambda(1),'','',gmmTEDU.tlambda(2),'','','','',gmmTEDU.hjp);
fprintf(fid,'%.2f & %s & %s & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %s & %s & %.2f & %.2f & %.2f \\\\\n',gmmTEDDispU.b(1),'','',gmmTEDDispU.b(2),'','',gmmTEDDispU.b(3),gmmTEDDispU.lambda(1),'','',gmmTEDDispU.lambda(2),'','',gmmTEDDispU.lambda(3),gmmTEDDispU.r2,gmmTEDDispU.hj);
fprintf(fid,'(%.2f) & %s & %s & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & (%.2f) & %s & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmTEDDispU.btstat(1),'','',gmmTEDDispU.btstat(2),'','',gmmTEDDispU.btstat(3),gmmTEDDispU.tlambda(1),'','',gmmTEDDispU.tlambda(2),'','',gmmTEDDispU.tlambda(3),'',gmmTEDDispU.hjp);
fprintf(fid,'%.2f & %s & %s & %s & %.2f & %s & %s & %.2f & %s & %s & %s & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmmNVIXU.b(1),'','','',gmmNVIXU.b(2),'','',gmmNVIXU.lambda(1),'','','',gmmNVIXU.lambda(2),'','',gmmNVIXU.r2,gmmNVIXU.hj);
fprintf(fid,'(%.2f) & %s & %s & %s & (%.2f) & %s & %s & (%.2f) & %s & %s & %s & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmmNVIXU.btstat(1),'','','',gmmNVIXU.btstat(2),'','',gmmNVIXU.tlambda(1),'','','',gmmNVIXU.tlambda(2),'','','',gmmNVIXU.hjp);
fprintf(fid,'%.2f & %s & %s & %s & %.2f & %s & %.2f & %.2f & %s & %s & %s & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmmNVIXDispU.b(1),'','','',gmmNVIXDispU.b(2),'',gmmNVIXDispU.b(3),gmmNVIXDispU.lambda(1),'','','',gmmNVIXDispU.lambda(2),'',gmmNVIXDispU.lambda(3),gmmNVIXDispU.r2,gmmNVIXDispU.hj);
fprintf(fid,'(%.2f) & %s & %s & %s & (%.2f) & %s & (%.2f) & (%.2f) & %s & %s & %s & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmNVIXDispU.btstat(1),'','','',gmmNVIXDispU.btstat(2),'',gmmNVIXDispU.btstat(3),gmmNVIXDispU.tlambda(1),'','','',gmmNVIXDispU.tlambda(2),'',gmmNVIXDispU.tlambda(3),'',gmmNVIXDispU.hjp);


fprintf(fid,'%.2f & %s & %s & %s & %s & %.2f & %s & %.2f & %s & %s & %s & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmmVIXU.b(1),'','','','',gmmVIXU.b(2),'',gmmVIXU.lambda(1),'','','','',gmmVIXU.lambda(2),'',gmmVIXU.r2,gmmVIXU.hj);
fprintf(fid,'(%.2f) & %s & %s & %s & %s & (%.2f) & %s & (%.2f) & %s & %s & %s & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmmVIXU.btstat(1),'','','','',gmmVIXU.btstat(2),'',gmmVIXU.tlambda(1),'','','','',gmmVIXU.tlambda(2),'','',gmmVIXU.hjp);


fprintf(fid,'%.2f & %s & %s & %s & %s & %.2f & %.2f & %.2f & %s & %s & %s & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmmVIXDispU.b(1),'','','','',gmmVIXDispU.b(2),gmmVIXDispU.b(3),gmmVIXDispU.lambda(1),'','','','',gmmVIXDispU.lambda(2),gmmVIXDispU.lambda(3),gmmVIXDispU.r2,gmmVIXDispU.hj);
fprintf(fid,'(%.2f) & %s & %s & %s & %s & (%.2f) & (%.2f) & (%.2f) & %s & %s & %s & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\n',gmmVIXDispU.btstat(1),'','','','',gmmVIXDispU.btstat(2),gmmVIXDispU.btstat(3),gmmVIXDispU.tlambda(1),'','','','',gmmVIXDispU.tlambda(2),gmmVIXDispU.tlambda(3),'',gmmVIXDispU.hjp);

fprintf(fid,'\n');
fclose(fid);







%% Cross-sectional asset pricing: Carry, VOL, CROWD, and SKEW
% ########################################################################### %
%{
    In this part, we conduct simple cross-sectional asset pricing tests. 
    We consider linear GMM estimation following Hansen (1982) and supplement 
    with Fama-MacBeth two-pass cross-sectional regressions.   
%}
% ########################################################################### %

disp('<*> Cross-sectional asset pricing: Carry, VOL, CROWD, and SKEW');

% GMM asset pricing horse races: All countries, one-month momentum
testAssets      = momPortfoliosBA(3:end,:);
gmmCRWD         = linearGMM(testAssets,[dolFactor(3:end,1) crowdFXinno(3:end,1)]);
gmmSKEW         = linearGMM(testAssets,[dolFactor(3:end,1) skewFXinno(3:end,1)]);
gmmCRWDDISP     = linearGMM(testAssets,[dolFactor(3:end,1) crowdFXinno(3:end,1) dispFXinno(3:end,1)]);
gmmSKEWDISP     = linearGMM(testAssets,[dolFactor(3:end,1) skewFXinno(3:end,1) dispFXinno(3:end,1)]);
gmmSKEWCRWD     = linearGMM(testAssets,[dolFactor(3:end,1) crowdFXinno(3:end,1) skewFXinno(3:end,1) dispFXinno(3:end,1)]);

testAssets       = [momPortfoliosBA(3:end,:) carPortfoliosBA(3:end,:) gapPortfoliosBA(3:end,:)];
gmmCRWD1         = linearGMM(testAssets,[dolFactor(3:end,1) crowdFXinno(3:end,1)]);
gmmSKEW1         = linearGMM(testAssets,[dolFactor(3:end,1) skewFXinno(3:end,1)]);
gmmCRWDDISP1     = linearGMM(testAssets,[dolFactor(3:end,1) crowdFXinno(3:end,1) dispFXinno(3:end,1)]);
gmmSKEWDISP1     = linearGMM(testAssets,[dolFactor(3:end,1) skewFXinno(3:end,1) dispFXinno(3:end,1)]);
gmmSKEWCRWD1     = linearGMM(testAssets,[dolFactor(3:end,1) crowdFXinno(3:end,1) skewFXinno(3:end,1) dispFXinno(3:end,1)]);


headerFile = {'b$_{DOL}$','b$_{\Delta CRWD}$','b$_{\Delta SKEW}$','b$_{\Delta CSV}$',...
    '$\lambda_{DOL}$','$\lambda_{\Delta CRWD}$','$\lambda_{\Delta SKEW}$',...
    '$\lambda_{\Delta CSV}$','R$^{2}$','HJ'};
fid = fopen('output/table7.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n',headerFile{:,:});
fprintf(fid,'%s \\\\\\midrule\n','\multicolumn{10}{c}{Panel A: Five momentum portfolios}');
fprintf(fid,'%.2f & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmmCRWD.b,'','',gmmCRWD.lambda,'','',gmmCRWD.r2,gmmCRWD.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmmCRWD.btstat,'','',gmmCRWD.tlambda,'','','',gmmCRWD.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %.2f & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmmSKEW.b(1),'',gmmSKEW.b(2),'',gmmSKEW.lambda(1),'',gmmSKEW.lambda(2),'',gmmSKEW.r2,gmmSKEW.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & (%.2f) & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmmSKEW.btstat(1),'',gmmSKEW.btstat(2),'',gmmSKEW.tlambda(1),'',gmmSKEW.tlambda(2),'','',gmmSKEW.hjp);
fprintf(fid,'%.2f & %.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmmCRWDDISP.b(1:2),'',gmmCRWDDISP.b(3),gmmCRWDDISP.lambda(1:2),'',gmmCRWDDISP.lambda(3),gmmCRWDDISP.r2,gmmCRWDDISP.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmCRWDDISP.btstat(1:2),'',gmmCRWDDISP.btstat(3),gmmCRWDDISP.tlambda(1:2),'',gmmCRWDDISP.tlambda(3),'',gmmCRWDDISP.hjp);
fprintf(fid,'%.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmmSKEWDISP.b(1),'',gmmSKEWDISP.b(2:3),gmmSKEWDISP.lambda(1),'',gmmSKEWDISP.lambda(2:3),gmmSKEWDISP.r2,gmmSKEWDISP.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\\midrule\n',gmmSKEWDISP.btstat(1),'',gmmSKEWDISP.btstat(2:3),gmmSKEWDISP.tlambda(1),'',gmmSKEWDISP.tlambda(2:3),'',gmmSKEWDISP.hjp);


fprintf(fid,'%s \\\\\\midrule\n','\multicolumn{10}{c}{Panel B: Five momentum, five carry, and five output gap portfolios}');
fprintf(fid,'%.2f & %.2f & %s & %s & %.2f & %.2f & %s & %s & %.2f & %.2f \\\\\n',gmmCRWD1.b,'','',gmmCRWD1.lambda,'','',gmmCRWD1.r2,gmmCRWD1.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & %s & (%.2f) & (%.2f) & %s & %s & %s & [%.2f] \\\\\n',gmmCRWD1.btstat,'','',gmmCRWD1.tlambda,'','','',gmmCRWD1.hjp);
fprintf(fid,'%.2f & %s & %.2f & %s & %.2f & %s & %.2f & %s & %.2f & %.2f \\\\\n',gmmSKEW1.b(1),'',gmmSKEW1.b(2),'',gmmSKEW1.lambda(1),'',gmmSKEW1.lambda(2),'',gmmSKEW1.r2,gmmSKEW1.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & %s & (%.2f) & %s & (%.2f) & %s & %s & [%.2f] \\\\\n',gmmSKEW1.btstat(1),'',gmmSKEW1.btstat(2),'',gmmSKEW1.tlambda(1),'',gmmSKEW1.tlambda(2),'','',gmmSKEW1.hjp);
fprintf(fid,'%.2f & %.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',gmmCRWDDISP1.b(1:2),'',gmmCRWDDISP1.b(3),gmmCRWDDISP1.lambda(1:2),'',gmmCRWDDISP1.lambda(3),gmmCRWDDISP1.r2,gmmCRWDDISP1.hj);
fprintf(fid,'(%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & %s & [%.2f] \\\\\n',gmmCRWDDISP1.btstat(1:2),'',gmmCRWDDISP1.btstat(3),gmmCRWDDISP1.tlambda(1:2),'',gmmCRWDDISP1.tlambda(3),'',gmmCRWDDISP1.hjp);
fprintf(fid,'%.2f & %s & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',gmmSKEWDISP1.b(1),'',gmmSKEWDISP1.b(2:3),gmmSKEWDISP1.lambda(1),'',gmmSKEWDISP1.lambda(2:3),gmmSKEWDISP1.r2,gmmSKEWDISP1.hj);
fprintf(fid,'(%.2f) & %s & (%.2f) & (%.2f) & (%.2f) & %s & (%.2f) & (%.2f) & %s & [%.2f] \\\\\n',gmmSKEWDISP1.btstat(1),'',gmmSKEWDISP1.btstat(2:3),gmmSKEWDISP1.tlambda(1),'',gmmSKEWDISP1.tlambda(2:3),'',gmmSKEWDISP1.hjp);



fprintf(fid,'\n');
fclose(fid);

















fundCond1 = fundCond1;

tmp         = momFactorBA(2:end,1).*12;
resMom1     = linRegNWA(tmp(fundCond1(1:end-1,1)==1,:),ones(sum(fundCond1(1:end-1,1)==1),1),0);
resMom2     = linRegNWA(tmp(fundCond1(1:end-1,1)==2,:),ones(sum(fundCond1(1:end-1,1)==2),1),0);
resMom3     = linRegNWA(tmp(fundCond1(1:end-1,1)==3,:),ones(sum(fundCond1(1:end-1,1)==3),1),0);
resMom4     = linRegNWA(tmp(fundCond1(1:end-1,1)==4,:),ones(sum(fundCond1(1:end-1,1)==4),1),0);

tmp         = momPortfoliosBA(2:end,1).*12;
resLos1     = linRegNWA(tmp(fundCond1(1:end-1,1)==1,:),ones(sum(fundCond1(1:end-1,1)==1),1),0);
resLos2     = linRegNWA(tmp(fundCond1(1:end-1,1)==2,:),ones(sum(fundCond1(1:end-1,1)==2),1),0);
resLos3     = linRegNWA(tmp(fundCond1(1:end-1,1)==3,:),ones(sum(fundCond1(1:end-1,1)==3),1),0);
resLos4     = linRegNWA(tmp(fundCond1(1:end-1,1)==4,:),ones(sum(fundCond1(1:end-1,1)==4),1),0);

tmp         = momPortfoliosBA(2:end,5).*12;
resWin1     = linRegNWA(tmp(fundCond1(1:end-1,1)==1,:),ones(sum(fundCond1(1:end-1,1)==1),1),0);
resWin2     = linRegNWA(tmp(fundCond1(1:end-1,1)==2,:),ones(sum(fundCond1(1:end-1,1)==2),1),0);
resWin3     = linRegNWA(tmp(fundCond1(1:end-1,1)==3,:),ones(sum(fundCond1(1:end-1,1)==3),1),0);
resWin4     = linRegNWA(tmp(fundCond1(1:end-1,1)==4,:),ones(sum(fundCond1(1:end-1,1)==4),1),0);


tmp         = dispFXinno(3:end,1);
resDisp1    = linRegNWA(tmp(fundCond1(2:end-1,1)==1,:),ones(sum(fundCond1(2:end-1,1)==1),1),0);
resDisp2    = linRegNWA(tmp(fundCond1(2:end-1,1)==2,:),ones(sum(fundCond1(2:end-1,1)==2),1),0);
resDisp3    = linRegNWA(tmp(fundCond1(2:end-1,1)==3,:),ones(sum(fundCond1(2:end-1,1)==3),1),0);
resDisp4    = linRegNWA(tmp(fundCond1(2:end-1,1)==4,:),ones(sum(fundCond1(2:end-1,1)==4),1),0);


figure;
subplot(2,2,1);
hold on
data  = diag([resDisp1.b(1) resDisp2.b(1) resDisp3.b(1) resDisp4.b(1)]);
data(data==0) = NaN;
b1 = bar(data(1,:));
b2 = bar(data(2,:));
b3 = bar(data(3,:));
b4 = bar(data(4,:));
% e1 = errorbar(1:4,[resDisp1.b(1) resDisp2.b(1) resDisp3.b(1) resDisp4.b(1)],[resDisp1.se(1) resDisp2.se(1) resDisp3.se(1) resDisp4.se(1)].*2,'.');
hold off
b1.FaceColor = colorBrewer(1);
b1.EdgeColor = colorBrewer(1);
b2.FaceColor = colorBrewer(2);
b2.EdgeColor = colorBrewer(2);
b3.FaceColor = colorBrewer(3);
b3.EdgeColor = colorBrewer(3);
b4.FaceColor = colorBrewer(4);
b4.EdgeColor = colorBrewer(4);
% e1.Color = 'k';
xticks(1:4);
xticklabels({'EE','ER','RR','RE'});
ylabel('Average innovations');
xlabel('Monetary conditions');
title('Dispersion innovations','FontWeight','normal');
box on
axis([0 5 -0.1 0.15]);
subplot(2,2,2);
hold on
data  = diag([resMom1.b(1) resMom2.b(1) resMom3.b(1) resMom4.b(1)]);
data(data==0) = NaN;
b1 = bar(data(1,:));
b2 = bar(data(2,:));
b3 = bar(data(3,:));
b4 = bar(data(4,:));
% e1 = errorbar(1:4,[resMom1.b(1) resMom2.b(1) resMom3.b(1) resMom4.b(1)],[resMom1.se(1) resMom2.se(1) resMom3.se(1) resMom4.se(1)].*2,'.');
hold off
b1.FaceColor = colorBrewer(1);
b1.EdgeColor = colorBrewer(1);
b2.FaceColor = colorBrewer(2);
b2.EdgeColor = colorBrewer(2);
b3.FaceColor = colorBrewer(3);
b3.EdgeColor = colorBrewer(3);
b4.FaceColor = colorBrewer(4);
b4.EdgeColor = colorBrewer(4);
% e1.Color = 'k';
xticks(1:4);
xticklabels({'EE','ER','RR','RE'});
ylabel('Average returns (in % p.a.)');
xlabel('Monetary conditions');
title('Momentum strategy','FontWeight','normal');
box on
axis([0 5 -5 10]);
subplot(2,2,3);
hold on
data  = diag([resLos1.b(1) resLos2.b(1) resLos3.b(1) resLos4.b(1)]);
data(data==0) = NaN;
b1 = bar(data(1,:));
b2 = bar(data(2,:));
b3 = bar(data(3,:));
b4 = bar(data(4,:));
% e1 = errorbar(1:4,[resLos1.b(1) resLos2.b(1) resLos3.b(1) resLos4.b(1)],[resLos1.se(1) resLos2.se(1) resLos3.se(1) resLos4.se(1)].*2,'.');
hold off
b1.FaceColor = colorBrewer(1);
b1.EdgeColor = colorBrewer(1);
b2.FaceColor = colorBrewer(2);
b2.EdgeColor = colorBrewer(2);
b3.FaceColor = colorBrewer(3);
b3.EdgeColor = colorBrewer(3);
b4.FaceColor = colorBrewer(4);
b4.EdgeColor = colorBrewer(4);
% e1.Color = 'k';
xticks(1:4);
xticklabels({'EE','ER','RR','RE'});
ylabel('Average returns (in % p.a.)');
xlabel('Monetary conditions');
title('Loser currencies','FontWeight','normal');
box on
axis([0 5 -5 10]);
subplot(2,2,4);
hold on
data  = diag([resWin1.b(1) resWin2.b(1) resWin3.b(1) resWin4.b(1)]);
data(data==0) = NaN;
b1 = bar(data(1,:));
b2 = bar(data(2,:));
b3 = bar(data(3,:));
b4 = bar(data(4,:));
% e1 = errorbar(1:4,[resWin1.b(1) resWin2.b(1) resWin3.b(1) resWin4.b(1)],[resWin1.se(1) resWin2.se(1) resWin3.se(1) resWin4.se(1)].*2,'.');
hold off
b1.FaceColor = colorBrewer(1);
b1.EdgeColor = colorBrewer(1);
b2.FaceColor = colorBrewer(2);
b2.EdgeColor = colorBrewer(2);
b3.FaceColor = colorBrewer(3);
b3.EdgeColor = colorBrewer(3);
b4.FaceColor = colorBrewer(4);
b4.EdgeColor = colorBrewer(4);
% e1.Color = 'k';
xticks(1:4);
xticklabels({'EE','ER','RR','RE'});
ylabel('Average returns (in % p.a.)');
xlabel('Monetary conditions');
box on
axis([0 5 -5 10]);
title('Winner currencies','FontWeight','normal');
set(gcf, 'PaperUnits', 'Centimeters','PaperSize',[21 19],'PaperPosition',[0.5 0.5 20 18]);
print(gcf,'-depsc','output/figure4.eps');







% ########################################################################### %
%% Computing code run time
% ########################################################################### %

tEnd = toc(tStart);
fprintf('<*> Runtime: %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
disp('<*> Routine Completed');

% ########################################################################### %
% [EOS]
% ########################################################################### %