%%
clear all
close all

%%
run(['..' filesep 'localdef_WIM_HB'])
addpath('plots')
addpath(path_ftrip) 
ft_defaults


%% [1] Gather data in pair-wise matrices

load('WIM_HB_HEP.mat')

ar = cell2mat(allRes);

data_MS_ERP=[];
group_MS_ERP=[];
data_MS_ECG=[];
times=t;
for nGr=unique(ar(:,32))'
    for nS=1:size(allERP,1)
        temp_dat = mean(allERP{nS}(:,:,allRes{nS}(:,32)==nGr),3);
        data_MS_ERP=cat(3,data_MS_ERP,temp_dat(1:end-1,times>=-0.3));
        data_MS_ECG=cat(1,data_MS_ECG,temp_dat(end,times>=-0.3));
        group_MS_ERP=cat(1,group_MS_ERP,[nGr nS]);
    end
end

data_VG_ERP=[];
group_VG_ERP=[];
data_VG_ECG=[];
for nGr=unique(ar(~isnan(ar(:,38)),38))'
    for nS=1:size(allERP,1)
        temp_dat = mean(allERP{nS}(:,:,allRes{nS}(:,38)==nGr),3);
        data_VG_ERP=cat(3,data_VG_ERP,temp_dat(1:end-1,times>=-0.3));
        data_VG_ECG=cat(1,data_VG_ECG,temp_dat(end,times>=-0.3));
        group_VG_ERP=cat(1,group_VG_ERP,[nGr nS]);
    end
end
times=times(times>=-0.3);

data_MS_GFP=[];
for nGr=1:size(group_MS_ERP,1)/size(allERP,1)
    temp_erp=squeeze(rms(data_MS_ERP(:,:,group_MS_ERP(:,1)==nGr), 1))';
    data_MS_GFP=cat(1,data_MS_GFP,temp_erp);
end

data_VG_GFP=[];
for nGr=1:size(group_VG_ERP,1)/size(allERP,1)
    temp_erp=squeeze(rms(data_VG_ERP(:,:,group_VG_ERP(:,1)==nGr), 1))';
    data_VG_GFP=cat(1,data_VG_GFP,temp_erp);
end

%%
totperm=250;

    [realpos_VG_GFP, realneg_VG_GFP, rd_VG_GFP, rdpv_VG_GFP, pd_VG_GFP, pdpv_VG_GFP, rdaov_VG_GFP, rdaovpv_VG_GFP, pdaov_VG_GFP, pdaovpv_VG_GFP]=get_cluster_permutation_lme(data_VG_GFP,group_VG_ERP,{'data','vg','subj'},[0 0 0],'data~1+vg+(1|subj)',0.05,0.05,totperm,times);
    [realpos_MS_GFP, realneg_MS_GFP, rd_MS_GFP, rdpv_MS_GFP, pd_MS_GFP, pdpv_MS_GFP, rdaov_MS_GFP, rdaovpv_MS_GFP, pdaov_MS_GFP, pdaovpv_MS_GFP]=get_cluster_permutation_lme(data_MS_GFP,group_MS_ERP,{'data','ms','subj'},[0 1 0],'data~1+ms+(1|subj)',0.05,0.05,totperm,times);
    [realpos_VG_ECG, realneg_VG_ECG, rd_VG_ECG, rdpv_VG_ECG, pd_VG_ECG, pdpv_VG_ECG, rdaov_VG_ECG, rdaovpv_VG_ECG, pdaov_VG_ECG, pdaovpv_VG_ECG]=get_cluster_permutation_lme(data_VG_ECG,group_VG_ERP,{'data','vg','subj'},[0 0 0],'data~1+vg+(1|subj)',0.05,0.05,totperm,times);
    [realpos_MS_ECG, realneg_MS_ECG_MS_ECG, rd_MS_ECG, rdpv_MS_ECG, pd_MS_ECG, pdpv_MS_ECG, rdaov_MS_ECG, rdaovpv_MS_ECG, pdaov_MS_ECG, pdaovpv_MS_ECG]=get_cluster_permutation_lme(data_MS_ECG,group_MS_ERP,{'data','ms','subj'},[0 1 0],'data~1+ms+(1|subj)',0.05,0.05,totperm,times);
    save([pwd filesep 'clusterperm_res_ERP_GFP_smooth'],'realpos_VG_GFP','realneg_VG_GFP','rd_VG_GFP','rdpv_VG_GFP','pd_VG_GFP','pdpv_VG_GFP','rdaov_VG_GFP','rdaovpv_VG_GFP','pdaov_VG_GFP','pdaovpv_VG_GFP',...
        'realpos_MS_GFP','realneg_MS_GFP','rd_MS_GFP','rdpv_MS_GFP','pd_MS_GFP','pdpv_MS_GFP','rdaov_MS_GFP','rdaovpv_MS_GFP','pdaov_MS_GFP','pdaovpv_MS_GFP',...
        'realpos_VG_ECG','realneg_VG_ECG','rd_VG_ECG','rdpv_VG_ECG','pd_VG_ECG','pdpv_VG_ECG','rdaov_VG_ECG','rdaovpv_VG_ECG','pdaov_VG_ECG','pdaovpv_VG_ECG',...
    'realpos_MS_ECG','realneg_MS_ECG_MS_ECG','rd_MS_ECG','rdpv_MS_ECG','pd_MS_ECG','pdpv_MS_ECG','rdaov_MS_ECG','rdaovpv_MS_ECG','pdaov_MS_ECG','pdaovpv_MS_ECG','totperm')

%%
clusteralpha=0.1;
montecarloalpha=0.05;
maxTime=0.9;
[realpos_VG_GFP,realneg_VG_GFP]=recompute_clusters_lme(rdaov_VG_GFP(times>-0.25 & times<maxTime), rdaovpv_VG_GFP(times>-0.25 & times<maxTime), pdaov_VG_GFP(:,times>-0.25 & times<maxTime), pdaovpv_VG_GFP(:,times>-0.25 & times<maxTime),...
    montecarloalpha,clusteralpha,totperm,times(times>-0.25 & times<maxTime));

[realpos_VG_ECG,realneg_VG_ECG]=recompute_clusters_lme(rdaov_VG_ECG(times>-0.25 & times<maxTime), rdaovpv_VG_ECG(times>-0.25 & times<maxTime), pdaov_VG_ECG(:,times>-0.25 & times<maxTime), pdaovpv_VG_ECG(:,times>-0.25 & times<maxTime),...
    montecarloalpha,clusteralpha,totperm,times(times>-0.25 & times<maxTime));

[realpos_MS_GFP,realneg_MS_GFP]=recompute_clusters_lme(rdaov_MS_GFP(times>-0.25 & times<maxTime,:), rdaovpv_MS_GFP(times>-0.25 & times<maxTime,:), pdaov_MS_GFP(:,times>-0.25 & times<maxTime,:), pdaovpv_MS_GFP(:,times>-0.25 & times<maxTime,:),...
    montecarloalpha,clusteralpha,totperm,times(times>-0.25 & times<maxTime));

[realpos_MS_ECG,realneg_MS_ECG]=recompute_clusters_lme(rdaov_MS_ECG(times>-0.25 & times<maxTime,:), rdaovpv_MS_ECG(times>-0.25 & times<maxTime,:), pdaov_MS_ECG(:,times>-0.25 & times<maxTime,:), pdaovpv_MS_ECG(:,times>-0.25 & times<maxTime,:),...
    montecarloalpha,clusteralpha,totperm,times(times>-0.25 & times<maxTime));

%% plots

Colors=[25.5000  153.0000   51.0000;... %On
    255.0000  127.5000   25.5000;...    %MW
    102.0000  165.7500  204.0000]./255; %MB

times2=times(times>-0.25 & times<maxTime);

h=figure('Position',[1441         229         965         689]);
subplot(2,2,1);
hplot=[];
for nGr=1:size(group_MS_ERP,1)/size(allERP,1)
    temp_erp=squeeze(data_MS_ECG(group_MS_ERP(:,1)==nGr,times>-0.25 & times<maxTime));
    [~,hplot(end+1)]=simpleTplot(times2,temp_erp,0,Colors(nGr,:),0,'-',0.2,1,0,1,0);
end
sig_samples=find(ismember(realpos_MS_ECG{1}.clusters,find(realpos_MS_ECG{1}.pmonte<montecarloalpha)));
scatter(times2(sig_samples),23*ones(1,length(sig_samples)),'Marker','s','MarkerFaceColor',Colors(3,:))

format_fig;
ylim([-25 60])
xlim([-0.25 maxTime])
ylabel('ECG ERP (\muV)')
xlabel('Time from R-peak (s)')
legend(hplot,{'ON','MW','MB'})

cmap2=cbrewer('seq','OrRd',6); cmap2(1:2,:)=[];
subplot(2,2,2);
ylim([-25 60])
xlim([-0.25 maxTime])
Cluster0=[0.268 0.388];
jbfill(Cluster0,[60 60],[-25 -25],[255 223 23]/255,[255 223 23]/255,1,0.2);

hplot=[];
for nGr=1:size(group_VG_ERP,1)/size(allERP,1)
    temp_erp=squeeze(data_VG_ECG(group_VG_ERP(:,1)==nGr,times>-0.25 & times<maxTime));
    [~,hplot(end+1)]=simpleTplot(times2,temp_erp,0,cmap2(nGr,:),0,'-',0.2,1,0,1,0);
end
sig_samples=find(ismember(realpos_VG_ECG{1}.clusters,find(realpos_VG_ECG{1}.pmonte<montecarloalpha)));
scatter(times2(sig_samples),22*ones(1,length(sig_samples)),'sk','filled')

format_fig;
xlabel('Time from R-peak (s)')
legend(hplot,{'V1','V2','V3','V4'})


subplot(2,2,3);
ylim([0 4])
xlim([-0.25 maxTime])
Cluster1=[0.524 0.728];
jbfill(Cluster1,[4 4],[0 0],[255 223 23]/255,[255 223 23]/255,1,0.2);

data_MS_GFP=[];
for nGr=1:size(group_MS_ERP,1)/size(allERP,1)
    temp_erp=squeeze(rms(data_MS_ERP(:,times>-0.25 & times<maxTime,group_MS_ERP(:,1)==nGr), 1))';
    simpleTplot(times2,temp_erp,0,Colors(nGr,:),0,'-',0.2,1,0,1,0);
    data_MS_GFP=cat(1,data_MS_GFP,temp_erp);
end
format_fig;
sig_samples=find(ismember(realpos_MS_GFP{1}.clusters,find(realpos_MS_GFP{1}.pmonte<montecarloalpha)));
scatter(times2(sig_samples),ones(1,length(sig_samples)),'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k')

ylabel('HEP GFP (\muV)')
xlabel('Time from R-peak (s)')

cmap2=cbrewer('seq','OrRd',6); cmap2(1:2,:)=[];
subplot(2,2,4);

data_VG_GFP=[];
for nGr=1:size(group_VG_ERP,1)/size(allERP,1)
    temp_erp=squeeze(rms(data_VG_ERP(:,times>-0.25 & times<maxTime,group_VG_ERP(:,1)==nGr), 1))';
    simpleTplot(times2,temp_erp,0,cmap2(nGr,:),0,'-',0.2,1,0,1,0);
    data_VG_GFP=cat(1,data_VG_GFP,temp_erp);
end
format_fig;
sig_samples=find(ismember(realpos_VG_GFP{1}.clusters,find(realpos_VG_GFP{1}.pmonte<montecarloalpha)));
scatter(times2(sig_samples),2*ones(1,length(sig_samples)),'sk','filled')


ylim([0 4])
xlim([-0.25 maxTime])

xlabel('Time from R-peak (s)')

exportgraphics(h, ['plots' filesep 'ERP_HEP_EEGandECG_GFP_VSandMS.png'], 'Resolution', 300)

%%

load('ft_layout_eeg_ecg_avref.mat'); 

layout.pos(ismember(layout.label,{'TP9','TP10','ECG'}),:)=[];
layout.width(ismember(layout.label,{'TP9','TP10','ECG'}),:)=[];
layout.height(ismember(layout.label,{'TP9','TP10','ECG'}),:)=[];
layout.label(ismember(layout.label,{'TP9','TP10','ECG'}),:)=[];

cmaptest=cbrewer('div','RdBu',64); cmaptest=flipud(cmaptest);
fprintf('%3.0f/%3.0f\n',0,length(layout.label)-2)
clear tVal_* pVal_* FVal_*
for nE=1:length(layout.label)-2
    fprintf('\b\b\b\b\b\b\b\b%3.0f/%3.0f\n',nE,length(layout.label)-2)
    for nCl=1 
        if nCl==1
            temp_table=array2table([squeeze(mean(data_MS_ERP(match_str({chanlocs.labels},layout.label{nE}),times>=Cluster1(1) & times<=Cluster1(2),:),2)) group_MS_ERP],'VariableNames',{'data','ms','subj'});
        else
            temp_table=array2table([squeeze(mean(data_MS_ERP(match_str({chanlocs.labels},layout.label{nE}),times>=Cluster2(1) & times<=Cluster2(2),:),2)) group_MS_ERP],'VariableNames',{'data','ms','subj'});
        end
        temp_table.subj=categorical(temp_table.subj);
        temp_table.ms=categorical(temp_table.ms);

        mdl=fitlme(temp_table,'data~1+ms+(1|subj)');
        tVal_topo_MW{nCl}(nE)=mdl.Coefficients.tStat(2);
        tVal_topo_MB{nCl}(nE)=mdl.Coefficients.tStat(3);
        pVal_topo_MW{nCl}(nE)=mdl.Coefficients.pValue(2);
        pVal_topo_MB{nCl}(nE)=mdl.Coefficients.pValue(3);

        aov_mdl=anova(mdl);
        FVal_topo_MS{nCl}(nE)=aov_mdl.FStat(2);
        pVal_topo_MS{nCl}(nE)=aov_mdl.pValue(2);
    end
end

cmap3=cbrewer('seq','BuPu',64); cmap3(cmap3<0)=0; cmap3(cmap3>1)=1;
for nCl=1
    figure('Position',[1502         553         305         256]);
    simpleTopoPlot_ft(FVal_topo_MS{nCl}', layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',find(pVal_topo_MS{nCl}<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
    caxis([0 1]*8)
    hmap=colorbar;
    colormap(cmap3);
    title({'MS effect',['Cluster ' num2str(nCl)]})
    format_fig;
    exportgraphics(gca, ['plots' filesep 'ERP_HEP_Topo_MSeffect_C' num2str(nCl) '.png'], 'Resolution', 300)
end

for nCl=1
    figure('Position',[1502         553         305         256]);
    simpleTopoPlot_ft(tVal_topo_MW{nCl}', layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',find(pVal_topo_MW{nCl}<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
    caxis([-1 1]*4)
    hmap=colorbar;
    colormap(cmaptest);
    title({'MW vs ON',['Cluster ' num2str(nCl)]})
    format_fig;
    exportgraphics(gca, ['plots' filesep 'ERP_HEP_Topo_MWvsON_C' num2str(nCl) '.png'], 'Resolution', 300)
end

for nCl=1
    figure('Position',[1502         553         305         256]);
    simpleTopoPlot_ft(tVal_topo_MB{nCl}', layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',find(pVal_topo_MB{nCl}<0.05),'pointsymbol','o','pointcolor','k','pointsize',64,'box','no','label','no')
    caxis([-1 1]*4)
    hmap=colorbar;
    colormap(cmaptest);
    title({'MB vs ON',['Cluster ' num2str(nCl)]})
    format_fig;
    exportgraphics(gca, ['plots' filesep 'ERP_HEP_Topo_MBvsON_C' num2str(nCl) '.png'], 'Resolution', 300)
end

