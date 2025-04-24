%%
clear
close

%% PLOT RESULTS OF HEP
load cmaps
load(['..' filesep 'WIM_HB_HEP_lmeEEG.mat'])

% tmaps & clusters
mT1 = Results(1).x_MW.Obs;
mT1topo = mT1;
mT1(not(Results(1).x_MW.Mask))=0;

t_pos = sort(Results(1).t_surr(:,:,twin,6), 1);
tp = squeeze( prctile(t_pos,97.5,1) );
surr_mask = mT1 > tp;
mT1(not(surr_mask))=0;

tclus = find(sum(mT1)>0);

mT2 = Results(1).x_MB.Obs;
mT2topo = mT2;
mT2(not(Results(1).x_MB.Mask))=0;

% define time points
tims = t(twin);
tc1=46:76;
tc2=86:100;


%% ERPs

% index channels
ch_ecg = 63;
chm = find(sum(abs(mT1(:,tc1))>0,2)>1 & sum(abs(mT1(:,tc2))>0,2)>1);

% index cols from design matrix
mw = 6; mb = 7;
v3 = 6; v2 = 7; v1 = 8;

figure; set(gcf,'units','centimeters', 'Position', [1 1 24 24])      

subplot(2,2,1)
a = fill([.2, .2, .6, .6], [-2, 3, 3, -2], [1 1 0]);
a.FaceAlpha = 0.1;
a.LineStyle = "none";
hold on
simpleTplot(t, Results(1).betas(ch_ecg,:,1), 0, msCols(1,:));
jbfill(t,Results(1).betas(ch_ecg,:,1)+Results(1).se(ch_ecg,:,1),Results(1).betas(ch_ecg,:,1)-Results(1).se(ch_ecg,:,1),msCols(1,:),msCols(1,:),1,0.5);
simpleTplot(t, Results(1).betas(ch_ecg,:,1)+Results(1).betas(ch_ecg,:,mw), 0, msCols(2,:));
jbfill(t,Results(1).betas(ch_ecg,:,1)+Results(1).betas(ch_ecg,:,mw)+Results(1).se(ch_ecg,:,mw),Results(1).betas(ch_ecg,:,1)+Results(1).betas(ch_ecg,:,mw)-Results(1).se(ch_ecg,:,mw),msCols(2,:),msCols(2,:),1,0.5);
simpleTplot(t, Results(1).betas(ch_ecg,:,1)+Results(1).betas(ch_ecg,:,mb), 0, msCols(3,:));
jbfill(t,Results(1).betas(ch_ecg,:,1)+Results(1).betas(ch_ecg,:,mb)+Results(1).se(ch_ecg,:,mb),Results(1).betas(ch_ecg,:,1)+Results(1).betas(ch_ecg,:,mb)-Results(1).se(ch_ecg,:,mb),msCols(3,:),msCols(3,:),1,0.5);
ylim([-1.2, 2.2])
xlim([-.3, .9])
ylabel("ECG (z.u.)")
xlabel('Time from R-peak (s)')
legend("", "ON", "", "MW", "", "MB", Location="northeast");
format_fig

subplot(2,2,2)
a = fill([.2, .2, .6, .6], [-2, 3, 3, -2], [1 1 0]);
a.FaceAlpha = 0.1;
a.LineStyle = "none";
hold on
simpleTplot(t, Results(2).betas(ch_ecg,:,1), 0, cmap_rdylbu(4,:));
jbfill(t,Results(2).betas(ch_ecg,:,1)+Results(2).se(ch_ecg,:,1),Results(2).betas(ch_ecg,:,1)-Results(2).se(ch_ecg,:,1),cmap_rdylbu(4,:),cmap_rdylbu(4,:),1,0.5);
simpleTplot(t, Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v3), 0, cmap_rdylbu(3,:));
jbfill(t,Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v3)+Results(2).se(ch_ecg,:,v3),Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v3)-Results(2).se(ch_ecg,:,v3),cmap_rdylbu(3,:),cmap_rdylbu(3,:),1,0.5);
simpleTplot(t, Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v2), 0, cmap_rdylbu(2,:));
jbfill(t,Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v2)+Results(2).se(ch_ecg,:,v2),Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v2)-Results(2).se(ch_ecg,:,v2),cmap_rdylbu(2,:),cmap_rdylbu(2,:),1,0.5);
simpleTplot(t, Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v1), 0, cmap_rdylbu(1,:));
jbfill(t,Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v1)+Results(2).se(ch_ecg,:,v1),Results(2).betas(ch_ecg,:,1)+Results(2).betas(ch_ecg,:,v1)-Results(2).se(ch_ecg,:,v1),cmap_rdylbu(1,:),cmap_rdylbu(1,:),1,0.5);
ylim([-1.2, 2.2])
xlim([-.3, .9])
xlabel('Time from R-peak (s)')
legend("", "V4", "", "V3", "", "V2", "", "V1", Location="northeast");
format_fig

subplot(2,2,3)
a = fill([.2, .2, .6, .6], [-2, 2, 2, -2], [1 1 0]);
a.FaceAlpha = 0.1;
a.LineStyle = "none";
hold on
simpleTplot(t, Results(1).betas(chm,:,1), 0, msCols(1,:));
simpleTplot(t, Results(1).betas(chm,:,1)+Results(1).betas(chm,:,mw), 0, msCols(2,:));
simpleTplot(t, Results(1).betas(chm,:,1)+Results(1).betas(chm,:,mb), 0, msCols(3,:));
plot(tims(tc1), repmat(-1, numel(tc1),1), 'k', 'LineWidth', 2.5)
plot(tims(tc2), repmat(-1, numel(tc2),1), 'k', 'LineWidth', 2.5)
ylim([-1.6, 1.6])
xlim([-.3, .9])
ylabel("HEP (\muV)")
xlabel('Time from R-peak (s)')
format_fig

subplot(2,2,4)
a = fill([.2, .2, .6, .6], [-2, 2, 2, -2], [1 1 0]);
a.FaceAlpha = 0.1;
a.LineStyle = "none";
hold on
simpleTplot(t, Results(2).betas(chm,:,1), 0, cmap_rdylbu(4,:));
simpleTplot(t, Results(2).betas(chm,:,1)+Results(2).betas(chm,:,v3), 0, cmap_rdylbu(3,:));
simpleTplot(t, Results(2).betas(chm,:,1)+Results(2).betas(chm,:,v2), 0, cmap_rdylbu(2,:));
simpleTplot(t, Results(2).betas(chm,:,1)+Results(2).betas(chm,:,v1), 0, cmap_rdylbu(1,:));
ylim([-1.6, 1.6])
xlim([-.3, .9])
xlabel('Time from R-peak (s)')
format_fig


% add annotations
annotation("textbox", str="A", FontSize = 24, FontWeight = 'bold', LineStyle = 'none', Position = [0 1 0 0] );
annotation("textbox", str="B", FontSize = 24, FontWeight = 'bold', LineStyle = 'none', Position = [0 .5 0 0] );

% export
try
    export_fig( 'figure4AB', '-png' )
catch
    hgexport(gcf, 'figure4AB', hgexport('factorystyle'), 'Format', 'png'); 
end



%% Topographies

% exclude extraneous channels (avoid topo-skirt)
ch = ~matches({chanlocs.labels}, {'FT9', 'FT10', 'ECG'});

figure; set(gcf,'units','centimeters', 'Position', [1 1 24 6])      
ax = tight_subplot(1,4,[.01 .01]);

axes(ax(1))
topoplot(mean(mT1topo(ch,tc1),2), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
    'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(mT1(ch,tc1)>0,2)),'o','k',4,1} );
    format_fig
    title(['\color[rgb]{',num2str(msCols(2,:)),'}', 'MW', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], 'FontSize', 20, 'FontWeight', 'Bold', 'interpreter', 'tex')

    axes(ax(2))
topoplot(mean(mT2topo(ch,tc1),2), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
    'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1} );
    format_fig
    title(['\color[rgb]{',num2str(msCols(3,:)),'}', 'MB', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], 'FontSize', 20, 'FontWeight', 'Bold', 'interpreter', 'tex')

axes(ax(3))
topoplot(mean(mT1topo(ch,tc2),2), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
    'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(mT1(ch,tc2)>0,2)),'o','k',4,1} );
    format_fig
    title(['\color[rgb]{',num2str(msCols(2,:)),'}', 'MW', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], 'FontSize', 20, 'FontWeight', 'Bold', 'interpreter', 'tex')

axes(ax(4))
topoplot(mean(mT2topo(ch,tc2),2), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
    'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1} );
    format_fig
    title(['\color[rgb]{',num2str(msCols(3,:)),'}', 'MB', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], 'FontSize', 20, 'FontWeight', 'Bold', 'interpreter', 'tex')


% add cbar
h = colorbar('LineWidth', 1.5);
ylabel(h, 't-value', 'rotation',270, 'VerticalAlignment','bottom', 'FontSize', 18);

% adjust panel spacing
posa = { [.05, .15, .2, .9]; [.26, .15, .2, .9]; [.52, .15, .2, .9]; [.73, .15, .2, .9] };
for sx=1:numel(posa)
    set(ax(sx),'position',posa{sx})
end

% Set everything to units pixels (avoids dynamic reposition)
set(ax,'units','pix')

% adjust figure size for labelling
posf = get(gcf,'position');
set(gcf,'position',[posf(1:2) posf(3)*1.05 posf(4)*1.25])

% add annotations
annotation("textbox", str="C", FontSize = 24, FontWeight = 'bold', LineStyle = 'none', Position = [0 1.04 0 0] );
annotation("textbox", str= sprintf('[%g, %g] ms', tims(tc1(1))*1000, tims(tc1(end))*1000), ...
    FontSize = 20, FontWeight = 'bold', LineStyle = 'none', Position = [.15 .15 .5 0] );
annotation("textbox", str= sprintf('[%g, %g] ms', tims(tc2(1))*1000, tims(tc2(end))*1000), ...
    FontSize = 20, FontWeight = 'bold', LineStyle = 'none', Position = [.6 .15 .5 0] );

% export
try
    export_fig( 'figure4C', '-png' )
catch
    hgexport(gcf, 'figure4C', hgexport('factorystyle'), 'Format', 'png'); 
end

