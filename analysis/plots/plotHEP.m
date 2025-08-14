%%
clear
close

%% PLOT RESULTS OF HEP
load cmaps
load(['..' filesep 'WIM_HB_HEP_lmeEEG_state.mat'])
load(['..' filesep 'WIM_HB_HEP_lmeEEG_surrog.mat'])

% index cols from design matrix
mw = find(contains(Results(1).cnams, 'x_MW')); mb = find(contains(Results(1).cnams, 'x_MB'));
v3 = find(contains(Results(2).cnams, 'x_V3')); v2 = find(contains(Results(2).cnams, 'x_V2')); v1 = find(contains(Results(2).cnams, 'x_V1'));

% tmaps & clusters
mT1 = Results(1).x_MW.Obs;
mT1topo = mT1;
mT1(not(Results(1).x_MW.Mask))=0;

mT2 = Results(1).x_MB.Obs;
mT2topo = mT2;
mT2(not(Results(1).x_MB.Mask))=0;

% surrogate correction
t_pos = sort(t_surr(:,:,twin,mw), 1);
tp = squeeze( prctile(t_pos,97.5,1) );
surr_mask = mT1 > tp;
mT1(not(surr_mask))=0;

% define time points
tims = t(twin);
tcpos = find(sum(mT1)>0);
tc1=48:70;
tc2=87:92;


%% Topographies

% exclude extraneous channels (avoid topo-skirt)
ch = ~matches({chanlocs.labels}, {'FT9', 'FT10', 'ECG'});

figure; set(gcf,'units','centimeters', 'Position', [1 1 24 30])      
ax = tight_subplot(3,4,[.01 .01]);

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


%% ERPs

% index channels
ch_ecg = 63;
chm = find(sum(abs(mT1(:,tc1))>0,2)>1 | sum(abs(mT1(:,tc2))>0,2)>1);

axes(ax(5))
a = fill([.2, .2, .6, .6], [-2, 2, 2, -2], [1 1 0]);
a.FaceAlpha = 0.1;
a.LineStyle = "none";
hold on
simpleTplot(t, Results(1).betas(chm,:,1), 0, msCols(1,:));
simpleTplot(t, Results(1).betas(chm,:,1)+Results(1).betas(chm,:,mw), 0, msCols(2,:));
simpleTplot(t, Results(1).betas(chm,:,1)+Results(1).betas(chm,:,mb), 0, msCols(3,:));
plot(tims(tc1), repmat(-1, numel(tc1),1), 'k', 'LineWidth', 2.5)
plot(tims(tc2), repmat(-1, numel(tc2),1), 'k', 'LineWidth', 2.5)
xline(0,'--')
ylim([-1.6, 1.6])
xlim([-.3, .9])
ylabel("HEP (\muV)")
title('Attentional State')
%legend("", "", "ON", "", "MW", "", "MB", Location="southeast");
format_fig

axes(ax(7))
a = fill([.2, .2, .6, .6], [-2, 2, 2, -2], [1 1 0]);
a.FaceAlpha = 0.1;
a.LineStyle = "none";
hold on
simpleTplot(t, Results(2).betas(chm,:,1), 0, cmap_rdylbu(4,:));
simpleTplot(t, Results(2).betas(chm,:,1)+Results(2).betas(chm,:,v3), 0, cmap_rdylbu(3,:));
simpleTplot(t, Results(2).betas(chm,:,1)+Results(2).betas(chm,:,v2), 0, cmap_rdylbu(2,:));
simpleTplot(t, Results(2).betas(chm,:,1)+Results(2).betas(chm,:,v1), 0, cmap_rdylbu(1,:));
xline(0,'--')
ylim([-1.6, 1.6])
xlim([-.3, .9])
title('Vigilance Level')
%legend("", "", "V4", "", "V3", "", "V2", "", "V1", Location="southeast");
format_fig

axes(ax(9))
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
xline(0,'--')
ylim([-1.2, 2.2])
xlim([-.3, .9])
ylabel("ECG (z.u.)")
xlabel('Time from R-peak (s)')
legend("", "ON", "", "MW", "", "MB", Location="northeast");
format_fig

axes(ax(11))
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
xline(0,'--')
ylim([-1.2, 2.2])
xlim([-.3, .9])
xlabel('Time from R-peak (s)')
legend("", "V4", "", "V3", "", "V2", "", "V1", Location="northeast");
format_fig


%% adjust positioning

% adjust panel spacing
posa = { [.07, .71, .2, .3]; [.26, .71, .2, .3]; [.52, .71, .2, .3]; [.71, .71, .2, .3];
        [.08, .39, .37, .26]; [.26, .39, 0, 0]; [.53, .39, .37, .26]; [.73, .39, 0, 0];
        [.08, .07, .37, .26]; [.26, .07, 0, 0]; [.53, .07, .37, .26]; [.73, .07, 0, 0] };

for sx=1:numel(posa)
    set(ax(sx),'position',posa{sx})
end


%% labels

% Set everything to units pixels (avoids dynamic reposition)
set(ax,'units','pix')

% adjust figure size for labelling
posf = get(gcf,'position');
%set(gcf,'position',[posf(1:2) posf(3)*1.05 posf(4)*1.1])

% add annotations
annotation("textbox", str="A", FontSize = 24, FontWeight = 'bold', LineStyle = 'none', Position = [0 1.01 0 0] );
annotation("textbox", str= sprintf('[%g, %g] ms', tims(tc1(1))*1000, tims(tc1(end))*1000), FontSize = 18, FontWeight = 'bold', LineStyle = 'none', Position = [.175 .77 .5 0] );
annotation("textbox", str= sprintf('[%g, %g] ms', tims(tc2(1))*1000, tims(tc2(end))*1000), FontSize = 18, FontWeight = 'bold', LineStyle = 'none', Position = [.625 .77 .5 0] );
annotation("textbox", str="B", FontSize = 24, FontWeight = 'bold', LineStyle = 'none', Position = [0 .71 0 0] );
annotation("textbox", str="C", FontSize = 24, FontWeight = 'bold', LineStyle = 'none', Position = [.49 .71 0 0] );


%% export

try
    export_fig( 'figure4', '-png', '-transparent' )
catch
    hgexport(gcf, 'figure4', hgexport('factorystyle'), 'Format', 'png'); 
end



