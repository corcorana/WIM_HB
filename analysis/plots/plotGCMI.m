%%
clear
close

%% PLOT RESULTS 
load cmaps
load(['..' filesep 'WIM_HB_GCMI_lmeEEG.mat'])

%% MI profiles

xlabs = {'1.0', '1.5', '2.2', '3.3', '4.8', '7.2', '10.7', '15.8'};

figure; set(gcf,'units','centimeters', 'Position', [1 1 24 10])    

subplot(1,2,1) % MS
hold on;

x = [(1:8)-.25; 1:8; (1:8)+.25]';

y = [mean( [Results(1,1).betas(:,:,1); Results(1,2).betas(:,:,1)] ); ...
    mean( [Results(1,1).betas(:,:,1)+Results(1,1).betas(:,:,4); Results(1,2).betas(:,:,1)+Results(1,2).betas(:,:,4)] ); ...
    mean( [Results(1,1).betas(:,:,1)+Results(1,1).betas(:,:,5); Results(1,2).betas(:,:,1)+Results(1,2).betas(:,:,5)] ) ]';

sd = [std( [Results(1,1).betas(:,:,1); Results(1,2).betas(:,:,1)] ); ...
    std( [Results(1,1).betas(:,:,1)+Results(1,1).betas(:,:,4); Results(1,2).betas(:,:,1)+Results(1,2).betas(:,:,4)] ); ...
    std( [Results(1,1).betas(:,:,1)+Results(1,1).betas(:,:,5); Results(1,2).betas(:,:,1)+Results(1,2).betas(:,:,5)] ) ]';

errorbar( x, y, sd, 'LineStyle', 'none', 'Capsize', 0, 'LineWidth', 1.5, 'Color', 'k');
scatter( x, y, 75, msCols, 'filled', 'd' );

xlim([.5, 8.5])
xticklabels(xlabs)
xlabel('EEG centre frequency (Hz)')

ylim([0, 0.12])
ytickformat('%.2f')
ylabel('MI (bits)')
title('Attentional State', 'FontWeight', 'bold', 'FontSize', 20)
format_fig
    
axes('Position',[.31 .47 .15 .35])
xlim([4.6, 8.4])
xticklabels(xlabs(5:8))
box on
hold on
errorbar( x(5:end,:), y(5:end,:), sd(5:end,:), 'LineStyle', 'none', 'Capsize', 0, 'LineWidth', 1.5, 'Color', 'k');
scatter( x(5:end,:), y(5:end,:), 75, msCols, 'filled', 'd' );
ylim([0, 0.0045])
h1 = legend("", "", "", "ON", "MW", "MB");
h1pos = get(h1, 'position');
format_fig

subplot(1,2,2) % VIG
hold on;

x = [(1:8)-.3; (1:8)-.1; (1:8)+.1; (1:8)+.3]';

y = [mean( [Results(2,1).betas(:,:,1); Results(2,2).betas(:,:,1)] ); ... % V4
    mean( [Results(2,1).betas(:,:,1)+Results(2,1).betas(:,:,4); Results(2,2).betas(:,:,1)+Results(2,2).betas(:,:,4)] ); ... % V3
    mean( [Results(2,1).betas(:,:,1)+Results(2,1).betas(:,:,5); Results(2,2).betas(:,:,1)+Results(2,2).betas(:,:,5)] ); ... % V2
    mean( [Results(2,1).betas(:,:,1)+Results(2,1).betas(:,:,6); Results(2,2).betas(:,:,1)+Results(2,2).betas(:,:,6)] ) ]'; % V1

sd = [std( [Results(2,1).betas(:,:,1); Results(2,2).betas(:,:,1)] ); ...
    std( [Results(2,1).betas(:,:,1)+Results(2,1).betas(:,:,4); Results(2,2).betas(:,:,1)+Results(2,2).betas(:,:,4)] ); ...
    std( [Results(2,1).betas(:,:,1)+Results(2,1).betas(:,:,5); Results(2,2).betas(:,:,1)+Results(2,2).betas(:,:,5)] ); ...
    std( [Results(2,1).betas(:,:,1)+Results(2,1).betas(:,:,6); Results(2,2).betas(:,:,1)+Results(2,2).betas(:,:,6)] ) ]';

errorbar( x, y, sd, 'LineStyle', 'none', 'Capsize', 0, 'LineWidth', 1.5, 'Color', 'k');
scatter( x, y, 75, flipud(cmap_rdylbu(1:4,:)), 'filled', 'd' );

xlim([.5, 8.5])
xticklabels(xlabs)
ytickformat('%.2f')
ylim([0, 0.12])
ylabel('')
xlabel('EEG centre frequency (Hz)')
title('Vigilance Level', 'FontWeight', 'bold', 'FontSize', 20)
format_fig


axes('Position',[.75 .47 .15 .35])
xlim([4.5, 8.5])
xticklabels(xlabs(5:8))
box on
hold on
errorbar( x(5:end,:), y(5:end,:), sd(5:end,:), 'LineStyle', 'none', 'Capsize', 0, 'LineWidth', 1.5, 'Color', 'k');
scatter( x(5:end,:), y(5:end,:), 75, flipud(cmap_rdylbu(1:4,:)), 'filled', 'd' );
ylim([0, 0.0045])
h2 = legend("", "", "", "", "V4", "V3", "V2", "V1");
h2pos = get(h2, 'position');
format_fig

% add annotation
annotation("textbox", str="A", FontSize = 22, FontWeight = 'bold', LineStyle = 'none', Position = [0 1.03 0 0] );


set(h1, 'position',[0.13 0.32 h1pos(3:4)], 'box', 'off')
set(h2, 'position',[0.57 0.335 h2pos(3:4)], 'box', 'off')


% export
try
    export_fig( 'figure5_MI', '-png', '-transparent' )
catch
    hgexport(gcf, 'figure5_MI', hgexport('factorystyle'), 'Format', 'png'); 
end



%% MI topos

ch = ~matches({chanlocs.labels}, {'FT9', 'FT10'}); % remove skirt from topos

% MS

bh_mw = Results(1,2).x_MW.Obs;
bh_mwTopo = bh_mw;
bh_mw(not(Results(1,2).x_MW.Mask))=0;

bh_mb = Results(1,2).x_MB.Obs;
bh_mbTopo = bh_mb;
bh_mb(not(Results(1,2).x_MB.Mask))=0;

hb_mw = Results(1,1).x_MW.Obs;
hb_mwTopo = hb_mw;
hb_mw(not(Results(1,1).x_MW.Mask))=0;

hb_mb = Results(1,1).x_MB.Obs;
hb_mbTopo = hb_mb;
hb_mb(not(Results(1,1).x_MB.Mask))=0;


figure; set(gcf,'units','centimeters', 'Position', [1 1 24 11])        
ncol = size(bh_mwTopo,2)+1;
ax = tight_subplot(4,ncol,[.01 .01]);
posa = get(ax,'position');
for p = 1:ncol

    axes(ax(p))
    try
    topoplot(bh_mwTopo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
        'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(bh_mw(ch,p-1)<0,2)),'o','y',4,1} );
    format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(msCols(2,:)),'}', 'MW', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], ...
                'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex');
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05])
    end
    
    axes(ax(p+ncol))
    try
    topoplot(bh_mbTopo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
        'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(bh_mb(ch,p-1)<0,2)),'o','y',4,1} );
    format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(msCols(3,:)),'}', 'MB', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], ...
            'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex')
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05]')
    end
    
    axes(ax(p+ncol*2))
    try
    topoplot(hb_mwTopo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
        'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(hb_mw(ch,p-1)<0,2)),'o','y',4,1} );
    format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(msCols(2,:)),'}', 'MW', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], ...
                'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex');
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05])
    end
    
    axes(ax(p+ncol*3))
    try
    topoplot(hb_mbTopo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_rdbu, ...
        'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(hb_mb(ch,p-1)<0,2)),'o','y',4,1} );
    format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(msCols(3,:)),'}', 'MB', '\color{black}', ' vs ', '\color[rgb]{',num2str(msCols(1,:)),'}', 'ON' ], ...
            'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex')
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05]')
    end

    if p>1
        subtitle( sprintf('[%.1f, %.1f]', fbands(p-1,1), fbands(p-1,2)), 'FontSize', 12, 'FontWeight', 'Bold')
        set(get(gca,'subtitle'), 'rotation', 0, 'Position', [0, -.8])
    end

end

% add labels
sgtitle( 'Attentional State', 'FontWeight', 'bold', 'FontSize', 20 )
annotation("textbox", str="B", FontSize = 22, FontWeight = 'bold', LineStyle = 'none', Position = [0 1.02 0 0] );

% add cbar
h = colorbar('LineWidth', 1.5);
ylabel(h, 't-value', 'rotation',270, 'VerticalAlignment','bottom', 'FontSize', 16);
% reset size of final topoplot
set(ax(p+ncol),'position',posa{p+ncol})
% adjust cbar positioning
set(h, 'Position', [.92, .385, .017, .25])

% adjust panel spacing
for sx=1:numel(posa)
    posb = posa{sx};
    posb = [posb(1) , 0.07+posb(2)*.9, posb(3:4)];
    set(ax(sx),'position',posb)
end


% Set everything to units pixels (avoids dynamic reposition)
set(ax,'units','pix')
% adjust figure size
posf = get(gcf,'position');
set(gcf,'position',[posf(1:2) posf(3)*1.05 posf(4)*1.05])

% add X axis label
title('EEG frequency band (Hz)', 'FontSize', 16, 'FontWeight', 'Bold')
set(get(gca,'title'), 'Position', [-4.6, -1.2])


% export
try
    export_fig( 'figure5_topoAS', '-png', '-transparent' )
catch
    hgexport(gcf, 'figure5_topoAS', hgexport('factorystyle'), 'Format', 'png'); 
end



%% VIG

bh_v3 = Results(2,2).x_V3.Obs;
bh_v3Topo = bh_v3;
bh_v3(not(Results(2,2).x_V3.Mask))=0;

bh_v2 = Results(2,2).x_V2.Obs;
bh_v2Topo = bh_v2;
bh_v2(not(Results(2,2).x_V2.Mask))=0;

bh_v1 = Results(2,2).x_V1.Obs;
bh_v1Topo = bh_v1;
bh_v1(not(Results(2,2).x_V1.Mask))=0;

hb_v3 = Results(2,1).x_V3.Obs;
hb_v3Topo = hb_v3;
hb_v3(not(Results(2,1).x_V3.Mask))=0;

hb_v2 = Results(2,1).x_V2.Obs;
hb_v2Topo = hb_v2;
hb_v2(not(Results(2,1).x_V2.Mask))=0;

hb_v1 = Results(2,1).x_V1.Obs;
hb_v1Topo = hb_v1;
hb_v1(not(Results(2,1).x_V1.Mask))=0;


figure; set(gcf,'units','centimeters', 'Position', [1 1 24 15])        
ncol = size(bh_v1Topo,2)+1;
ax = tight_subplot(6,ncol,[.01 .01]);
posa = get(ax,'position');
for p = 1:ncol

    axes(ax(p))
    try
        topoplot(bh_v3Topo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_pror, ...
            'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(bh_v3(ch,p-1)<0,2)),'o','y',4,1} );
        format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(cmap_rdylbu(3,:)),'}', 'V3', '\color{black}', ' vs ', '\color[rgb]{',num2str(cmap_rdylbu(4,:)),'}', 'V4' ], ...
                'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex');
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05])
    end
    
    axes(ax(p+ncol))
    try
        topoplot(bh_v2Topo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_pror, ...
            'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(bh_v2(ch,p-1)<0,2)),'o','y',4,1} );
        format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(cmap_rdylbu(2,:)),'}', 'V2', '\color{black}', ' vs ', '\color[rgb]{',num2str(cmap_rdylbu(4,:)),'}', 'V4' ], ...
            'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex')
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05]')
    end
    
    axes(ax(p+ncol*2))
    try
        topoplot(bh_v1Topo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_pror, ...
            'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(bh_v1(ch,p-1)<0,2)),'o','y',4,1} );
        format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(cmap_rdylbu(1,:)),'}', 'V1', '\color{black}', ' vs ', '\color[rgb]{',num2str(cmap_rdylbu(4,:)),'}', 'V4' ], ...
            'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex')
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05]')
    end

    axes(ax(p+ncol*3))
    try
        topoplot(hb_v3Topo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_pror, ...
            'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(hb_v3(ch,p-1)<0,2)),'o','y',4,1} );
        format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(cmap_rdylbu(3,:)),'}', 'V3', '\color{black}', ' vs ', '\color[rgb]{',num2str(cmap_rdylbu(4,:)),'}', 'V4' ], ...
                'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex');
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05])
    end

    axes(ax(p+ncol*4))
    try
        topoplot(hb_v2Topo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_pror, ...
            'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(hb_v2(ch,p-1)<0,2)),'o','y',4,1} );
        format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(cmap_rdylbu(2,:)),'}', 'V2', '\color{black}', ' vs ', '\color[rgb]{',num2str(cmap_rdylbu(4,:)),'}', 'V4' ], ...
            'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex')
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05]')
    end

    axes(ax(p+ncol*5))
    try
        topoplot(hb_v1Topo(ch,p-1), chanlocs(ch), 'style', 'map', 'gridscale', 300, 'colormap', cmap_pror, ...
            'maplimits', [-4 4], 'emarker',{'.',[.5 .5 .5],[],1},'emarker2',{find(sum(hb_v1(ch,p-1)<0,2)),'o','y',4,1} );
        format_fig
    catch
        set(gca,'Visible','off')
    end
    if p == 2
        title(['\color[rgb]{',num2str(cmap_rdylbu(1,:)),'}', 'V1', '\color{black}', ' vs ', '\color[rgb]{',num2str(cmap_rdylbu(4,:)),'}', 'V4' ], ...
            'FontSize', 12, 'FontWeight', 'Bold', 'interpreter', 'tex')
        set(get(gca,'title'), 'rotation', 90, 'Position', [-.6, 0.05]')
    end

    if p>1
        subtitle( sprintf('[%.1f, %.1f]', fbands(p-1,1), fbands(p-1,2)), 'FontSize', 12, 'FontWeight', 'Bold')
        set(get(gca,'subtitle'), 'rotation', 0, 'Position', [0, -.8])
    end
end
    
    

% add labels
sgtitle( 'Vigilance Level', 'FontWeight', 'bold', 'FontSize', 20 )
annotation("textbox", str="C", FontSize = 22, FontWeight = 'bold', LineStyle = 'none', Position = [0 1.01 0 0] );

% add cbar
h = colorbar('LineWidth', 1.5);
ylabel(h, 't-value', 'rotation',270, 'VerticalAlignment','bottom', 'FontSize', 16);
% reset size of final topoplot
set(ax(p+ncol),'position',posa{p+ncol})
% adjust cbar positioning
set(h, 'Position', [.92, .41, .017, .19])

% adjust panel spacing
for sx=1:numel(posa)
    posb = posa{sx};
    posb = [posb(1) , 0.05+posb(2)*.9, posb(3:4)];
    set(ax(sx),'position',posb)
end


% Set everything to units pixels (avoids dynamic reposition)
set(ax,'units','pix')
% adjust figure size
posf = get(gcf,'position');
set(gcf,'position',[posf(1:2) posf(3)*1.05 posf(4)])

% add X axis label
title('EEG frequency band (Hz)', 'FontSize', 16, 'FontWeight', 'Bold')
set(get(gca,'title'), 'Position', [-4.6, -1.2])


% export
try
    export_fig( 'figure5_topoVL', '-png', '-transparent' )
catch
    hgexport(gcf, 'figure5_topoVL', hgexport('factorystyle'), 'Format', 'png'); 
end



