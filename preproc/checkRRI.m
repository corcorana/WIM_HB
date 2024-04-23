function checkRRI(EEG, mkr, unit)
% EEG = EEGlab structure
% mkr = char vector specifying EEG event marker code indicating R-peak event
% unit = optional input for histrograms: z = z-scored intervals; ms = millisecond intervals

if nargin <3
    unit = 'z';
end

% open figure
figure('units','normalized','outerposition',[0.2 0.2 .7 .5])
sgtitle([EEG.filename(1:end-4), ': Summary of R-R intervals' ], 'FontSize', 22)

% estimate successive interbeat (R-R) intervals (seconds)
rpk_latencies = [EEG.event(strcmp({EEG.event.type}, mkr)).latency]/EEG.srate;
rri = diff(rpk_latencies)*1000;

% histogram of IBIs
subplot(1,4,1)
if strcmp(unit, 'z')
    histogram(normalize(rri))
    line([0 0], get(gca,'YLim'), 'Color', 'r', 'LineWidth', 1.5)
    line([-3 -3], get(gca,'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
    line([3 3], get(gca,'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
    xlabel('R-R Interval (z.u.)', 'FontSize', 18)
elseif strcmp(unit, 'ms')
    histogram(rri)
    line([mean(rri) mean(rri)], get(gca,'YLim'), 'Color', 'r', 'LineWidth', 1.5)
    line([mean(rri)-std(rri)*3 mean(rri)-std(rri)*3], get(gca,'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
    line([mean(rri)+std(rri)*3 mean(rri)+std(rri)*3], get(gca,'YLim'), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
    xlabel('R-R Interval (m.s.)', 'FontSize', 18)
end

% timeseries of IBIs (tachogram)
subplot(1,4,2:4)
plot(rpk_latencies(2:end)*1000, rri);
xlabel('Session Time (s)', 'FontSize', 18)
ylabel('R-R Interval (m.s.)', 'FontSize', 18)

end