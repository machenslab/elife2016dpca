function dpca_plot_romo(data, time, yspan, explVar, compNum, events, signif, marg)

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;

if strcmp(data, 'legend')
    hold on
    
    for f=1:6
        plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
    end
    
    plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
    plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
    
    text(1.2, -2, 'F1 < F2')
    text(1.2, -3, 'F1 > F2')
    text(1.2, 1, 'F1 = 10 Hz')
    text(1.2, 2, 'F1 = 14 Hz')
    text(1.2, 3, 'F1 = 18 Hz')
    text(1.2, 4, 'F1 = 24 Hz')
    text(1.2, 5, 'F1 = 30 Hz')
    text(1.2, 6, 'F1 = 34 Hz')

    axis([0 3 -4.5 7.5])
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    set(gca,'Visible','off')
        
    return
end

if isempty(time)
    time = 1:size(data,4);
end

axis([time(1) time(end) yspan])
hold on
   
for f=1:size(data,2)
    plot(time, squeeze(data(1, f, 1, :)), '--', 'color', colors(f,:), 'LineWidth', 2)
    plot(time, squeeze(data(1, f, 2, :)), 'color', colors(f,:), 'LineWidth', 2)
end

if ~isempty(explVar)
    title(['Component #' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'])
else
    title(['Component #' num2str(compNum)])
end

if ~isempty(events)
    plot([events', events'], yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3)
end
