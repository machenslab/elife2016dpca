function dpca_plot_constantinidis(data, time, yspan, explVar, compNum, events, signif, marg)

colors = jet(7);
colors = colors([7 6 4 2 1],:);
colors(3,:) = colors(3,:)*0.6;

if strcmp(data, 'legend')
    hold on
    
    for f=1:5
        plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
    end
    
    plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
    plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
    
    text(1.2, -2, 'Match trial')
    text(1.2, -3, 'Non-match trial')
    text(1.2, 1, 'Last stimulus at preferred direction')
    text(1.2, 2, 'Last stimulus at 45^\circ')
    text(1.2, 3, 'Last stimulus at 90^\circ')
    text(1.2, 4, 'Last stimulus at 135^\circ')
    text(1.2, 5, 'Last stimulus at 180^\circ')

    axis([0 3 -4.5 5.5])
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
    plot(time(time<2), squeeze(data(1, f, 1, time<2)), '--', 'color', colors(f,:), 'LineWidth', 2)
    plot(time(time>=2), squeeze(data(1, f, 1, time>=2)), '--', 'color', colors(5-f+1,:), 'LineWidth', 2)
    plot(time, squeeze(data(1, f, 2, :)),  'color', colors(f,:), 'LineWidth', 2)
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