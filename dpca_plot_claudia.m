function dpca_plot_claudia(data, time, yspan, explVar, compNum, events, signif, marg)

colors = [0 0 1; 1 0 0];

if strcmp(data, 'legend')
    hold on
    
    for f=1:2
        plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
    end
    
    plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
    plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
    
    plot([0.5 1], [-6 -6], 'Color', [.7 .7 .7], 'LineWidth', 4)
    
    text(1.2, 1, 'Left odour')
    text(1.2, 2, 'Right odour')
    text(1.2, -2, 'Rat went to the right')
    text(1.2, -3, 'Rat went to the left')
    text(1.2, -6, 'Rewarded condition')

    axis([0 3 -7.5 2+1.5])
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

if marg == 4
    rewardThickness = 4;
else
    rewardThickness = 2;
end

plot(time, squeeze(data(1, 1, 1, :)), '--', 'color', colors(1,:), 'LineWidth', rewardThickness)
plot(time, squeeze(data(1, 1, 2, :)), 'color', colors(1,:), 'LineWidth', 2)
plot(time, squeeze(data(1, 2, 1, :)), '--', 'color', colors(2,:), 'LineWidth', 2)
plot(time, squeeze(data(1, 2, 2, :)), 'color', colors(2,:), 'LineWidth', rewardThickness)

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
