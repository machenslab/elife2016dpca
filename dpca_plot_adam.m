function dpca_plot_adam(data, time, yspan, explVar, compNum, events, signif, marg)

colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
colors = flipud(colors);

if strcmp(data, 'legend')
    hold on
    
    for f=1:6
        plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
    end
    
    plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
    plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
    
    plot([0.5 1], [-6 -6], 'Color', [.7 .7 .7], 'LineWidth', 4)

    text(1.2, -2, 'Rat went to the right')
    text(1.2, -3, 'Rat went to the left')
    text(1.2, 1, 'Right odour: 0%')
    text(1.2, 2, 'Right odour: 32%')
    text(1.2, 3, 'Right odour: 44%')
    text(1.2, 4, 'Left odour:  56%')
    text(1.2, 5, 'Left odour:  68%')
    text(1.2, 6, 'Left odour:  100%')
    text(1.2, -6, 'Rewarded condition')
    
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    set(gca,'Visible','off')
    axis([0 3 -7.5 6+1.5])
        
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
   
for f=1:size(data,2)
    if f<=3
        plot(time, squeeze(data(1, f, 1, :)), '--', 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, :)), 'color', colors(f,:), 'LineWidth', rewardThickness)
    else
        plot(time, squeeze(data(1, f, 1, :)), '--', 'color', colors(f,:), 'LineWidth', rewardThickness)
        plot(time, squeeze(data(1, f, 2, :)), 'color', colors(f,:), 'LineWidth', 2)
    end
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