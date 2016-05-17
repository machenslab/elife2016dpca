%% CLAU

clear all
load data_claudia.mat
% load data_claudia_piecewise.mat
% time = (time-46)*0.01;
% timeEvents = time([45:91:455]);

firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);

D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 2 & meanFiringRate < 50);
t = 1:length(time);

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
decodingClasses = {[1 1; 2 2], [1 2; 1 2], [], [1 2; 3 4]};
timeSplits = find(time>timeEvents(3),1);

load('cv_claudia_withnoise.mat', 'optimalLambda') 
% optimalLambda = 1.2975e-05 % for piecewise

load 'classification_claudia_noTimeSplits_withnoise.mat'

[~,~,Cnoise] = dpca_getNoiseCovariance(firingRatesAverage(n,:,:,t), ...
    firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:));

[W,V,whichMarg] = dpca(firingRatesAverage(n,:,:,t), 50, ...
    'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

% [W,V,whichMarg] = dpca_ldaTest(firingRatesAverage(n,:,:,t), 15, combinedParams, 0, Cnoise);%, 'manova');
% componentsSignif = [];
 
% W = dpca_NIPS(firingRatesAverage(n,:,:,t), 50, [], []);
% V = W;
% componentsSignif = [];
% explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, V, 'combinedParams', combinedParams);
% [~, whichMarg] = max(explVar.margVar);
% [~, ind] = sort(explVar.componentVar, 'descend');
% W = W(:,ind);
% V = V(:,ind);
% whichMarg = whichMarg(:,ind);

explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, V, ...
    'combinedParams', combinedParams, 'Cnoise', Cnoise, ...
    'numOfTrials', numOfTrials(n,:,:));

%demixingIndices = max(explVar.margVar(:,1:15))./explVar.componentVar(1:15);
%noiseVar = diag(W'*Cnoise*W) / explVar.totalVar * 100;
%eee = [explVar.margVar(:,1:5); noiseVar(1:5)'];
%max(eee)./sum(eee)

dpca_plot(firingRatesAverage(n,:,:,t), W, V, @dpca_plot_claudia, ... %_piecewise, ...
    'whichMarg', whichMarg,                 ...
    'time', time(t),                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,               ...
    'ylims', [20 50 150 70],                ...
    'componentsSignif', componentsSignif,   ...
    'legendSubplot', 8,                     ...
    'marginalizationNames', margNames,      ...
    'explainedVar', explVar,                ...
    'marginalizationColours', margColours);


% ------------- PCA in marginalizations -------------
% dpca_perMarginalization(firingRatesAverage(n,:,:,t), @dpca_plot_claudia, 'combinedParams', combinedParams);
% ---------------------------------------------------

% ------------- PCA -------------
% X = firingRatesAverage(n,:,:,t);
% X = X(:,:);
% X = bsxfun(@minus, X, mean(X,2));
% [W,~,~] = svd(X);
% explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, W, 'combinedParams', combinedParams);
% dpca_plot(firingRatesAverage(n,:,:,t), W, W, @dpca_plot_claudia, 'explainedVar', explVar)
% -------------------------------

% -------------- "Rerouting" analysis ---------------
% dpca_test_rerouting(firingRatesAverage(n,:,:,t), time(t), timeEvents, 1);
% ---------------------------------------------------

%% ADAM

clear all
load data_adam.mat

firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
clear firingRatesPerTrial_sparse

D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,2:5,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage(:,2:5,:,:), D, []), 2);
n = find(minN >= 2 & meanFiringRate < 50);
t = 1:length(time);
%t = find(time>=timeEvents(3) & time<=timeEvents(4));

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
timeSplits = find(time>timeEvents(3),1);
decodingClasses = {[1 1; 2 2; 3 3; 4 4], [1 2; 1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6; 7 8]};

load('cv_adam_withnoise.mat', 'optimalLambda')
load classification_adam_noTimeSplits_withnoise.mat

[~,~,Cnoise] = dpca_getNoiseCovariance(firingRatesAverage(n,2:5,:,t), ...
    firingRatesPerTrial(n,2:5,:,t,:), numOfTrials(n,2:5,:));

[W,V,whichMarg] = dpca(firingRatesAverage(n,2:5,:,t), 50, ...
    'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

% W = dpca_NIPS(firingRatesAverage(n,2:5,:,t), 50, [], []);
% V = W;
% componentsSignif = [];
% explVar = dpca_explainedVariance(firingRatesAverage(n,2:5,:,t), W, V, 'combinedParams', combinedParams);
% [~, whichMarg] = max(explVar.margVar);
% [~, ind] = sort(explVar.componentVar, 'descend');
% W = W(:,ind);
% V = V(:,ind);
% whichMarg = whichMarg(:,ind);

explVar = dpca_explainedVariance(firingRatesAverage(n,2:5,:,t), W, V, ...
    'combinedParams', combinedParams, 'Cnoise', Cnoise, ...
    'numOfTrials', numOfTrials(n,2:5,:));

% demixingIndices = max(explVar.margVar(:,1:15))./explVar.componentVar(1:15);

dpca_plot(firingRatesAverage(n,2:5,:,t), W, V, @dpca_plot_adam, ...
    'X_extra', firingRatesAverage(n,:,:,t), ...
    'whichMarg', whichMarg,                 ...
    'time', time(t),                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,               ...
    'ylims', [20 20 100 50],                 ...
    'componentsSignif', componentsSignif,   ...
    'legendSubplot', 8,                     ...
    'marginalizationNames', margNames,      ...
    'explainedVar', explVar,                ...
    'marginalizationColours', margColours);

X = firingRatesAverage(n,2:5,:,t);
XF = firingRatesAverage(n,:,:,t);
X = X(:,:);
XF = XF(:,:);
XF = bsxfun(@minus, XF, mean(X,2));
ZF = W' * XF ;
Zfull = reshape(ZF(find(whichMarg==4,3),:), [3 6 2 length(t)]);

figure
for c=1:3
    subplot(1,3,c)
    conf = squeeze(mean(Zfull(c,:,:,time(t)>timeEvents(3)&time(t)<timeEvents(4)),4));
    ratios = [0    32    44    56    68    100];
    hold on
    plot(ratios, [conf(1:3,1)' conf(4:6, 2)'], '-ro', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    plot(ratios, [conf(1:3,2)' conf(4:6, 1)'], '-o', 'LineWidth', 2, 'Color', [0 0.6 0], 'MarkerSize', 5, 'MarkerFaceColor', [0 0.6 0])
    xlabel('Odour ratio')
    axis square
end


%% ROMO

clear all
load data_romo_14_15.mat
addpath dpca/matlab

firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);

D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 5 & meanFiringRate < 50 & ismember(areaMask, [1 3 4])');
t = 1:length(time);

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'Interaction'};
timeSplits = [find(time>timeEvents(2),1) find(time>timeEvents(3),1)];% find(time>timeEvents(4),1)];
decodingClasses = {[1 1; 2 2; 3 3; 4 4; 5 5; 6 6], [1 2; 1 2; 1 2; 1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6; 7 8; 9 10; 11 12]};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% linear model
% drop = [2 1; 3 1; 4 1; 5 1; 2 2; 3 2; 4 2; 5 2];
% for i=1:size(drop,1)
%     firingRatesAverage(n,drop(i,1),drop(i,2),t) = nan;
%     firingRatesPerTrial(n,drop(i,1),drop(i,2),t,:) = 0;
%     numOfTrials(n,drop(i,1),drop(i,2)) = 0;
% end
% st = [10 14 18 24 30 34; 10 14 18 24 30 34]';
% dec = [1 1 1 1 1 1; 2 2 2 2 2 2]';
% intercept = ones(size(st));
% for i = 1:length(n)
%     d = squeeze(firingRatesAverage(n(i),:,:,t));
%     Y = reshape(d,length(st(:)),[]);
%     X = [intercept(:) st(:) dec(:)];
%     indnnan = find(~isnan(Y(:,1)));
%     beta = (X(indnnan,:)'*X(indnnan,:))\X(indnnan,:)'*Y(indnnan,:);
%     dhat = X*beta;
%     d(isnan(d(:))) = dhat(isnan(d(:)));
%     firingRatesAverage(n(i),:,:,t) = d;
% end

% end of linear model part

load('cv_romo_withnoise.mat', 'optimalLambda')
load classification_romo_noTimeSplits_withnoise.mat

[~,~,Cnoise] = dpca_getNoiseCovariance(firingRatesAverage(n,:,:,t), ...
     firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:));
  
[W,V,whichMarg] = dpca(firingRatesAverage(n,:,:,t), 50, ...
    'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);%, 'scale', 'yes');

% componentsSignif = [];
componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

% [W,V,whichMarg] = dpca_NIPS_dk(firingRatesAverage(n,:,:,t), 50, combinedParams);

% [W,V,whichMarg] = dpca_ldaTest(firingRatesAverage(n,:,:,t), 15, combinedParams, 1, Cnoise);
% componentsSignif = [];

% explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, V, ...
%     'combinedParams', combinedParams);

% explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, V, ...
%      'combinedParams', combinedParams, 'X_trial', firingRatesPerTrial(n,:,:,t,:), ...
%      'numOfTrials', numOfTrials(n,:,:));

explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, V, ...
     'combinedParams', combinedParams, ...
     'Cnoise', Cnoise, 'numOfTrials', numOfTrials(n,:,:));

% demixingIndices = max(explVar.margVar(:,1:15))./explVar.componentVar(1:15);
% display(['Demixing Indices of the first 15 components: ' ...
%     num2str(mean(demixingIndices),2) '+-' num2str(std(demixingIndices),2) '% (mean+-SD)'])

dpca_plot(firingRatesAverage(n,:,:,t), W, V, @dpca_plot_romo, ...
    'whichMarg', whichMarg,                 ...
    'time', time(t),                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,               ...
    'ylims', [150 150 400 150],             ...
    'componentsSignif', componentsSignif,   ...
    'legendSubplot', 16,                    ...
    'marginalizationNames', margNames,      ...
    'explainedVar', explVar,                ...
    'marginalizationColours', margColours);%, 'numCompToShow', 50);

% dpca_perMarginalization(firingRatesAverage(n,:,:,t), @dpca_plot_romo, 'combinedParams', combinedParams,...
%     'marginalizationNames', margNames)

% ------------- 3D visualization ------------
X = firingRatesAverage(n,:,:,t);
X = X(:,:);
X = bsxfun(@minus, X, mean(X,2));
Z = W(:,whichMarg==1)'*X;
Zfull = reshape(Z(1:3,:), [3 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
colors = jet(7);
colors = colors([1 2 4 5 6 7],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
ZZ = mean(Zfull, 3);
t1 = [41:110];
t2 = [391:460];

figure('Position', [2000 0 1300 1300])
axes('Units', 'Pixels', 'Position', [1 1 1300 1300], 'Color', [0 0 0])
hold on
for i=1:6
    plot3(squeeze(ZZ(1,i,:)), squeeze(ZZ(2,i,:)), squeeze(ZZ(3,i,:)), 'LineWidth', 2, 'Color', colors(i,:))
    plot3(squeeze(ZZ(1,i,t1)), squeeze(ZZ(2,i,t1)), squeeze(ZZ(3,i,t1)), 'LineWidth', 4, 'Color', colors(i,:))
    plot3(squeeze(ZZ(1,i,t2)), squeeze(ZZ(2,i,t2)), squeeze(ZZ(3,i,t2)), 'LineWidth', 4, 'Color', colors(i,:))
    plot3(squeeze(ZZ(1,i,t2(1:2:end))), squeeze(ZZ(2,i,t2(1:2:end))), squeeze(ZZ(3,i,t2(1:2:end))), '.', ...
        'MarkerSize', 25, 'Color', colors(i,:))
end
axis normal
axis(110*[-1 1 -1 1 -1 1])
axis vis3d

vidObj = VideoWriter('movie_romo.avi');
open(vidObj);
thetas = 1:1:360;
for i = 1:length(thetas)
    view([cosd(-135+thetas(i)) sind(-135+thetas(i)) 0.2])
    pause(0.01)
    thisFrame = getframe(gca);
    thisFrame.cdata = thisFrame.cdata(1:1300,1:1300,:);
    writeVideo(vidObj, thisFrame);
end
close(vidObj);

figure('Position', [2000 0 1300 1300])
axes('Units', 'Pixels', 'Position', [1 1 1300 1300], 'Color', [0 0 0])
hold on
axis(110*[-1 1 -1 1])
cols = [];
for i=1:6
    cols = [cols; linspace(colors(i,1),0,30)' linspace(colors(i,2),0,30)' linspace(colors(i,3),0,30)'];
end
colormap(cols)
vidObj = VideoWriter('movie_romo4.avi');
open(vidObj);
for t = 1:size(ZZ,4)
    if t == 50
        plot([-100 100], [-40 40], 'Color', [.5 .5 .5])
        text(80, 45, 'F1 onset', 'Color', [.9 .9 .9], 'FontSize', 13);
    elseif t == 400
        plot([-100 100], [40 -40], 'Color', [.5 .5 .5])
        text(80, -45, 'F2 onset', 'Color', [.9 .9 .9], 'FontSize', 13);
    elseif t == 125
        ht = text(-10, 25, 'Delay period', 'Color', [.9 .9 .9], 'FontSize', 13);
    elseif t == 375
        delete(ht)
    end
    
    for i=1:6
        ind = max(1,t-29):t;
        hh(i) = surface('XData', [squeeze(ZZ(1,i,ind)) squeeze(ZZ(1,i,ind))],             ... % N.B.  XYZC Data must have at least 2 cols
            'YData', [squeeze(ZZ(2,i,ind)) squeeze(ZZ(2,i,ind))],             ...
            'ZData', zeros(length(ind),2)-1, ...
            'CData', [(i-1)*30+(length(ind):-1:1)' (i-1)*30+(length(ind):-1:1)'],             ...
            'FaceColor', 'none',        ...
            'EdgeColor', 'interp',      ...
            'Marker', 'none',           ...
            'LineWidth', 3);
%         hh(i) = plot(squeeze(ZZ(1,i,ind)), squeeze(ZZ(2,i,ind)), ...
%             'LineWidth', 2, 'Color', colors(i,:));
        h(i) = plot(squeeze(ZZ(1,i,t)), squeeze(ZZ(2,i,t)), 'o', ...
            'MarkerSize', 10, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    end
    pause(0.01)
    thisFrame = getframe(gca);
    thisFrame.cdata = thisFrame.cdata(1:1296,1:1296,:);
    writeVideo(vidObj, thisFrame);
    
    delete(h)
    delete(hh)
end
close(vidObj);

% ------------- DEMIXING Indices ------------
% iDPCA = max(explVar.margVar(:,1:15))./explVar.componentVar(1:15);
% X = firingRatesAverage(n,:,:,t);
% X = X(:,:);
% X = bsxfun(@minus, X, mean(X,2));
% [U,~,~] = svd(X);
% explVarPCA = dpca_explainedVariance(firingRatesAverage(n,:,:,t), U(:,1:15), U(:,1:15), ...
%     'combinedParams', combinedParams);
% iPCA = max(explVarPCA.margVar(:,1:15))./explVarPCA.componentVar(1:15);
% Xm = dpca_marginalize(firingRatesAverage(n,:,:,t), 'combinedParams', combinedParams, 'ifFlat', 'yes');
% vars = [var(Xm{1},[],2)'; var(Xm{2},[],2)'; var(Xm{3},[],2)'; var(Xm{4},[],2)'];
% iN = max(vars) ./ sum(vars);
% boxplot([iDPCA iPCA iN], [ones(1,15) ones(1,15)*2 ones(1,length(iN))*3])
% -------------------------------------------

% -------- EXAMPLE NEURONS -----------
% nn = find(minN >= 5 & meanFiringRate < 40 & meanFiringRate > 20 & ismember(areaMask, [1 3 4])');
% ind = nn(randperm(length(nn),100));
% ind(1:6) = [270 250 2061 261 1949 223];
% figure
% for i=1:12%ind
%     %subplot(2,3,i)
%     %subplot(10,10,i)
%     subplot(3,4,i)
%     dpca_plot_romo(firingRatesAverage(ind(i),:,:,:), time(t), [0 100], 0, 0, timeEvents, [], 0)
%     title('')
%     set(gca, 'XTick', [])
%     set(gca, 'YTick', [])
% end
% -----------------------------------

% ----------- PCA -------------
% X = firingRatesAverage(n,:,:,t);
% X = X(:,:);
% X = bsxfun(@minus, X, mean(X,2));
% [U,S,~] = svd(X);
% explVar_pca = diag(S.^2)/sum(X(:).^2) * 100;
% Z = U(:,1:20)'*X;
% Zfull = reshape(Z, [20 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
% figure
% for i=1:6%19
%     subplot(2,3,i)
%     %subplot(4,5,i)
%     if i==1
%         yl = [-450 450];
%     else
%         yl = [-200 200];
%     end
%     dpca_plot_romo(Zfull(i,:,:,:), time(t), yl, explVar_pca(i), i, timeEvents, [], 0)
% end
% subplot(4,5,20)
% dpca_plot_romo('legend')
% ----------------------------


%% CONSTANTINIDIS

clear all
preorpost = 'post';

load(['data_constantinidis_' preorpost '_merged.mat'])

firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);

D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 5 & (monkeyMask==0 | monkeyMask==3) & meanFiringRate < 50);
t = 1:length(time);

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'Interaction'};
timeSplits = [find(time>timeEvents(2),1) find(time>timeEvents(3),1) find(time>timeEvents(4),1)];
decodingClasses = {[1 1; 2 2; 3 3; 4 4; 5 5], [1 2; 1 2; 1 2; 1 2; 1 2], [], [1 2; 3 4; 5 6; 7 8; 9 10]};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

[~,~,Cnoise] = dpca_getNoiseCovariance(firingRatesAverage(n,:,:,t), ...
     firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:));

load(['cv_constantinidis_' preorpost '_merged_withnoise.mat'], 'optimalLambda')
load(['classification_constantinidis_' preorpost '_merged_noTimeSplits_withnoise.mat'])

[W,V,whichMarg] = dpca(firingRatesAverage(n,:,:,t), 50, ...
    'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

% [W,V,whichMarg] = dpca_NIPS_dk(firingRatesAverage(n,:,:,t), 50, combinedParams);
% componentsSignif = [];
% 
% [W,V,whichMarg] = dpca_ldaTest(firingRatesAverage(n,:,:,t), 15, combinedParams, 1, Cnoise);
% componentsSignif = [];

explVar = dpca_explainedVariance(firingRatesAverage(n,:,:,t), W, V, ...
    'combinedParams', combinedParams, 'Cnoise', Cnoise, ...
    'numOfTrials', numOfTrials(n,:,:));

dpca_plot(firingRatesAverage(n,:,:,t), W, V, @dpca_plot_constantinidis, ...
    'whichMarg', whichMarg,                 ...
    'time', time(t),                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,               ...
    'ylims', [100 30 150 100],               ...
    'componentsSignif', componentsSignif,   ...
    'legendSubplot', 16,                    ...
    'marginalizationNames', margNames,      ...
    'explainedVar', explVar,                ...
    'marginalizationColours', margColours);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUPPLEMENTARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-validation for optimal regularization parameter, results

clear all

colors = [23 100 171; 187 20 25; 150 150 150; 114 97 171; 0 0 0]/256;

filenames = {'cv_claudia_withnoise.mat', 'cv_adam_withnoise.mat', 'cv_romo_withnoise.mat', ...
    'cv_constantinidis_post_merged_withnoise.mat', 'cv_constantinidis_pre_merged_withnoise.mat'};

titles = {'Feierstein dataset', 'Kepecs dataset', 'Romo dataset', 'Constantinidis dataset', 'Constantinidis PRE dataset'};

xticks = [1e-07:1e-07:1e-06 2e-06:1e-06:1e-05 2e-05:1e-05:1e-04 2e-04:1e-04:1e-03];
xtickLabels = num2cell(xticks);
for i=setdiff(1:length(xticks), [1 10 19 28])
    xtickLabels{i} = '';
end

% figure
% hold on
% title({'Relative cross-validation errors (1-R^2) with 10*4 demixed components,', 'mean over 10 repetitions'})
% xlabel('Regularization parameter, lambda')
% ylabel('Residual variance over total test variance')
%
% for i=1:5
%     load(filenames{i})
%     meanError = mean(errors,2);
%     [~, ind] = min(meanError);
%     
%     hh = patch([log(lambdas) fliplr(log(lambdas))], [min(errors,[],2)' fliplr(max(errors,[],2)')], colors(i,:)*0.2+[1 1 1]*0.8);
% %    set(hh, 'FaceAlpha', 0.2)
%     set(hh, 'EdgeColor', 'none')
%     
%     h(i) = plot(log(lambdas), meanError, 'Color', colors(i,:), 'LineWidth', 2);
%     plot(log(lambdas(ind)), meanError(ind), '.', 'Color', colors(i,:), 'MarkerSize', 30)
% end
%
% set(gca,'XTick', log(xticks))
% set(gca,'XTickLabel', xtickLabels)
% 
% legend(h, 'Feierstein dataset', 'Kepecs dataset', 'Romo dataset', 'Constantinidis dataset', ...
%     'Constantinidis PRE dataset', 'Location', 'SouthWest')
% axis([log(lambdas(1)) log(lambdas(end)) 0 1.2])
% plot(xlim, [1 1], 'k')

figure('Position', [100 100 1800 400])
for dataset=1:4
    subplot(1,4,dataset)
    hold on
    xlabel('Regularization parameter, lambda')
    ylabel('Residual variance over total test variance')
    title(titles{dataset})
    
    load(filenames{dataset})
    
    errorsMarg = cat(1, errorsMarg, shiftdim(errors,-1));
    meanError = mean(errorsMarg,3);
    [~, ind] = min(meanError, [], 2);
    
    h = [];
    for i=1:5
        err = squeeze(errorsMarg(i,:,:));
        hh = patch([log(lambdas) fliplr(log(lambdas))], [min(err,[],2)' fliplr(max(err,[],2)')], colors(i,:)*0.2+[1 1 1]*0.8);
        %set(hh, 'FaceAlpha', 0.2)
        set(hh, 'EdgeColor', 'none')
    end
    for i=1:5
        h(i) = plot(log(lambdas), meanError(i,:), 'Color', colors(i,:), 'LineWidth', 2);
        plot(log(lambdas(ind(i))), meanError(i,ind(i)), '.', 'Color', colors(i,:), 'MarkerSize', 30)
    end
    
    set(gca,'XTick', log(xticks))
    set(gca,'XTickLabel', xtickLabels)
    if dataset==4
        legend(h, {'Stimulus', 'Decision', 'Cond.-independent', 'Interaction', 'Together'}, 'Location', 'SouthEast')
        legend boxoff
    end
    axis([log(lambdas(1)) log(lambdas(end)) 0 1.2])
    plot(xlim, [1 1], 'k')
end

%% Accuracies -- summary figure for the supplementaries

clear all
figure
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [18 27]/2.54);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 18 27]/2.54);

load('data_romo_14_15.mat', 'time')
load classification_romo_noTimeSplits_withnoise.mat
dpca_accuracyPlot(time, accuracy, brier, accuracyShuffle, brierShuffle, [6 6 6 2 2 2 12 12 12], 99, 1)

preorpost = 'post';
load(['data_constantinidis_' preorpost '_merged.mat'], 'time')
load(['classification_constantinidis_' preorpost '_merged_noTimeSplits_withnoise.mat'])
dpca_accuracyPlot(time', accuracy, brier, accuracyShuffle, brierShuffle, [5 5 5 2 2 2 10 10 10], 99, 2)

% preorpost = 'pre';
% load(['data_constantinidis_' preorpost '_merged.mat'], 'time')
% load(['classification_constantinidis_' preorpost '_merged_noTimeSplits_withnoise.mat'])
% dpca_accuracyPlot(time', accuracy, brier, accuracyShuffle, brierShuffle, [5 5 5 2 2 2 10 10 10], 99, 3)

load('data_claudia.mat', 'time')
load 'classification_claudia_noTimeSplits_withnoise.mat'
dpca_accuracyPlot(time, accuracy, brier, accuracyShuffle, brierShuffle, [2 2 2 2 2 2 4 4 4], 99, 3)

load('data_adam.mat', 'time')
load 'classification_adam_noTimeSplits_withnoise.mat'
dpca_accuracyPlot(time, accuracy, brier, accuracyShuffle, brierShuffle, [4 4 4 2 2 2 8 8 8], 99, 4)

% set(gcf, 'renderer', 'painters');
%print(gcf, '-dpdf', 'my-figure.pdf');
%print(gcf, '-dpng', 'my-figure.png');
% print(gcf, '-depsc2', 'decoding.eps');

%% Mante approach - summary figure for the supplementaries

clear all
figure('Position', [100 100 1800 800])

load data_romo_14_15.mat
firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 5 & meanFiringRate < 50 & ismember(areaMask, [1 3 4])');
t = 1:length(time);

[Wmante_romo, betas_romo] = dpca_mante(firingRatesAverage(n,:,:,t), firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:), [10 14 18 24 30 34], [1 2 4 3]);
% save('W_mante_new.mat', 'Wmante_romo', 'betas_romo')
% load('W_mante_new.mat', 'Wmante_romo');
Wmante = Wmante_romo;

X = firingRatesAverage(n,:,:,t);
X = X(:,:); 
X = bsxfun(@minus, X, mean(X,2));
Z = Wmante'*X;
Z = bsxfun(@times, Z, 1./std(Z,[],2));
Zfull = reshape(Z, [4 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
ylims = [6 6 6 6];%[300 300 300 300];
for i=1:4
    subplot(4,5,(i-1)*5 + 1)
    dpca_plot_romo(Zfull(i,:,:,:), time(t), ylims(i)*[-1 1], [], i, timeEvents, [], 0)
    title('')
end

Y = X;
X = Z([1 2 4],:);
B = Y*X'/(X*X');
explVar = (1 - sum(sum((Y-B*X).^2)) / sum(Y(:).^2)) * 100;
display(['Expl var: ' num2str(explVar)])

preorpost = 'post';
load(['data_constantinidis_' preorpost '_merged.mat'])
firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 5 & (monkeyMask==0 | monkeyMask==3) & meanFiringRate < 50);
t = 1:length(time);

[Wmante_const, betas_const] = dpca_mante(firingRatesAverage(n,:,:,t), firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:), [0 45 90 135 180], [1 4 2 3]);
% save('W_mante_new.mat', 'Wmante_const', 'betas_const', '-append')
% load('W_mante_new.mat', 'Wmante_const');
Wmante = Wmante_const;

X = firingRatesAverage(n,:,:,t);
X = X(:,:);
X = bsxfun(@minus, X, mean(X,2));
Z = Wmante'*X;
Z = bsxfun(@times, Z, 1./std(Z,[],2));
Zfull = reshape(Z, [4 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
ylims = [6 6 6 6];
for i=1:4
    subplot(4,5,(i-1)*5 + 2)
    dpca_plot_constantinidis(Zfull(i,:,:,:), time(t), ylims(i)*[-1 1], [], i, timeEvents, [], 0)
    title('')
end

preorpost = 'pre';
load(['data_constantinidis_' preorpost '_merged.mat'])
firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 5 & (monkeyMask==0 | monkeyMask==3) & meanFiringRate < 50);
t = 1:length(time);

[Wmante_constPRE, betas_constPRE] = dpca_mante(firingRatesAverage(n,:,:,t), firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:), [0 45 90 135 180], [1 4 2 3]);
% save('W_mante_new.mat', 'Wmante_constPRE', 'betas_constPRE', '-append')
% load('W_mante_new.mat', 'Wmante_constPRE');
Wmante = Wmante_constPRE;

X = firingRatesAverage(n,:,:,t);
X = X(:,:);
X = bsxfun(@minus, X, mean(X,2));
Z = Wmante'*X;
Z = bsxfun(@times, Z, 1./std(Z,[],2));
Zfull = reshape(Z, [4 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
ylims = [6 6 6 6];
for i=1:4
    subplot(4,5,(i-1)*5 + 3)
    dpca_plot_constantinidis(Zfull(i,:,:,:), time(t), ylims(i)*[-1 1], [], i, timeEvents, [], 0)
    title('')
end

load data_claudia.mat
firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,:,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage, D, []), 2);
n = find(minN >= 2 & meanFiringRate < 50);
t = 1:length(time);

[Wmante_clau, betas_clau] = dpca_mante(firingRatesAverage(n,:,:,t), firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:), [0 1], [4 2 1 3]);
% save('W_mante_new.mat', 'Wmante_clau', 'betas_clau', '-append')
% load('W_mante_new.mat', 'Wmante_clau');
Wmante = Wmante_clau;

X = firingRatesAverage(n,:,:,t);
X = X(:,:);
X = bsxfun(@minus, X, mean(X,2));
Z = Wmante'*X;
Z = bsxfun(@times, Z, 1./std(Z,[],2));
Zfull = reshape(Z, [4 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
ylims = [6 6 6 6];
for i=1:4
    subplot(4,5,(i-1)*5 + 4)
    dpca_plot_claudia(Zfull(i,:,:,:), time(t), ylims(i)*[-1 1], [], i, timeEvents, [], 0)
    title('')
end

load data_adam.mat
firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
clear firingRatesPerTrial_sparse
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,2:5,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage(:,2:5,:,:), D, []), 2);
n = find(minN >= 2 & meanFiringRate < 50);
t = 1:length(time);

[Wmante_adam, betas_adam] = dpca_mante(firingRatesAverage(n,:,:,t), firingRatesPerTrial(n,:,:,t,:), numOfTrials(n,:,:), [0 32 44 56 68 100]/100, [4 2 1 3]);
% save('W_mante_new.mat', 'Wmante_adam', 'betas_adam', '-append')
% load('W_mante_new.mat', 'Wmante_adam');
Wmante = Wmante_adam;

X = firingRatesAverage(n,:,:,t);
X = X(:,:);
X = bsxfun(@minus, X, nanmean(X,2));
Z = Wmante'*X;
Z = bsxfun(@times, Z, 1./nanstd(Z,[],2));
Zfull = reshape(Z, [4 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
ylims = [6 6 6 6];
for i=1:4
    subplot(4,5,(i-1)*5 + 5)
    dpca_plot_adam(Zfull(i,:,:,:), time(t), ylims(i)*[-1 1], [], i, timeEvents, [], 0)
    title('')
end

subplot(451)
title('ROMO')
subplot(452)
title('CONSTANTINIDIS')
subplot(453)
title('Constantinidis PRE')
subplot(454)
title('CLAUDIA')
subplot(455)
title('KEPECS')

subplot(451)
ylabel('STIMULUS')
subplot(456)
ylabel('DECISION')
subplot(4,5,11)
ylabel('TIME')
subplot(4,5,16)
ylabel('INTERACTION')

%% Clustering the neurons based on encoders, all datasets
% And also histograms of encoder weights, all datasets
% And also counting pairs of axes from different marginalizations

clear all

datasetsDataFiles = {'data_romo_14_15.mat', 'data_constantinidis_post_merged.mat', ...
    'data_constantinidis_pre_merged.mat', 'data_claudia.mat', 'data_adam.mat'};
datasetsLambdaFiles = {'cv_romo_withnoise.mat', ...
    'cv_constantinidis_post_merged_withnoise.mat', 'cv_constantinidis_pre_merged_withnoise.mat', ...
    'cv_claudia_withnoise.mat', 'cv_adam_withnoise.mat'};
datasetsNames = {'Romo', 'Const', 'Const PRE', 'Claudia', 'Adam'};
datasetsMinN = [5 5 5 2 2];
datasetsStimuliSubset = {1:6, 1:5, 1:5, 1:2, 2:5};

mode = 'recompute';
%mode = 'load';

datasetsToUse = [1 2 4 5];

if strcmp(mode, 'recompute')
    for dd = 1:length(datasetsToUse);
        dataset = datasetsToUse(dd);
        
        display(['Running dPCA on dataset ' num2str(dataset) ' [' datasetsNames{dataset} ']...'])
        pause(0.001)
        
        load(datasetsDataFiles{dataset})
        firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
        clear firingRatesPerTrial_sparse
        D = size(numOfTrials,1);
        minN = min(reshape(numOfTrials(:,datasetsStimuliSubset{dataset},:), D, []), [], 2);
        meanFiringRate = mean(reshape(firingRatesAverage(:,datasetsStimuliSubset{dataset},:,:), D, []), 2);
        n = find(minN >= datasetsMinN(dataset) & meanFiringRate < 50);
        if dataset == 1
             n = find(minN >= datasetsMinN(dataset) & meanFiringRate < 50 & ismember(areaMask, [1 3 4])');
        elseif dataset == 2 || dataset == 3
             n = find(minN >= datasetsMinN(dataset) & meanFiringRate < 50 & (monkeyMask==0 | monkeyMask==3));
        end
        t = 1:length(time);
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        load(datasetsLambdaFiles{dataset}, 'optimalLambda')
        
        [~,~,Cnoise] = dpca_getNoiseCovariance(firingRatesAverage(n,datasetsStimuliSubset{dataset},:,t), ...
            firingRatesPerTrial(n,datasetsStimuliSubset{dataset},:,t,:), numOfTrials(n,datasetsStimuliSubset{dataset},:));
        [W,V,whichMarg] = dpca(firingRatesAverage(n,datasetsStimuliSubset{dataset},:,t), 50, ...
            'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);
        encoders{dataset} = V;
        decoders{dataset} = W;
        whichMargs{dataset} = whichMarg;
        save('encoders-for-clustering.mat', 'encoders', 'decoders', 'whichMargs');
    end
else
    load('encoders-for-clustering.mat')
end

figure('Position', [100 100 1800 400])
for dd = 1:length(datasetsToUse);
    dataset = datasetsToUse(dd);
    subplot(1,4,dd)
    title(datasetsNames{dataset})
    hold on
    densityClustering(encoders{dataset}, 0.25)
    ylim([0 2])
end

figure('Position', [100 100 1800 400])
for dd = 1:length(datasetsToUse);
    dataset = datasetsToUse(dd);
    subplot(1,4,dd)
    title(datasetsNames{dataset})
    hold on
    
    bin = 0.005;
    side = 0.1;
    e = -side:bin:side;
    c = e(1:end-1)+bin/2;
    hh = [];
    for i=1:15
        hh(i,:) = histc(encoders{dataset}(:,i), e);
    end
    hh = hh(:,1:end-1)/size(encoders{dataset},1)/bin;
    plot(c, hh, 'k');
    plot(c, hh(1,:), 'r', 'LineWidth', 2);
    xlim([-side side])
    ylim([0 60])
end

dots = [];
dotsNoTime = [];
for dataset = datasetsToUse
    for i=1:15
        for j=i+1:15
            if whichMargs{dataset}(i) ~= whichMargs{dataset}(j)
                dots = [dots abs(encoders{dataset}(:,i)' * encoders{dataset}(:,j))];
            end
            if whichMargs{dataset}(i) ~= whichMargs{dataset}(j) && ...
                    whichMargs{dataset}(i) ~= 3 && whichMargs{dataset}(j) ~= 3
                dotsNoTime = [dotsNoTime abs(encoders{dataset}(:,i)' * encoders{dataset}(:,j))];
            end
        end
    end
end
length(dots)
mean(dots)
length(dotsNoTime)
mean(dotsNoTime)

%% Comparing dPCA-2015, dPCA-2011, MANOVA/LDA

clear all

datasetsDataFiles = {'data_romo_14_15.mat', 'data_constantinidis_post_merged.mat', ...
    'data_constantinidis_pre_merged.mat', 'data_claudia.mat', 'data_adam.mat'};
datasetsLambdaFiles = {'cv_romo_withnoise.mat', ...
    'cv_constantinidis_post_merged_withnoise.mat', 'cv_constantinidis_pre_merged_withnoise.mat', ...
    'cv_claudia_withnoise.mat', 'cv_adam_withnoise.mat'};
datasetsNames = {'Romo', 'Const', 'Const PRE', 'Claudia', 'Adam'};
datasetsMinN = [5 5 5 2 2];
datasetsStimuliSubset = {1:6, 1:5, 1:5, 1:2, 2:5};
plotFunctions = {@dpca_plot_romo, @dpca_plot_constantinidis, @dpca_plot_constantinidis, ...
    @dpca_plot_claudia, @dpca_plot_adam};
WmanteTable = {'Wmante_romo', 'Wmante_const', 'Wmante_constPRE', 'Wmante_clau', 'Wmante_adam'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

datasetsManteOrders = {[1 2 4 3],[1 4 2 3],[1 4 2 3],[4 2 1 3],[4 2 1 3]};
datasetsManteStimX = {[10 14 18 24 30 34],[0 45 90 135 180],[0 45 90 135 180],[0 1],[0 32 44 56 68 100]/100};

mode = 'recompute';
%mode = 'load';

componentsToShow = [1 3 1; 1 1 1; 1 2 1; 2 3 1; 2 2 1; 4 3 2; 4 4 1]; % dataset / marginalization / component number
%componentsToShow = [4 3 2; 4 4 1];
%ylims = [500 200 200 200 200 200 200];

figure

if strcmp(mode, 'recompute')
    for dataset = unique(componentsToShow(:,1))'
        display(['Running dPCA on dataset ' num2str(dataset) ' [' datasetsNames{dataset} ']...'])
        pause(0.001)
        
        load(datasetsDataFiles{dataset})
        firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);
        clear firingRatesPerTrial_sparse
        D = size(numOfTrials,1);
        minN = min(reshape(numOfTrials(:,datasetsStimuliSubset{dataset},:), D, []), [], 2);
        meanFiringRate = mean(reshape(firingRatesAverage(:,datasetsStimuliSubset{dataset},:,:), D, []), 2);
        n = find(minN >= datasetsMinN(dataset) & meanFiringRate < 50);
        if dataset == 1
             n = find(minN >= datasetsMinN(dataset) & meanFiringRate < 50 & ismember(areaMask, [1 3 4])');
        elseif dataset == 2 || dataset == 3
             n = find(minN >= datasetsMinN(dataset) & meanFiringRate < 50 & (monkeyMask==0 | monkeyMask==3));
        end
        t = 1:length(time);
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        load(datasetsLambdaFiles{dataset}, 'optimalLambda')
        
        [~,~,Cnoise] = dpca_getNoiseCovariance(firingRatesAverage(n,datasetsStimuliSubset{dataset},:,t), ...
            firingRatesPerTrial(n,datasetsStimuliSubset{dataset},:,t,:), numOfTrials(n,datasetsStimuliSubset{dataset},:));
        [W,V,whichMarg] = dpca(firingRatesAverage(n,datasetsStimuliSubset{dataset},:,t), 50, ...
            'combinedParams', combinedParams, 'lambda', optimalLambda, 'Cnoise', Cnoise);
        
        [W2011,V2011,whichMarg2011] = dpca_NIPS_dk(firingRatesAverage(n,:,:,t), 50, combinedParams);
        [Wmanova,Vmanova,whichMargManova] = dpca_ldaTest(firingRatesAverage(n,:,:,t), 15, ...
            combinedParams, 1, Cnoise, 'manova'); 
        [Wnaive,Vnaive,whichMargNaive] = dpca_ldaTest(firingRatesAverage(n,:,:,t), 15, ...
            combinedParams, 1, Cnoise, 'naive'); 
        [Wlda,Vlda,whichMargLda] = dpca_ldaTest(firingRatesAverage(n,:,:,t), 15, ...
            combinedParams, 1, Cnoise); 
        Wmante = dpca_mante(firingRatesAverage(n,:,:,t), firingRatesPerTrial(n,:,:,t,:), ...
            numOfTrials(n,:,:), datasetsManteStimX{dataset}, datasetsManteOrders{dataset});
        
        X = firingRatesAverage(n,:,:,t);
        X = X(:,:); 
        X = bsxfun(@minus, X, mean(X,2));
        
        Xm = dpca_marginalize(firingRatesAverage(n,:,:,t), 'combinedParams', combinedParams, 'ifFlat', 'yes');
        
        for comp = 1:size(componentsToShow,1)
            if componentsToShow(comp,1) ~= dataset
                continue
            end
            
            ind = find(whichMarg == componentsToShow(comp,2));
            ind = ind(componentsToShow(comp,3));
            decoders = W(:,ind);
            ind = find(whichMarg2011 == componentsToShow(comp,2));
            if ~isempty(ind)
                ind = ind(componentsToShow(comp,3));
                if corr(W2011(:,ind), decoders(:,1)) > 0
                    decoders = [decoders W2011(:,ind)];
                else
                    decoders = [decoders -W2011(:,ind)];
                end
            else
                decoders = [decoders nan(size(decoders,1),1)];
            end
            ind = find(whichMargManova == componentsToShow(comp,2));
            if ~isempty(ind)
                ind = ind(componentsToShow(comp,3));
                if corr(Wmanova(:,ind), decoders(:,1)) > 0
                    decoders = [decoders Wmanova(:,ind)];
                else
                    decoders = [decoders -Wmanova(:,ind)];
                end
            else
                decoders = [decoders nan(size(decoders,1),1)];
            end
            ind = find(whichMargNaive == componentsToShow(comp,2));
            if ~isempty(ind)
                ind = ind(componentsToShow(comp,3));
                if corr(Wnaive(:,ind), decoders(:,1)) > 0
                    decoders = [decoders Wnaive(:,ind)];
                else
                    decoders = [decoders -Wnaive(:,ind)];
                end
            else
                decoders = [decoders nan(size(decoders,1),1)];
            end
            ind = find(whichMargLda == componentsToShow(comp,2));
            if ~isempty(ind)
                ind = ind(componentsToShow(comp,3));
                if corr(Wlda(:,ind), decoders(:,1)) > 0
                    decoders = [decoders Wlda(:,ind)];
                else
                    decoders = [decoders -Wlda(:,ind)];
                end
            else
                decoders = [decoders nan(size(decoders,1),1)];
            end
            
%            if componentsToShow(comp,2) ~= 3
                ind = componentsToShow(comp,2);
%                if ind == 4
%                    ind = 3;
%               elseif ind == 3
%                Wmante = load('W_mante.mat', WmanteTable{dataset});
%                Wmante = struct2array(Wmante);
                if corr(Wmante(:,ind), decoders(:,1)) > 0
                    decoders = [decoders Wmante(:,ind)];
                else
                    decoders = [decoders -Wmante(:,ind)];
                end
%            else
%                decoders = [decoders nan(size(decoders,1),1)];
%            end
            
            %save('decoders-for-comparison.mat', decoders)
                        
            for d = find(~isnan(mean(decoders)))
                subplot(size(componentsToShow,1)+1, 7, (comp-1)*7+d)
                
                Z = decoders(:,d)'*X;
                Z = Z/std(Z);
                Zfull = reshape(Z, [1 size(firingRatesAverage,2) size(firingRatesAverage,3) length(t)]);
                plotFunctions{dataset}(Zfull(1,:,:,:), time(t), 5*[-1 1], [], [], [], [], 0)
                title('')
                
                for mm=1:length(Xm)
                    v(mm) = sum((decoders(:,d)'*Xm{mm}).^2);
                end
                vv(comp,d,:) = v/sum(v);
            end
            
            subplot(size(componentsToShow,1)+1, 7, (comp-1)*7+7)
            bar(squeeze(vv(comp,:,:)), 'stacked')
            axis([0 7 0 1])
        end    
    end
    
    for c = 1:size(vv,2)
        subplot(size(componentsToShow,1)+1, 7, size(componentsToShow,1)*7+c)
        bar(squeeze(vv(:,c,:)), 'stacked')
        axis([0 size(vv,1)+1 0 1])
    end    
end

subplot(size(componentsToShow,1)+1, 7, 1)
title('dPCA-2015')
subplot(size(componentsToShow,1)+1, 7, 2)
title('dPCA-2011')
subplot(size(componentsToShow,1)+1, 7, 3)
title('MANOVA-based demixing')
subplot(size(componentsToShow,1)+1, 7, 4)
title('Naive demixing')
subplot(size(componentsToShow,1)+1, 7, 5)
title('LDA-based demixing')
subplot(size(componentsToShow,1)+1, 7, 6)
title('Mante et al.')

colormap(margColours)
