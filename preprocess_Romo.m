clear all

maindir = '/home/dmitry/Dropbox/Machens Lab Sharing/Data/Ranulfo Romo/';

outputFileName = 'data_romo_15.mat';

% dirs = {'rr014/prefront', 'rr014/s2', 'rr015/prefront/left', 'rr015/prefront/right'};
% monkeys = [1 1 2 2];    % monkey mask for each folder
% areas   = [1 2 3 4];    % area mask for each folder
% switchAreaAt = NaN;     % only applies to monkeys 24-25 (where recordings from two areas are stored in the same data-file)
% f1stimulationset = [10 14 18 24 30 34]; % set of F1 stimulation frequencies
% f1onlyLower =  [];      % f1 which were always lower than f2 (none for monkeys 14-15)
% f1onlyHigher = [];      % f1 which were always higher than f2 (none for monkeys 14-15)
% electrodeNum = 8;       % number of electrodes

dirs = {'rr015/prefront/left', 'rr015/prefront/right'};
monkeys = [1 1];    % monkey mask for each folder
areas   = [1 2];    % area mask for each folder
switchAreaAt = NaN;     % only applies to monkeys 24-25 (where recordings from two areas are stored in the same data-file)
f1stimulationset = [10 14 18 24 30 34]; % set of F1 stimulation frequencies
f1onlyLower =  [];      % f1 which were always lower than f2 (none for monkeys 14-15)
f1onlyHigher = [];      % f1 which were always higher than f2 (none for monkeys 14-15)
electrodeNum = 8;       % number of electrodes

% dirs = {'rr013/m1', 'rr013/prefront/left', 'rr013/prefront/right', 'rr013/s2cont', 'rr013/sma'};
% monkeys = [1 1 1 1 1];
% areas   = [1 2 3 4 5];
% switchAreaAt = NaN;
% f1stimulationset = [10 14 18 22 26 30 34];
% f1onlyLower =  [1 2]; 
% f1onlyHigher = [6 7]; 
% electrodeNum = 8;     

% dirs = {'R24/', 'R25/'};
% monkeys = [1 2];
% areas   = [1 2]; % NB: units 1-42 are from S1 and 43-84 from S2
% switchAreaAt = 42;
% f1stimulationset = [10 14 16 18 20 22 24 26 28 30 34];
% f1onlyLower =  [1 2 3 5];
% f1onlyHigher = [7 9 10 11];
% electrodeNum = 84;

%%

f1Both = setdiff(1:length(f1stimulationset), [f1onlyLower f1onlyHigher]);

unitCounter = 0;
rejectedUnits = 0;
usedUnits = 0;
rejectedSessions = 0;
okSessions = 0;
unstableNeurons = 0;

x = -1000:1000;
gaussKernel = 1/sqrt(2*pi)/50 * exp(-x.^2/50^2/2);

time = -0.5:0.01:7.5;

monkeyMask = [];
areaMask = [];
stableMask = [];

totalMistakes = [];

% Loop over all folders
for d = 1:length(dirs)
    display(['Processing folder ' dirs{d} ' (' num2str(d) ' out of ' num2str(length(dirs)) ')'])
    
    files = what([maindir dirs{d}]);
    
    % Loop over all files in the folder
    for f = 1:length(files.mat)
        fprintf(['Processing session ' num2str(f) ' out of ' num2str(length(files.mat))])
        
        load ([maindir dirs{d} '/' files.mat{f}])
        
        freqs = cell2mat(result(2:end, 4));
        [freqset, ~, freqs] = unique(freqs);
        freqsToUse = find(ismember(freqset, f1stimulationset));
        
        % Checking that this session contains all frequencies
        if length(freqsToUse) ~= length(f1stimulationset)
            fprintf([' ...skipping, not all frequencies present: ' num2str(freqset') '\n'])
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        % sometimes some info is missing - omit these trials
        for tt = 2:size(result,1)
            if isempty(result{tt,9})
                omit(tt-1) = 1;
                fprintf(['omitting some trials in ' files.mat{f}])
            else
                omit(tt-1) = 0;
            end
        end
        
        decisions = double(cell2mat(result(2:end, 4)) > cell2mat(result(2:end,5)));
        mistakes = ~cell2mat(result(2:end, 3));
        decisions(mistakes) = ~decisions(mistakes);
                
        decisions(omit==1) = nan;
        
        % Checking that this session contains all decisions for those
        % frequencies for which f2>f1 and f2<f1 were both present
        freqsToUseBothDec = find(ismember(freqset, f1stimulationset(f1Both)));
        subset = find(ismember(freqs,freqsToUseBothDec) & mistakes==0);
        if(size(unique([freqs(subset) decisions(subset)], 'rows'), 1) ~= length(f1Both)*2)
            fprintf(' ...skipping, not enough decisions for "main" frequencies\n')
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        % Checking that this session contains expected decisions for the
        % frequencies for which always f2<f1 or always f2<f1
        ifReject = false;
        for ii = f1onlyLower
            if(isempty(find(freqs==freqsToUse(ii) & decisions==0 & mistakes==0)))
                ifReject = true;
                break
            end
        end
        for ii = f1onlyHigher
            if(isempty(find(freqs==freqsToUse(ii) & decisions==1 & mistakes==0)))
                ifReject = true;
                break
            end
        end
        if ifReject
            fprintf([' ...skipping: no right decision for frequency ' num2str(f1stimulationset(ii)) '\n'])
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        % If this point is reached, the session is fine
        okSessions = okSessions + 1;
        
        % Counting mistakes
        totalMistakes = [totalMistakes; mistakes];
        %continue
        
        % Selecting all units that have at least 1 spike in at least 1 trial
        % Nota bene: electrode #9 is recording stimulus! (for monkeys 13,14,15)
        unitsToUse = nan(1,electrodeNum);       %length(result{2,6})
        for tr=2:size(result,1)
            for u=1:electrodeNum 
                if ~isempty(result{tr,6}{u})
                    unitsToUse(u) = u;
                end
            end
        end
        rejectedUnits = rejectedUnits + length(isnan(unitsToUse));
        unitsToUse = find(~isnan(unitsToUse));
        
        % Loop over all responding units
        for unit = unitsToUse
            fprintf('.')
                        
            % Checking stability
            baselineSpikes = zeros(1, length(decisions));
            for tr = find(~isnan(decisions))
                sp = result{tr+1, 6};
                so1 = result{1+tr, 9};
                sp = sp{unit} - so1;
                baselineSpikes(tr) = length(find(sp>-500 & sp<0));
            end
            B = 10;
            baselineSpikesBatches = reshape(baselineSpikes(1:floor(length(decisions)/B)*B), B, []);
            p = anova1(baselineSpikesBatches, [], 'off');
            if p < 0.001
                unstableNeurons = unstableNeurons + 1;
                                
%                 pp = ttest(baselineSpikesBatches);
%                 pp(isnan(pp))=0;
%                 pp1 = find(diff([0 pp])==1);
%                 pp2 = find(diff([pp 0])==-1);
%                 [~, pppos] = max(pp2-pp1);
%                 ppStart = pp1(pppos);
%                 ppEnd = pp2(pppos);
%                 if isempty(ppStart)
%                     continue
%                 end
%                 ppStart = (ppStart-1) * B + 1;
%                 ppEnd = ppEnd * B;
%                 if length(decisions) - ppEnd <= B
%                     ppEnd = length(decisions);
%                 end
                
                unitCounter = unitCounter + 1;
                stableMask(unitCounter) = 0;
%                 stableRange = ppStart:ppEnd;
                stableRange = 1:length(decisions); 
            else
                unitCounter = unitCounter + 1;  
                stableMask(unitCounter) = 1;
                stableRange = 1:length(decisions); 
            end
            
            % Filling the masks
            monkeyMask(unitCounter) = monkeys(d);
            if isnan(switchAreaAt)
                areaMask(unitCounter) = areas(d);
            else
                if unit > switchAreaAt
                    areaMask(unitCounter) = areas(2);
                else
                    areaMask(unitCounter) = areas(1);
                end
            end
            
            % Loop over all frequencies
            for fr = 1:length(freqsToUse)
                
                % Loop over all decisions
                for dec = 0:1
                    
                    % For impossible/very rare frequency-decision
                    % combinations we fill the data with NaNs
                    if(ismember(fr, f1onlyLower) && dec==1 || ismember(fr, f1onlyHigher) && dec==0)
                        rate(unitCounter, fr, dec+1, :) = nan(501,1);
                        continue
                    end                    
                    
                    % For possible combinations -- Loop over trials
                    repetitions = [];
                    for tr = intersect(find(freqs == freqsToUse(fr) & decisions == dec & mistakes == 0)', stableRange)
                        sp = result{1+tr,6};
                        so1 = round(result{1+tr,9});
                        spiketrain = round(sp{unit});
                        
                        if isempty(spiketrain)
                            % If there are no spikes, taking zero firing rate
                            repetitions = [repetitions; zeros(1,801)];
                        else
                            % Otherwise convolving the spike train with the
                            % Gaussian filter
                            if spiketrain(1)==0
                                spiketrain = spiketrain(2:end);
                            end
                            
                            spiketrainfull = zeros(1,30000);
                            spiketrainfull(round(spiketrain)) = 1;
                            psth = conv(spiketrainfull, gaussKernel, 'same');
                            
                            psthChunk = psth(so1-500:so1+7500);
                            repetitions = [repetitions; psthChunk(1:10:end) * 1000];
                        end
                    end
                    
                    if isempty(repetitions)
                        rateN(unitCounter, fr, dec+1, :) = 0;
                        continue
                    end
                    
                    % Averaging over trials
                    rate(unitCounter, fr, dec+1, :) = mean(repetitions, 1);
                    N = size(repetitions, 1);
                    rateN(unitCounter, fr, dec+1, :) = N;
                    
                    % Filling in the huge array with all individual trials
                    for trtr = 1:size(repetitions,1)
                        rateAllTrials(unitCounter, fr, dec+1, :, trtr) = repetitions(trtr,:);
                    end
                    
%                     if N > 1
%                         rateSTD(unitCounter, fr, dec+1, :) = std(repetitions, [], 1);
%                         
%                         % we want 10 samples
%                         rateNoise(unitCounter, fr, dec+1, :, :) = nan(size(repetitions,2), 10);
%                         if N*(N-1)/2 >= 10
%                             % following 3 lines generate 10 random different pairs
%                             % from integers 1..N
%                             k = randperm(N*(N-1)/2, 10);
%                             q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);
%                             p = k - (q-1).*(q-2)/2;
%                             nn = 10;
%                         else
%                             k = randperm(N*(N-1)/2, N*(N-1)/2);
%                             q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);
%                             p = k - (q-1).*(q-2)/2;
%                             nn = N*(N-1)/2;
%                         end
%                         for ii = 1:nn
%                             rateNoise(unitCounter, fr, dec+1, :, ii) = ...
%                                 (repetitions(p(ii),:) - repetitions(q(ii),:)) /sqrt(2*N);
%                         end
%                     else
%                         rateSTD(unitCounter, fr, dec+1, :) = nan(1, size(repetitions,2));
%                         for ii = 1:10
%                             rateNoise(unitCounter, fr, dec+1, :, ii) = nan(1, size(repetitions,2));
%                         end
                    %end
                end
            end
        end
        fprintf('\n')
    end
end

display(['Number of rejected sessions: ' num2str(rejectedSessions)])
display(['Number of processed sessions: ' num2str(okSessions)])
display(['Number of rejected (empty) units: ' num2str(rejectedUnits)])
display(['Number of unstable units (not rejected!): ' num2str(unstableNeurons)])
display(['Number of processed units: ' num2str(unitCounter)])

% sparsifying and renaming the variables
firingRatesPerTrial_size = size(rateAllTrials);
firingRatesPerTrial_sparse = sparse(rateAllTrials(:));
firingRatesAverage = rate;
numOfTrials = rateN;
unstableNeuronsMask = ~stableMask;
timeEvents = [0 0.5 3.5 4.0 7.0];
timeEventsNames = {'F1start', 'F1stop', 'F2start', 'F2stop', 'Cue'};

display('Saving...')
save(outputFileName, 'firingRatesPerTrial_sparse', 'firingRatesPerTrial_size', ...
                     'firingRatesAverage', 'numOfTrials', ...
                     'time', 'timeEvents', 'timeEventsNames', 'unstableNeuronsMask', 'monkeyMask', 'areaMask')
                 
display('Done')

%% HDF5 export

% h5name = 'romo_monkeys_13.h5';
% h5create(h5name, '/data', size(rate))
% h5write( h5name, '/data', rate)
% h5create(h5name, '/sem', size(rateSEM))
% h5write( h5name, '/sem', rateSEM)
% h5create(h5name, '/filter/subject', size(monkeyMask))
% h5write( h5name, '/filter/subject', monkeyMask)
% h5create(h5name, '/filter/area', size(areaMask))
% h5write( h5name, '/filter/area', areaMask)
% h5create(h5name, '/meta/t', 501)
% h5write( h5name, '/meta/t', -500:10:4500)
% h5create(h5name, '/meta/sampling_in_ms', 1)
% h5write( h5name, '/meta/sampling_in_ms', 10)
% h5create(h5name, '/meta/gauss_filter_std_in_ms', 1)
% h5write( h5name, '/meta/gauss_filter_std_in_ms', 50)
% h5disp(h5name, '/', 'min') 

%% TO DISPLAY USED FREQUENCIES (exploring the datasets...)

% figure
% hold on
% xlabel('Sessions')
% ylabel('Frequencies')
% axis([0 3000 0 50])
% x=1;
% for d = 1:length(dirs)
%     files = what([maindir dirs{d}]);
%     
%     for f = 1:length(files.mat)
%         load ([maindir dirs{d} '/' files.mat{f}])
%         freqs = cell2mat(result(2:end, 4));
%         [freqset, ~, freqs] = unique(freqs);
%         plot(x, freqset, 'b*')
%         x = x+1;
%     end
%     
%     if(d>1 && monkeys(d)>monkeys(d-1))
%         plot([x x], [0 50], 'k', 'LineWidth', 4)
%     else
%         plot([x x], [0 50], 'k')
%     end
% end
% for y=[10 14 18 22 26 30 34]
%     plot([0 3000], [y y], 'r')
% end
