clear all

filename = '/home/dmitry/Dropbox/Machens Lab Sharing/Dmitry_Wieland/data extraction/Christos Constantinidis/2012/Constantinidis2012_spatial.hdf5';
postorpre = 'post';
outputFileName = 'data_constantinidis_post.mat';

%% Combine Wieland's HDF5 files and put them together in one MAT

% load all trials from Wieland's HDF5 file and join them together
firingRatesPerTrial = [];
for i=1:10
    display(['Loading trial #' num2str(i) ' out of 10'])
    A = h5read(filename, ['/data/' postorpre '/trial ' num2str(i) '/StimulusOnset']);
    A = double(permute(A, [4 3 2 1]));
    firingRatesPerTrial = cat(5, firingRatesPerTrial, A);
end

numOfTrials = h5read(filename, ['/meta/' postorpre '/trial_numbers']);
numOfTrials = double(permute(numOfTrials, [3 2 1]));

areaMask = h5read(filename, ['/filter/' postorpre '/area']);
monkeyMask = h5read(filename, ['/filter/' postorpre '/monkey']);
time = h5read(filename, ['/meta/' postorpre '/t']) / 1000;

firingRatesPerTrial = firingRatesPerTrial(:,:,:,time<=4.5,:);
time = time(time<=4.5);

s = size(firingRatesPerTrial);
firingRatesAverage = zeros(s(1:4));

% averaging over trials in each condition
display('Averaging over trials...')
for n = 1:size(firingRatesPerTrial,1)
    for s = 1:size(firingRatesPerTrial,2)
        for d = 1:size(firingRatesPerTrial,3)
            firingRatesAverage(n,s,d,:) = squeeze(mean(firingRatesPerTrial(n,s,d,:,1:numOfTrials(n,s,d)),5));
        end
    end
end

% sparsifying the full trial matrix
display('Sparsifying the full matrix...')
firingRatesPerTrial(isnan(firingRatesPerTrial)) = 0;
firingRatesPerTrial_size = size(firingRatesPerTrial);
firingRatesPerTrial_sparse = sparse(firingRatesPerTrial(:));

timeEvents = [0 0.5 2 2.5 4];
timeEventsNames = {'Stimulus1 on', 'Stimulus1 off', 'Stimulus2 on', 'Stimulus2 off', 'Response cue on'};

display('Saving the file...')
save(outputFileName, 'firingRatesPerTrial_sparse', 'firingRatesPerTrial_size', ...
                     'firingRatesAverage', 'numOfTrials', ...
                     'areaMask', 'monkeyMask', 'time', 'timeEvents', 'timeEventsNames');

display('Done')

%% Merge the directions

clear all
inputFileName = 'data_constantinidis_pre.mat';
outputFileName = 'data_constantinidis_pre_merged.mat';

load(inputFileName)
firingRatesPerTrial = reshape(full(firingRatesPerTrial_sparse), firingRatesPerTrial_size);

% kick out the centre target
firingRatesAverage = firingRatesAverage(:,1:8,:,:);
numOfTrials = numOfTrials(:,1:8,:);
firingRatesPerTrial = firingRatesPerTrial(:,1:8,:,:,:);

% rotate each neuron to preferred direction
tuning = mean(reshape(firingRatesAverage(:, :, :, time>=0 & time<0.5), ...
    size(firingRatesAverage,1), size(firingRatesAverage,2), []), 3);

for i = 1:length(tuning)
    tun = tuning(i,:);
    [~, j] = max(tun);
    ind = [j:8 1:j-1];
    tuning(i,:) = tun(ind);
    
    firingRatesAverage(i,:,:,:) = firingRatesAverage(i,ind,:,:);
    numOfTrials(i,:,:) = numOfTrials(i,ind,:);
    firingRatesPerTrial(i,:,:,:,:) = firingRatesPerTrial(i,ind,:,:,:);
end

% pool together targets equally distant from the preferred direction
merge = [2 8; 3 7; 4 6];
N = size(firingRatesPerTrial, 5);
firingRatesPerTrial = cat(5, firingRatesPerTrial, nan(size(firingRatesPerTrial)));

for i=1:3
    a = merge(i,1);
    b = merge(i,2);
    n = numOfTrials(:,a,:);
    m = numOfTrials(:,b,:);
    firingRatesAverage(:,a,:,:) = bsxfun(@times, n./(n+m), firingRatesAverage(:,a,:,:)) + ...
                                  bsxfun(@times, m./(n+m), firingRatesAverage(:,b,:,:));
    
    for neur = 1:size(firingRatesPerTrial, 1)
        for d = 1:2
            nn = numOfTrials(neur, a, d);
            mm = numOfTrials(neur, b, d);
            firingRatesPerTrial(neur,a,d,:,nn+1:nn+mm) = firingRatesPerTrial(neur,b,d,:,1:mm);
        end
    end
                 
    numOfTrials(:,a,:) = n+m;
end

firingRatesAverage = firingRatesAverage(:,1:5,:,:);
numOfTrials = numOfTrials(:,1:5,:);
firingRatesPerTrial = firingRatesPerTrial(:,1:5,:,:,:);

% sparsifying the full trial matrix
display('Sparsifying the full matrix...')
firingRatesPerTrial(isnan(firingRatesPerTrial)) = 0;
firingRatesPerTrial_size = size(firingRatesPerTrial);
firingRatesPerTrial_sparse = sparse(firingRatesPerTrial(:));

display('Saving the file...')
save(outputFileName, 'firingRatesPerTrial_sparse', 'firingRatesPerTrial_size', ...
                     'firingRatesAverage', 'numOfTrials', ...
                     'areaMask', 'monkeyMask', 'time', 'timeEvents', 'timeEventsNames');

display('Done')
