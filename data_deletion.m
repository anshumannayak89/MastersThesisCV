%% Data manipulation function
%% Author: Anshuman
%% Originally Coded: 191022 (Anshuman)
%% Updated: 191024 (Musa)
%% Updated: 200109 (Musa / Anshuman)

function [new_weights,new_codebook] = data_deletion(weights,codebook,location,qImage,OFFSET)

%OFFSET = 50;
dbImage = qImage-OFFSET;
% index = find([codebook.frame_mappings]==image_no);

neuron_no = location(qImage).BMU; % Neuron after which you want to clear

new_codebook = struct('frame_mappings',[],'frequency',0);
new_weights = [];

COUNT = 0;

j = 1;

indices = [];

% range of neurons you want to manipulate data
for i = 1: neuron_no

    frame_mappings = codebook(i).frame_mappings;

    % update neuron's frequency and frame mappings
    arr = frame_mappings <= dbImage;

    if sum(arr) > 0

        new_codebook(j).frequency = sum(arr);
        new_codebook(j).frame_mappings = frame_mappings(arr);
        new_weights = [new_weights, weights(:,i)];
        j = j + 1;
    else

        indices = [indices, i];
        COUNT = COUNT + (neuron_no - i > 15);
        continue;
    end
end

if COUNT > 0

    warning('>>>  Neurons pruned from the middle indices!  <<<')
    indices
end
