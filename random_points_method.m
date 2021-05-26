%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In this function we compute KF/ Best matching frame
%% based on RANDOM NO OF INLIERS, i.e. we select equal no of inliers RANDOMLY
%% and try to compute the error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[bestNeuron,bestKeyFrame,bestMatchingFrame]     =   random_points_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)

error_1     = 1000;         % default error for neuron 1
error_2     = 1000;         % default error for neuron 2

% Analysis of 1st Closest Neuron / Closest Neuron (Assuming only 1 KF exist)
idx_1       = find( cellfun(@(x)isequal(x,nr1),{neuron_KeyFrame.neuron}) );

try
    keyFrame1_1         = neuron_KeyFrame(idx_1).key_frames(1);
    error_1             = image_feature_match_random_points_method(qImage,keyFrame1_1,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
    finalKeyFrame_1     = keyFrame1_1;
    
    %if a neuron has more than 1 key frame, compute the best key frame
    iNoOfKeyFrames_1    = numel(neuron_KeyFrame(idx_1).frequency);
    
    if iNoOfKeyFrames_1 > 1
        for i = 2:iNoOfKeyFrames_1
            KF          = neuron_KeyFrame(idx_1).key_frames(i);
            error_KF_n  = image_feature_match_random_points_method(qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
            
            if error_KF_n < error_1
                finalKeyFrame_1     = KF;           % Best Key Frame for Neuorn 1
                error_1             = error_KF_n;
            end
        end
    end
catch
    nr2 = nr2_temp; % If we do not find a Key Frame in closest Neuron we will reconsider 2nd closest neruron
end


% Analysis of 2nd Closest Neuron (Assuming only 1 KF exist)
if nr2 > 0      % Proceed if we have 2nd closest neuron
    idx_2   = find( cellfun(@(x)isequal(x,nr2),{neuron_KeyFrame.neuron}) );
    
    try
        keyFrame1_2         = neuron_KeyFrame(idx_2).key_frames(1);
        error_2             = image_feature_match_random_points_method(qImage,keyFrame1_2,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
        finalKeyFrame_2     = keyFrame1_2;
        
        % if a neuron has more than 1 key frame, compute the best key frame
        iNoOfKeyFrames_2    = numel(neuron_KeyFrame(idx_2).frequency);
        
        if iNoOfKeyFrames_2 > 1
            for i = 2:iNoOfKeyFrames_2
                KF          = neuron_KeyFrame(idx_1).key_frames(i);
                error_KF_n  = image_feature_match_random_points_method(qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
                
                if error_KF_n < error_2
                    finalKeyFrame_2     = KF;           % Best Key Frame for Neuorn 2
                    error_2             = error_KF_n;
                end
            end
        end
    catch
        %        fprintf('\n Warning: Not enough matching points, Key Frame from Neuron 2 (Next Closest Neuron skipped).');
    end
end

if error_2 > error_1
    bestKeyFrame  = finalKeyFrame_1;        % Winner Key Frame
    bestNeuron    = nr1;
else
    bestKeyFrame  = finalKeyFrame_2;        % Winner Key Frame
    bestNeuron    = nr2;
end


%%  Once we find Key Frames, we look through frames associated with the key frames
%%  and find the most similar/matching image
%%  -----------------------------------------------------------------------
% Winner Key Frame -- finalkeyFrame
length  = numel(keyFrame_Frame);

for iLength = 1:length
    % if we have multiple key frames, proceed further once we find the...
    % best Key Frame
    if keyFrame_Frame(iLength).key_frame == bestKeyFrame
        
        fprintf('\n\nComputing best matching frame...');
        frames_temp         = keyFrame_Frame(iLength).frames;        
        frameBMU            = keyFrame_Frame(iLength).frames(1);
        
        try
            error_1             = image_feature_match_random_points_method(qImage,frameBMU,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
        catch
            continue
        end
        bestMatchingFrame   = frameBMU;      % Best Matching Frame
        finalError          = error_1;
        
        % Find no of frames associated with key frame
        iNoOfFrames         = numel(frames_temp);
        
        % if multiple frames, find best frame among them 
        if iNoOfFrames > 1
            for i = 2:iNoOfFrames
                BF      = keyFrame_Frame(iLength).frames(i);
                try
                    error_2 = image_feature_match_random_points_method(qImage,BF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
                catch
                    continue
                end
                
                if error_2 > finalError
                    bestMatchingFrame   = frameBMU;         % Best Matching Frame
                    finalError          = error_1;
                else
                    bestMatchingFrame   = BF;               % Best Matching Frame
                    finalError          = error_2;
                end
            end
        end
    end
end
end%end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [finalError] = image_feature_match_random_points_method(qImage,dbImage,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)
% qImage    - Query Image
% dbImage   - Image we have in DataBase
% location_data - Location of images in system

image_Q     = image_name(qImage,location_data,dataset);
image_I     = image_name(dbImage,location_data,dataset);

% Read images.
% RGB images converted to Grayscale
if size(imread(image_Q),3) > 1              
    IQgray  = rgb2gray(imread(image_Q));
    Igray   = rgb2gray(imread(image_I));
else
    IQgray  = imread(image_Q);
    Igray   = imread(image_I);
end

% Detect feature points each image.
imagePointsQ    = extractFeaturesMethodType(IQgray);
imagePointsI    = extractFeaturesMethodType(Igray);

% Detected feature points converted to Descriptor points.
imagePointsQ    = featuresDescriptorMethodType(imagePointsQ.Location);
imagePointsI    = featuresDescriptorMethodType(imagePointsI.Location);

% Extract feature descriptors from each image.
featuresQ       = extractFeatures(IQgray,imagePointsQ,'Method','AUTO');
featuresI       = extractFeatures(Igray,imagePointsI,'Method','AUTO');

% Match features across the dbImage and Querry Image .
indexPairs_QI   = matchFeatures(featuresQ,featuresI);

matchedPointsQ  = imagePointsQ(indexPairs_QI(:,1),:);
matchedPointsI  = imagePointsI(indexPairs_QI(:,2),:);

rng('shuffle');
% Find INLIERS KF and Querry Image
try
    [fMAT, inliers] = estimateFundamentalMatrix(matchedPointsI,matchedPointsQ,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',0.01,'Confidence',99);
catch
    Warning ('Insufficient points to compute F-Matrix');
end
inliers_1           = find(inliers==1);
noOfInliers         = numel(inliers_1);

finalError          = 100;                      % Initial error value
randomNoPoints      = randi([10 noOfInliers]);  % Random number of points to compute error
minNoIteration      = 25;                       % Minimum number of iterations performed before we finalize error
errorMat            = zeros(1,50);              % MAtrix to hold error values after each iteration

for i = 1:50
    if randomNoPoints > noOfInliers
%        fprintf('\n No of inliers exceeded');
        return 
    end
    
    randomInliers   = inliers_1(randperm(noOfInliers, randomNoPoints));
    
    MPLocI          = matchedPointsI.Location(randomInliers,:); % Location of matched points in db Image
    MPLocQ          = matchedPointsQ.Location(randomInliers,:); % Location of matched points in Query Image
    inliers_MPLocI  = [MPLocI ones(size(MPLocI,1),1)];          % tarnsforming to the form [Xi Yi 1]
    inliers_MPLocQ  = [MPLocQ ones(size(MPLocI,1),1)];          % tarnsforming to the form [Xi Yi 1]
    
    j               = size(inliers_MPLocI,1);
    Total_Error     = 0;
    
    for i10=1:j
        LQ          = MPLocI(i10,:);
        Q11         = inliers_MPLocI(i10,:);
        Q12         = inliers_MPLocQ(i10,:);
        
        Error_1     = Q12*fMAT*Q11';
        Total_Error = Total_Error + abs(Error_1);
    end
    if Total_Error < finalError
        finalError = Total_Error;
    end
    errorMat(i) = finalError;      % update 
    
    if i >= minNoIteration && (sum(errorMat(1,i-9:i)))/10 == finalError
        % After certain minimum number of iterations are done and
        % If the error value has not changed for last 10 iterations
        % we consider that error as the minimum error
        %disp(finalError)
        break
    end    
end
end%end of function

