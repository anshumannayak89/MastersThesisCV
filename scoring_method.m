%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In this function we assign scores of keyframe/images based on :
%% 1. Error (Less is better)
%% 2. No of inlier (More is better)
%% 3. Maximum distance of inlier from mean (More is better)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestNeuron,bestKeyFrame,bestMatchingFrame] = scoring_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)

% initializing default values
finalInlierCount_1  = 0 ;
finalInlierCount_2  = 0 ;

finalMaxDistance_1  = 0 ;
finalMaxDistance_2  = 0 ;

finalError_1        = 1000 ;
finalError_2        = 1000 ;

% Analysis of 1st Closest Neuron / Closest Neuron
idx_1   = find( cellfun(@(x)isequal(x,nr1),{neuron_KeyFrame.neuron}) );   % Index of Neuron being analyzed
try
    
    keyFrame1_1     = neuron_KeyFrame(idx_1).key_frames(1);                     % Retrieving Key Frame 1
    
    % Computing Error, inliers, distance of farthest inlier from mean of all inliers
    [iErrorKeyFrame1_1,iInlierCount_1,maxDistance_1]    = image_feature_match_score(qImage,keyFrame1_1,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
    
    
    finalKeyFrame_1     = keyFrame1_1;      % Best Key Frame for Neuorn 1
    finalError_1        = iErrorKeyFrame1_1;
    finalInlierCount_1  = iInlierCount_1;
    finalMaxDistance_1  = maxDistance_1;
    
    % if a neuron has more than 1 key frame, compute the best key frame
    iNoOfKeyFrames_1    = numel(neuron_KeyFrame(idx_1).frequency);
    if iNoOfKeyFrames_1 > 1
        for i = 2:iNoOfKeyFrames_1
            KF                              = neuron_KeyFrame(idx_1).key_frames(i);
            [error,inlierCount,maxDistance] = image_feature_match_score(qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
            Score_1 = 0;
            Score_2 = 0;
            
            if finalInlierCount_1 < inlierCount
                Score_2     = Score_2 +1 ;
            else
                Score_1     = Score_1 +1 ;
            end
            
            if finalError_1 > error
                Score_2     = Score_2 +1 ;
            else
                Score_1     = Score_1 +1 ;
            end
            
            if finalMaxDistance_1 < maxDistance
                Score_2     = Score_2 +1 ;
            else
                Score_1     = Score_1 +1 ;
            end
            
            if Score_1 < Score_2
                finalKeyFrame_1     = KF;      % Best Key Frame for Neuorn 1
                finalError_1        = error;
                finalInlierCount_1  = inlierCount;
                finalMaxDistance_1  = maxDistance;
            end
        end
    end
catch
    nr2 = nr2_temp; % If we do not find a Key Frame in closest Neuron we will reconsider 2nd closest neruron
end

% Analysis of 2nd Closest Neuron
if nr2  > 0      % Proceed if we have 2nd closest neuron
    idx_2   = find( cellfun(@(x)isequal(x,nr2),{neuron_KeyFrame.neuron}) );
    
    try        
        keyFrame1_2                                 = neuron_KeyFrame(idx_2).key_frames(1);
        [iErrorBMU_1,iInlierCount_2,maxDistanc_2]   = image_feature_match_score(qImage,keyFrame1_2,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
        
        finalKeyFrame_2     = keyFrame1_2;      % Best Key Frame for Neuorn 2
        finalError_2        = iErrorBMU_1;
        finalInlierCount_2  = iInlierCount_2;
        finalMaxDistance_2  = maxDistanc_2;
        
        % if a neuron has more than 1 key frame, compute the best key frame
        iNoOfKeyFrames_2    = numel(neuron_KeyFrame(idx_2).frequency);
        if iNoOfKeyFrames_2 > 1
            for i = 2:iNoOfKeyFrames_2
                KF                              = neuron_KeyFrame(idx_2).key_frames(i);
                [error,inlierCount,maxDistance] = image_feature_match_score(qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
                Score_1 = 0;
                Score_2 = 0;
                
                if finalInlierCount_2 < inlierCount
                    Score_2     = Score_2 +1 ;
                else
                    Score_1     = Score_1 +1 ;
                end
                
                if finalError_2 > error
                    Score_2     = Score_2 +1 ;
                else
                    Score_1     = Score_1 +1 ;
                end
                
                if finalMaxDistance_2 < maxDistance
                    Score_2     = Score_2 +1 ;
                else
                    Score_1     = Score_1 +1 ;
                end
                
                if Score_1 < Score_2
                    finalKeyFrame_2     = KF;      % Best Key Frame for Neuorn 2
                    finalError_2        = error;
                    finalInlierCount_2  = inlierCount;
                    finalMaxDistance_2  = maxDistance;
                end
            end
        end
    catch
        %      fprintf('\n Warning: Not enough matching points, Key Frame from Neuron 2 (Next Closest Neuron skipped).');
    end
end

%Comparision between KeyFrames of Neuron 1 and Neuron 2
Score_1     = 0;
Score_2     = 0;
if finalInlierCount_2 > finalInlierCount_1
    Score_2     = Score_2 +1 ;
else
    Score_1     = Score_1 +1 ;
end

if finalError_2 < finalError_1
    Score_2     = Score_2 +1 ;
else
    Score_1     = Score_1 +1 ;
end

if finalMaxDistance_2 > finalMaxDistance_1
    Score_2     = Score_2 +1 ;
else
    Score_1     = Score_1 +1 ;
end


if Score_1 > Score_2
    bestKeyFrame  = finalKeyFrame_1;        % Winner Key Frame
    bestNeuron    = nr1;
else
    bestKeyFrame  = finalKeyFrame_2;        % Winner Key Frame
    bestNeuron    = nr2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MODEULE : Frame Association / Geometric Validation
%%  Once we find Key Frames, we look through frames associated with the key frames
%%  and find the most similar/matching image
%%  -----------------------------------------------------------------------
% Winner Key Frame -- finalkeyKrame
length  = numel(keyFrame_Frame);

for i   = 1:length
    if keyFrame_Frame(i).key_frame == bestKeyFrame
        
        fprintf('\n\nComputing best matching frame...');
        frames_temp     = keyFrame_Frame(i).frames;        
        frameBMU        = keyFrame_Frame(i).frames(1);
        
        try
            [iErrorBMU_1,iInlierCount_1,maxDistance_1]  = image_feature_match_score(qImage,frameBMU,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
        catch
            continue
        end
        
        bestMatchingFrame   = frameBMU;      % Best Matching Frame
        finalError          = iErrorBMU_1;
        finalInlierCount    = iInlierCount_1;
        finalMaxDistance    = maxDistance_1;
        
        % if a Key Frame has more than 1 frame, compute the best key frame
        iNoOfFrames = numel(frames_temp);
        if iNoOfFrames > 1
            for i2 = 2:iNoOfFrames
                
                BF  = keyFrame_Frame(i).frames(i2);
                
                try
                    [iErrorBMU_2,iInlierCount_2,maxDistanc_2]   = image_feature_match_score(qImage,BF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
                catch
                    continue
                end
                
                finalError_2        = iErrorBMU_2;
                finalInlierCount_2  = iInlierCount_2;
                finalMaxDistance_2  = maxDistanc_2;
                
                Score_1 = 0;
                Score_2 = 0;
                if finalInlierCount_2 > finalInlierCount
                    Score_2     = Score_2 +1 ;
                else
                    Score_1     = Score_1 +1 ;
                end
                
                if finalError_2 < finalError
                    Score_2     = Score_2 +1 ;
                else
                    Score_1     = Score_1 +1 ;
                end
                
                if finalMaxDistance_2 > finalMaxDistance
                    Score_2     = Score_2 +1 ;
                else
                    Score_1     = Score_1 +1 ;
                end
                
                if Score_2 > Score_1
                    bestMatchingFrame = BF;       % Best Frame
                end
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [error,inlier,max_distance] = image_feature_match_score(qImage,dbImage,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)
% qImage    - Query Image
% dbImage   - Image we have in DataBase
% locationKITTI - Location of images in system

image_Q = image_name(qImage,location_data,dataset);
image_I = image_name(dbImage,location_data,dataset);

% Read images.
% RGB image to Grayscale image
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
[fMAT, inliers]     = estimateFundamentalMatrix(matchedPointsI,matchedPointsQ,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',0.01,'Confidence',99);
inliers_1           = find(inliers==1);
noOfInliers         = numel(inliers_1);

MPLocI  = matchedPointsI.Location(inliers_1,:); % Location of matched points in db Image
MPLocQ  = matchedPointsQ.Location(inliers_1,:); % Location of matched points in Query Image

mu      = mean(MPLocI);      % Mean of all the inliers

inliers_MPLocI  = [MPLocI ones(size(MPLocI,1),1)]; % tarnsforming to the form [Xi Yi 1]
inliers_MPLocQ  = [MPLocQ ones(size(MPLocI,1),1)]; % tarnsforming to the form [Xi Yi 1]
j               = size(inliers_MPLocI,1);
Total_Error_1   = 0;
max_distance    = 0;

for i10=1:j
    LQ      = MPLocI(i10,:);
    Q11     = inliers_MPLocI(i10,:);
    Q12     = inliers_MPLocQ(i10,:);
    
    w       = norm(mu-LQ);        % Distance between mean and the key points
    Error_1 = 1/w * Q12*fMAT*Q11';
    
    if(abs(w) > max_distance)
        max_distance    = abs(w);
    end
    Total_Error_1       = Total_Error_1 + abs(Error_1);
end

error   = Total_Error_1;
inlier  = noOfInliers;
end %end of function