%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In this function we compute KF/ Best matching frame
%% based on MINIMUM NO OF INLIERS, i.e. we select equal no of inliers as
%% per the frame having the least no of inliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestNeuron,bestKeyFrame,bestMatchingFrame] = mean_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)
error_11    = 1000;
error_22    = 1000;

% Analysis of 1st Closest Neuron / Closest Neuron (Assuming only 1 KF exist)
idx_1 = find( cellfun(@(x)isequal(x,nr1),{neuron_KeyFrame.neuron}) );

try
    keyFrame1_1                 = neuron_KeyFrame(idx_1).key_frames(1);
    [iInlierCount_1,error_11]   = image_feature_match(qImage,keyFrame1_1,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
    
    finalKeyFrame_1     = keyFrame1_1;      % Best Key Frame for Neuorn 1
    finalInlierCount_1  = iInlierCount_1;
    
    % if a neuron has more than 1 key frame, compute the best key frame
    iNoOfKeyFrames_1    = numel(neuron_KeyFrame(idx_1).frequency);
    
    if iNoOfKeyFrames_1 > 1
        for i = 2:iNoOfKeyFrames_1
            KF                  = neuron_KeyFrame(idx_1).key_frames(i);
            finalInlierCount_n  = image_feature_match(qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
            
            %Comparision between KeyFrames
            if finalInlierCount_n > finalInlierCount_1
                min_inlier  = finalInlierCount_1;
            else
                min_inlier  = finalInlierCount_n;
            end
            
            error_11        = image_feature_match_mean(min_inlier,qImage,finalKeyFrame_1,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
            error_KF_n      = image_feature_match_mean(min_inlier,qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
            
            if error_KF_n < error_11
                finalKeyFrame_1     = KF;      % Best Key Frame for Neuorn 1
                error_11            = error_KF_n;
            end
        end
    end
catch
    nr2 = nr2_temp; % If we do not find a Key Frame in closest Neuron we will reconsider 2nd closest neruron
end

% Analysis of 2nd Closest Neuron (Assuming only 1 KF exist)
iNoOfKeyFrames_2 = 0;
if nr2 > 0      % Proceed if we have 2nd closest neuron
    idx_2 = find( cellfun(@(x)isequal(x,nr2),{neuron_KeyFrame.neuron}) );
    
    try
        keyFrame1_2                 = neuron_KeyFrame(idx_2).key_frames(1);
        [iInlierCount_2,error_22]   = image_feature_match(qImage,keyFrame1_2,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
        
        finalKeyFrame_2             = keyFrame1_2;      % Best Key Frame for Neuorn 2
        finalInlierCount_2          = iInlierCount_2;
        
        % if a neuron has more than 1 key frame, compute the best key frame
        iNoOfKeyFrames_2 = numel(neuron_KeyFrame(idx_2).frequency);
        if iNoOfKeyFrames_2 > 1
            for i   = 2:iNoOfKeyFrames_2
                KF                      = neuron_KeyFrame(idx_2).key_frames(i);
                finalInlierCount_n      = image_feature_match(qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
                
                %Comparision between KeyFrames
                if finalInlierCount_n > finalInlierCount_2
                    min_inlier  = finalInlierCount_2;
                else
                    min_inlier  = finalInlierCount_n;
                end
                
                error_22        = image_feature_match_mean(min_inlier,qImage,finalKeyFrame_2,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
                error_KF_n      = image_feature_match_mean(min_inlier,qImage,KF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
                
                if error_KF_n < error_22
                    finalKeyFrame_2     = KF;      % Best Key Frame for Neuorn 2
                    error_22            = error_KF_n;
                end
            end
        end
    catch
        %        fprintf('\n Warning: Not enough matching points, Key Frame from Neuron 2 (Next Closest Neuron skipped).');
    end
end

if iNoOfKeyFrames_2 > 0 && iNoOfKeyFrames_2 < 2
    %Comparision between KeyFrames of Neuron 1 and Neuron 2
    if finalInlierCount_2 > finalInlierCount_1
        min_inlier  = finalInlierCount_1;
    else
        min_inlier  = finalInlierCount_2;
    end
    error_11        = image_feature_match_mean(min_inlier,qImage,finalKeyFrame_1,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
    error_22        = image_feature_match_mean(min_inlier,qImage,finalKeyFrame_2,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
end

if error_22 > error_11
    bestKeyFrame  = finalKeyFrame_1;        % Winner Key Frame
    bestNeuron    = nr1;
else
    bestKeyFrame  = finalKeyFrame_2;        % Winner Key Frame
    bestNeuron    = nr2;
end

%%  Once we find Key Frames, we look through frames associated with the key frames
%%  and find the most similar/matching image
%%  -----------------------------------------------------------------------
% Winner Key Frame -- finalkeyKrame
length  = numel(keyFrame_Frame);
for i12 = 1:length
    if keyFrame_Frame(i12).key_frame == bestKeyFrame
        fprintf('\n\nComputing best matching frame...');
        frames_temp         = keyFrame_Frame(i12).frames;        
        frameBMU            = keyFrame_Frame(i12).frames(1);
        try
            iInlierCount_1      = image_feature_match(qImage,frameBMU,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
        catch
            continue
        end
        
        bestMatchingFrame       = frameBMU;      % Best Matching Frame
        finalInlierCount        = iInlierCount_1;
        
        iNoOfFrames             = numel(frames_temp);
        
        if iNoOfFrames > 1
            for i2 = 2:iNoOfFrames
                BF              = keyFrame_Frame(i12).frames(i2);
                try
                    iInlierCount_2  = image_feature_match(qImage,BF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType);
                catch
                    continue
                end
                %                finalInlierCount_2 = iInlierCount_2;
                
                if iInlierCount_2 > finalInlierCount
                    min_inlier  = finalInlierCount;
                else
                    min_inlier  = iInlierCount_2;
                end
                
                try
                    error_11        = image_feature_match_mean(min_inlier,qImage,bestMatchingFrame,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
                    error_22        = image_feature_match_mean(min_inlier,qImage,BF,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType); % calling the new function computing using mean
                catch
                    continue
                end
                
                if error_22 > error_11
                    bestMatchingFrame       = frameBMU;        % BMU
                    finalInlierCount    = iInlierCount_1;
                else
                    bestMatchingFrame       = BF;        % BMU
                    finalInlierCount    = iInlierCount_2;
                end
            end
        end
    end
end
end %end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODULE: image_feature_match_mean()
%% Function: Function to find error to find closest neuron
%% ------------------------------------------------------------------------
function [error] = image_feature_match_mean(min_inlier,qImage,dbImage,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)
%min_inlier - Count of minumum number of inliers as determined from the
%images
% qImage    - Query Image
% dbImage   - Image we have in DataBase
% locationKITTI - Location of images in system
min_inlier_count    = min_inlier;
image_Q             = image_name(qImage,location_data,dataset);
image_I             = image_name(dbImage,location_data,dataset);

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
[fMAT, inliers]     = estimateFundamentalMatrix(matchedPointsI,matchedPointsQ,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',0.01,'Confidence',99);
inliers_1           = find(inliers==1);
no_inliers          = numel(inliers_1);

if no_inliers < min_inlier_count
    min_inlier_count = no_inliers;
end


MPLocI          = matchedPointsI.Location(inliers_1,:); % Location of matched points in db Image
MPLocI_t        = MPLocI';                              % Transposing to find distance from the mean
mu              = mean(MPLocI);                         % Finding Mean of all the points

edist           = sqrt(sum((repmat(mu',1,size(MPLocI_t,2))-MPLocI_t).^2,1));  % finding distance from all the mean
[dist,dist_ind] = sort(edist);

new_inlier_ind  = dist_ind(1:min_inlier_count)'; % Finding indexes of matched points which are closest to mean as per min_inlier_count

MPLocI          = matchedPointsI.Location(new_inlier_ind,:);
MPLocQ          = matchedPointsQ.Location(new_inlier_ind,:);

inliers_MPLocI  = [MPLocI ones(size(MPLocI,1),1)]; % tarnsforming to the form [Xi Yi 1]
inliers_MPLocQ  = [MPLocQ ones(size(MPLocI,1),1)]; % tarnsforming to the form [Xi Yi 1]
j               = size(inliers_MPLocI,1);
Total_Error_1   = 0;

for i10=1:j
    %         LQ = MPLocI(i10,:);
    Q11             = inliers_MPLocI(i10,:);
    Q12             = inliers_MPLocQ(i10,:);    
    Error_1         = Q12*fMAT*Q11';    
    Total_Error_1   = Total_Error_1 + abs(Error_1);
end
error   = Total_Error_1;
end %end of function
