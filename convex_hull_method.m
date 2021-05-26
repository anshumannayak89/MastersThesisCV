%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In this function we compute KF/ Best matching frame
%% based on AREA enclosed by inliers
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bestNeuron,bestKeyFrame,bestMatchingFrame] = convex_hull_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,featureDetectorMethodType,featuresDescriptorMethodType)
area_KF_Nr_1 = 0;
area_KF_Nr_2 = 0;

% Analysis of 1st Closest Neuron / Closest Neuron
idx_1   = find( cellfun(@(x)isequal(x,nr1),{neuron_KeyFrame.neuron}) );

try
    keyFrame1_1         = neuron_KeyFrame(idx_1).key_frames(1);
    [area_KF_1]         = image_feature_match_convex_hull(qImage,keyFrame1_1,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType);    
    finalKeyFrame_1     = keyFrame1_1;      % Best Key Frame for Neuorn 1
    area_KF_Nr_1        = area_KF_1;
    
    iNoOfKeyFrames_1    = numel(neuron_KeyFrame(idx_1).frequency);
    if iNoOfKeyFrames_1 > 1
        for i = 2:iNoOfKeyFrames_1
            KF          = neuron_KeyFrame(idx_1).key_frames(i);
            [area_KF_n] = image_feature_match_convex_hull(qImage,keyFrame1_1,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType);
            
            if area_KF_1 < area_KF_n
                finalKeyFrame_1     = KF;      % Best Key Frame for Neuorn 1
                area_KF_Nr_1        = area_KF_n;
            end
        end
    end
catch
    nr2 = nr2_temp;     % If we do not find a Key Frame in closest Neuron we will reconsider 2nd closest neruron
end

% Analysis of 2nd Closest Neuron
if nr2 > 0
    idx_2   = find( cellfun(@(x)isequal(x,nr2),{neuron_KeyFrame.neuron}) );
    try
        keyFrame1_2         = neuron_KeyFrame(idx_2).key_frames(1);
        [area_KF_2]         = image_feature_match_convex_hull(qImage,keyFrame1_2,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType);        
        finalKeyFrame_2     = keyFrame1_2;      % Best Key Frame for Neuorn 2
        area_KF_Nr_2        = area_KF_2;
        
        iNoOfKeyFrames_2    = numel(neuron_KeyFrame(idx_2).frequency);
        if iNoOfKeyFrames_2 > 1
            for i = 2:iNoOfKeyFrames_2
                KF              = neuron_KeyFrame(idx_2).key_frames(i);
                [area_KF_n]     = image_feature_match_convex_hull(qImage,KF,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType);
                
                if area_KF_Nr_2 < area_KF_n
                    finalKeyFrame_2     = KF;      % Best Key Frame for Neuorn 2
                    area_KF_Nr_2        = area_KF_n;
                end
            end
        end
    catch
        %             fprintf('\n Warning: Not enough matching points, Key Frame from Neuron 2 (Next Closest Neuron skipped).');
    end
end

if area_KF_Nr_1 < 5000 && area_KF_Nr_2 < 5000
    %        NO MATCH FOUND
    return
else
    if area_KF_Nr_1 > area_KF_Nr_2
        bestKeyFrame  = finalKeyFrame_1;        % Winner Key Frame
        bestNeuron    = nr1;
    else
        bestKeyFrame  = finalKeyFrame_2;        % Winner Key Frame
        bestNeuron    = nr2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MODEULE : Frame Association / Geometric Validation
%%  Once we find Key Frames, we look through frames associated with the key frames
%%  and find the most similar/matching image
%%  -----------------------------------------------------------------------
% Winner Key Frame -- finalkeyKrame
length = numel(keyFrame_Frame);

for i = 1:length
    if keyFrame_Frame(i).key_frame == bestKeyFrame
        %            fprintf('\n\nComputing best matching frame...');
        frames_temp         = keyFrame_Frame(i).frames;        
        frameBMU            = keyFrame_Frame(i).frames(1);
        try
            [area_BMU_1]        = image_feature_match_convex_hull(qImage,frameBMU,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType);        
        catch
        end
        bestMatchingFrame   = frameBMU;      % Best Matching Frame
        final_area_BMU      = area_BMU_1;
        
        iNoOfFrames     = numel(frames_temp);
        if iNoOfFrames > 1
            for i2 = 2:iNoOfFrames
                BF              = keyFrame_Frame(i).frames(i2);
                try
                    [area_BMU_n]    = image_feature_match_convex_hull(qImage,BF,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType);
                %    fprintf('\n Area: %d ',area_BMU_n);
                catch
                    area_BMU_n = 0;
                end
                
                if area_BMU_n > final_area_BMU
                    bestMatchingFrame   = BF;           % Best Frame
                    final_area_BMU  = area_BMU_n;
                end
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [area] = image_feature_match_convex_hull(qImage,dbImage,location_data,dataset,featureDetectorMethodType,featuresDescriptorMethodType)
% qImage    - Query Image
% dbImage   - Image we have in DataBase
% locationKITTI - Location of images in system

image_Q     = image_name(qImage,location_data,dataset);
image_I     = image_name(dbImage,location_data,dataset);

% Read images.
if size(imread(image_Q),3) > 1              % RGB images converted to Grayscale
    IQgray  = rgb2gray(imread(image_Q));
    Igray   = rgb2gray(imread(image_I));
else
    IQgray  = imread(image_Q);
    Igray   = imread(image_I);
end

% Detect feature points each image.
imagePointsQ    = featureDetectorMethodType(IQgray);
imagePointsI    = featureDetectorMethodType(Igray);

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
[fMAT, inliers]     = estimateFundamentalMatrix(matchedPointsI,matchedPointsQ,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',0.1,'Confidence',99);
inliers_1           = find(inliers==1);
noOfInliers         = numel(inliers_1);
%disp(noOfInliers);
MPLocI              = matchedPointsI.Location(inliers_1,:); % Location of matched points in db Image

data                = [MPLocI(:,1) MPLocI(:,2)];
data                = double(data);
[k,av]              = convhull(data);
area                = av;
end %end of function