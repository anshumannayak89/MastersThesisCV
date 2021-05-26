%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODULE: image_feature_match()  
%% Function: Function to find inliers to find closest neuron
%% ------------------------------------------------------------------------
function [inlier,error] = image_feature_match(qImage,dbImage,location_data,dataset,extractFeaturesMethodType,featuresDescriptorMethodType)
% qImage    - Query Image
% dbImage   - Image we have in DataBase
% locationKITTI - Location of images in system
    
    image_Q = image_name(qImage,location_data,dataset); 
    image_I = image_name(dbImage,location_data,dataset);

    % Read images.
    if size(imread(image_Q),3) > 1
        IQgray = rgb2gray(imread(image_Q));
        Igray = rgb2gray(imread(image_I));
    else
        IQgray = imread(image_Q);
        Igray = imread(image_I);
    end
    
    % Detect feature points each image.
    imagePointsQ = extractFeaturesMethodType(IQgray);
    imagePointsI = extractFeaturesMethodType(Igray);
    
    % Detected feature points converted to Descriptor points.
    imagePointsQ    = featuresDescriptorMethodType(imagePointsQ.Location);
    imagePointsI    = featuresDescriptorMethodType(imagePointsI.Location);
    
    % Extract feature descriptors from each image.
    featuresQ = extractFeatures(IQgray,imagePointsQ,'Method','AUTO');
    featuresI = extractFeatures(Igray,imagePointsI,'Method','AUTO');

    % Match features across the dbImage and Querry Image .
    indexPairs_QI = matchFeatures(featuresQ,featuresI);
    
%  %   size(indexPairs_QI,1);
%     
%     if size(indexPairs_QI,1) >= 8
% %        disp([': ',num2str(qImage),' ',num2str(dbImage)]);
%     else
% %        disp(['ERR: ',num2str(qImage),' ',num2str(dbImage)]);
%     end
        
    matchedPointsQ = imagePointsQ(indexPairs_QI(:,1),:);
    matchedPointsI = imagePointsI(indexPairs_QI(:,2),:);
    
    rng('shuffle');
    % Find INLIERS KF and Querry Image
    [fMAT, inliers] = estimateFundamentalMatrix(matchedPointsI,matchedPointsQ,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',0.01,'Confidence',99);
    inliers_1= find(inliers==1);
    noOfInliers = numel(inliers_1);
    
    MPLocI = matchedPointsI.Location(inliers_1,:); % Location of matched points in db Image
    MPLocQ = matchedPointsQ.Location(inliers_1,:); % Location of matched points in Query Image

    mu = mean(MPLocI);
   
    inliers_MPLocI = [MPLocI ones(size(MPLocI,1),1)]; % tarnsforming to the form [Xi Yi 1]
    inliers_MPLocQ = [MPLocQ ones(size(MPLocI,1),1)]; % tarnsforming to the form [Xi Yi 1]
    j=size(inliers_MPLocI,1); 
    Total_Error_1 = 0;
    max_distance = 0;

    for i10=1:j        
        LQ = MPLocI(i10,:);
        Q11=inliers_MPLocI(i10,:);
        Q12=inliers_MPLocQ(i10,:);
        
        w = norm(mu-LQ);
        Error_1 = 1/w * Q12*fMAT*Q11';
        if(abs(w) > max_distance)
            max_distance = abs(w);
        end
        Total_Error_1 = Total_Error_1 + abs(Error_1);
    end

    error = Total_Error_1;
    inlier = noOfInliers;
end %end of function