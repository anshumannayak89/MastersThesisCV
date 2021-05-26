%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME        : main_script.m                                            %%
% DESC        : Method Approach                                          %%
% Date Created: August 10, 2019                                          %%                                    
%             :                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clc
clear
fprintf('\nPreparing Trained Data...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VARIABLES
beta            = 64;                    % Controls Spatial Resolution of Network
keyFrame_Frame  = {};                    % Cell holding Key Frame and frames associated with it
bad_image       = [];
neuron_KeyFrame = struct('neuron',[],'key_frames',[]);
status          = 000 ;
count_RESULT    = 1;                     % Counter For Result Table
RESULT          = {};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   USER INPUTS
%   -----------
%   Select desirable DATASET (KITTI, St.Lucia...)
%%  -----------------------------------------------------------------------
fprintf('\nAvailable Datasets: ');
fprintf('\n 1. KITTI');
fprintf('\n 2. St. Lucia');
dataset = input('\n Please select the Dataset: ');
switch dataset
    case 1
        disp            (' KITTI selected');
        fprintf         ('\nLoading KITTI Data...');
        load            ('/home/nayak/Desktop/NewFolder/WORKSPACE/DATASETS/kitti00_splr.mat');
        location_data   = '/eda/media/visual_odometry/00/image_0/';         % Location of KITTI images in system
        OFFSET          = 50;                                               % No of images to ignore preceeding Query Image
        distLimit       = 5;                                                % Max distance between Query Image and Ground Truth (in Meter)
    case 2
        disp            (' St.Lucia selected');
        fprintf         ('\nLoading St.Lucia Data...');
        load            ('/home/nayak/Desktop/NewFolder/WORKSPACE/DATASETS/stlucia_190809_0845_splr.mat');
        location_data   = '/eda/datasets/stlucia_190809_0845/images/';      % Location of St.Lucia images in system
        OFFSET          = 100;                                              % No of images to ignore preceeding Query Image
        distLimit       = 10.5;                                               % Max distance between Query Image and Ground Truth (in Meter)
    otherwise
        disp            (' Invalid choice! KITTI selected as default.');
        fprintf         ('\nLoading KITTI Data...');
        load            ('/home/nayak/Desktop/NewFolder/WORKSPACE/DATASETS/kitti00_splr.mat');
        location_data   = '/eda/media/visual_odometry/00/image_0/';         % Location of KITTI images in system
        OFFSET          = 50;                                               % No of images to ignore preceeding Query Image
        distLimit       = 5;                                                % Max distance between Query Image and Ground Truth (in Meter)
end

%% Select desirable method to compute
% -------------------------------------------------------------
fprintf('\nAvailable methods: ');
fprintf('\n 1. Score method');
fprintf('\n 2. Mean method: Select Min-Inliers w.r.t Mean');
fprintf('\n 3. Random points method: Select K-random Inliers');
fprintf('\n 4. Convex Hull method: Total Area enclosed by Inliers \n');

method_type     = input('\n Please select the desired method: ');           % List of avilable methods
methods_list    = {'score','mean','randompts','convex'};                    % default: convex

switch method_type
    case 1
        disp(' Score method selected');
    case 2
        disp(' Mean method selected');
    case 3
        disp(' Random points method selected');
    case 4
        disp(' Convex hull method selected');
    otherwise
        disp(' Invalid choice! Convex hull method selected as default.');
end

%% Ask user to select desirable method to EXTRACT FEATURES (CORNOR DETECTORS) of KEY FRAME
% (DATE: 15-Apr-2020)
% ---------------------------------------------------------------------------------------
fprintf('\nAvailable methods for Key Frames: ');
fprintf('\n 1. BRISK');
fprintf('\n 2. FAST');
fprintf('\n 3. HARRIS');
fprintf('\n 4. Shi-Tomasi -"Good Features to Track"');
fprintf('\n 5. ORB');
detectFeatures_KF_methodType = input('\n Please select the desired method: ');
switch detectFeatures_KF_methodType
    case 1
        disp(' BRISK selected');
        detectFeatures_KF_Method = 'detectBRISKFeatures';
    case 2
        disp(' FAST selected');
        detectFeatures_KF_Method = 'detectFASTFeatures';
        
    case 3
        disp(' HARRIS selected');
        detectFeatures_KF_Method = 'detectHarrisFeatures';
    case 4
        disp(' Shi-Tomasi -"Good Features to Track" selected');
        detectFeatures_KF_Method = 'detectMinEigenFeatures';
        
    case 5
        disp(' ORB selected');
        detectFeatures_KF_Method = 'detectORBFeatures';        
    otherwise
        disp(' Invalid choice! FAST selected as default.');
        detectFeatures_KF_Method = 'detectFASTFeatures';
end
detectFeatures_KF_MethodType   = str2func(detectFeatures_KF_Method);

%% Select desirable method to DETECT FEATURES
%  --------------------------------------------------------------------
fprintf('\nAvailable methods for FEATURE DETECTOR: ');
fprintf('\n 1. ORB');
fprintf('\n 2. BRISK');
fprintf('\n 3. SURF');
fprintf('\n 4. Shi-Tomasi -"Good Features to Track"');
fprintf('\n 5. MSER');
fprintf('\n 6. HARRIS');
fprintf('\n 7. FAST');
fprintf('\n 8. KAZE');

featureDetector_methodType      = input('\n Please select the desired method: ');
switch featureDetector_methodType                                           %%  Feature Type   %% Scale Independent
    case 1
        disp(' ORB selected');                                      
        detectFeaturesMethod   = 'detectORBFeatures';                       %   Corner          %%  NO
    case 2
        disp(' BRISK selected');
        detectFeaturesMethod   = 'detectBRISKFeatures';                     %   Corner          %%  YES
    
    case 3
        disp(' SURF selected');
        detectFeaturesMethod   = 'detectSURFFeatures';                      %   Blob            %%  YES
    case 4
        disp(' Shi-Tomasi -"Good Features to Track" selected');
        detectFeaturesMethod   = 'detectMinEigenFeatures';                  %   Corner          %%  NO
    
    case 5
        disp(' MSER selected');                                             %   Region with     %%  YES 
        detectFeaturesMethod   = 'detectMSERFeatures';                      %   Uniform Intensity
    case 6
        disp(' HARRIS selected');
        detectFeaturesMethod   = 'detectHarrisFeatures';                    %   Corner          %%  NO

    case 7
        disp(' FAST selected');
        detectFeaturesMethod   = 'detectFASTFeatures';                      %   Corner          %%  NO
    case 8
        disp(' KAZE selected');
        detectFeaturesMethod   = 'detectKAZEFeatures';                      %   Blob            %%  NO
    
    otherwise
        disp(' Invalid choice! ORB selected as default.');
        detectFeaturesMethod   = 'detectORBFeatures';
end
featureDetectorMethodType = str2func(detectFeaturesMethod);

%% This section asks user to select desirable method to EXTRACT / DECRIPTOR FEATURES 
%  ---------------------------------------------------------------------
fprintf('\nAvailable methods for FEATURE DESCRIPTOR: ');
fprintf('\n 1. ORB');
fprintf('\n 2. BRISK');
fprintf('\n 3. SURF');
fprintf('\n 4. FREAK / corner Points selected');
fprintf('\n 5. MSER');
fprintf('\n 6. KAZE');
featureDescriptor_methodType      = input('\n Please select the desired method: ');
switch featureDescriptor_methodType
    case 1
        disp(' ORB selected');
        featuresDescriptorMethod   = 'ORBPoints';
    case 2
        disp(' BRISK selected');
        featuresDescriptorMethod   = 'BRISKPoints';
    case 3
        disp(' SURF selected');
        featuresDescriptorMethod   = 'SURFPoints';
    case 4
        disp(' FREAK / corner Points selected');
        featuresDescriptorMethod   = 'cornerPoints';
    case 5
        disp(' MSER selected');
        featuresDescriptorMethod   = 'MSERRegions';
    case 6
        disp(' KAZE selected');
        featuresDescriptorMethod   = 'KAZEPoints';
    otherwise
        disp(' Invalid choice! ORB selected as default.');
        featuresDescriptorMethod   = 'ORBPoints';
end
featureDescriptorMethodType = str2func(featuresDescriptorMethod);

%% This section asks user to enter the similarity threshold / distance
% between 2 neurons
% -------------------------------------------------------------------
distRatioThreshold  = input('\nEnter desired distance ratio. (0.5 - 0.95): ');
if distRatioThreshold >0.95 || distRatioThreshold < 0.5
    distRatioThreshold = 0.9;                                               % default value = 0.9 
    fprintf(' Improper distance ratio entered, defaulting to 0.9\n');
end

%% This sections asks user to enter one QUERY IMAGE or multiple QUERY IMAGES
% -------------------------------------------------------------------------
fprintf('\nSelect 1. for a Single image');
fprintf('\nSelect 2. for a Range of images');

imgRange = input('\n Please enter your choice: ');                          % Variable to hold user choice for single or range of image
switch imgRange
    case 1
        qImageStart =   input('\nEnter query image : ');                    % Accepts a single query image 
        qImageEnd   =   qImageStart;                                        % Sets the Start Image as last Image for loop
        endRange    =   1 ;                                                 % Keeps track of images to be evaluated
    case 2
        qImageStart =   input('\nEnter start image : ');                    % Accepts the beginning Image of series 
        endRange    =   input('\nEnter Range (10,20,50 etc...): ');         % Keeps track of total images to be evaluated
        qImageEnd   =   qImageStart + endRange;                             % Sets the end of the loop
end

%% Evaluation
%  ----------
%  Actual evaluation begins
%  ---------------------------------
for qImage = qImageStart:qImageEnd

    if rem(qImageEnd - qImage,10) == 0                                      % Informs user of remaining images to be evaluated
        fprintf('\nRemaining %d',qImageEnd - qImage);
    end
    
    status = 000;                                                           % Resetting status to 000 for a new evaluation
    
    % Loading appropriate DATASET as per user choice
    if dataset      == 1
        load('kitti00_splr.mat');
    elseif dataset  == 2
        load('stlucia_190809_0845_splr.mat')
    else
        load('kitti00_splr.mat');
    end
    
    % Uniforming structure of GPS data
    if size(GPS,1) > 2
        GPS = GPS';
    end
    
    % Preparing "weights" & "codebook" for evalutaion  
    [weights,codebook] = data_deletion(weights,codebook,location,qImage,OFFSET);

    %answer = questdlg("Do you want to visualize the image sequence?");
    answer  = 'no';
    isDraw  = strcmp(lower(answer),'yes');
    
    % Finding closest neurons (nr1 - closest)(nr2 - 2nd closest))
    % Function returns closest Nr, 2nd closest neuron as ratio of distance
    % --------------------------------------------------------------------
    [nr1,nr2,distRatio]     = nearest_neighbour(X_mN,qImage,weights,beta);
    nr2_temp                = nr2;                                          % keeping nr2 for possible later use 
    
    if distRatio > distRatioThreshold    
        neuron  = [nr1 nr2];
    else
        % fprintf('\n  Distance exceeds threshold value of %f, ignoring Neuron No: %d \n',distRatioThreshold,nr2);
        nr2     = 0;    % setting nr2 = 0 when its too far
        neuron  = nr1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Module: Consecutive frames and keyframes
    %% ------------------------------------------------------------------------
    % fprintf('\nComputuing consecutive frames and forming constellations ...');
    % fprintf('\nComputuing KEY FRAME...\n');
    count_KF_F = 0;                                                         % To keep track of keyFrame_Frame matrix/structure
    
    iNoOfNeurons = numel(neuron);                                           % Number of Neuron
    for i1  = 1: iNoOfNeurons
        i2  = neuron(i1);
        
        % Temporal segregation of frames (constellation) in neurons as per time stamp
        constellation   = block_frame(i2,codebook);
        noOfCell        = numel(constellation);                             % no of constellations
        
        % Creating a structurefor Neurons and its associated Keyframes
        % Contains 3 fields: neurons, key_frames & frequency (of key_frames)
        % ------------------------------------------------------------------
        neuron_KeyFrame(i1).neuron      = i2;                               % entering neuron numbers to neuron_KeyFrame
        neuron_KeyFrame(i1).key_frames  = [];
        
        for i3 = 1:noOfCell
            iLengthOfCell = numel(constellation{i3});                       % length of individual constellation
            
            % If 1 frame by default its the Key Frame
            % ---------------------------------------
            if iLengthOfCell < 2                
                count_KF_F                      = count_KF_F + 1;
                iKeyFrame                       = constellation{i3}(1);
                neuron_KeyFrame(i1).key_frames  = [neuron_KeyFrame(i1).key_frames iKeyFrame];
                
                % Forming structure keyFrame_Frame
                % Contains 4 fields:
                % --------------------------------
                keyFrame_Frame(count_KF_F).key_frame     = iKeyFrame;       % key_frame
                keyFrame_Frame(count_KF_F).neuron        = i2;              % neuron (associated)
                keyFrame_Frame(count_KF_F).frames        = iKeyFrame;       % frames (associated)
                keyFrame_Frame(count_KF_F).frequency     = numel(keyFrame_Frame(count_KF_F).frames);  % Number of associated frames
                
            % If multiple frames then we COMPUTE the KEY FRAME
            % -------------------------------------------------
            else
                count_KF_F = count_KF_F + 1;
                try
                    [iKeyFrame,indexKeyFrame,bad_images_1]  = block_keyframe(constellation{i3},location_data,isDraw,dataset,detectFeatures_KF_MethodType); % returns the key frame and key_frame index
                    bad_image                               = [bad_image bad_images_1];     % stores images inelegible to be Key Frames
                    
                    if isempty(iKeyFrame)
                        errID   = 'MATLAB:badKF';
                        msgtext = 'No eligible KFs in Neuron %d: ';
                        ME      = MException(errID,msgtext,i2);
                        throw(ME)
                    end
                    
                catch ME
                    disp('Error Message:')
                    disp(ME.message)
                end
                
                neuron_KeyFrame(i1).key_frames          = [neuron_KeyFrame(i1).key_frames iKeyFrame];
                
                % creating connection between Key Frame and images
                keyFrame_Frame(count_KF_F).key_frame    = iKeyFrame;
                keyFrame_Frame(count_KF_F).neuron       = i2;
                keyFrame_Frame(count_KF_F).frames       = [];
                
                for i5 = 1:iLengthOfCell
                    keyFrame_Frame(count_KF_F).frames   = [keyFrame_Frame(count_KF_F).frames constellation{i3}(i5)];
                    strKeyFrameName                     = strcat("FID-",num2str(iKeyFrame));
                    strFrameName                        = strcat("FID-",num2str(constellation{i3}(i5)));     %Frame_IDs to string and Naming frames as FId-1, FId-2
                end
                keyFrame_Frame(count_KF_F).frequency    = numel(keyFrame_Frame(count_KF_F).frames);
            end
            neuron_KeyFrame(i1).frequency               = numel(neuron_KeyFrame(i1).key_frames);            % Updates the no of Key Frames in a Neuron
        end
    end
    % By this stage we found the Key Frames of Neurons 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Evaluate Matching Frame for Query Image
% ---------------------------------------    
    % Execute following function as per user choice
    % ---------------------------------------------
    try
        switch (method_type)
            
            case 1
                % score between: 1. No if Inliers, 2. Error, and 3. Distance of farthest point from mean
                [bestNeuron,bestKeyFrame,bestMatchingFrame]     = scoring_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,featureDetectorMethodType,featureDescriptorMethodType);
                
            case 2
                [bestNeuron,bestKeyFrame,bestMatchingFrame]     = mean_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,featureDetectorMethodType,featureDescriptorMethodType);
                
            case 3
                [bestNeuron,bestKeyFrame,bestMatchingFrame]     = random_points_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,featureDetectorMethodType,featureDescriptorMethodType);
                
            otherwise
                [bestNeuron,bestKeyFrame,bestMatchingFrame]     = convex_hull_method(qImage,location_data,neuron_KeyFrame,keyFrame_Frame,nr1,nr2,nr2_temp,dataset,featureDetectorMethodType,featureDescriptorMethodType);
        end
        
    catch ME
        if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
            fprintf('\nNo Match Found for Frame : %d',qImage);
            
        elseif strcmp(ME.identifier,'vision:points:notEnoughMatchedPts')
            %matchedPoints1 and matchedPoints2 do not have enough points.
            %The number of points in each set must be at least 8.
            fprintf('\nNot enough Matched Points.');
            fprintf('\nNo Match Found for Frame : %d',qImage);
            
        elseif strcmp(ME.identifier,'vision:points:notEnoughInlierMatches')
            %Could not find enough inliers in matchedPoints1 and matchedPoints2.            
            fprintf('\nNot enough Inliers.');
            fprintf('\nNo Match Found for Frame : %d',qImage);
            
        elseif strcmp(ME.identifier,'MATLAB:unassignedOutputs')
            fprintf('\nNo Winner Neuron');
            fprintf('\nNo Match Found for Frame : %d',qImage);
            
        else
            fprintf('\nERROR in Frame : %d\n',qImage);
            disp('Error Message:')
            disp(ME.message)
        end
        nr2     = nr2_temp;     % In case of an exception skipped NR 2 is again considered
        status  = 404 ;         % If no match found, status is changes 404 aka Match not found
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate Ground Truth
% ----------------------
    % Get index of closest image as per GPS data (ground truth)
    % ---------------------------------------------------------
       
    globalTruth = compute_ground_truth_global(X_mN,GPS,qImage,0,OFFSET,beta,distLimit); % 0 as we are sending individual images
    gTruth      = globalTruth(1).Match;
    
    dist_Q_G    = norm (GPS(:,qImage) - GPS(:,gTruth));                   % compute distance between Q Image and Ground truth
    
    % if we find a match / Visited Place
    % ----------------------------------
    if status ~= 404
        % Distance between Query Image and Winner frame
        dist_Q_W = norm (GPS(:,qImage) - GPS(:,bestMatchingFrame));
        
        % Distance between Ground Truth and Winner frame2
        dist_G_W = norm (GPS(:,gTruth) - GPS(:,bestMatchingFrame));
        
    % If no match found/ Unvisited Place
    % ----------------------------------
    else
        bestMatchingFrame = 0000 ;  % We may not need these values
        dist_G_W = 9999;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classifying Positive and Negatives 
% ----------------------------------
    if status ~= 404 && dist_Q_W <= distLimit       % True Positive     :Match (Yes) & distace < 5m /10.5m   
        visited    = 'YES';
        match      = bestMatchingFrame;
        distance   = dist_Q_W;
        apr        = 'TP';
        
    elseif status ~=404 && dist_Q_W > distLimit    % False Positive    :Match (Yes) & distance > 5m/10.5m
        visited    = 'YES';
        match      = bestMatchingFrame; % bestMatchingFrame
        distance   = dist_Q_W;
        apr        = 'FP';
       
    elseif status == 404 && dist_Q_G <= distLimit   % False Negative    :Match (NO) & distance < 5m/10.5m
        visited    = 'NO';
        match      = 'NONE';
        distance   = dist_Q_G;
        apr        = 'FN';
        
    else
        %if status == 404 && dist_Q_G > distLimit   % True Negative     :Match (NO) & distance > 5m /10.5m
        visited    = 'NO';
        match      = 'NONE';
        distance   = dist_Q_G;
        apr        = 'TN';
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TABULATING RESULT
% --------------------------
    RESULT(count_RESULT).Q_Image    = qImage;
    RESULT(count_RESULT).Visited    = visited;
    RESULT(count_RESULT).Match      = match;
    RESULT(count_RESULT).G_Truth    = gTruth;
    RESULT(count_RESULT).Distance   = distance;
    RESULT(count_RESULT).Status     = status;
    RESULT(count_RESULT).APR        = apr;
   
    count_RESULT                    = count_RESULT +1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output DISPLAY
%  Display output image only if 1 image is selected for evaluation
% --------------------------
    if imgRange < 2             
        if status < 1
            imageQuery      = image_name(qImage,location_data,dataset);
            imageFrameBMU   = image_name(bestMatchingFrame,location_data,dataset);
            imageGTruth     = image_name(gTruth,location_data,dataset);
            
            subplot(2,2,1);     %Top Left
            imshow(imageQuery);
            QI              = strcat('Query Image: ',num2str(qImage));
            title(QI);
            
            subplot(2,2,3);     %Bottom Left
            BMF             = strcat('Best Matching Frame: ',num2str(bestMatchingFrame));
            imshow(imageFrameBMU);
            title(BMF);
            
            subplot(2,2,4);     %Bottom Right
            GT              = strcat('Ground Truth: ',num2str(gTruth));
            imshow(imageGTruth);
            title(GT);
            
            %%building tree / heirarchichal structure
            build_tree(qImage,nr1,nr2,bestNeuron,bestKeyFrame,bestMatchingFrame,gTruth,neuron_KeyFrame,keyFrame_Frame);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %End of Actual Evaluation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Accuracy Precision and Recall
% -----------------------------------------
true_Positive   = numel(find(strcmp({RESULT.APR}, 'TP')==1));
true_Negative   = numel(find(strcmp({RESULT.APR}, 'TN')==1));
false_Positve   = numel(find(strcmp({RESULT.APR}, 'FP')==1));
false_Negative  = numel(find(strcmp({RESULT.APR}, 'FN')==1));

accuracy        = (true_Positive + true_Negative) / (true_Positive + true_Negative + false_Positve + false_Negative);
recall          = true_Positive / (true_Positive + false_Negative);
precision       = true_Positive / (true_Positive + false_Positve);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output DISPLAY
%% --------------------------
fprintf('\nStart                    : %d',  qImageStart);
fprintf('\nTotal Images             : %f',  endRange);
fprintf('\nMethod                   : %d',  method_type);
fprintf('\nDist Ratio               : %f',  distRatioThreshold);
fprintf('\nKF-Feature DETECTOR      : ');   disp(detectFeatures_KF_Method);
fprintf('IP-Feature DETECTOR      : ');   disp(detectFeaturesMethod);
fprintf('IP-Feature EXTRACTOR     : ');   disp(featuresDescriptorMethod);

fprintf('Re-Visited Loc           : %d\n',numel(find(strcmp({RESULT.Visited}, 'YES')==1)));
fprintf('Accuracy                 : %d',  accuracy);
fprintf('\nRecall                   : %d',  recall);
fprintf('\nPrecision                : %d\n',  precision);

toc