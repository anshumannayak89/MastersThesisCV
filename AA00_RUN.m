%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME        : AA00_RUN.m                                               %%
% DESC        : To run a batch of Simulation                             %%
%             :                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
tic
%count_RESULT_ALL    = 0;        % Counter for RESULT_ALL Table
RESULT_ALL          = {}; 

%% Description of Variables used / passed on to AA00_main_script_func.m
% We pass the following variables to run a batch of simulation
% dataset                         % number 1,2
% method_type                     % number 1 2 3 4
% extractFeatures_methodType      % number 1 2 3 4
% extractFeatures_KF_methodType   % number 1 2 3 4 
% distRatioThreshold              % 0.9,0.8,0.75
% imgRange                        % My case always 2 (Range of Images)
% qImageStart                     % Starting Image
% endRange                        % Number of Images 


%% USER INPUTS
%% -----------
% This section asks user to select desirable DATASET (KITTI, St.Lucia...)
% -----------------------------------------------------------------------
fprintf('\n Available Datasets: ');
fprintf('\n 1. KITTI');
fprintf('\n 2. St. Lucia');
dataset = input('\n Please select the Dataset: ');

% This sections asks user to select desirable method to compute
% -------------------------------------------------------------
fprintf('\n Available methods: ');
fprintf('\n 1. Score method');
fprintf('\n 2. Mean method: Select Min-Inliers w.r.t Mean');
fprintf('\n 3. Random points method: Select K-random Inliers');
fprintf('\n 4. Convex Hull method: Total Area enclosed by Inliers \n');
method_type = input('\n Please select the desired method: ');

%distLimit = input('\nEnter GPS accuracy (5-20 meters) : ');

distRatioThreshold = input('\nEnter desired distance ratio. (0.5 - 0.95): ');
imgRange    =   2 ;                                                         % as we always use this to compute multiple images multiple times
qImageStart =   input('\nEnter start image : ');
endRange    =   input('\nEnter Range (10,20,50 etc...): ');
%distRatio_mat   =   [0.9 0.8 0.75];
distLimit_mat   =   [10.5];%12.5 15 25 50];

parfor i = 1:5
    RESULT_ALL_test  = {};
    count_RESULT_ALL = 0;
    %Selects which cornor detectors to use to compute KEY FRAMES
    %BRISK, FAST, HARRIS, MIN EIGEN /Shi-Tomasi
    detectFeatures_KF_methodType = i;
   
    for j = 7:7
        %Selects which feature Detector to use to MATCH IMAGES
        %ORB, BRISK, SURF, MIN EIGEN / Shi-Tomasi 
        featureDetector_methodType = j;
        
        for k = 1:6
            %Selects which feature extractor to use to MATCH IMAGES
            %ORB, BRISK, SURF
            featureDescriptor_methodType = k;
        
            for l = 1:numel(distLimit_mat)
                count_RESULT_ALL    = count_RESULT_ALL + 1;
    %           distRatioThreshold  = distRatio_mat(l);    % Passing  
                distLimit           = distLimit_mat(l);     %

                % Calling function AA00_main_script_func 
                [reVisitedLoc,accuracy,recall,precision] = AA00_main_script_func(dataset,method_type,detectFeatures_KF_methodType,featureDetector_methodType,featureDescriptor_methodType,distRatioThreshold,imgRange,qImageStart,endRange,distLimit);

                if detectFeatures_KF_methodType       == 1
                    KF_methodType   =   'BRISK';
                elseif detectFeatures_KF_methodType   == 2
                    KF_methodType   =   'FAST';

                elseif detectFeatures_KF_methodType   == 3
                    KF_methodType   =   'HARRIS';
                elseif detectFeatures_KF_methodType   == 4
                    KF_methodType   =   'Shi-Tomasi';  
                    
                elseif detectFeatures_KF_methodType   == 5
                    KF_methodType   =   'ORB';  
                end
                
                
                if featureDetector_methodType       == 1
                    fDetector_methodType   =   'ORB';
                elseif featureDetector_methodType   == 2
                    fDetector_methodType   =   'BRISK';

                elseif featureDetector_methodType   == 3
                    fDetector_methodType   =   'SURF';
                elseif featureDetector_methodType   == 4
                    fDetector_methodType   =   'Shi-Tomasi';
                    
                elseif featureDetector_methodType   == 5
                    fDetector_methodType   =   'MSER';
                elseif featureDetector_methodType   == 6
                    fDetector_methodType   =   'HARRIS';
                    
                elseif featureDetector_methodType   == 7
                    fDetector_methodType   =   'FAST';
                elseif featureDetector_methodType   == 8
                    fDetector_methodType   =   'KAZE';  
                end


                if featureDescriptor_methodType        == 1
                    fDescriptor_methodType   =   'ORB';
                elseif featureDescriptor_methodType    == 2
                    fDescriptor_methodType   =   'BRISK';

                elseif featureDescriptor_methodType    == 3
                    fDescriptor_methodType   =   'SURF';
                elseif featureDescriptor_methodType    == 4
                    fDescriptor_methodType   =   'FREAK/CornerPoints';
                
                elseif featureDescriptor_methodType    == 5
                    fDescriptor_methodType   =   'MSER';
                elseif featureDescriptor_methodType    == 6
                    fDescriptor_methodType   =   'KAZE';    
                    
                end


                RESULT_ALL_test(count_RESULT_ALL).Method        = method_type;
                RESULT_ALL_test(count_RESULT_ALL).DistRatio     = distRatioThreshold;
                RESULT_ALL_test(count_RESULT_ALL).GPS_distance  = distLimit;
                RESULT_ALL_test(count_RESULT_ALL).KF            = KF_methodType;
                RESULT_ALL_test(count_RESULT_ALL).Detector      = fDetector_methodType;
                RESULT_ALL_test(count_RESULT_ALL).Descriptor    = fDescriptor_methodType;
                RESULT_ALL_test(count_RESULT_ALL).Revisit       = reVisitedLoc;
                RESULT_ALL_test(count_RESULT_ALL).Accuracy      = accuracy;
                RESULT_ALL_test(count_RESULT_ALL).Recall        = recall;
                RESULT_ALL_test(count_RESULT_ALL).Precision     = precision;
            end
        end
    end
    save_results_func(RESULT_ALL_test,i);
end
toc

function save_results_func(RESULT_ALL_test,i)
    eval(['save RESULT_ALL_test',int2str(i),'.mat']) 
end
        