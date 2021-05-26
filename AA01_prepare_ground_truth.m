%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              COMPUT GROUND TRUTH FROM GPS DATA                         %
%%  Date     :  03 Mar 2020 (Fixed)                                       %
%%  Objective:  Computes GPS ground truth for a single or multiple images % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%clc
%clear
tic

%% variables
beta = 64;  % a constant

%% USER INPUTS
%% ------------------------------------------------------------------------
% Ask user to select desirable DATASET
fprintf('\n Available Datasets: ');
fprintf('\n 1. KITTI');
fprintf('\n 2. St. Lucia');

dataset = input('\n Please select the Dataset: ');
switch dataset
    case 1
        disp(' KITTI selected');
        fprintf('\nLoading KITTI Data...');
        load('/home/nayak/Desktop/NewFolder/WORKSPACE/DATASETS/kitti00_splr.mat');
        location_data = '/eda/media/visual_odometry/00/image_0/';          % Location of KITTI images in system
        OFFSET = 50;
        distLimit = 5; % in M                                          % Max distance between Query Image and Ground Truth
    case 2
        disp(' St.Lucia selected');
        fprintf('\nLoading St.Lucia Data...');
        load('/home/nayak/Desktop/NewFolder/WORKSPACE/DATASETS/stlucia_190809_0845_splr.mat');
        location_data = '/eda/datasets/stlucia_190809_0845/images/';       % Location of St.Lucia images in system
        OFFSET = 100;
        distLimit = 10.5; % in M                                         % Max distance between Query Image and Ground Truth
    otherwise
        disp(' Invalid choice! KITTI selected as default.');
        fprintf('\nLoading KITTI Data...');
        load('/home/nayak/Desktop/NewFolder/WORKSPACE/DATASETS/kitti00_splr.mat');
        location_data = '/eda/media/visual_odometry/00/image_0/';          % Location of KITTI images in system
        OFFSET = 50;
        distLimit = 5; % in M                                          % Max distance between Query Image and Ground Truth
end

% Uniforming structure of GPS data
if size(GPS,1) > 2 
    GPS = GPS';
end

%% QUERY IMAGE/S
%% ------------------------------------------------------------------------
% Asking user to enter query image
fprintf('\n Select 1. for a Single image');
fprintf('\n Select 2. for a Range of images');

imgRange = input('\n Please enter your choice: ');                         % Variable to hold user choice for single or range of image 
switch imgRange
    case 1
        qImage=input('\nEnter query image : ');
        qImageEnd = 0;                                                     % 0 as we are computing a single image
        globalTruth = compute_ground_truth_global(X_mN,GPS,qImage,qImageEnd,OFFSET,beta,distLimit);
        fprintf('\nQ-Image: %d \nGround Truth: %d \nDistance: %d\n',globalTruth(1).Q_Image,globalTruth(1).Match,globalTruth(1).Distance);        
    case 2
        qImageStart = input('\nEnter start image : ');
        qImageEnd = input('\nEnter Range (10,20,50 etc...): ');
        globalTruth = compute_ground_truth_global(X_mN,GPS,qImageStart,qImageEnd,OFFSET,beta,distLimit);
    otherwise
        fprintf('\n Single image selected by default');
        qImage=input('\nEnter query image : ');
        qImageEnd = 0;                                                     % 0 as we are computing a single image
        globalTruth = compute_ground_truth_global(X_mN,GPS,qImage,qImageEnd,OFFSET,beta,distLimit);
        fprintf('\nQ-Image: %d \nGround Truth: %d\n',globalTruth(1).Q_Image,globalTruth(1).Match);
end
toc