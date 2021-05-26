%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MODEULE: Key Frame
%%  Function: BLOCK_KEYFRAME(); Function to compute key frames and key_frame index
%% ------------------------------------------------------------------------
function [iKeyFrame,indexKeyFrame,bad_images] = block_keyframe(cell_block,location_data,isDraw,dataset,extractFeatures_KF_MethodType)
bad_images          = [];
feature_final       = [];
temp_index_tracking = [];

iLengthOfCell   = numel(cell_block);                                        % length of individual constellation (cell_block = constellation{i3})
for i4 = 1:iLengthOfCell
    image   = cell_block(i4);
    image_1 = image_name(image,location_data,dataset);
    
    I1 = imread(image_1);
    
    if isDraw
        imshow(I1);
    end
    
    % Before we begin resizing
    % Step 1 : convert RGB 2 GRAY
    if size(I1,3) > 1
        I1  = rgb2gray(I1);
    end
    
    % Step 2 : Resize while maintaining aspect ratio with 256
    % Step 2a: Converting portratit images to landscape
    if size(I1,1) > size(I1,2)
        I1  = I1';
    end
    % Step 2b: Determine aspect ratio
    aspRatio = size(I1,2)/size(I1,1);
    
    % Step 2c: Resizing image to maintian aspect ratio of original image w.r.t 256
    I1_resize = imresize(I1,[256 aspRatio*256],'method','bilinear');
    
    % Step 3 : Determine Crop size
    cropSize    = mod(size(I1_resize,2),256);
    Lcrop       = round(cropSize/2);
    Rcrop       = round(cropSize/2);
    
    %
    cropLimit = fix(size(I1_resize,2)/256);
    if cropLimit    == 1
        cropWindow  = 256;
    else
        cropWindow  = 256 * (cropLimit - 1);
    end
    
    feature_Points  = extractFeatures_KF_MethodType(I1_resize);      % Detecting features to determine good candidates for Keyframes
    location_Points = feature_Points.Count;                          % gives actual number of features detected
    
    if location_Points > 300 %% FIND A SUITABLE THRESHOLD
        temp_index_tracking = [temp_index_tracking cell_block(i4)];
        I1_cropped          =  imcrop(I1_resize,[Lcrop,0,size(I1,2)-Rcrop,256]); % cropping image to aspect ratio 768x256
        feature             = [];
        
        for x1 = 0:256:cropWindow            
            I1_crop_temp            = imcrop(I1_cropped,[x1,0,256,256]);    % HOG detection window 256x256
            features_I1_crop_temp   = extractHOGFeatures(I1_crop_temp,'CellSize',[64 64],'NumBins',12,'BlockSize',[1 1],'UseSignedOrientation',true);
            feature                 = [feature features_I1_crop_temp];
        end
    else
        %            fprintf('\n NOT ENOUGH FAST FEATURES for Frame id : %d ', image);
        bad_images = [bad_images image-1];
        continue
    end
    feature_T       = transpose(feature);
    feature_final   = [feature_final feature_T];
end

%% changes as per discussion with Musa
mu                  = mean(feature_final,2);
dist2               = sum((feature_final - mu).^2,1);
[M,indexKeyFrame]   = min(dist2);
iKeyFrame           = temp_index_tracking(indexKeyFrame);
%end %end of function