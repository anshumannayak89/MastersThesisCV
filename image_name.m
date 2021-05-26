%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODULE   : Image Naming
%% Function : IMAGE_NAME(); Function to resolve image number to image name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [strImageName] = image_name(image_no,location_data,dataset)

% Kitti Data
if dataset == 1 
    image_no = image_no - 1 ;
    zeros_str = '00000';
    if image_no > 999
        zeros_str = '00';
    elseif image_no > 99
        zeros_str = '000';
    elseif image_no > 9
        zeros_str = '0000';
    end    
    strImageName =  strcat(location_data,zeros_str,num2str(image_no),'.png');   % Generating full image name as per KITTI data set

% St.Lucia Data
elseif dataset == 2 
    image_no = image_no ;
    zeros_str = '0000';
    if image_no > 9999
        zeros_str = '';        
    elseif image_no > 999
        zeros_str = '0';
    elseif image_no > 99
        zeros_str = '00';
    elseif image_no > 9
        zeros_str = '000';
    end    
    strImageName =  strcat(location_data,'frame_',zeros_str,num2str(image_no),'.png');   % Generating full image name as per KITTI data set
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
