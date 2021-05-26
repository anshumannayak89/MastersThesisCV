%%  Function: BLOCK_FRAME() , Fuction to create constellation of continuion frames
%%  ------------------------------------------------------------------------
function constellation = block_frame(N,codebook)
    frameMat_mat = codebook(N).frame_mappings;                              % Matrix containing Frames 
    
    frameMat_diff_mat = diff(frameMat_mat);                                 % Calculate the difference between frames forrming a ... 
                                                                            % ...matrix like [1 1 1 3 1 1 6] etc                             
    frameMat_startIndex_temp_mat = find(frameMat_diff_mat > 1 ) +1 ;        % temp matrix to hold values 
    frameMat_startIndex_mat = [ones(1), frameMat_startIndex_temp_mat];      % Starting Index of consecutive frames
    
    frameMat_endIndex_temp_mat = find(frameMat_diff_mat > 1 );              % temp matrix to hold values
    frameMat_endIndex_mat = [frameMat_endIndex_temp_mat length(frameMat_mat)]; % End Index of consecutive Frames
    
    size_mat = (frameMat_endIndex_mat - frameMat_startIndex_mat) +1 ;       % Matrix to hold size of each constellation
    
    constellation = mat2cell(frameMat_mat,1,size_mat) ;                     % Converting frames into smaller constellations 