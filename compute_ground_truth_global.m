%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module to compute GROUND TRUTH using GLOBAL features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function globalTruth = compute_ground_truth_global(X_mN,GPS,qImageStart,qImageEnd,OFFSET,beta,distLimit)
globalTruth = {};
count = 1;

if qImageEnd == 0
    % For a Single Image
    img_Feature     = X_mN(:,qImageStart);                  % gives GIST vector
    X_mN_new        = X_mN(:,[1:qImageStart-OFFSET]);
    similarity      = exp(-beta*sum(bsxfun(@minus,X_mN_new,img_Feature).^2,1));
    [sortSim,index] = sort(similarity,'descend');
    ground_truth    = index(1);                             % closest image as per global features

    
    dist_Q_G        = norm (GPS(:,qImageStart) - GPS(:,ground_truth));   % compute distance between Q Image and Ground truth
    globalTruth(count).Q_Image = qImageStart;
    
%     if dist_Q_G <= distLimit
        globalTruth(count).Match    = ground_truth;
        globalTruth(count).Distance = dist_Q_G;

   
else
    % For a Range of Image
    for i = qImageStart:(qImageStart+qImageEnd-1)           % Start Index : End Index
        img_Feature     = X_mN(:,i);                        % gives GIST vector
        X_mN_new        = X_mN(:,[1:i-OFFSET]);
        similarity      = exp(-beta*sum(bsxfun(@minus,X_mN_new,img_Feature).^2,1));
        [sortSim,index] = sort(similarity,'descend');
        ground_truth    = index(1);                         % closest image as per global features
        
        dist_Q_G        = norm (GPS(:,i) - GPS(:,ground_truth));   % compute distance between Q Image and Ground truth
        
        globalTruth(count).Q_Image = i;
        
        if dist_Q_G <= distLimit
            globalTruth(count).Match    = ground_truth;
            globalTruth(count).Distance = dist_Q_G;

% We are finding ground truth using global feature , i.e. what is the
% closest image as per global feature
        else
            globalTruth(count).Match    = 'NONE';
            globalTruth(count).Distance = dist_Q_G;
        end
        count = count +1;
    end
end
%fprintf('\nComputed ground truth as per Global features. Result stored in globalTruth \n');
end