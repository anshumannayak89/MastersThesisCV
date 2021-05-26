function [nr1,nr2,distRatio] = nearest_neighbour(X_mN,qImage,weights,beta)
%fprintf('\nFetching GIST features from X_mN');
qImage_feature = X_mN(:,qImage);  % gives GIST vector

% Determining winner Neuron [c=argmin||g-Wi||^²] Eq-4
% Simultaneous Place Learning and Recognition (S.M. Ali Musa Kazmi)
%beta = 64 ;                 % Controls Spatial Resolution of Network
distance = exp(-beta*sum(bsxfun(@minus,weights,qImage_feature).^2,1));

[sortDist,index] = sort(distance,'descend');

nr1 = index(1);         % closest neuron
nr2 = index(2);         % 2nd closest neuron

% Ratio between 2nd closed neuron to closest neuron
distRatio = sortDist(2)/sortDist(1);