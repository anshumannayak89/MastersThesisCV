%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULE: Build Tree 
% function to build heiarchichal structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function build_tree(qImage,nr1,nr2,winnerNeuron,winnerKeyFrame,finalFrameBMU,gTruth,neuron_KeyFrame,keyFrame_Frame)
closestNeuron_1 = find( cellfun(@(x)isequal(x,nr1),{neuron_KeyFrame.neuron}) );
closestNeuron_2 = find( cellfun(@(x)isequal(x,nr2),{neuron_KeyFrame.neuron}) );

% Building the base of tree by adding Query Image and 2 closest neuron
if nr2>0
    tree_P = [strcat("QFID: ",num2str(qImage)) strcat("QFID: ",num2str(qImage)) ];
    tree_C = [strcat("NID: ",num2str(nr1)) strcat("NID: ",num2str(nr2))];
else
    tree_P = [strcat("QFID: ",num2str(qImage))];% strcat("QFID: ",num2str(qImage)) ];
    tree_C = [strcat("NID: ",num2str(nr1))];% strcat("NID: ",num2str(nr2))];
end

% to define matche from Global features
tree_Path_global = [strcat("QFID: ",num2str(qImage)) strcat("NID: ",num2str(nr1))];     % Global fetaures as computed using weights & X_mN
% to define matche from local features
if winnerKeyFrame == finalFrameBMU
    tree_Path_local = [strcat("NID: ",num2str(winnerNeuron)) strcat("KFID: ", ...
                                  num2str(winnerKeyFrame))];
else
    tree_Path_local = [strcat("NID: ",num2str(winnerNeuron)) strcat("KFID: ", ...
                                  num2str(winnerKeyFrame)) strcat("KFID: ", ...
                                  num2str(winnerKeyFrame)) strcat("FID: ", ...
                                  num2str(finalFrameBMU))];
end

tree_KF = [];

%building branch by adding Key frame/s of closest neuron
for z1=1:neuron_KeyFrame(closestNeuron_1).frequency
    KF_temp = neuron_KeyFrame(closestNeuron_1).key_frames(z1);
    index_KF = find([keyFrame_Frame.key_frame]==KF_temp);

    tree_P = [tree_P strcat("NID: ",num2str(nr1))];
    tree_C = [tree_C strcat("KFID: ",num2str(KF_temp))];
    tree_KF = [tree_KF strcat("KFID: ",num2str(KF_temp))];
    
    % adding frames to respective key frames. be aware if there are multiple    
    for z2=1:keyFrame_Frame(index_KF).frequency
        if keyFrame_Frame(index_KF).frames(z2) ~= KF_temp
            tree_P = [tree_P strcat("KFID: ",num2str(KF_temp))];
            tree_C = [tree_C strcat("FID: ",num2str(keyFrame_Frame(index_KF).frames(z2)))];
        end
    end

end

if nr2>0
    %building branch by adding key frame/s of 2nd closest neuron
    for z3=1:neuron_KeyFrame(closestNeuron_2).frequency
        KF_temp = neuron_KeyFrame(closestNeuron_2).key_frames(z3);
        index_KF = find([keyFrame_Frame.key_frame]==KF_temp);

        tree_P = [tree_P strcat("NID: ",num2str(nr2))];
        tree_C = [tree_C strcat("KFID: ",num2str(KF_temp))];
        tree_KF = [tree_KF strcat("KFID: ",num2str(KF_temp))];

        % adding frames to respective key frames. be aware if there are multiple   
        for z2=1:keyFrame_Frame(index_KF).frequency
            if keyFrame_Frame(index_KF).frames(z2) ~= KF_temp
                tree_P = [tree_P strcat("KFID: ",num2str(KF_temp))];
                tree_C = [tree_C strcat("FID: ",num2str(keyFrame_Frame(index_KF).frames(z2)))];
            end
        end

    end
end

subplot(2,2,2);
G = graph(tree_P,tree_C);
h = plot(G,'o','EdgeColor','k', 'MarkerSize',4);
highlight(h,tree_Path_global,'EdgeColor','r','LineWidth',2);    % Highlighting Global path in RED
highlight(h,tree_Path_local,'EdgeColor','g','LineWidth',2);     % Highlighting Local feature path in Green
highlight(h,strcat("QFID: ",num2str(qImage)),'Marker','^','MarkerSize',12,'NodeColor','b'); % Colouring the Query Image Blue
if nr2>0
    highlight(h,[strcat("NID: ",num2str(nr1)) strcat("NID: ",num2str(nr2))],...
                                                             'Marker','s','MarkerSize',12,'NodeColor','r'); % Colouring Neurons as Red
else
    highlight(h,strcat("NID: ",num2str(nr1)),'Marker','s','MarkerSize',12,'NodeColor','r'); % Colouring Neurons as Red
end
if winnerKeyFrame == finalFrameBMU 
    highlight(h,tree_KF,'Marker','s','MarkerSize',12,'NodeColor','g'); %Colouring KEY FRAMES Black
else
    highlight(h,tree_KF,'Marker','o','MarkerSize',12,'NodeColor','k'); %Colouring KEY FRAMES Black
    highlight(h,strcat("FID: ",num2str(finalFrameBMU)),'Marker','s','MarkerSize',12,'NodeColor','g'); %Colouring/Highlighting winner frame with a green '*'
end
hold on;

% Dummy plots for Legend
LH(1) = plot(nan, nan, '-r');
L{1} = 'Closest Neuron (Gist)';
LH(2) = plot(nan, nan, '-g');
L{2} = 'Geometric Validation (Local Feature)';
LH(3) = plot(nan,nan,'^','MarkerSize',10,'MarkerFaceColor','b');
L{3} = 'Query Image';
LH(4) = plot(nan,nan,'sr','MarkerSize',10,'MarkerFaceColor','r');
L{4} = 'Neurons';
if winnerKeyFrame == finalFrameBMU 
    LH(5) = plot(nan,nan,'s','MarkerSize',10','MarkerFaceColor','g');
    L{5} = 'Key Frame/Best Matching Frame';
    LH(6) = plot(nan,nan,'ob','MarkerSize',4','MarkerFaceColor','b');
    L{6} = 'Frames';
    hold on;
    LH(7) = plot(nan,nan,'linestyle','none');
    str = strcat('Ground Truth: FID: ',num2str(gTruth));        % Obtained from Ground Truth
    L{7} = str;
    legend(LH,L);
else
    LH(5) = plot(nan,nan,'ok','MarkerSize',10,'MarkerFaceColor','k');
    L{5} = 'Key Frames';
    LH(6) = plot(nan,nan,'ob','MarkerSize',4','MarkerFaceColor','b');
    L{6} = 'Frames';
    LH(7) = plot(nan,nan,'s','MarkerSize',10','MarkerFaceColor','g');
    L{7} = 'Best Matching Frame';
    LH(8) = plot(nan,nan,'linestyle','none');
    str = strcat('Ground Truth: FID: ',num2str(gTruth));        % Obtained from Ground Truth
    L{8} = str;
    legend(LH,L);
end