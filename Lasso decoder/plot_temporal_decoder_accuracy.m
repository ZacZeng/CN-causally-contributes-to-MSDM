% Visualize predicted accuracy of cross temporal decoder
% SVM and Lasso now
% ZZ @20220424

function plot_temporal_decoder_accuracy(accuracy, decoder_type, st, decoder_tcenters)

% accuracy, cell structure, number of the way of temporal alignment * bootstrapN

% decoder_type, 1,SVM    2, Lasso

% st, stimulus type, 1, Vestibular   2, Visual    3, Combined

% decoder_tcenters, center time of each decoder
%% Pre-definition
% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')

stimtype_name = {'Vestibular', 'Visual', 'Combined'};
stimtype_color = [41 89 204; 248 28 83; 14 153 46]/255;

% image_colormap = [0:0.05:1; 0:0.05:1; 0:0.05:1]' ;

decoder_type_name = {'SVM', 'LASSO'};

%% Average N-bootstrap accuracy
figure('pos', [70 50 1200 700], 'name', [decoder_type_name{decoder_type} ', test st=' num2str(st), ', Prediction accuracy']);
% ha = tight_subplot(1, size(accuracy,1), [],[],[],[], [size(accuracy{1,1},2), size(accuracy{2,1},2)]);
ha = tight_subplot(1, length(accuracy), [],[],[],[], [size(accuracy{1},1), size(accuracy{2},1)]);

for j = 1:length(accuracy)  %size(accuracy, 1)
%     accuracy_all_boots = reshape(cell2mat(accuracy(j,:)), size(accuracy{j,1},1), size(accuracy{j,1},1), size(accuracy,2));
    accuracy_all_boots = accuracy{j};
    accuracy_mean = mean(accuracy_all_boots,3);
    accuracy_sem = std(accuracy_all_boots,0,3) / sqrt(size(accuracy,2));
    
    pp = []; 
    % permutation test
    progressbar('Cross-modal accuracy')

        for tt= 1:size(accuracy{j},1)
            parfor ttt = 1:size(accuracy{j},1)
%                 pp(tt,ttt) = sum(squeeze(accuracy_all_boots(tt,ttt,:)) <= 0.5) / size(accuracy_all_boots(tt,ttt,:),3);
                pp(tt,ttt) = permutationTest(squeeze(accuracy_all_boots(tt,ttt,:)),0.5, 1000, 'exact',1);
            end
            progressbar(tt / size(accuracy{j},1));
        end
        
        set(gcf, 'CurrentAxes', ha(j));
        imagesc(decoder_tcenters{j}, decoder_tcenters{j}, accuracy_mean, [0.4 1]);
%         colormap jet;
        
        hold on;
        contour(decoder_tcenters{j}, decoder_tcenters{j}, pp<0.01, [1 1], 'm-',  'linewidth',2);
        
        plot([0 0], ylim, 'k--');  
        plot(xlim, [0 0], 'k--');
        
        axis tight;  axis square;
        
        if j == 2
            colorbar('southoutside')
        else
            xlim([0-50 1500+50]); ylim([0-50 1500+50]); 
            xticks(0:500:1500); yticks(0:500:1500);
            xlabel('Testing'); ylabel('Training')
        end
        
        set(gca, 'YDir','normal');

end
SetFigure();

end