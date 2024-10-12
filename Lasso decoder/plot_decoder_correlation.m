% Visualize the correlation of each temporal decoder
% Averaged Correlation image or Angle image
% SVM or Lasso now
% ZZ @20220424

function plot_decoder_correlation(model, decoder_type, st, decoder_tcenters)

% model: Raw data from script 'decoder_train', cells, each cell has many
% repetitions of the same decoder

% decoder_type, 1, SVM    2, Lasso

%st, stimulus type, 1, vestibular   2, visual,   3, combined

% decoder_tcenters, center time of each time window

%% Pre-definition
% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')

stimtype_name = {'Vestibular', 'Visual', 'Combined'};
stimtype_color = [41 89 204; 248 28 83; 14 153 46]/255;

image_colormap = [0:0.05:1; 0:0.05:1; 0:0.05:1]' ;

decoder_type_name = {'SVM', 'LASSO'};
%% Extract decoding plane
if decoder_type == 1
    for j = 1:size(model,1)       % Which period
        for nn = 1: size(model, 2) % Which repetition
            for tt = 1:sum(cellfun(@(x) ~isempty(x), model(j,nn,:)))  % Which time window
                
                % By default, 'autoscale' = true, so here i should rescale back the weight
                weight{j,nn}(:,tt) = model{j,nn,tt}.SupportVectors' * model{j,nn,tt}.Alpha;
            end
        end
    end
    
elseif decoder_type == 2
    for j = 1:size(model.B,1)       % Which period
        for nn = 1: size(model.B, 2) % Which repetition
%             for tt = 1:sum(cellfun(@(x) ~isempty(x), model.B(j,nn,:)))  % Which time window
                
                weight{j,nn} = [model.B{j,nn,:}];
%                 min_mse_ind = model.FitInfo{j,nn,tt}.IndexMinMSE;
%                 weight{j,nn}(:,tt) = model.B{j,nn,tt}(:,min_mse_ind);
%             end
        end
    end
end


%% ================================== Plotting =================================
% 1.) Correlation
figure('pos', [50 50 1200 700], 'name', [decoder_type_name{decoder_type} ', st=' num2str(st), ', Weight Correlation']);
ha = tight_subplot(1, size(weight,1), [],[],[],[], [size(weight{1,1},2), size(weight{2,1},2)]);
for j = 1:size(weight,1)
    set(gcf, 'CurrentAxes',  ha(j)); hold on;
    
    weight_all_boot = reshape(cell2mat(weight(j,:)),size(weight{j,1},1),size(weight{j,1},2),size(weight,2));
    [cor, pp] = corr(mean(weight_all_boot,3), 'type','Spearman');
    
    imagesc(decoder_tcenters{j},decoder_tcenters{j}, cor,[0,1]);
    colormap(image_colormap);
    if j == 2
        colorbar('southoutside')
    end
    
    axis tight;
    axis square;
    
end
SetFigure();

% 2.) Angle
figure('pos', [70 50 1200 700], 'name', [decoder_type_name{decoder_type} ', st=' num2str(st), ', Decoder Angle']);
ha = tight_subplot(1, size(weight,1), [],[],[],[], [size(weight{1,1},2), size(weight{2,1},2)]);

for j = 1:size(weight,1)
    for nn = 1:size(weight,2)
        angle_beta = [];
        
        for tt = 1:size(weight{j,nn},2)
            for ttt = tt : size(weight{j,nn},2)
                
                % Angle between different time-window decoders
                angle_beta(tt,ttt) = abs(acos(weight{j,nn}(:,tt)' * weight{j,nn}(:,ttt) ./ sqrt(sum(weight{j,nn}(:,tt) .^2)) ./ sqrt(sum(weight{j,nn}(:,ttt) .^2))) / pi * 180);
            end
        end
        angle_all_boot{j,nn} = angle_beta + triu(angle_beta,1)';   % Copy values in upper triangle to the lower
    end
    angle{j} = reshape(cell2mat(angle_all_boot(j,:)), size(angle_beta,1), size(angle_beta,1),nn);
    
    angle_mean{j} = nanmean(angle{j}, 3);
    angle_sem{j} = nanstd(angle{j},0,3) ./ sqrt(sum(~isnan(angle{j}),3));
    
%     % Permutation test
%     for tt = 1:size(weight{j,nn},2)
%         for ttt = tt : size(weight{j,nn},2)
%             pp_angle{j}(tt,ttt) = permutationTest(squeeze(angle{j}(tt,ttt,:)), 90, 1000, 'exact', 1);
%         end
%     end

    set(gcf, 'CurrentAxes', ha(j)); hold on;
    
    imagesc(decoder_tcenters{j},decoder_tcenters{j},angle_mean{j}, [0 100]);
    colormap (image_colormap)
    
%     % Significance indicating 
%     contour(decoder_tcenters{j}, decoder_tcenters{j}, pp_angle{j}<0.05, [1 1], 'm-', 'linewidth',2)
    
    if j == 2
        colorbar('southoutside')
    end
    
    axis tight; axis square;
end
SetFigure();


% How many non-zero beta for Lasso decoder
if decoder_type ==2    % only for Lasso, for regularization
    figure('pos', [70 50 800 700], 'name', [decoder_type_name{decoder_type} ', st=' num2str(st), ', Non-zero weight number']);
    ha = tight_subplot(1, size(weight,1), [],[],[],[], [size(weight{1,1},2), size(weight{2,1},2)]);
    
    non_zero_betaN_temp = cellfun(@(x) sum(x~=0), weight, 'uniform',0);
    
    for j = 1:size(non_zero_betaN_temp, 1)
        non_zero_betaN{j} = reshape(cell2mat(non_zero_betaN_temp(j,:)),[],size(non_zero_betaN_temp,2));
        
        non_0_betaN_mean{j} = mean(non_zero_betaN{j},2);
        non_0_betaN_std{j} = std(non_zero_betaN{j},0,2);
        
        set(gcf, 'CurrentAxes', ha(j));
        shadedErrorBar(decoder_tcenters{j}, non_0_betaN_mean{j}, non_0_betaN_std{j})
        
        if j ==2
            set(gca,'ytick',[],'ycolor','w')
        end
        
        axis tight;
    end
    max_mean = max(cell2mat(non_0_betaN_mean'));
    set(findall(gcf, 'type', 'axes'), 'YLim', [-10 max_mean+20]);
    
    SetFigure;
end

end