function h = LinearCorrelation(rawx,rawy,varargin) % (figN,logx,logy,seperate)
% Linear correlation, scatter points, etc.
% 
% result = LinearCorrelation(x,y,varargin)
%   If x and y are vectors, direct linear correlation
%   If x and y are cell arrays, linear correlation for each element pairs +
%         combinations in 'CombinedIndex'
%
%   Optional properties
%       'CombinedIndex': binary coded combinations of different elements
%                       (e.g., [3 7] = 1 and 2 and 3 and 1&2 and 1&2&3)
%       'figN': 2499(default)
%       'logx' & 'logy': 0 (default) or 1
%       'Method': 'Pearson' (default) or 'Spearman'
%       'LineStyles': {'k-','k:','r-','b-','g-','c-','m-'}
%       'Markers': {'o','o','s','^','v','<','>'};
%       'MarkerSize': 15
%       'FaceColors': {'k','none',[0.8 0.8 0.8]};
%       'Axes': axes
%
% @HH20141130

% ------ Parse input parameters -------
paras = inputParser;

addOptional(paras,'figN',2499);
addOptional(paras,'logx',0);
addOptional(paras,'logy',0);

addOptional(paras,'MethodOfCorr','Pearson');
addOptional(paras,'FittingMethod',1); % Type 1: Normal squared error (LSM); Type 2: Perpendicular error (solved by PCA)

addOptional(paras,'CombinedIndex',[]);
addOptional(paras,'PlotCombinedOnly',-1);   %  If true, then only plot those in CombinedIndex, no single ones.
addOptional(paras,'Style','ko--');
addOptional(paras,'LineStyles',{'k-','k:','r-','b-','g-','c-','m-'});
addOptional(paras,'Markers',{'o','o','s','^','v','<','>'});
addOptional(paras,'MarkerSize',10);
% addOptional(paras,'FaceColors',{'k','none',[0.8 0.8 0.8]});
addOptional(paras,'FaceColors',{'none'});  % Changed by ZZ 20220117
addOptional(paras,'EdgeColors',{'none'});  % Added by ZZ 20230404
addOptional(paras,'Xlabel','X');
addOptional(paras,'Ylabel','Y');

% Add same scale. HH20141217
addOptional(paras,'SameScale',0);
addOptional(paras, 'AxisSquare', 0);   % Added by ZZ 20220117
addOptional(paras, 'Diagonal', 0);   % Plot diagonal or not

% Add histogram.  HH20141217
addOptional(paras,'XHist',0);
addOptional(paras,'XHistStyle','grouped');
addOptional(paras,'YHist',0);
addOptional(paras,'YHistStyle','grouped');
addOptional(paras,'Axes',[]);

parse(paras,varargin{:});

for ii = 1:length(paras.Results.EdgeColors)
    EdgeColors{ii} = paras.Results.EdgeColors{ii};
%     EdgeColors{ii} = paras.Results.LineStyles{ii}(1);
end
% --------- End input parser ----------

% -------- Reorganize data (all to cell format) --------
duplicateX = 0;
if ~iscell(rawx) % raw x is matrix
    if size(rawx,2) ~=  size(rawy,2)
        if size(rawx,2)==1
            rawx = repmat(rawx,1,size(rawy,2));
            duplicateX = 1;
        else
            disp('Size error...');
            return;
        end
    end
    % Transform to cell
    x = mat2cell(rawx,size(rawx,1),ones(1,size(rawx,2)));
    y = mat2cell(rawy,size(rawy,1),ones(1,size(rawy,2)));
else % x is cell
    if length(rawx)~=length(rawy)
        disp('Size error...');
        return;
    end
    x = rawx;
    y = rawy;
end

for j = 1:length(x)
    if length(x{j}) ~= length(y{j})
        disp('Size error...');
        return;
    end;        
end

% -------- End reorganizing data --------

if isempty(paras.Results.Axes)  % Generate a new figure
    
    set(figure(paras.Results.figN),'position',[  36 81 1098 729]); % clf;
    ymult = 1098/ 729;
    
    if paras.Results.XHist
        h.ax_xhist = axes('position',[0.1 0.75 0.4 0.206]);
        set(h.ax_xhist,'xtick',[]);
        ylabel('Cases');  hold on;
    end
    if paras.Results.YHist
        h.ax_yhist = axes('position',[0.52 0.11 0.206/ymult 0.4*ymult]);
        set(h.ax_yhist,'xtick',[]);
        ylabel('Cases');  hold on;
        view([90 -90])
    end
    
    h.ax_raw = axes('Position',[0.1 0.11 0.4 0.4*ymult]);
    
else % Use the requested axis
    h.ax_raw = paras.Results.Axes;
    set(h.ax_raw,'color','none');
    pos = get(h.ax_raw,'position');
    
    if paras.Results.XHist
        h.ax_xhist = axes('position',[pos(1) pos(2)+pos(4) pos(3) pos(4)/5]);
        set(h.ax_xhist,'xtick',[]); hold on; axis off;
%         ylabel('Cases'); 
    end
    if paras.Results.YHist
        h.ax_yhist = axes('position',[pos(1)+pos(3) pos(2) pos(3)/5 pos(4)]);
        set(h.ax_yhist,'xtick',[]); hold on; axis off;
%         ylabel('Cases');  
        view([90 -90])
    end
    
end

PlotCombinedOnly = paras.Results.PlotCombinedOnly;
if PlotCombinedOnly == -1  % Not specified, automatically define. Backward compatibility
    if isempty(paras.Results.CombinedIndex) % Clearly, the user wanted to plot fittings for each raw data
        PlotCombinedOnly = 0;
    else
        PlotCombinedOnly = 1;
    end        
end

corr_inds = [2.^((1:length(x))-1) paras.Results.CombinedIndex];

xx_all = [];
yy_all = [];
non_empty_dot = [];
non_empty_line = [];

for i = 1:length(corr_inds)
    
    % -- Combine data   
    xx = [];
    yy = [];
    combs = find(fliplr(dec2bin(corr_inds(i)))=='1');
    for cc = combs % Parse binary code
        xx = [xx ; x{cc}(:)];
        yy = [yy ; y{cc}(:)];
    end
    
    nanInd = or(isnan(xx),isnan(yy));
    xx(nanInd) = [];
    yy(nanInd) = [];
    
    if paras.Results.logx
        xx = log(xx);
        x_text = ['Log(' paras.Results.Xlabel ')'];
    else
        x_text = paras.Results.Xlabel;
    end
    
    if paras.Results.logy
        yy = log(yy);
        y_text = ['Log(' paras.Results.Ylabel ')'];
    else
        y_text = paras.Results.Ylabel;
    end
    
    % -- Linear correlation   
    r = nan; p = nan; k = nan;
    
    if ~isempty(xx)
       
        if ~ PlotCombinedOnly || i > length(x) 
        
            [r,p] = corr(xx,yy,'type',paras.Results.MethodOfCorr);
            
            if paras.Results.FittingMethod == 2 % Minimize perpendicular distance (PCA)
                
                %             coeff = pca([xx yy]);
                %             if ~isempty(coeff)
                %                 linPara(1) = coeff(2) / coeff(1);
                %                 linPara(2) = mean(yy)- linPara(1) * mean(xx);
                %             else
                %                 linPara = [nan nan];
                %             end
                
                % New version of regressing perpendicular error HH20180613
                fitType1 = regress_perp(xx,yy,1,0.05);
                
                if ~isempty(fitType1)
                    linPara(1) = fitType1.k;
                    linPara(2) = fitType1.b;
                else
                    linPara = [nan nan];
                end
                linParaSE = [std(fitType1.bootK) std(fitType1.bootB)];
                
                h.group(i).fitType1 = fitType1;
                
                % -- Plotting
                xxx = linspace(min(xx),max(xx),1000);
                Y = linPara(1) * xxx + linPara(2);
                
                xxx = xxx(min(yy) <= Y & Y <= max(yy));
                Y = Y(min(yy) <= Y & Y <= max(yy));
                
            elseif paras.Results.FittingMethod == 1 % LMS method
                
                % [linPara,S] = polyfit(xx,yy,1);
                
                [b,~,stat]= glmfit(xx,yy); % glmfit can directly give us std of coeffs (se)
                linPara = fliplr(b'); % Note the definitions different
                linParaSE = fliplr(stat.se');
                
                % -- Plotting
                xxx = linspace(min(xx),max(xx),1000);
                Y = polyval(linPara,xxx);
                
            end
            
            % -- Saving
            h.group(i).line = plot(h.ax_raw, xxx,Y,paras.Results.LineStyles{1+mod(i-1,length(paras.Results.LineStyles))},'linewidth',2);
            non_empty_line = [non_empty_line i];
            
            %     h.group(i).text = text(xxx(end),Y(end),sprintf('\\itr\\rm = %3.3f, \\itp\\rm = %3.3f',r,p),'color',paras.Results.LineStyles{i}(1));
            line_leg{i} = ['[' num2str(combs) '] ' sprintf('\n\\itr^2\\rm = %3.3g,\\itp\\rm = %2.2g, \\itk\\rm = %3.3g\\pm%.1g',r,p,linPara(1),linParaSE(1))];

%             h.group(i).r_square = r^2;
            h.group(i).r_square = r;
            h.group(i).p = p;
            h.group(i).para = linPara;
            h.group(i).paraSE = linParaSE;
            % h.group(i).S = S;

        end
        
        axes(h.ax_raw);
        
        if i <= length(x)   % Stuffs for raw data (no combinations)
            hold on;
            h.group(i).dots = plot(xx,yy,paras.Results.Markers{1+mod(i-1,length(paras.Results.Markers))},...
                'Color',paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))},...
                'MarkerFaceColor',paras.Results.FaceColors{1+mod(i-1,length(paras.Results.FaceColors))},...
                'MarkerSize',paras.Results.MarkerSize);
%             h.group(i).dots = plot(xx,yy,paras.Results.Markers{1+mod(i-1,length(paras.Results.Markers))},...
%                 'Color',paras.Results.LineStyles{1+mod(i-1,length(paras.Results.LineStyles))}(1),...
%                 'MarkerFaceColor',paras.Results.FaceColors{1+mod(i-1,length(paras.Results.FaceColors))},...
%                 'MarkerSize',paras.Results.MarkerSize);
            
            dot_leg{i} = num2str(i);
            non_empty_dot = [non_empty_dot i];        
            
            % Cache raw data for histogram later
            xx_for_hist{i} = xx;
            xx_all = [xx_all; xx];
            yy_for_hist{i} = yy;
            yy_all = [yy_all; yy];
        end
        
    else
        xx_for_hist{i} = xx;
        yy_for_hist{i} = yy;

    end

end

axes(h.ax_raw);

xlabel(x_text);
ylabel(y_text);
axis tight;

% SetFigure(20);

xlims = xlim;
ylims = ylim;

if paras.Results.SameScale  % Same scale comparison
    
    xlims_new = [min(xlims(1),ylims(1)) max(xlims(2),ylims(2))];
    xlims_new = [xlims_new(1)-range(xlims_new)/20 xlims_new(2)+range(xlims_new)/20];
    
    axis([xlims_new xlims_new]);
    h.diag = plot([xlims_new(1) xlims_new(2)],[xlims_new(1) xlims_new(2)],'k:');
    xticks = get(gca,'xtick');
    set(gca,'ytick',xticks);
    
else 
    xlim([xlims(1)-range(xlims)/20 xlims(2)+range(xlims)/20]);
    ylim([ylims(1)-range(ylims)/20 ylims(2)+range(ylims)/20]);
    
end


try
    leg = legend([h.group(non_empty_dot).dots h.group(non_empty_line).line],...
        {dot_leg{non_empty_dot} line_leg{non_empty_line}},'AutoUpdate','off'); % ,'Location', 'EastOutside');

catch
    keyboard
end

% set(leg,'Fontsize',15,'Position',[0.65 0.70 0.346 0.284]);

h.leg = leg;

% Added by ZZ 20220117
if paras.Results.AxisSquare
%     axis equal;
    axis square;
end

xlims = xlim;
ylims = ylim;

if paras.Results.Diagonal
    plot([0 0], ylims, 'k--');
    plot(xlims, [0 0], 'k--');
    plot([xlims(1) xlims(2)], [ylims(1) ylims(2)], 'k-', 'linewidth',1.5);
end

%% Histogram goes here
        
if paras.Results.XHist
    % Find unique xcenters
    x_centers = linspace(min(xx_all),max(xx_all),abs(paras.Results.XHist));
    
    for i = 1:length(x) % Raw data
        hist_x(:,i) = hist(xx_for_hist{i},x_centers);
        if strcmp(paras.Results.XHistStyle,'stacked')
            mean_x(i) = median(cell2mat(xx_for_hist(:)));
        else
            mean_x(i) = median(xx_for_hist{i});
        end
    end
    
    if range(x_centers) > 0
        
        axes(h.ax_xhist);
        if paras.Results.XHist > 0  % Hist
            
            if duplicateX  % Raw data only have one Xs
                h_b = bar(x_centers,hist_x(:,1),1,paras.Results.XHistStyle);
                set(h_b(1),'EdgeColor','k',...
                    'FaceColor','k');
                plot([mean_x(1) mean_x(1)],ylim,'k--','linew',2);
                
            else
                h_b = bar(x_centers,hist_x,1,paras.Results.XHistStyle);
                
                for i = 1:length(x) % Raw data
                    set(h_b(i),'EdgeColor',EdgeColors{1+mod(i-1,length(paras.Results.LineStyles))},...
                        'FaceColor',paras.Results.FaceColors{1+mod(i-1,length(paras.Results.FaceColors))});
                    plot([mean_x(i) mean_x(i)],ylim,[paras.Results.LineStyles{1+mod(i-1,length(paras.Results.LineStyles))}],'linew',2);
                end
            end
            
            
        else % Cumulative hist
            for i = 1:length(x)
                plot(x_centers,cumsum(hist_x),paras.Results.LineStyles{1+mod(i-1,length(paras.Results.LineStyles))},'linew',2);
            end
        end
        
        axis tight;
        % xlim(xlims);
        linkaxes([h.ax_xhist h.ax_raw],'x');
    end
end

if paras.Results.YHist
    % Find unique xcenters
    y_centers = linspace(min(yy_all),max(yy_all),abs(paras.Results.YHist));
    
    for i = 1:length(x) % Raw data
        hist_y(:,i) = hist(yy_for_hist{i},y_centers);
        if strcmp(paras.Results.YHistStyle,'stacked')
            mean_y(i) = median(cell2mat(yy_for_hist(:)));
        else
            mean_y(i) = median(yy_for_hist{i});
        end
    end
    
    axes(h.ax_yhist);
    
    if range(y_centers) > 0
        if paras.Results.YHist > 0  % Hist
            
            h_b = bar(y_centers,hist_y,1,paras.Results.YHistStyle);
            
            for i = 1:length(x) % Raw data
                set(h_b(i),'EdgeColor',EdgeColors{1+mod(i-1,length(paras.Results.LineStyles))},...
                    'FaceColor',paras.Results.FaceColors{1+mod(i-1,length(paras.Results.FaceColors))});
                plot([mean_y(i) mean_y(i)],ylim,[paras.Results.LineStyles{1+mod(i-1,length(paras.Results.LineStyles))}],'linew',2);
            end
            
        else % Cumulative hist
            for i = 1:length(x)
                plot(y_centers,cumsum(hist_y(:,i)),paras.Results.LineStyles{1+mod(i-1,length(paras.Results.LineStyles))},'linew',2);
            end
        end
        
        axis tight;
        xlim(ylims);
    end
end

axes(h.ax_raw);


