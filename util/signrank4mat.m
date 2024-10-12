function [p,h,stats] = signrank4mat(x,varargin)

% Signrank for matrix
% ZZ @20221228
% Input:
% x: vector or matrix
% varargin:
% y: if existed, should have the same rank with x
% 1 or 2: row-wise or column-wise, default is column's comparison

switch length(varargin)
    case 0     % Only input x
        disp('Default column-wise')
        column = 1;
        
    case 1
        if varargin{1} == 1   % Input x and 1
            disp('Column-wise')
            column = varargin{1};
            
        elseif varargin{1} == 2  % Input x and 2
            disp('Row-wise')
            column = varargin{1};
            
        elseif length(varargin{1}) > 1      % the other input is y for paired test
            disp('Default column-wise! Paired test...')
            column = 1;
            
            if size(x,2) ~= size(varargin{1},2)
                error('X and Y have different numbers of columns')
            else
                yy = varargin{1};
            end
            
        end
        
    case 2
        column = varargin{2};
        if length(varargin{2}) >1
            error('The Third input sholud be 1 or 2')
        else
            if column == 1
                disp('Column-wise')
                if size(x,2) ~= size(varargin{1}, 2)
                    error('X and Y have different numbers of columns')
                end
                
            elseif column == 2
                disp('Row-wise')
                if size(x,1) ~= size(varargin{1}, 1)
                    error('X and Y have different numbers of rows')
                end
            else
                error('The Third input should be 1 or 2')
            end
                           
            yy = varargin{1};
        end
        
    otherwise
        error('Number of inputs is too many!');
end


% -- Signrank test
if exist('yy')
    if column == 1
        for c = 1:size(x,2)
            [p(:,c),h(:,c),stats(:,c)] = signrank(x(:,c), yy(:,c));
        end
        
    elseif column == 2
        for r = 1:size(x,1)
            [p(r,:),h(r,:),stats(r,:)] = signrank(x(r,:), yy(r,:));
        end
    end
    
else
    if column == 1
        for c = 1:size(x,2)
            [p(:,c),h(:,c),stats(:,c)] = signrank(x(:,c));
        end
        
    elseif column == 2
        for r = 1:size(x,1)
            [p(r,:),h(r,:),stats(r,:)] = signrank(x(r,:));
        end
    end
end

end