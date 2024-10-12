function plot_optim_params
%% plots parameters for the optim_powcomb model, per subject

%% settings
data_folder = 'fit_data';
data_file_powcomb = @(subj_id) [data_folder filesep 'mf_optim_powcomb_' subj_id '.mat'];
burnin_powcomb = 22000;
subj_ids = {'231139','231140','548203','548214','950091','950107','967716','jason','eliana','kalpana'};
subj_num = length(subj_ids);
subj_colors = [255 0 0; 0  255 0; 0 0 255; 255 0 255; 0 255  255; ...
    255 255 0; 0 0 0; 112 219 147; 181 166 66; 95 159 159; ...
    184 115 51; 47 79 47; 153 50 205; 135 31 120; 133 94 66; ...
    84 84 84; 142 35 35; 245 204 176; 35 142 35; 205 127 50; ...
    219 219 112; 192 192 192; 82 127 118; 159 159 95; 142 35 107; ...
    47 47 79; 235 199 158; 207 181 59; 255 127 0; 219 112 219; ...
    217 217 243; 89 89 171; 140 23 23; 35 142 104; 107 66 38]*(200/255/255);



%% load different parameters & statistics (mode, var)
k_vest = NaN(subj_num, 2);
bound_vest = NaN(subj_num, 2);
tnds = NaN(subj_num, 3, 2);
lp = NaN(subj_num, 2);
biases = cell(subj_num, 1);
for subj_idx = 1:subj_num
    d = load(data_file_powcomb(subj_ids{subj_idx}));
    s = d.s((burnin_powcomb+1):end,:);
    k_vest(subj_idx,:) = [d.best_p(5) sqrt(var(s(:,5)))];
    bound_vest(subj_idx,:) = [d.best_p(6) sqrt(var(s(:,6)))];
    tnds(subj_idx,:,:) = [d.best_p(10) sqrt(var(s(:,10)));
        d.best_p(11) sqrt(var(s(:,11)));
        d.best_p(12) sqrt(var(s(:,12)))];
    lp(subj_idx,:) = [d.best_p(end) sqrt(var(s(:,end)))];
    c_num = length(d.cohs);
    biases{subj_idx} = [d.best_p(end-(2*c_num+1):(end-1)); ...
        sqrt(var(s(:,end-(2*c_num+1):(end-1)),[],1))];
end


%% plot parameters
pbar = 1;
% kvest
figure('Color','white');  hold on;
ylim([0.5 (0.5+subj_num)]);  xlim([0 100]);
for subj_idx = 1:subj_num
    rectangle('Position',[0 (subj_num-subj_idx+0.6) k_vest(subj_idx,1) 0.8],...
        'EdgeColor','none','FaceColor',subj_colors(subj_idx,:));
    plot(k_vest(subj_idx,1)+k_vest(subj_idx,2)*[-1 1],...
        (subj_num-subj_idx+1)*[1 1],'k-');
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);
xlabel('k_{vest}');  ylabel('subject');
% boundvest
figure('Color','white');  hold on;
ylim([0.5 (0.5+subj_num)]);  xlim([0 0.45]);
for subj_idx = 1:subj_num
    rectangle('Position',[0 (subj_num-subj_idx+0.6) bound_vest(subj_idx,1) 0.8],...
        'EdgeColor','none','FaceColor',subj_colors(subj_idx,:));
    plot(bound_vest(subj_idx,1)+bound_vest(subj_idx,2)*[-1 1],...
        (subj_num-subj_idx+1)*[1 1],'k-');
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);
xlabel('\theta_{vest}');  ylabel('subject');
% tnd
figure('Color','white');  hold on;
ylim([0.5 (0.5+subj_num)]);  xlim([0 1]);
for subj_idx = 1:subj_num
    for i = 1:3
        rectangle('Position',[0 (subj_num-subj_idx+0.6+(0.7/3+0.05)*(i-1))...
            tnds(subj_idx,i,1) (0.7/3)],'EdgeColor','none','FaceColor',...
            subj_colors(subj_idx,:)*(1/i)+[1 1 1]*0.5*(1-1/i));
        plot(tnds(subj_idx,i,1)+tnds(subj_idx,i,2)*[-1 1],...
            (subj_num-subj_idx+0.6+(0.7/3/2)+(0.7/3+0.05)*(i-1))*[1 1],'k-');
    end
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);
xlabel('tnd');  ylabel('subject');
% biases
figure('Color','white');  hold on;
ylim([0.5 (0.5+subj_num)]);
for subj_idx = 1:subj_num
    subj_biases = biases{subj_idx};
    c_num = size(subj_biases,2);
    for i = 1:c_num
        y = 0;  w = subj_biases(1,i);
        if w < 0;  y = w;  w = -w;  end;
        rectangle('Position',[y (subj_num-subj_idx+0.65+0.7*1.3*(i-1)/(1.3*c_num-0.3))...
            w (0.7/(1.3*c_num-0.3))],'EdgeColor','none','FaceColor',...
            subj_colors(subj_idx,:)*(1/i)+[1 1 1]*0.5*(1-1/i));
        plot(subj_biases(1,i)+subj_biases(2,i)*[-1 1],...
            (subj_num-subj_idx+0.65+0.7*1.3*(i-1)/(1.3*c_num-0.3)+0.35/(1.3*c_num-0.3))*[1 1],...
            'k-');
    end
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);
xlabel('biases');  ylabel('subject');
% lapses
figure('Color','white');  hold on;
ylim([0.5 (0.5+subj_num)]);  xlim([0 0.07]);
for subj_idx = 1:subj_num
    rectangle('Position',[0 (subj_num-subj_idx+0.6) lp(subj_idx,1) 0.8],...
        'EdgeColor','none','FaceColor',subj_colors(subj_idx,:));
    plot(lp(subj_idx,1)+lp(subj_idx,2)*[-1 1],...
        (subj_num-subj_idx+1)*[1 1],'k-');
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);
xlabel('p_l');  ylabel('subject');



%figure; bar(tnds(:,:,1));