% constructing the 5-D matrix to do dPCA
% the way below is very rude, need to be refined


PSTH_ves_expand = cell(2,1);
PSTH_ves_expand_tomat = cell(2,1);
for j = 1:length(PSTH_ves)
    for reps = 1:20
        for c = 3: size(PSTH_ves{j}.raw,1)
            tobeappended = nan((20-size(PSTH_ves{j}.raw{c},1)), size(PSTH_ves{j}.raw{c},2));
            PSTH_ves_expand{j}.raw{c,1} = ([PSTH_ves{j}.raw{c};tobeappended])';
            PSTH_ves_expand_tomat{j} =  [PSTH_ves_expand_tomat{j} PSTH_ves_expand{j}.raw{c,1}(:,reps)];
        end
    end
%     PSTH_ves_expand_tomat{j} = (PSTH_ves_expand_tomat{j})';
end

PSTH_ves_expand_tomat_reshape{1} = reshape(PSTH_ves_expand_tomat{1}, 200,8,20);

PSTH_ves_expand_tomat_reshape_permute{1} = permute(PSTH_ves_expand_tomat_reshape{1},[2,1,3]);

PSTH_ves_expand_tomat_reshape_permute_3cond = repmat(PSTH_ves_expand_tomat_reshape_permute{1},3,1);

PSTH_3Cond = reshape(PSTH_ves_expand_tomat_reshape_permute_3cond, 8,3,200,20);

% for i=1:20
%     PSTH_ves_expand_tomat_reshape{1}(:,:,i) = PSTH_ves_expand_tomat_reshape{1}(:,:,i)';
% end
