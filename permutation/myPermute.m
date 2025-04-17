function [posinvpmap,neginvpmap,voxstatmap] = myPermute(brain,mask,initp,covar,stattype,correctiontype,corrtype)
    if nargin < 7; corrtype = 'Spearman'; end % by default spearman correlation
    assert(~any(isnan(brain(:))),'there is nan in brain data')
    assert(ismember(corrtype,{'Spearman','Pearson'}),'Correlation type should either be ''Pearson'' or ''Spearman''')
    assert(ismember(correctiontype,{'FDR','FWER','uncorrected'}),'Correction type should be one of ''FDR'', ''FWER'', or ''uncorrected''')
    assert(ismember(stattype,{'mass','size'}),'Stat type should be one of ''mass'' or ''size''')

    nrep = 5000; brain = double(brain); mask = double(mask);

    if ~isempty(covar); covar = double(covar);
        if strcmp(corrtype,'Spearman'); brain = tiedrank(brain); covar = tiedrank(covar); end % prepare for spearman by pre-ranking the data
    end

    % get original test p-value and clusters formed
    [pos_cstat,neg_cstat,posCC,negCC,voxstatmap] = calcinvp(brain,covar,mask,0,1-initp,stattype);
    posinvpmap = zeros(size(mask)); neginvpmap = zeros(size(mask));

    if posCC.NumObjects > 0 || negCC.NumObjects > 0 % there's at least one cluster at initial threshold
        % do permutation to calculate null
        posnull = cell(nrep,1); negnull = cell(nrep,1);
        parfor i = 1:size(posnull,1)
            i
            [posnull{i},negnull{i}] = calcinvp(brain,covar,mask,1,1-initp,stattype);
        end

        if strcmp(correctiontype,'FWER') % max cluster stat within perm
            posnull = cellfun(@max,posnull); negnull = cellfun(@max,negnull);
        else % random selection of cluster within perm
            posnull = cellfun(@(x) x(randi(length(x))),posnull); negnull = cellfun(@(x) x(randi(length(x))),negnull);
        end
        pospval = sum(posnull>pos_cstat')./nrep; negpval = sum(negnull>neg_cstat')./nrep;
        
        if strcmp(correctiontype,'FDR') % FDR adjustment
            [~, ~, ~, pospval]=fdr_bh(pospval); [~, ~, ~, negpval]=fdr_bh(negpval);
        end

        % assigning inverse p-values to clusters
        for i = 1:posCC.NumObjects; posinvpmap(posCC.PixelIdxList{i}) = 1-pospval(i); end
        for i = 1:negCC.NumObjects; neginvpmap(negCC.PixelIdxList{i}) = 1-negpval(i); end
    end
end

function [pos_cstat,neg_cstat,pos_CC,neg_CC,statmap] = calcinvp(brain,covar,mask,shuffle,init_thresh,stattype) % calculate cluster statistics
    if isempty(covar) % t-test
        if shuffle == 1; brain = brain .* (randsample([-1, 1],size(brain,1),true)'); end % sign flipping shuffle
        stat = mean(brain) ./ (std(brain) ./ sqrt(size(brain,1))); p = tcdf(-stat, size(brain,1) - 1);
    else % correlation test
        if shuffle == 1; covar = covar( randsample(length(covar),length(covar)) ); end % random shuffling of pairs
        [stat,p] = corr(brain,covar,'tail','right');
    end
    statmap = mask; statmap(mask(:)==1) = stat; voxinvpmap = mask;

    pos_cstat = 0; voxinvpmap(mask(:)==1) = 1-p; pos_CC = bwconncomp(voxinvpmap>init_thresh);
    if ~isempty(pos_CC.PixelIdxList)
        if strcmp(stattype,'mass') % cluster mass
            pos_cstat = cellfun(@(x) sum(statmap(x)), pos_CC.PixelIdxList)';
        else % cluster size
            pos_cstat = cellfun(@length,pos_CC.PixelIdxList)';
        end
    end
    neg_cstat = 0; voxinvpmap(mask(:)==1) = p; neg_CC = bwconncomp(voxinvpmap>init_thresh);
    if ~isempty(neg_CC.PixelIdxList)
        if strcmp(stattype,'mass') % cluster mass
            neg_cstat = cellfun(@(x) -sum(statmap(x)), neg_CC.PixelIdxList)';
        else % cluster size
            neg_cstat = cellfun(@length,neg_CC.PixelIdxList)';
        end
    end
end