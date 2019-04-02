function chist = chrodist(rgb,Nbin)
    mask = rgb<0.01;
    mask = ~logical(sum(mask,1)) & mean(rgb,1)>0.05;

    C = [1,0,0;0,1,0;1,1,1]; % base converse matrix
    pC = C*rgb(:,mask); % chromaticity array
    hpC = bsxfun(@rdivide,pC,pC(3,:)); % homogenous
    
    hr = linspace(0,1,Nbin+1);
    sz = zeros(1,2);
    X = hpC(1:2,:)';
    
    loc = zeros(size(X));
    for i = 1:2
        [~,loc(:,i)] = histc(X(:,i), hr , 1);
        sz(i) = length(hr)-1;
    end
    chist = (accumarray(loc(:,:), 1, sz))';
    chist = chist/max(chist(:));
end