function [H ,D]= H_from_feat(oim,rim,res,Nl)
    co = chrodist(reshape(oim,[],3)',res);
    cr = chrodist(reshape(rim,[],3)',res);
    Nl_max = Nl+1;
    npts = 0;
    while (npts<4 && Nl<Nl_max)
        Nl = Nl+1;
        Io = co.^(1/Nl);
        Ir = cr.^(1/Nl);
        [p1,p2] = distmatchbyfeat(Io,Ir);
        npts = size(p1,2);
    end
    if npts==0
        H = [];
    else
        [inliers, H,D] = ransacfithomography(p1,p2,0.01);
    end
    C = [1,0,0;0,1,0;1,1,1];
    disp('Estimated H');
    disp(H);

    fprintf('Matches: %d, Nl: %d\n',npts,Nl);
    if ~isempty(H)
        poim = reshape(oim,[],3)';
        poC = C*poim;
        size(poC)
        peC = H*poC;
        pe = C\peC;
        Ie = chrodist(pe,res).^(1/Nl);
        eim = reshape(pe',size(oim));
        sc = mean(oim(:))/mean(eim(:));
        eim = eim*sc;    
    end

    figure;
    imshow(1-[Io,Ir]);
    if ~isempty(H)
        x = round([p1(1,:); p2(1,:) + size(Io,2)/res]*res);
        y = round([p1(2,:); p2(2,:)]*res);
        line(x(:,inliers), y(:,inliers),'LineWidth',2);

        figure;
        subplot(1,3,1);
        imshow(oim);
        title('original');
        subplot(1,3,2);
        imshow(rim);
        title('reference');
        subplot(1,3,3);
        imshow(max(eim,0));
        title('estimated');

        figure;
        subplot(1,2,1);
        imshowpair(Io,Ir);
        title('Inital Gamuts');
        subplot(1,2,2);
        imshowpair(Ie,Ir);
        title('Aligned Gamuts');
    end

end

