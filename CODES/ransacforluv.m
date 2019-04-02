function [M, inliers] = ransacforluv(x, fittingfn, s, t,maxDataTrials, maxTrials,white)
    p = 0.99;
    [~, npts] = size(x);
    bestM = [];     
    trialcount = 0;
    bestscore =  0;
    N = 1;            
    while N > trialcount
        degenerate = 1;
        count = 1;
        while degenerate
            ind = randsample(npts, s);
            degenerate = isdegenerate( x(:,ind));
            if ~degenerate
                M = feval(fittingfn, x(:,ind));
                if isempty(M)
                    degenerate = 1;
                end
            end
            count = count + 1;
            if count > maxDataTrials
                break
            end
        end

        if ~degenerate
            [inliers, M] = distfn( M, x, t,white);
        else
            inliers = [];
        end
        ninliers = length(inliers);
        if ninliers > bestscore   
            bestscore = ninliers;  
            bestinliers = inliers;
            bestM = M;
            fracinliers =  ninliers/npts;
            pNoOutliers = 1 - fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  
            pNoOutliers = min(1-eps, pNoOutliers);
            N = log(1-p)/log(pNoOutliers);
        end
        
        trialcount = trialcount+1;
        if trialcount > maxTrials
            warning( sprintf('ransac reached the maximum number of %d trials',maxTrials));
            break
        end
    end
      
    if ~isempty(bestM)  
        M = bestM;
        inliers = bestinliers;
    else
        M = [];
        inliers = [];
        warning('ransac was unable to find a useful solution');
    end

    function r = isdegenerate(x)
        nd=3;
        xcomb = combnk(1:s,3);
        ncomb = size(xcomb,1);
        x1 = x(1:nd,:);   % Extract x1 and x2 from x
        x2 = x((nd+1):end,:);
        ir1=zeros(ncomb,1);
        ir2=zeros(ncomb,1);
        for i=1:ncomb
            %ir1(i)=iscolinear_n(x1(:,xcomb(i,:)));
            P=x1(:,xcomb(i,:));
            ir1(i)=norm(cross(P(:,2)-P(:,1), P(:,3)-P(:,1))) < eps ;
            P=x2(:,xcomb(i,:));
            ir2(i)=norm(cross(P(:,2)-P(:,1), P(:,3)-P(:,1))) < eps ;

        end

        r = any([ir1,ir2]);
    end

    function [inliers, H] = distfn(H, x, t,white)
        nd=3;
        lx1 = x(1:nd,:);   % Extract x1 and x2 from x
        lx2 = x((nd+1):end,:);

        % Calculate, in both directions, the transfered points    
        Hx1 = H*lx1;

        % Calculate lab distance
        luv_ref = HGxyz2luv(lx2',white)'; % reference LUV
        luv_est = HGxyz2luv(Hx1',white)'; % reference LUV

        uv_ref = bsxfun(@rdivide,luv_ref(2:3,:),max(luv_ref(1,:),eps));
        uv_est = bsxfun(@rdivide,luv_est(2:3,:),max(luv_est(1,:),eps));

        d = sqrt(sum((uv_ref-uv_est).^2,1));
        inliers = find(d<t);
    end

end