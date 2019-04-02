function [M, inliers] = ransac(x, fittingfn, s, t,maxDataTrials, maxTrials)
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
            [inliers, M] = distfn( M, x, t);
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

    function [inliers, H] = distfn(H, x, t)
        nd=3;

        x1 = x(1:nd,:);   % Extract x1 and x2 from x
        x2 = x((nd+1):end,:);     
        Hx1    = H*x1;
        invHx2 = H\x2;    
        %x1     = (x1);
        %x2     = (x2);
        Hx1    = hnormalise(Hx1);
        invHx2 = hnormalise(invHx2); 

        d2 = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
        inliers = find(abs(d2) < t);
    end

    function nx = hnormalise(x)

        [rows,npts] = size(x);
        nx = x;
        finiteind = find(abs(x(rows,:)) > eps);
        for r = 1:rows-1
        nx(r,finiteind) = x(r,finiteind)./x(rows,finiteind);
        end
        nx(rows,finiteind) = 1;
    end  
end