function [H,inliers] = alsransac_luv(x1, x2, white, t)

    mit = 500;
    [nd,~] = size(x1);
    s = nd+1;
    
    fittingfn = @wrap_als;
    [~,inliers] = ransacforluv([x1;x2], fittingfn, s, t, 100, mit,white);
    
    if numel(inliers)>=4
        H = H_als(x1(:,inliers),x2(:,inliers));
    else
        H = H_als(x1,x2);
    end
    
    function H = wrap_als(x)
        H = H_als(x(1:nd,:),x((nd+1):end,:));
    end

end
