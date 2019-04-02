function [inliers,H,D] = ransacfithomography(x1, x2, t)

    mit = 5000;

    [nd,~] = size(x1);
    s = nd+1;
    fittingfn = @wrap_vgg_homographynd;
    [~,inliers] = ransac([x1;x2], fittingfn,s,t, 100, mit);
    if nargout > 1
        if numel(inliers)>=4
            [H,~,D] = H_als(x1(:,inliers),x2(:,inliers));
        else
            [H,~,D] = H_als(x1,x2);
        end
    else
        H = [];
    end

    function H = wrap_vgg_homographynd(x)
         [H] = H_als(x(1:nd,:),x((nd+1):end,:));
    end

end
