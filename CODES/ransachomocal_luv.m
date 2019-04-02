function M = ransachomocal_luv(rgb,xyz,white,rgb_u)

    b_err = Inf;
    n_trail = 100;
    for i = 1:n_trail
        t_M = alsransac_luv(rgb',xyz',white',0.2);
        %xyz_est = homocvt(rgb_u,t_M);
        H=t_M;
        p=rgb_u;
        C = eye(size(H));
        pr = reshape(p,[],size(H,1))';
        prC = C*pr; % chromaticity array
        qC = H*prC; % apply inverse homography
        xyz_est = (C\qC)'; % convert back to rgb

        t_err = mean(luv_err(xyz_est,xyz));

        if t_err<b_err
            M = t_M; b_err = t_err;
        end
    end
end

function err = luv_err(xyz_est,xyz_std)
    XYZ_est = xyz_est./xyz_est(4,2);
    
    % LUV error
    luv_est = HGxyz2luv(XYZ_est,xyz_std(4,:));
    luv_ref = HGxyz2luv(xyz_std,xyz_std(4,:)); % reference LUV
    err = sqrt(sum((luv_ref - luv_est).^2,2));
end

