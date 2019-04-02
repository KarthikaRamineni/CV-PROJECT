% configrationfff
dbpath = 'HG_ColourChecker/'; % path of rawdata
fmethod = {@alshomocal;@ransachomocal_luv;@lscal;@rpcal};

method=[1,2,3,4];


dirName=[dbpath,'patch_real'];
dirData = dir(dirName);     
dirIndex = [dirData.isdir]; 
idx = ~dirIndex;
if isempty(idx)
    list = {};
else
    list = {dirData(idx).name}'; 
end
fn=list;
Npic = numel(fn);
Npatch = 24;
Nmethod = numel(fmethod);

% non-uniform shading errors

errorlab_n = zeros(Npatch,Npic,Nmethod);
errorluv_n = zeros(Npatch,Npic,Nmethod);
errorrgb_n = zeros(Npatch,Npic,Nmethod);

% uniform shading errors
errorlab_u = zeros(Npatch,Npic,Nmethod);
errorluv_u = zeros(Npatch,Npic,Nmethod);
errorrgb_u = zeros(Npatch,Npic,Nmethod);

% relative difference of correction matrix
md = zeros(1,Npic,Nmethod);

for i = 1:Npic
    
    % ref cat
    cat = regexp(fn{i},'^[^_]+','match');
    % load data
    load([dbpath,'patch_real/',fn{i}]);
    % load reference data
    load([dbpath,'ref_real-',cat{1},'.mat']);

    xyz_std = ref.XYZ./ref.XYZ(4,2); % refernece XYZ
    lab_ref = HGxyz2lab(xyz_std,xyz_std(4,:)); % reference LAB
    luv_ref = HGxyz2luv(xyz_std,xyz_std(4,:)); % reference LUV
    rgb_ref = xyz2rgb(xyz_std,'WhitePoint',xyz_std(4,:)); % reference LUV

    % compute colour correction matrix
    fsv = reshape(cap.sv,[],3);
    fsv_uniform = reshape(cap.sv_uniform,[],3);
    for m = 1:Nmethod
        % compute the color correction transform
        %switch Method(m)
        
        if (method(m)==2)
            M_n = fmethod{m}(fsv,xyz_std,xyz_std(4,:),fsv_uniform);
            M_u = fmethod{m}(fsv_uniform,xyz_std,xyz_std(4,:),fsv_uniform);
        else
            M_n = fmethod{m}(fsv,xyz_std);
            M_u = fmethod{m}(fsv_uniform,xyz_std);
        end

        % compute xyz using the ground truth RGBs
        %switch Method{m}
        if (method(m)==1 || method(m) ==2)
            xyz_est_n = homocvt(fsv_uniform,M_n);
            xyz_est_u = homocvt(fsv_uniform,M_u);
        elseif(method(m)==4)
            xyz_est_n = M_n.cfun(fsv_uniform,M_n.matrix,M_n.terms);
            xyz_est_u = M_u.cfun(fsv_uniform,M_u.matrix,M_u.terms);
        else
            xyz_est_n = fsv_uniform*M_n;
            xyz_est_u = fsv_uniform*M_u;
        end

        % normalize by a white patch's green intensity
        %XYZ_est_n = xyz_est_n./xyz_est_u(4,2);
        %XYZ_est_u = xyz_est_u./xyz_est_u(4,2);
        XYZ_est_n = xyz_est_n./xyz_est_n(4,2);
        XYZ_est_u = xyz_est_u./xyz_est_u(4,2);
        
        % Evaluation
        % non-uniform test DE LAB
        lab_est_n = HGxyz2lab(XYZ_est_n,xyz_std(4,:));
        errorlab_n(:,i,m) = sqrt(sum((lab_ref - lab_est_n).^2,2)); 
        % uniform test DE LAB
        lab_est_u = HGxyz2lab(XYZ_est_u,xyz_std(4,:));
        errorlab_u(:,i,m) = sqrt(sum((lab_ref - lab_est_u).^2,2));

        % LUV error
        luv_est_n = HGxyz2luv(XYZ_est_n,xyz_std(4,:));
        luv_est_u = HGxyz2luv(XYZ_est_u,xyz_std(4,:));
        errorluv_n(:,i,m) =sqrt(sum((luv_ref - luv_est_n).^2,2)); 
        errorluv_u(:,i,m) = sqrt(sum((luv_ref - luv_est_u).^2,2));
        
        % RGB error
        rgb_est_n = xyz2rgb(XYZ_est_n,'WhitePoint',xyz_std(4,:));
        rgb_est_u = xyz2rgb(XYZ_est_u,'WhitePoint',xyz_std(4,:));
        errorrgb_n(:,i,m) = sqrt(sum((rgb_ref - rgb_est_n).^2,2));
        errorrgb_u(:,i,m) = sqrt(sum((rgb_ref - rgb_est_n).^2,2));
    end
end

% print evaluation results (non-uniform shading)
trgb_n = eval_table(errorrgb_n,'RGB (Non-Uniform)');
t76_n = eval_table(errorlab_n,'DeltaE LAB 1976 (Non-Uniform)');
tluv_n = eval_table(errorluv_n,'DeltaE LUV (Non-Uniform)');


trgb_u = eval_table(errorrgb_u,'RGB (Uniform)');
t76_u = eval_table(errorlab_u,'DeltaE LAB 1976 (Uniform)');
tluv_u = eval_table(errorluv_u,'DeltaE LUV (Uniform)');

