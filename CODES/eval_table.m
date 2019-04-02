function tab = eval_table(de,de_str)
    %   Copyright 2018 Han Gong <gong@fedoraproject.org>
    %   University of East Anglia
    disp_flag = true;
    Method = {'ALS';'ALS_RANSAC';'LS';'RP'};
    summary = @(f,x) squeeze(mean(f(x,1),2));

    de_mean = summary(@mean,de);
    de_median = summary(@median,de);
    de_95 = summary(@(x,d) quantile(x,.95,d),de);
    de_max = summary(@(x,d) max(x,[],d),de);

    Variable = {'Mean','Median','pct95','Max'};
    tab = table(de_mean,de_median,de_95,de_max,'RowNames',Method,'VariableNames',Variable);

    if disp_flag
        disp(de_str);
        disp(tab);
    end
end
