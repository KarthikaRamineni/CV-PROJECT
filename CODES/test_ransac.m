clear all;
tpath = 'H_test_im/';
gamma = 1;
oim = imread([tpath,'1_1im.png']);
rim = imread([tpath,'1_2im.png']);
oim = im2double(oim);
rim = im2double(rim);
H = H_from_feat(oim,rim,320,4);