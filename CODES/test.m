im_path = 'H_test_im/';

% query image
a = imread([im_path,'1_1im.png']);
a = im2double(a);
A = reshape(a,[],3)';

% reference image
b = (imread([im_path,'1_2im.png']));
b = im2double(b);
B = reshape(b,[],3)';

% ALS
[H,~,D] = H_als(A,B);
als1 = H*A*D;
ALS_shade = reshape(als1',size(b));
als2 = H*A;
ALS = reshape(als2',size(b));

% LS
M = B/A;
ls = M*A;
LS = reshape(ls',size(b));

%comparision
figure;
ch1 = chrodist(A,128);
ch2 = chrodist(B,128);
subplot(3,1,1); imshowpair(ch1,ch2);title('original chromaticity diag of both images')

ch1 = chrodist(als2,128);
ch2 = chrodist(B,128);
subplot(3,1,2);imshowpair(ch1,ch2);title('after using H chromaticity diag of both images')

ch1 = chrodist(als1,128);
ch2 = chrodist(B,128);
subplot(3,1,3);imshowpair(ch1,ch2);title('after using both H,D chromaticity diag of both images')


figure;
subplot(2,3,1); imshow(a); title('A');
subplot(2,3,2); imshow(b); title('B');
subplot(2,3,3); imshow(ALS_shade); title('A to B with shading correction');
subplot(2,3,4); imshow(ALS); title('A to B without shading correction');
subplot(2,3,5); imshow(LS); title('A to B by Least-Squares');
