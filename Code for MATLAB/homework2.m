%% INNIT
clc
clear all
close all

blockSize = 16;
searchRange = 8;

frame_1 = rgb2gray(imread('stefan/01.bmp'));
frame_2 = rgb2gray(imread('stefan/05.bmp'));

[rows, cols] = size(frame_1);

%% STEP 1
[time_dfd, min_d_dfd, motion_vectors_dfd] = EBMA(frame_1, frame_2, rows, cols, blockSize, searchRange, 'dfd');
[mp_error_dfd, rec_dfd] = prediction_error(min_d_dfd, frame_1, frame_2, searchRange);

fig = figure(1);
subplot(131); imshow(rec_dfd(1:(rows-blockSize), 1:(cols-blockSize)));title('rec image')

err_img = im2double(frame_2) - im2double(rec_dfd); 
fig = figure(1);
subplot(132); imshow(err_img(1:(rows-blockSize), 1:(cols-blockSize)));title('error image')


%% STEP 2
T = dctmtx(8);
dct = @(block_struct) T * block_struct.data * T';
B = blockproc(err_img,[8 8],dct);
coeff = zigzag(T);
mask = zeros(8,8);
mask(1,1) = 1;
mask(1,2) = 1;
mask(2,1) = 1;
mask(3,1) = 1;
B2 = blockproc(B,[8 8],@(block_struct) mask .* block_struct.data);
% perform the 8*8 inverse DCT
invdct = @(block_struct) T' * block_struct.data * T;
rec_err_img1 = blockproc(B2,[8 8],invdct);
% B = zeros(rows,cols);
for a=1:8:rows
    for b=1:8:cols
        B(a:a+7,b:b+7) = dct_2d(err_img(a:a+7,b:b+7));
    end
end
B2 = blockproc(B,[8 8],@(block_struct) mask .* block_struct.data);
invdct = @(block_struct) T' * block_struct.data * T;
rec_err_img2 = blockproc(B2,[8 8],invdct);  % 2-1
% test image
figure();
subplot(1,3,1);
imshow(rec_err_img2(1:(rows-blockSize), 1:(cols-blockSize)));title('2')
subplot(1,3,2);
imshow(B(1:(rows-blockSize), 1:(cols-blockSize)));title('2')
% Compute the PSNR of the reconstructed error image.
k = 8;
% k is the number of binary digits used by the image to represent the pixels, that is, the bit depth.
fmax = 2.^k - 1;
a = fmax.^2;
MSE =(double(im2uint8(rec_err_img2)) -double( im2uint8(err_img))).^2;
b = mean(mean(MSE));
PSNR2 = 10*log10(a/b);
fprintf('PSNR_s2');
disp(PSNR2);  % 12.058  11.9912

%% STEP 3
q = quantizer();
qB2 = quantize(q, B2);
qrec_err_img = blockproc(qB2,[8 8],invdct);

subplot(1,3,2);
imshow(rec_err_img1(1:(rows-blockSize), 1:(cols-blockSize)));
subplot(1,3,3);
imshow(qrec_err_img(1:(rows-blockSize), 1:(cols-blockSize)));

reclim_frame = im2double(rec_dfd) + rec_err_img1;
rec_frame = im2double(rec_dfd) + rec_err_img2;
fig = figure();
subplot(1,2,1);
imshow(rec_frame(1:(rows-blockSize), 1:(cols-blockSize)))
subplot(1,2,2);
imshow(reclim_frame(1:(rows-blockSize), 1:(cols-blockSize)))

fig = figure();
qreclim_frame = im2double(rec_dfd) + qrec_err_img;
subplot(1,3,1);
imshow(rec_frame(1:(rows-blockSize), 1:(cols-blockSize)));title('3')
subplot(1,3,2);
imshow(reclim_frame(1:(rows-blockSize), 1:(cols-blockSize)));title('3')
subplot(1,3,3);
imshow(qreclim_frame(1:(rows-blockSize), 1:(cols-blockSize)))


% fprintf('PSNR_s3');
% psnr(rec_frame, frame_2)

% Compute the PSNR of the reconstructed error image.
k = 8;
% k is the number of binary digits used by the image to represent the pixels, that is, the bit depth.
fmax = 2.^k - 1;
a = fmax.^2;
MSE =(double(im2uint8(rec_frame)) -double( im2uint8(reclim_frame))).^2;
b = mean(mean(MSE));
PSNR3 = 10*log10(a/b);
fprintf('PSNR_s3');
disp(PSNR3);  % 12.058

%% STEP 4

qrec_err_img = round(rescale(qrec_err_img,-255,255));

symbol = (-255:255);

p2 = histogram(qrec_err_img, 511, 'Normalization', 'probability');
imgEntropy = CalcImgEntropy(im2uint8(qrec_err_img), 2 );  % 0.9951  entropy step1 errimg 4.45
%[symbol,p] = hist(qreclim_frame(:)', double(unique(qreclim_frame)));

%p = softmax(p.Values);

p = p2.Values';

[dict,avglen] = huffmandict(symbol, p);  % dict : histogram  avglen : codelength

comp = huffmanenco(qrec_err_img(:)', dict);

qreclim_frame2 = huffmandeco(comp, dict); % vector

qreclim_frame2 = rescale(reshape(qreclim_frame2,[rows,cols]), -1, 1);
qrec_err_img = rescale(qrec_err_img, -1, 1);
% Compute the compression ratio. (the ratio between the data size before and after compression)
% quality = 100;
% [ output_image, compressed_vector, ratio,entropy_out ] = jpeg_computing(qrec_err_img, quality);
% fprintf('ratio',ratio);  % 23.7394


fig = figure();
subplot(1,2,1);
imshow(qrec_err_img);
subplot(1,2,2);
imshow(qreclim_frame2);

% [M N]=size(qrec_err_img);
% fun1=zigzag(qrec_err_img);
% [count1,x] = imhist(fun1); %// Change
% p1 = count1/ numel(fun1);
% [dict1,avglen1]=huffmandict(x,p1); %// Change
% comp1= huffmanenco(fun1,dict1); 
% Im1 = huffmandeco(comp1,dict1);


%% step 5
err_img2 = qrec_err_img;
[err_rows, err_cols] = size(qreclim_frame2);
qreclim_frame2 = mat2gray(qreclim_frame2);
[time_dfd5, min_d_dfd5, motion_vectors_dfd5] = EBMA(qreclim_frame2,err_img2, err_rows, err_cols, blockSize, searchRange, 'dfd');
[mp_error_dfd_5, rec_dfd_5] = prediction_error(min_d_dfd5, qreclim_frame2,err_img2, searchRange);

fig = figure();
subplot(131); imshow(rec_dfd_5(1:(err_rows-blockSize), 1:(err_cols-blockSize)));title('rec image')

final_err_img = im2double(err_img2) - im2double(rec_dfd_5); 
% final_err_img = mat2gray(final_err_img, [-255, 255]);
fig = figure();
subplot(132); imshow(final_err_img(1:(err_rows-blockSize), 1:(err_cols-blockSize)));title('final error image')

% Compute the PSNR of the reconstructed error image.
k = 8;
% k is the number of binary digits used by the image to represent the pixels, that is, the bit depth.
fmax = 2.^k - 1;
a = fmax.^2;
MSE =(double(im2uint8(err_img2)) -double( im2uint8(final_err_img))).^2;
b = mean(mean(MSE));
PSNR5 = 10*log10(a/b);
fprintf('PSNR_s5');
disp(PSNR5);  % 12.058
