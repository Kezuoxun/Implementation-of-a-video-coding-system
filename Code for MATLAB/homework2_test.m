%%
clear all;
clc;
%% Set up
Img1_Orig = imread('stefan/01.png');
IMG1_Gray = rgb2gray(Img1_Orig);
figure(1);
subplot(2,1,1);
imshow(IMG1_Gray);
title("01.png");
IMG2_Orig = imread('stefan/05.png');
IMG2_Gray = rgb2gray(IMG2_Orig);
figure(1);
subplot(2,1,2);
imshow(IMG2_Gray);
title("04.png");
mbsize = 16; %block size
p = 8; % searchRange
CountError_Algorithm = "MAD";
%%  step1 EBMA & Three step method  
% calling the motion estimation function to get the motion vectors
timer=tic;
[y] = motionEst_ES(IMG1_Gray,IMG2_Gray,mbsize, p,CountError_Algorithm); %EBMA
fprintf('time: %f\n',toc(timer));%timer
%% step1 EBMA output
[imgComp] = motionComp(IMG1_Gray, double(y), mbsize);
% calculating vertical entropy
Entropy_Y = entropy(y(1,:));
fprintf('Entropy Vertical: %f\n', Entropy_Y);
z=10;
% calculating horizontal entropy
Entropy_X = entropy(y(2,:));
fprintf('Entropy Horizontal: %f\n', Entropy_X);
[row,col] = size(IMG1_Gray);

% Calculating motion Compensated Frame Difference, its MSE and its entropy value
%MCF
figure(4);
imshow(imgComp);
title("imgComp");
%imwrite(uint8(imgComp),"matrix.png");

% Calculating MCFD MSE value and entropy value below
Z = immse(imgComp,IMG2_Gray);
%Z = Z/z;
fprintf('Mean Square Error Frame Difference(imgComp-IMG2_Gray): %f\n', Z);
peaksnr = psnr(imgComp,IMG2_Gray);%PSNR
fprintf('PSNR:%f\n',peaksnr);
fprintf('Entropy Frame Difference: %f\n', entropy(imgComp-IMG2_Gray));
%% costransform  step2-1
errorimg = double(IMG2_Gray) - double(imgComp);
cosmatrix = costrans(errorimg);
% mask(1,1) = 1;
% mask(1,2) = 1;
% mask(2,1) = 1;
% mask(3,1) = 1;
z = zig(cosmatrix);
I = invcos(z);  % IDCT

figure(5);
imshow(I);
peaksnr = psnr(double(I),errorimg,255);  %PSNR 須從double資料型態轉回unit8才可運
fprintf('PSNRerrorimg(2-1):%f\n',peaksnr);
%% quentize  step3
q = quentize(z);
dequimg = invcos(q);
I3 = rescale(dequimg);
peaksnr = psnr(dequimg,errorimg,511);  %PSNR 須從double資料型態轉回unit8才可運算
fprintf('PSNRerrorimg(3-1):%f\n',peaksnr);
figure(6);
imshow(I3)  % output 2.1
%% recontruction Step5
I5 = mat2gray(I);
% a = uint8(I5) + uint8(dequimg);  
a = imgComp + uint8(dequimg);  % ori is : imgComp
peaksnr = psnr(a,IMG2_Gray,511); %PSNR 須從double資料型態轉回unit8才可運算
fprintf('PSNRrecon(5):%f\n',peaksnr);
figure();
imshow(a);
%% Step4

%% Step 4-a
%計算Error影像中每個像素出現的次數
Error = errorimg;
for i=-255:1:255
    Error_Pixel_Amount(i+256) = length(find(Error==i));                 %統計每個pixel出現的次數 
    Error_Pixel_Frequency(i+256) = (Error_Pixel_Amount(i+256))/101376;  %統計每個pixel出現的機率 
    Record_Pixel_Frequency(i+256,1) = i;
    Record_Pixel_Frequency(i+256,2) = Error_Pixel_Frequency(i+256);
end

x=-15:1:15;   %決定histogram橫軸
figure; histogram(Error,x,'Normalization','probability');

Entropy = entropy(Error);  %計算entropy

%% Step 4-b 
%全部個機率做排序，並且pixel值要對應到機率值!!
for i = 511:-1:1
    for j =1:1:i-1
            if (Record_Pixel_Frequency(j,2) < Record_Pixel_Frequency(j+1,2))
                temp=Record_Pixel_Frequency(j,2);
                Record_Pixel_Frequency(j,2)=Record_Pixel_Frequency(j+1,2);
                Record_Pixel_Frequency(j+1,2)=temp;
                temp2=Record_Pixel_Frequency(j,1);
                Record_Pixel_Frequency(j,1)=Record_Pixel_Frequency(j+1,1);
                Record_Pixel_Frequency(j+1,1)=temp2;
            end
     end
end

S=sum(Record_Pixel_Frequency(16:511,2));           %機率第16大~511全部相加當作第16個
Huffman_Sequence=Record_Pixel_Frequency(1:15,2);     %把機率前15大的值存到矩陣
Huffman_Sequence(16,1)=S;                            %把機率16~511總和存到矩陣  
Huffman_Sequence=sort(Huffman_Sequence);       %Huffman機率排序(由小到大)

HT=hufftree(Huffman_Sequence);  %呼叫function 
prhcode(HT,[]);  %呼叫function

%% ESCAPE symbols
for i=1:1:486
    length_coding(i,1)=i+15;
    s=dec2bin(i,9);   %轉換成binary，並且至少9bit
    ESCAPE(i,1)=strcat('000001',cellstr(s));
end

for i=1:1:511
    Record_Pixel_Frequency(i,3)=Record_Pixel_Frequency(i,2)*101376;  %計算每個pixel出現次數
end

for i=1:1:1 %前五大機率的bit數分配為3
    Record_Pixel_Frequency(i,4)=1;
end

for i=2:1:3 %第6~8大機率的bit數分配為4
    Record_Pixel_Frequency(i,4)=3;
end

for i=4:1:5 %第9~12大機率的bit數分配為5
    Record_Pixel_Frequency(i,4)=4;
end

for i=6:1:8 %第13~15大機率的bit數分配為6
    Record_Pixel_Frequency(i,4)=5;
end

for i=9:1:12 %第13~15大機率的bit數分配為6
    Record_Pixel_Frequency(i,4)=6;
end

for i=13:1:15 %第13~15大機率的bit數分配為6
    Record_Pixel_Frequency(i,4)=7;
end

for i=12:1:15 %第13~15大機率的bit數分配為6
    Record_Pixel_Frequency(i,4)=7;
end

for i=16:1:511 %第16~511大機率的bit數分配為12
    Record_Pixel_Frequency(i,4)=12;
end

%% Step 4-c
[height, width] = size(errorimg);
% % step 4-3
% B=length(comp);  %原始图像比特长度
% sumcode=length(qreclim_frame_decode);  %编码后比特长度
% CR=sumcode/B;  %计算压缩率

code_bits=sum(Record_Pixel_Frequency(1:511,4));
Origin_bits=8*height*width;

ratio = Origin_bits/code_bits;


%% Step 5
Reconstru_image = inp_1'+Quantify_IDCT;
figure;imshow(Reconstru_image,[]);title('重建影像');
diff=abs(inp_1'-(Reconstru_image)); 
MSE=(sum(sum(diff.*diff)))/(height*width);
PSNR3=10*log10((255*255)/MSE);





%% 副程式 
function HT = hufftree(Huffman_Sequence)
  tree = cell(length(Huffman_Sequence),1);   %產生與陣列p相同長度的cell
for i=1:length(Huffman_Sequence)
    tree{i}=i;                             %每個節點的代號初始值依序為1,2,3...
end

  while numel(tree)>2
     [Huffman_Sequence, pos] = sort(Huffman_Sequence); %從小排到大   
     Huffman_Sequence(2)=Huffman_Sequence(1)+Huffman_Sequence(2);                %將最小及次小的值相加,存到次小位置
     Huffman_Sequence(1)=[];        %刪除最小值的元素,因此陣列p的個數少一個
     tree = tree(pos);              %節點依照p值大小排序
     tree{2}={tree{1},tree{2}};     %最小及次小的節點合併，存至次小節點內
     tree(1)=[];                          %刪除最小的節點
  end
  HT = tree;
end  

function prhcode(HT,code)
   if isa(HT,'cell')  %判斷tree是否為cell
      prhcode(HT{1},[code '0']);
      prhcode(HT{2},[code '1']);
   else
         mystr = strcat(num2str(HT) , '=' , code);
         disp(mystr);
   end
end

% [rows, cols] = size(errorimg);  % use step 1 error image
% qrec_err_img = round(rescale(errorimg,-255,255));
% img_entropy = entropy(qrec_err_img) ;
% 
% symbol = (-255:255);
% 
% % step 2
% p2 = histogram(qrec_err_img, 511, 'Normalization', 'probability');
% % step 4-2 test
% p = p2.Values';  % get probability of symbol
% symbols=find(p~=0);
% p=p(symbols);
% [p,sortindex]=sort(p);   
% 
% %将符号按照出现的概率大小排序
% 
% symbols=symbols(sortindex);
% 
% len=length(symbols);
% 
% symbols_index=num2cell(1:len);
% 
% codeword_tmp=cell(len,1);
% 
% while length(p)>1     %生产Huffman树，得到码字编码表
% 
%     index1=symbols_index{1};
% 
%     index2=symbols_index{2};
% 
%     codeword_tmp(index1)=addnode(codeword_tmp(index1),uint8(0));
% 
%     codeword_tmp(index2)=addnode(codeword_tmp(index2),uint8(1));
% 
%     f=[sum(p(1:2)) f(3:end)];
% 
%     symbols_index=[{[index1,index2]} symbols_index(3:end)];
% 
%     [f,sortindex]=sort(f);
% 
%     symbols_index=symbols_index(sortindex);
% 
% end
% 
% %
% 
% [dict,avglen] = huffmandict(symbol, p); % dict : histogram  avglen : codelength
% 
% comp = huffmanenco(qrec_err_img(:)', dict);  % encode  use 1D image
% 
% qreclim_frame_decode = huffmandeco(comp, dict);  % decode 
% 
% qreclim_frame2 = rescale(reshape(qreclim_frame_decode,[rows,cols]), -1, 1);  % vector to image
% qrec_err_img = rescale(qrec_err_img, -1, 1);
% 
% % step 4-3
% B=length(comp);  %原始图像比特长度
% sumcode=length(qreclim_frame_decode);  %编码后比特长度
% CR=sumcode/B;  %计算压缩率
% 
% fig = figure();
% subplot(1,2,1);
% imshow(qrec_err_img);
% subplot(1,2,2);
% imshow(qreclim_frame2);





