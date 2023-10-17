function [] = tb_qrcodeDec_v2(inputIMG)
%% QR Code decoder
% - This example shows how to decode QR Code encoded in version2 format.
% *QR Code is a registered trademark of DENSO WAVE INCORPORATED

% �o�[�W����2�A���������x��L��QR�R�[�h�𕜍����闬����m�F���邽�߂�
% �T���v���R�[�h�ł�
% *QR�R�[�h��(��)�f���\�[�E�F�[�u�̓o�^���W�ł�
%clear all, close all, clc;

%% �摜�f�[�^�ǂݍ���
% - Read Image

% �\�ߍ쐬���Ă������摜�f�[�^��ǂݍ��݂܂��BwebCam�Ŏ΂߂���B�e���Ă���A
% �c��ł��܂��B�܂��AQR�R�[�h��ɂ�MathWorks�̃��S���`���ꂽ�����u����A
% �R�[�h�̈ꕔ���B��Ă����Ԃł�
%img = imread('mwprinciple_webcam.png');
%img = imrotate(img, 180);
img=inputIMG;
imshow(img)

%% �摜�f�[�^��2�l��
% - Convert RGB to grayscale image and binarize
bw = imbinarize(rgb2gray(img));
% imshow(bw)

%% QR�R�[�h�̈ʒu���o�p�p�^��(�t�@�C���_)���o
% - Find where the finder patterns are
[Centroid, bw2, flag, width, bbox] = detectFinder(bw);

%% �t�@�C���_�̒��S3���W�𒊏o
% - Extract XY coordinates of center of finder patterns
markers = Centroid(logical(flag),:);

%% �c�ݕ␳�p�A���C�����g�p�^�[�����o
% - Find where the alignment pattern is
bwa = imclearborder(bw);
bwa = bwpropfilt(bwa, 'EulerNumber', [0, 0]);
bwa = imfill(bwa, 'holes');
stats = regionprops(bwa, 'Area', 'Centroid');
area = [stats.Area];
[mv, i] = min(area);
alp = stats(i).Centroid;

%% �t�@�C���_&�A���C�����g�p�^�[���̈ʒu������
% - Visualize finder and alignment patterns

% �t�@�C���_�̈ʒu�����p�}�X�N
bwf = imfill(bw2, 'holes');

% �A���C�����g�p�^�[���̈ʒu�����p�}�X�N
bwl = bwlabel(bwa);
bwl = bwl == i;
bwlf = bwf | bwl;
mask3 = repmat(bwlf, [1 1 3]);
imgf = img;
% �t�@�C���_&�A���C�����g�p�^�[���ȊO�̋P�x�l�𗎂Ƃ�
imgf(~mask3) = imgf(~mask3) * 0.2;
% imshow(imgf)

%% 3�̃t�@�C���_�̈ʒu�֌W���Z�o���A�ǂ̂悤�ȕ�������QR�R�[�h���B�e���ꂽ�̂���c�����܂�
% - Understand positional relation of each patterns

[idxC, idxX, idxY] = finderPos(markers);
% 3�R�[�i�[�ɔz�u�����t�@�C���_�̒��ŁA���㕔�ɔz�u�����t�@�C���_�̈ʒu��
% �ԐF�}�[�J��z�u
imgf = insertShape(imgf, 'FilledCircle', [markers(idxC, :) 10], 'Color', 'red');
%imshow(imgf)
% X��(����)�����ɑ΂���p�x���ŏ��ƂȂ�ӂ𐬂��t�@�C���_�̈ʒu��
% �F�}�[�J��z�u
imgf = insertShape(imgf, 'FilledCircle', [markers(idxX, :) 10], 'Color', 'blue');
%imshow(imgf)
% Y��(����)�����ɑ΂���p�x���ŏ��ƂȂ�ӂ𐬂��t�@�C���_�̈ʒu��
% ���F�}�[�J��z�u
imgf = insertShape(imgf, 'FilledCircle', [markers(idxY, :) 10], 'Color', 'yellow');
%imshow(imgf)

%% �c��ł���QR�R�[�h���􉽊w�I�ϊ��ɂ��␳���܂�
% - Create x and y coordinates of control points for geometric transformation

[Cref, Xref, Yref, Aref] = createControlPts(idxC, idxX, idxY, width, width, markers);

post = [Cref; Xref; Yref; Aref];
org = [markers(idxC, :); markers(idxX, :); markers(idxY, :); alp];

%% �􉽊w�ϊ��p�s����Z�o
% - Fit geometric transformation to control point pairs
T = fitgeotrans(org, post,'projective'); % projective2d �N���X�i�ϊ��s��Ǝ����j
T.T             % �������ꂽ�s��̊m�F

%% �􉽊w�I�ϊ���A���ʂ�\��
% - Apply geometric transformation to image

bwt = imwarp(bw, T);
% imshow(bwt); % �f�t�H���g�͐��`���

%% �ׂ����m�C�Y���������܂�
% - Remove noises
bwti = ~bwt;
%Iinv = imclose(Iinv, strel(ones(2)));
bwti = imopen(bwti, strel(ones(3)));
bwti = imclearborder(bwti);
bwt = ~bwti;
%imshow(bwt)

%% �c�ݕ␳�����������̂ŁA���߂ăt�@�C���_�̈ʒu�����߂܂�
% - Find where the finder patterns are for transposed image
[Centroid, bw, flag, width, ~] = detectFinder(bwt);

%% 3���W�𒊏o
% - Extract XY coordinates of center of finder patterns
markers = Centroid(logical(flag),:);

%% 3�̃t�@�C���_�̈ʒu�֌W���Z�o���܂�
% - Understand positional relation of each patterns
[idxC, idxX, idxY] = finderPos(markers);

%% ���W�ʒu�ɉ����ĉ摜��]
% - Rotate image and coordinates to locate the corner of finder patterns to
% top left
[idxC, idxX, idxY, bwt] = rotPosition(idxC, idxX, idxY, bwt, markers);
%imshow(bwt)

%% �{���ł�25x25�Z���̃o�[�W����2 QR�R�[�h��p���Ă��܂����A�摜�ɂ�����1�Z��������̑傫�������߂܂�
% - Calculate the size of cell for version2 QR Code in an image

% 1�Z���̍���������
% �t�@�C���_�}�[�J�[�Ԃ̋�����18�Z���ł���̂ŁAregionprops�ŋ��߂����W��
% ������18�Ŋ���A1�Z���Ԃ̍��������߂܂�
%xy = round(markers(idxC,:));
xy = round(idxC);
%height = (markers(idxY, 2) - xy(2)) / 18;
height = (idxY(1, 2) - xy(2)) / 18;
% ���l�ɂ���1�Z���̕������߂܂�
width = (idxX(1, 1) - xy(1)) / 18;

% QR�R�[�h�ɖ��ߍ��܂�Ă���^�C�~���O�p�^�[���𗘗p����������萳�m�ł����A
% �ȗ����̂��߂ɖ{���ł͊������܂��B

%% ���ɁA�e�Z���̒��S���W��\���z����쐬���܂�
% - Create array that represents x and y coordinates of each cells

%QR�R�[�h : �o�[�W����2�͏c��25�Z��
xref = xy(1) - width*3;
yref = xy(2) - height*3;
xtable = (1:25) * width - width;
xtable = round(xtable + xref);
ytable = (1:25) * height - height;
ytable = round(ytable + yref);

%% ���߂����S���W������ɁA�t�H�[�}�b�g����ǂݎ��܂�
% - Read format information to understand error correction level and mask
% pattern used in QR Code

% �f�R�[�h���s�����߂ɕK�v�Ȍ��������x����A�{���ꂽ�}�X�N�p�^���̏��
% �L�q����Ă��܂��B
bwo = ~bwt;
fmt1 = bwo(ytable(9), xtable([1:6,8]));
fmt2 = bwo(ytable([1:6, 8:9]), xtable(9));

%fmt2 = bwo(xtable([1:6,8:9]), ytable(9));
%fmt1 = bwo(xtable(9), ytable([1:6, 8]));

fdata = [fmt1, fliplr(fmt2')];

%% �t�H�[�}�b�g���̕���
% - Decode format information using BCH decoder

% �t�H�[�}�b�g���ɂ̓}�X�N�p�^�����{����Ă��邽�߁Axor���Ƃ�܂��B
bitp = logical([1 0 1 0 1 0 0 0 0 0 1 0 0 1 0]);
fdata = xor(fdata, bitp);
% BCH�����ɂ�蕄��������Ă���̂ŁA�������܂�
finfo = bchdec(gf(fdata), 15, 5);
minfo = finfo(3:5); %�}�X�N�p�^�����

%% ���ɁA�f�[�^��ǂݎ�鏇�Ԃ������C���f�b�N�X�e�[�u�����쐬���܂�
% - Create index table that represents how to read data bits and error
% correction bits

% �ǂݏo���͉E���p���珇�ɁA�ׂ荇����ƌ��݂ɓǂݏo���Ă����A
% �R�[�h�̏�[�������͉��[�ɓ��B�����ꍇ�͓ǂݏo��������]�����܂��B

% ��e�[�u���쐬
xidx = 0:24;
xidx = repmat(xidx, [25, 1]);
yidx = 0:24;
yidx = repmat(yidx', [1, 25]);
idxtbl = cat(3, xidx, yidx);
% �f�[�^�ǂݎ��ӏ��C���f�b�N�X�Z�o
[xidx yidx reftable] = readtable(28+16);
% �ǂݏo����������Ȃ̂ŁA�f�[�^�ƒ����������i�[����Ă���̈��
% �������Ċm�F���܂��B
reftable = label2rgb(reftable);
%figure, imshow(reftable)
%title('QR code layout : Version2, Error correction level = L')

%% �쐬�����C���f�b�N�X�ɏ]���A�f�[�^��ǂݎ��܂�
% - Read data and error correction bits according to the table
data = bwo(ytable(yidx), xtable(xidx));
%data = bwo(xtable(xidx), ytable(yidx));

data = diag(data);
mx = idxtbl(:,:,1);
mx = mx(yidx, xidx);
mx = diag(mx);
my = idxtbl(:,:,2);
my = my(yidx, xidx);
my = diag(my);

% �t�H�[�}�b�g���ɋL�ڂ��ꂽ�}�X�N�p�^���̏�����ɁA
% �Y���ӏ��̃f�[�^�𔽓]���܂��B
mask = cat(3, mx, my);
invidx =  maskGen(mask, minfo);
data(invidx) = ~data(invidx);

%% �ǂݏo�����r�b�g��𕜍����Ă݂܂�
% - Try to decode data without error correction

% version2, ���������x��L�̏ꍇ�A�������͍ő�28�A
% ������������16�ƂȂ�܂��B

% �f�[�^�̐擪4bit��QR�R�[�h�̃��[�h�w���q�ɂȂ��Ă��܂��B
% �ȍ~�̃f�R�[�h�ɕK�v�Ȃ��߁A�܂����[�h���m�F���܂��B
% 1 : �������[�h
% 2 : �p�������[�h
% 3 : 8bit���[�h
% 4 : �������[�h
data = num2str(data)';
mode = data(1:4);
mode = bin2dec(mode);

% �����ł́A8bit���[�h�̂ݎ��グ�܂��B
if mode == 4 % 8bit���[�h
    %�������w���q�̕���
    length = data(5:12);
    length = bin2dec(length);
    D = data(13:end);
    text = cell(1, length);
    for i = 1:length
        tmp = D((i-1)*8+1:i*8);
        text{1, i} = char(bin2dec(tmp));
    end
end
% �����������ʂ�\��
text = strjoin(text)

% ���炩�̕����񂪊m�F�ł��܂������AQR�R�[�h��ɒu���ꂽ���ɂ����
% �ꕔ���B��Ă������߁A�����񂪐������ǂ߂Ă��Ȃ��\������ł��B

%% ���ʑ}��(����������)
% - Insert decoded data into image and visualize it
I1 = insertObjectAnnotation(img, 'rectangle', bbox.BoundingBox,...
    text, 'TextBoxOpacity',0.9,'FontSize',26,'LineWidth',3);
figure, imshow(I1)

%% �����������𗘗p���A������������������܂�
% - Decode data bits with error correction

% QR�R�[�h�ł̓u���b�N���ɋ���RS(���[�h�\������)�����p����Ă��܂��B
% �܂��A2�i�̃r�b�g���10�i�ɕϊ����A�����₷�����܂��B
data2 = zeros(1, 44);
for i = 1:44
    tmp = data((i-1)*8+1:i*8);
    data2(1, i) = bin2dec(tmp);
end
%�f�[�^�ƒ��������ɕ������܂��B
data = data2(1, 1:28);
tail = data2(1, 29:end);

%%
% RS(���[�h�\������)��������̃R���t�B�O���[�V����
m = 8; % 1�V���{��������̃r�b�g��
n = 2^m-1; % �����꒷
k = n-size(tail, 2); % ���b�Z�[�W��

hDec = comm.RSDecoder(n,k);
hDec.GeneratorPolynomialSource = 'Property'; % �����̐����������̓��[�U�[�w��
hDec.GeneratorPolynomial       = rsgenpoly(n,k,[],0); % ����������

% �Z�k�����ł��邽�߁A�O�����[���p�f�B���O
shortmsg = [zeros(1,n-size(data, 2)-size(tail,2)) data tail]';

decoded = step(hDec, shortmsg);
decoded = decoded(size(decoded,1)-size(data,2)+1:size(decoded,1),1);

% 10�i>2�i�ϊ�
tmp = 0;
for i = 1:28
    tmp = [tmp, dec2bin(decoded(i),8)];
end
databox = tmp(1,2:end);

mode = databox(1:4);
mode = bin2dec(mode);
% �����ł́A8bit���[�h�̂ݎ��グ�܂��B
if mode == 4 % 8bit���[�h
    %�������w���q�̕���
    length = databox(5:12);
    length = bin2dec(length);
    D = databox(13:end);
    text = cell(1, length);
    for i = 1:length
        tmp = D((i-1)*8+1:i*8);
        text{1, i} = char(bin2dec(tmp));
    end
end
% �����������ʂ�\��
text = strjoin(text)

%% ���ʑ}��(����������)
% - Insert decoded data into image and visualize it
I2 = insertObjectAnnotation(img, 'rectangle', bbox.BoundingBox,...
    text, 'TextBoxOpacity',0.9,'FontSize',26,'LineWidth',3,'Font','Calibri');
figure, imshow(I2)
end
%%
% Copyright 2019 The MathWorks, Inc.