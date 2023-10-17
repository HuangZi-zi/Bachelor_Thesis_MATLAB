clear;

% Read Image
img = imread('Resource\mwprinciple_webcam.png');
%img = imrotate(img, 180);
imshow(img)

% Convert RGB to grayscale image and binarize
bw = imbinarize(rgb2gray(img));
imshow(bw)

% Find where the finder patterns are
[Centroid, bw2, flag, width, bbox] = detectFinder(bw);


% Extract XY coordinates of center of finder patterns
markers = Centroid(logical(flag),:);

% Find where the alignment pattern is
bwa = imclearborder(bw);
bwa = bwpropfilt(bwa, 'EulerNumber', [0, 0]);%只保留开区域
bwa = imfill(bwa, 'holes');
stats = regionprops(bwa, 'Area', 'Centroid');
area = [stats.Area];
[mv, i] = min(area);
alp = stats(i).Centroid;

% Visualize finder and alignment patterns

% ファインダの位置可視化用マスク
bwf = imfill(bw2, 'holes');

% アライメントパターンの位置可視化用マスク
bwl = bwlabel(bwa);
bwl = bwl == i;
bwlf = bwf | bwl;
mask3 = repmat(bwlf, [1 1 3]);
imgf = img;
% ファインダ&アライメントパターン以外の輝度値を落とす
imgf(~mask3) = imgf(~mask3) * 0.2;
imshow(imgf)

% - Understand positional relation of each patterns

[idxC, idxX, idxY] = finderPos(markers);
% 3コーナーに配置されるファインダの中で、左上部に配置されるファインダの位置に
% 赤色マーカを配置
imgf = insertShape(imgf, 'FilledCircle', [markers(idxC, :) 10], 'Color', 'red');
imshow(imgf)
% X軸(水平)方向に対する角度が最小となる辺を成すファインダの位置に
% 青色マーカを配置
imgf = insertShape(imgf, 'FilledCircle', [markers(idxX, :) 10], 'Color', 'blue');
imshow(imgf)
% Y軸(垂直)方向に対する角度が最小となる辺を成すファインダの位置に
% 黄色マーカを配置
imgf = insertShape(imgf, 'FilledCircle', [markers(idxY, :) 10], 'Color', 'yellow');
imshow(imgf)


% - Create x and y coordinates of control points for geometric transformation

[Cref, Xref, Yref, Aref] = createControlPts(idxC, idxX, idxY, width, width, markers);

post = [Cref; Xref; Yref; Aref];
org = [markers(idxC, :); markers(idxX, :); markers(idxY, :); alp];


% - Fit geometric transformation to control point pairs
T = fitgeotrans(org, post,'projective'); % projective2d クラス（変換行列と次元）
T.T             % 生成された行列の確認


% - Apply geometric transformation to image

bwt = imwarp(bw, T); 
imshow(bwt); % デフォルトは線形補間


% - Remove noises
bwti = ~bwt;
%Iinv = imclose(Iinv, strel(ones(2)));
bwti = imopen(bwti, strel(ones(3)));
bwti = imclearborder(bwti);
bwt = ~bwti;
imshow(bwt)


% - Find where the finder patterns are for transposed image
[Centroid, bw, flag, width, ~] = detectFinder(bwt);


% - Extract XY coordinates of center of finder patterns
markers = Centroid(logical(flag),:);

% - Understand positional relation of each patterns
[idxC, idxX, idxY] = finderPos(markers);


% - Rotate image and coordinates to locate the corner of finder patterns to
% top left
[idxC, idxX, idxY, bwt] = rotPosition(idxC, idxX, idxY, bwt, markers);
imshow(bwt)


% - Calculate the size of cell for version2 QR Code in an image

% 1セルの高さを決定
% ファインダマーカー間の距離は18セルであるので、regionpropsで求めた座標の
% 距離を18で割り、1セル間の高さを求めます
%xy = round(markers(idxC,:));
xy = round(idxC);
%height = (markers(idxY, 2) - xy(2)) / 18;
height = (idxY(1, 2) - xy(2)) / 18;
% 同様にして1セルの幅を求めます
width = (idxX(1, 1) - xy(1)) / 18;

% QRコードに埋め込まれているタイミングパターンを利用した方がより正確ですが、
% 簡略化のために本例題では割愛します。


% - Create array that represents x and y coordinates of each cells

%QRコード : バージョン2は縦横25セル
xref = xy(1) - width*3;
yref = xy(2) - height*3;
xtable = (1:25) * width - width;
xtable = round(xtable + xref);
ytable = (1:25) * height - height;
ytable = round(ytable + yref);

% - Read format information to understand error correction level and mask
% pattern used in QR Code

% デコードを行うために必要な誤り訂正レベルや、施されたマスクパタンの情報が
% 記述されています。
bwo = ~bwt;
fmt1 = bwo(ytable(9), xtable([1:6,8]));
fmt2 = bwo(ytable([1:6, 8:9]), xtable(9));

%fmt2 = bwo(xtable([1:6,8:9]), ytable(9));
%fmt1 = bwo(xtable(9), ytable([1:6, 8]));

fdata = [fmt1, fliplr(fmt2')];


% - Decode format information using BCH decoder

% フォーマット情報にはマスクパタンが施されているため、xorをとります。
bitp = logical([1 0 1 0 1 0 0 0 0 0 1 0 0 1 0]);
fdata = xor(fdata, bitp);
% BCH符号により符号化されているので、復号します
finfo = bchdec(gf(fdata), 15, 5);
minfo = finfo(3:5); %マスクパタン情報

% - Create index table that represents how to read data bits and error
% correction bits

% 読み出しは右下角から順に、隣り合う列と交互に読み出していき、
% コードの上端もしくは下端に到達した場合は読み出す方向を転換します。

% 基準テーブル作成
xidx = 0:24;
xidx = repmat(xidx, [25, 1]);
yidx = 0:24;
yidx = repmat(yidx', [1, 25]);
idxtbl = cat(3, xidx, yidx);
% データ読み取り箇所インデックス算出
[xidx yidx reftable] = readtable(28+16);
% 読み出し方が特殊なので、データと訂正符号が格納されている領域を
% 可視化して確認します。
reftable = label2rgb(reftable);
figure, imshow(reftable)
title('QR code layout : Version2, Error correction level = L')


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

% フォーマット情報に記載されたマスクパタンの情報を基に、
% 該当箇所のデータを反転します。
mask = cat(3, mx, my);
invidx =  maskGen(mask, minfo);
data(invidx) = ~data(invidx);


% - Try to decode data without error correction

% version2, 誤り訂正レベルLの場合、文字数は最大28、
% 誤り訂正符号は16となります。

% データの先頭4bitはQRコードのモード指示子になっています。
% 以降のデコードに必要なため、まずモードを確認します。
% 1 : 数字モード
% 2 : 英数字モード 
% 3 : 8bitモード
% 4 : 漢字モード
data = num2str(data)';
mode = data(1:4);
mode = bin2dec(mode);

% ここでは、8bitモードのみ取り上げます。
if mode == 4 % 8bitモード
    %文字数指示子の復号
    length = data(5:12);
    length = bin2dec(length);
    D = data(13:end);
    text = cell(1, length);
    for i = 1:length
        tmp = D((i-1)*8+1:i*8);
        text{1, i} = char(bin2dec(tmp));
    end
end
% 復号した結果を表示
text = strjoin(text)

% 何らかの文字列が確認できましたが、QRコード上に置かれた紙によって
% 一部が隠れていたため、文字列が正しく読めていない可能性が大です。


% - Insert decoded data into image and visualize it
I1 = insertObjectAnnotation(img, 'rectangle', bbox.BoundingBox,...
    text, 'TextBoxOpacity',0.9,'FontSize',26,'LineWidth',3);
figure, imshow(I1)


% - Decode data bits with error correction

% QRコードではブロック誤りに強いRS(リードソロモン)が利用されています。
% まず、2進のビット列を10進に変換し、扱いやすくします。
data2 = zeros(1, 44);
for i = 1:44
    tmp = data((i-1)*8+1:i*8);
    data2(1, i) = bin2dec(tmp);
end
%データと訂正符号に分離します。
data = data2(1, 1:28);
tail = data2(1, 29:end);


% RS(リードソロモン)復号化器のコンフィグレーション
m = 8; % 1シンボルあたりのビット数
n = 2^m-1; % 符号語長
k = n-size(tail, 2); % メッセージ長

hDec = comm.RSDecoder(n,k);
hDec.GeneratorPolynomialSource = 'Property'; % 符号の生成多項式はユーザー指定
hDec.GeneratorPolynomial       = rsgenpoly(n,k,[],0); % 生成多項式

% 短縮符号であるため、前方をゼロパディング
shortmsg = [zeros(1,n-size(data, 2)-size(tail,2)) data tail]'; 

decoded = step(hDec, shortmsg);
decoded = decoded(size(decoded,1)-size(data,2)+1:size(decoded,1),1);

% 10進>2進変換
tmp = 0;
for i = 1:28
    tmp = [tmp, dec2bin(decoded(i),8)];
end
databox = tmp(1,2:end);

mode = databox(1:4);
mode = bin2dec(mode);
% ここでは、8bitモードのみ取り上げます。
if mode == 4 % 8bitモード
    %文字数指示子の復号
    length = databox(5:12);
    length = bin2dec(length);
    D = databox(13:end);
    text = cell(1, length);
    for i = 1:length
        tmp = D((i-1)*8+1:i*8);
        text{1, i} = char(bin2dec(tmp));
    end
end
% 復号した結果を表示
text = strjoin(text)


% - Insert decoded data into image and visualize it
I2 = insertObjectAnnotation(img, 'rectangle', bbox.BoundingBox,...
    text, 'TextBoxOpacity',0.9,'FontSize',26,'LineWidth',3,'Font','Calibri');
figure, imshow(I2)

%%
% Copyright 2019 The MathWorks, Inc.

%% detectFinder
function [Centroid, bw2, flag, wsum, bbox] = detectFinder(img)

% 颜色反转
bwo = ~img;

% 清理边界
bw = imclearborder(bwo);

% Bounding Box
bw2 = imopen(bw, ones(5));
bw2 = imclose(bw2, ones(25));
bbox = regionprops(bw2, 'BoundingBox');

% EulerNumber
bw = bwpropfilt(bw, 'EulerNumber', [0, 0]);
% 孔洞填充
bwa = imfill(bw, 'holes');

% 重心
stats = regionprops(bwa, 'Centroid');
Centroid = cat(1, stats.Centroid);

% 線幅が1:1:3:1:1の比となる領域(ファインダーの位置)を探索
sz = size(Centroid,1);
flag = zeros(sz(1),1);
bwl = bwlabel(bwa);% 将逻辑值矩阵中连通的"1"看成一块，并统一编号
bw2 = bw;
wsum = 0;
for i = 1:sz(1)
    roi = bwl == i;
    area = bwo;
    area(~roi) = 0;

    xy = round(Centroid(i,:));
    x = area(:, xy(1));
    
    pw = pulsewidth(double(x)); %パルス幅(明部)
    ps = pulsesep(double(x)); %パルス間隔(暗部)
    
    % 明部3,暗部2の組み合わせが検出されていない場合はBreak
    flagout = ~(size(ps, 1) == 2 && size(pw, 1) == 3);
    if flagout == 1
        bw2(bwl == i) = 0;
    else
        % 1セルの幅の平均値を求め、比率が1:1:3:1:1に近いかどうかを判定
        width = sum([pw;ps]) / 7;
        pat = [pw(1), ps(2), pw(2), ps(2), pw(3)] / width;
        if pdist2(pat, [1 1 3 1 1]) < 0.6
            flag(i) = 1;
        else
            bw2(bwl == i) = 0;
        end
        wsum = wsum + width;
    end
end
wsum = wsum/sz(1);
end

%% finderPos
function [idxC, idxX, idxY] = finderPos(markers)

% 2直線の成す角度が90度に近い点を探索
angle = zeros(3,1);
for i = 1:3
    tmp = markers;
    ref = markers(i,:);
    tmp(i,:) = [];
    % 線を表す式 : y = [x 1] * [a; b] のabを求める為に、polyfitを利用
    ab2 = polyfit([ref(1,2); tmp(1,2)], [ref(1,1); tmp(1,1)], 1);
    ab1 = polyfit([ref(1,2); tmp(2,2)], [ref(1,1); tmp(2,1)], 1);
    % ドット積を利用して角度を求める
    vect1 = [1 ab1(1)];
    vect2 = [1 ab2(1)];
    dp = dot(vect1, vect2);
    % ベクタ長の計算
    length1 = sqrt(sum(vect1.^2));
    length2 = sqrt(sum(vect2.^2));
    % 角度を算出(larger angle)
    angle(i) = 180-acos(dp/(length1*length2))*180/pi;
end
[~, idxC] = min(abs(abs(angle) - 90));

% X軸との角度を求める
tblidx = [1 2 3];
tmp = markers;
tmp(idxC,:) = [];
tblidx(idxC) = [];
ref = markers(idxC,:);
angle = zeros(2, 1);
for i = 1:2
    % 線を表す式 : y = [x 1] * [a; b] のabを求める為に、polyfitを利用
    ab2 = polyfit([ref(1,2); tmp(i,2)], [ref(1,1); tmp(i,1)], 1);
    %ab1 = polyfit([ref(1,2); ref(1,2)], [ref(1,1); ref(1,1) + 15], 1);
    ab1(1) = 10^8; %X軸に水平な直線なので、係数は10^8程度に設定しておく
    
    % ドット積を利用して角度を求める
    vect1 = [1 ab1(1)];
    vect2 = [1 ab2(1)];
    dp = dot(vect1, vect2);
    % ベクタ長の計算
    length1 = sqrt(sum(vect1.^2));
    length2 = sqrt(sum(vect2.^2));
    % 角度を算出(larger angle)
    angle(i) = acos(dp/(length1*length2))*180/pi;
end
idx = angle > 90;
angle(idx) = 180 - angle(idx);

[~, i] = min(abs(angle));
idxX = tblidx(i);
tblidx(i) = [];
idxY = tblidx;
end

%% createControlPts
function [Cref, Xref, Yref, Aref] = createControlPts(idxC, idxX, idxY, width, height, markers)

% ファインダの位置と対になる基準点を算出します。
ost = 7; %左上基準点からのマージン(セル数)
ratio = 2;
refw = (1+ost)*width;
refh = (1+ost)*height;
w = (18)*width;
wa = (15)*width;
h = (18)*height;
ha = (15)*height;

% 基準となるファインダが存在する象限を求める
cent(1, 1) = (markers(idxC, 1) + markers(idxX, 1) + markers(idxY, 1)) / 3;
cent(1, 2) = (markers(idxC, 2) + markers(idxX, 2) + markers(idxY, 2)) / 3;

%2 or 3
if markers(idxC, 1) < cent(1, 1)
    % 2 -- 基準
    if markers(idxC, 2) < cent(1, 2)
        mode = 2;
    else
        mode = 3;
    end
%1 or 4
else
    % 1
    if markers(idxC, 2) < cent(1, 2)
        mode = 1;
    % 4
    else
        mode = 4;
    end
end

% 基準位置に応じた座標指定
if mode == 3
    Cref = [refw, refh+h];
    Xref = [refw+w, refh+h]
    Yref = [refw, refh];
    Aref = [refw+wa, refh+(h-ha)];
elseif mode == 2
    Cref = [refw, refh];
    Xref = [refw+w, refh]
    Yref = [refw, refh+h];
    Aref = [refw+wa, refh+ha];
elseif mode == 1
    Cref = [refw+w, refh];
    Xref = [refw, refh]
    Yref = [refw+w, refh+h];
    Aref = [refw+(w-wa), refh+ha];
else
    Cref = [refw+w, refh+h];
    Xref = [refw, refh+h]
    Yref = [refw+w, refh];
    Aref = [refw+(w-wa), refh+(h-ha)];
end

end

%% rotPosition
function [idxC2, idxX2, idxY2, bwt]  = rotPosition(idxC, idxX, idxY, bwt, markers)

% 基準となるファインダが存在する象限を求める
cent(1, 1) = (markers(idxC, 1) + markers(idxX, 1) + markers(idxY, 1)) / 3;
cent(1, 2) = (markers(idxC, 2) + markers(idxX, 2) + markers(idxY, 2)) / 3;

%2 or 3
if markers(idxC, 1) < cent(1, 1)
    % 2 -- 基準
    if markers(idxC, 2) < cent(1, 2)
        mode = 2;
    else
        mode = 3;
    end
%1 or 4
else
    % 1
    if markers(idxC, 2) < cent(1, 2)
        mode = 1;
    % 4
    else
        mode = 4;
    end
end

sz = size(bwt);

% 画像&座標回転
if mode == 2
    %何もしない
    theta = 0;
    idxC2 = [markers(idxC, 1), markers(idxC, 2)];
    idxX2 = [markers(idxX, 1), markers(idxX, 2)];
    idxY2 = [markers(idxY, 1), markers(idxY, 2)];         
elseif mode == 1
    theta = 90;
    idxC2 = [markers(idxC, 2), sz(2)-markers(idxC, 1)];
    idxX2 = [markers(idxY, 2), sz(2)-markers(idxY, 1)];
    idxY2 = [markers(idxX, 2), sz(2)-markers(idxX, 1)];            
elseif mode == 3
    theta = -90;
    idxC2 = [sz(1) - markers(idxC, 2), markers(idxC, 1)];
    idxX2 = [sz(1) - markers(idxY, 2), markers(idxY, 1)];
    idxY2 = [sz(1) - markers(idxX, 2), markers(idxX, 1)]; 
else
    theta = 180;
    idxC2 = [sz(2) - markers(idxC, 1), sz(1) - markers(idxC, 2)];
    idxX2 = [sz(2) - markers(idxX, 1), sz(1) - markers(idxX, 2)];
    idxY2 = [sz(2) - markers(idxY, 1), sz(1) - markers(idxY, 2)];     
end

bwt = imrotate(bwt, theta);

end