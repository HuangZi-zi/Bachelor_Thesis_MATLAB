function [Centroid, bw2, flag, wsum, bbox] = detectFinder(img)

%% 輝度反転
bwo = ~img;

%% 境界に接しているオブジェクトをクリア
bw = imclearborder(bwo);

%% QRコードの領域(Bounding Box)を算出
bw2 = imopen(bw, ones(5));
bw2 = imclose(bw2, ones(25));
bbox = regionprops(bw2, 'BoundingBox');

%% EulerNumberでフィルタリング
bw = bwpropfilt(bw, 'EulerNumber', [0, 0]);

%% 穴埋め
bwa = imfill(bw, 'holes');

%% 重心の位置を算出
stats = regionprops(bwa, 'Centroid');
Centroid = cat(1, stats.Centroid);

%% 線幅が1:1:3:1:1の比となる領域(ファインダーの位置)を探索
sz = size(Centroid,1);
flag = zeros(sz(1),1);
bwl = bwlabel(bwa);
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

%%
% Copyright 2019 The MathWorks, Inc.