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

%%
% Copyright 2019 The MathWorks, Inc.