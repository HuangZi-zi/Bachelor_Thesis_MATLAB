function [idxC2, idxX2, idxY2, bwt]  = rotPosition(idxC, idxX, idxY, bwt, markers)

%% 基準となるファインダが存在する象限を求める
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

%% 画像&座標回転
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