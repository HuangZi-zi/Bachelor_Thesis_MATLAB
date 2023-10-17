function [x y table2] = readtable(num)

% データ・誤り訂正符号読み出し順インデックステーブル
table = zeros(25, 25);
table(1:9, 1:9) = 1; % 左上ファインダ
table(18:25, 1:9) = 2; % 左下ファインダ
table(1:9, 18:25) = 3; % 右上ファインダ
table(17:21, 17:21) = 4; % アライメントパターン
table(7, 10:17) = 5; % タイミングパターン1
table(10:17, 7) = 6; % タイミングパターン2
table2 = table;

ptx = 25;
pty = 25;
flag = true;
swap = false;
done = 0;

x = zeros(1, num*8);
y = zeros(1, num*8);
while (1)
    % 読み出し可能領域(ファインダパターン等以外)の場合
    if table(pty, ptx - swap) == 0
        done = done + 1;        
        x(1, done) = ptx - swap;
        y(1, done) = pty;
        table(pty, ptx - swap) = 6 + done + 1;
        table2(pty, ptx - swap) = fix((done-1)/8)*2 + 20;
    end
    % 所定の回数以上読み出した場合
    if done >= 8*num
        break;
    end
    % 右側 > 左側 > 右側のように左右交互に読み出すため、
    % swapフラグによりどちらを読むかを管理
    if swap == 1
        % 上方向に読み出していく場合
        if flag == 1
            pty = pty - 1;
        % 下方向に読み出していく場合    
        else
            pty = pty + 1;
        end
    end
    
    % 上方向に読み出し時、上限に達した場合
    if pty < 1
        pty = 1;
        if ptx == 9
            %例外処理
            ptx = 9 - 3;
        else
            ptx = ptx - 2;
        end
        flag = ~flag;
    end
    % 下方向に読み出し時、下限に達した場合
    if pty > 25
        pty = 25;
        ptx = ptx - 2;
        flag = ~flag;
    end    
    swap = ~swap;
end

end

%%
% Copyright 2019 The MathWorks, Inc.