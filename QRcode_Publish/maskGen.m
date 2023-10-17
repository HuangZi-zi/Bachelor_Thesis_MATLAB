function invidx = maskGen(mask, minfo)

info = minfo == [1 1 1];
sel = bin2dec(num2str(info));

switch sel
    case 0
        % "000" = (x + y) mod 2 = 0
        invidx = mod((mask(:,:,2) + mask(:,:,1)), 2) == 0;
    case 1
        % "001" = x mod 2 = 0
        invidx = mod(mask(:,:,2), 2) == 0;
    case 2
        % "010" = y mod 3 = 0
        invidx = mod(mask(:,:,1), 3) == 0;
    case 3
        % "011" = (x + y) mod 3 = 0
        invidx = mod((mask(:,:,2) + mask(:,:,1)), 3) == 0;
    case 4
        % "100" = ((y div 2)+(x div 3)) mod 2 = 0
        invidx = mod((fix(mask(:,:,1)/2) + fix(mask(:,:,2)/3)), 2) == 0;
    case 5
        % "101" = (xy) mod 2 + (xy) mod 3 = 0
        invidx = mod((mask(:,:,2) .* mask(:,:,1)), 2) + mod((mask(:,:,2) .* mask(:,:,1)), 3) == 0;
    case 6
        % "110" = ((xy) mod 2 + (xy) mod 3) mod 2 = 0
        invidx =  mod((mod((mask(:,:,2) .* mask(:,:,1)), 2) + mod((mask(:,:,2) .* mask(:,:,1)), 3)), 2) == 0;
    otherwise
        % "111" = ((xy) mod 3 + (x+y) mod 2) mod 2 = 0
        invidx = mod((mod((mask(:,:,2) .* mask(:,:,1)), 3) + mod((mask(:,:,2) + mask(:,:,1)), 2)), 2) == 0;
end

end

%%
% Copyright 2019 The MathWorks, Inc.