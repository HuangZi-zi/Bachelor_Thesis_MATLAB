function outputIMG = u_basic_process(inputIMG)


% 转化为灰度图处理
% IMGgray=rgb2gray(inputIMG);



kernel=strel('square',3);

%模糊
IMGblur = imgaussfilt(inputIMG,4);

% %膨胀
% IMGdil=imdilate(IMGblur,kernel);
% 
% %腐蚀
% IMGero=imerode(IMGdil,kernel);
% 
% outputIMG=gray2rgb(IMGero);

% 自适应二值化
outputIMG=imbinarize(IMGblur,'adaptive');


% % 将图像从空域转换为频域
% frequencyDomainImage = fft2(double(IMGgray));
% 
% % Define the size of the filter (should match the size of the frequency domain image)
% filterSize = size(frequencyDomainImage);
% 
% % 带通滤波
% % Define the center frequency of the passband (adjust as needed)
% centerFrequency = [filterSize(1)/2, filterSize(2)/2]; % Adjust the coordinates as needed
% 
% % Define the width of the passband (adjust as needed)
% bandwidth = 20; % Adjust the bandwidth as needed
% 
% % Create frequency coordinates
% [X, Y] = meshgrid(-centerFrequency(2):filterSize(2)-centerFrequency(2)-1, ...
%                    -centerFrequency(1):filterSize(1)-centerFrequency(1)-1);
% 
% % Create a bandpass filter in the frequency domain
% distances = sqrt(X.^2 + Y.^2);
% bandpassFilter = (distances >= (centerFrequency - bandwidth/2)) & ...
%                  (distances <= (centerFrequency + bandwidth/2));
% 
% 
% % 在频域中应用滤波器
% filteredFrequencyDomainImage = frequencyDomainImage .* bandpassFilter;
% 
% % 将图像从频域转换回空域
% filteredImage = ifft2(filteredFrequencyDomainImage);

% outputIMG=uint8(abs(filteredImage));

figure;
subplot(1, 2, 1);
imshow(outputIMG);
title("图像基本处理")

subplot(1, 2, 2);
imshow(inputIMG);
title('Origin');
end