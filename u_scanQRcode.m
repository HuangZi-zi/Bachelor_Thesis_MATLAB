% function [] = u_scanQRcode(inputIMG)



% end

image=imread("Resource\QRcode.jpg");

decodedData = decodeQRCodeWithZBar(image);

% Display the decoded data
disp(decodedData);