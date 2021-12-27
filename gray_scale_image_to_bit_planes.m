%% APP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the original color image
secret = imread('cameraman.png');
imshow(secret)
%title('Original - RGB')

grayImage = rgbToGrayScaleImage(secret);

% Read in original RGB image.
rgbImage = imread('lena256X256.jpeg');
imshow(rgbImage)

% Create Stego Image by Applying Bit Plane Coding Alg.

bitPlaneToCode = 5;
firstMsbChoiceOfSecret = 8;
secondMsbChoiceOfSecret = 7;
thirdMsbChoiceOfSecret = 6;

recombinedRGBImage = bitPlaneCodingToCreateStego(rgbImage, grayImage, bitPlaneToCode, firstMsbChoiceOfSecret,secondMsbChoiceOfSecret,thirdMsbChoiceOfSecret);


% K-N Secret Sharing
k = 5;
n = 6;

allShares = knGenerateShares(recombinedRGBImage,n,k);

% LSB Stegonography

imageArray = getImagesFromFolder('envelopeİmages/','*.jpeg');

envelopedImages = cell(1,length(imageArray));

for imageIndex = 1:length(imageArray)
  fprintf(1, 'Enveloping Share %d\n', imageIndex);
  envelopedImages{imageIndex} = (lsbReplacement(allShares{imageIndex},imageArray{imageIndex}));
end

% LSB Decryption

getSharesInEnvelope = cell(1,length(imageArray));

for imageIndex = 1:length(envelopedImages)
  fprintf(1, 'Get Share %d from envelope\n', imageIndex);
  getSharesInEnvelope{imageIndex} = lsbDecryption(envelopedImages{imageIndex});
end

% K-N Secret Sharing Decryption

knDecrptedImage = knDecryption(getSharesInEnvelope,k);

% Bit Plane DeCoding to get Secret

secretFromCover = bitPlaneDecoding(knDecrptedImage,bitPlaneToCode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIT PLANE DECODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function secretImage = bitPlaneDecoding(stegoImage,bitPlaneToCode)

    
    % Extract color channels.
    redChannel = (stegoImage(:,:,1)); % Red channel
    greenChannel = (stegoImage(:,:,2)); % Green channel
    blueChannel = (stegoImage(:,:,3)); % Blue channel

    redChannelBitPlanes = bitPlaneGenerator(redChannel);
    blueChannelBitPlanes = bitPlaneGenerator(blueChannel);
    greenChannelBitPlanes = bitPlaneGenerator(greenChannel);
    % 
    
    temp = zeros(256,256,8);
    
    temp(:,:,8) = redChannelBitPlanes(:,:,bitPlaneToCode);
    temp(:,:,7) = greenChannelBitPlanes(:,:,bitPlaneToCode);
    temp(:,:,6) = blueChannelBitPlanes(:,:,bitPlaneToCode);
    
    secretImage = cast(fetchOriginal(temp),'uint8');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECREATE IMAGE FROM ITS SEPARETED CHANNELS INTO BIT PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function recombinedImage = bitPlaneCodingToCreateStego(cover, secret, bitPlaneToCoding, firstMsb, secondMsb, thirdMsb)
   
    % Extract color channels.
    redChannel = (cover(:,:,1)); % Red channel
    greenChannel = (cover(:,:,2)); % Green channel
    blueChannel = (cover(:,:,3)); % Blue channel

    redChannelBitPlanes = bitPlaneGenerator(redChannel);
    blueChannelBitPlanes = bitPlaneGenerator(blueChannel);
    greenChannelBitPlanes = bitPlaneGenerator(greenChannel);
    % 
    secretBitPlanes = bitPlaneGenerator(secret);

    newBlueChannel = bitPlaneCoding(blueChannelBitPlanes,bitPlaneToCoding,secretBitPlanes,thirdMsb);
    newRedChannel = bitPlaneCoding(redChannelBitPlanes,bitPlaneToCoding,secretBitPlanes,firstMsb);
    newGreenChannel = bitPlaneCoding(greenChannelBitPlanes,bitPlaneToCoding,secretBitPlanes,secondMsb);
    % 
    % 
    recombinedImage = cat(3, newRedChannel, newGreenChannel, newBlueChannel);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMBINED RGBs AFTER STENOGRAPHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newChannelBitPlanes=bitPlaneCoding(coverChannelBitPlanes,coverIndex,secretChannelBitPlanes,secretIndex)

    coverChannelBitPlanes(:,:,coverIndex) = secretChannelBitPlanes(:,:,secretIndex);
    newChannelBitPlanes = cast(fetchOriginal(coverChannelBitPlanes),'uint8');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSB REPLACEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shareEmbeddedEnvelope=lsbReplacement(share,envelope)

    envelopeOf24BitStream = to24bitStream(envelope);
    
    shareOf24BitStream = to24bitStream(share);
    
    flattenShare = reshape(shareOf24BitStream,1,[]);
    for k=1:size(envelopeOf24BitStream,1)
        tempEnvelope24bit = envelopeOf24BitStream(k,:);
        if (k-1)*2 == size(flattenShare,1)-2
            break
        end
        tempEnvelope24bit(24) = flattenShare((k-1)*2 + 1);
        tempEnvelope24bit(16) = flattenShare((k-1)*2 + 2);
        envelopeOf24BitStream(k,:) = tempEnvelope24bit;
    end
    
    shareEmbeddedEnvelope = uint8(bitStreamOf24toImage(envelopeOf24BitStream, 1024, 768));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSB RETRIEVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function retrieveLsb = lsbDecryption(envelope)

    envelopeOf24BitStream = to24bitStream(envelope);
    
    flattenShare = zeros(1,256*256*24);
    
    for i=1:size(envelopeOf24BitStream,1)
        if i == size(envelopeOf24BitStream,1)-2
            break
        end
        flattenShare((i-1)*2 + 1) = envelopeOf24BitStream(i,24);
        flattenShare((i-1)*2 + 2) = envelopeOf24BitStream(i,16);
    end
    
    retrieveLsb = bitStreamOf24toImage(reshape(flattenShare,256*256,24),256,256);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET IMAGES AS AN ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function images=getImagesFromFolder(folder,pattern)
    
    if ~isfolder(folder)
      errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
      uiwait(warndlg(errorMessage));
      return;
    end
    
    filePattern = fullfile(folder, pattern);
    jpegFiles = dir(filePattern);
    noOfJpegFiles = length(jpegFiles);
    
    images = cell(1,noOfJpegFiles);
    
    for imageIndex = 1:noOfJpegFiles
      baseFileName = jpegFiles(imageIndex).name;
      fullFileName = fullfile(folder, baseFileName);
      fprintf(1, 'Now reading %s\n', fullFileName);
      images{imageIndex} = imread(fullFileName);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIT PLANES GENERATOR FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bitPlanes=bitPlaneGenerator(grayImage)

% Store the original image size
[N, D]=size(grayImage);

% converts a nonnegative decimal integer d to a binary row vector.
% pixel intensity to binary
% for ex: (uint8)
% array = [1 2 3 4]
% ----------------------------------------
% de2bi(array)
% ----------------------------------------
% 0    0    0    0    0    0    0    1
% 0    0    0    0    0    0    1    0
% 0    0    0    0    0    0    1    1
% 0    0    0    0    0    1    0    0
% -----------------------------------------
grayImageBinary = de2bi(double(grayImage));

% Since we have used uint8 , it returns 8
numberOfPlane = size(grayImageBinary,2);

% 256 x 256 x 8 Matrix to store each bit plane matrix
bitPlanes = zeros(N, D, numberOfPlane);

% Create and calculate each bit plane
for i = 1:numberOfPlane
    bitPlane = grayImageBinary(:,i);
    bitPlanes(:,:,i) = reshape(bitPlane, N, D);
end

% Show the each bit plane
% figure
% for i = 1:numberOfPlane
%     plane = bitPlanes(:,:,i);
%     subplot(2,4,i);
%     if i == 1
%         imshow(plane);title('LSB Bit Plane');
%     elseif i == 8
%         imshow(plane);title('MSB Bit Plane');
%     else
%         imshow(plane);title(sprintf('%d nd Bit Plane',i));
%     end
% end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RGB TO GRAYSCALE IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function grayImage=rgbToGrayScaleImage(RGB)
    if length(size(RGB)) == 3
        grayImage = rgb2gray(RGB);
        imshow(grayImage)
        title('Gray Scale Image of Original Image')
    else
        grayImage = RGB;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FETCHING ORIGINAL MESSAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cc=fetchOriginal(allBitPlanes)
cc = allBitPlanes(:,:,8);
for a = 1:7
    cc = 2 * cc + allBitPlanes(:,:,8-a);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE RANDOM VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rand=generateRandomVar(n)
    recons = 2;
    rand = [0 0];
    count = 1;
    while count < recons+1
        random = randi([1 n],1,1);
        while (~ismember(random, rand))
            rand(count) = random;
            count = count + 1;
        end
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 24 BIT STREAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bitStreamOf24 = to24bitStream(image)

    redChannelofStego = (image(:,:,1)); % Red channel
    greenChannelofStego = (image(:,:,2)); % Green channel
    blueChannelofStego = (image(:,:,3)); % Blue channel
    
    bitSteamOfRedChannelOfStego = double(dec2bin(double(redChannelofStego)) == '1');
    bitSteamOfGreenChannelOfStego = double(dec2bin(double(greenChannelofStego)) == '1');
    bitSteamOfBlueChannelOfStego = double(dec2bin(double(blueChannelofStego)) == '1');
    
    bitStreamOf24 = [bitSteamOfRedChannelOfStego bitSteamOfGreenChannelOfStego bitSteamOfBlueChannelOfStego];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 24 BIT STREAM to IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constructedImage = bitStreamOf24toImage(bitStreamOf24,N,D)

        shareR=bitStreamOf24(:,1:8);
        shareG=bitStreamOf24(:,9:16);
        shareB=bitStreamOf24(:,17:24);

        shareRConstracted = bin2dec(num2str(shareR));
        shareGConstracted = bin2dec(num2str(shareG));
        shareBConstracted = bin2dec(num2str(shareB));
        
        
        shareRReshaped = reshape(shareRConstracted, N, D);
        shareGReshaped = reshape(shareGConstracted, N, D);
        shareBReshaped = reshape(shareBConstracted, N, D);

        constructedImage = (cat(3, shareRReshaped, shareGReshaped, shareBReshaped));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K-N SECRET SHARING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nShares=knGenerateShares(image,n,k)

    
    [N, D, F]=size(image);

    stegoImageDividedInto24 = to24bitStream(image);
    
    %K-N Secret Sharing parameters
    
    nShares = cell(1,n);
    for cell_idx = 1 : n
        nShares{cell_idx} = image;
    end

    fullFilledZeroShares = zeros(size(stegoImageDividedInto24,1), size(stegoImageDividedInto24,2), n);

    zeroConstuctedShares = cell(1,n);
    for cell_idx = 1 : n
        zeroConstuctedShares{cell_idx} = fullFilledZeroShares(:,:,cell_idx);
    end    
    rowSize=size(stegoImageDividedInto24,1);
    columnSize=size(stegoImageDividedInto24,2);
    
    for row=1:rowSize
        for column=1:columnSize
            if stegoImageDividedInto24(row,column) == 1
                rand = generateRandomVar(n);
                for shareNo=1:(n-k+1)
                    zeroConstuctedShares{rand(shareNo)}(row,column) = 1;
                end
            end
        end
    end
    
    for shares=1:n
        
        shareRGB=zeroConstuctedShares{shares};
        
        nShares{shares} = bitStreamOf24toImage(shareRGB, N, D);
        fprintf(1, 'K-N For Share %d\n', shares);
    end
    
    figure
    for i = 1:n
        shareS = nShares{i};
        subplot(3,2,i);
        imshow(shareS);title(sprintf('%d nd Share',i));
    end
    
    disp("End of KN");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K-N SECRET SHARING DECRYPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function decryptedImage=knDecryption(shares,k)
    
    [N, D, F] = size(shares{1});

    temp = zeros(256*256,24);

    for cell_idx = 1:k

        share = shares{cell_idx};

        stegoImageDividedInto24 = to24bitStream(share);

        temp = temp | stegoImageDividedInto24;
    end

    decryptedImage = uint8(bitStreamOf24toImage(temp, N, D));

end