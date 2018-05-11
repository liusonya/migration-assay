clear
close all
img_name='20180430_G6_stopper-0025_GFP.TIF';
rng(0)

wrong_pixel_penalty=1;
%intensity_threshold=112.5; %255/2
%intensity_threshold=85; %255/3
%intensity_threshold=63.75; %255/4
intensity_threshold=51; %255/5
%intensity_threshold=25.5; %255/10

% false, no calibrate. true, calibrate
train_mask = false;

%Define "training" file here:
if train_mask


    kernelSize=3;
    img=imread(img_name);
    I2=imgradient(img);
    I1=stdfilt(I2);
    I1=medfilt2(I1,[kernelSize,kernelSize]);
    I=mat2gray(I1);

    se = strel('disk', 40);
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    mask = im2bw(Iobrcbr, graythresh(Iobrcbr));
    Lrgb = label2rgb(mask, 'jet', 'w', 'shuffle');

    % figure, imshow(img), imagesc, hold on
    % himage = imshow(Lrgb);  
    % set(himage, 'AlphaData', 0.3);


    % Create a logical image of a circle with specified
    % diameter, center, and image size.
    % First create the image.
    [imageSizeX, imageSizeY] = size(img);

    centerX = imageSizeX/2;
    centerY = imageSizeY/2;
    radius = imageSizeX*3/8;



    lb=[centerX-200, centerY-200, round(imageSizeX/8)]; 
    ub=[centerX+200, centerY+200, round(imageSizeX/2-10)];

    [columnsInImage, rowsInImage ] = meshgrid(1:imageSizeY, 1:imageSizeX);
    f = @(x)obj_func(x,~mask, wrong_pixel_penalty, columnsInImage, rowsInImage);

    % Compute the optimal radius
    %x= ga(f, 3, [],[], [],[], lb,ub, [], [1,2,3])
    init_guess = round((lb + ub)/2);
    tic
    x = CGA(f,init_guess, lb,ub)
    toc

    % circlePixels is the "mask" i.e. the circumference that we are interested on.
    circlePixels = (rowsInImage - x(1)).^2 + (columnsInImage - x(2)).^2 <= x(3).^2;

    Lrgb = label2rgb(circlePixels, 'jet', 'w', 'shuffle');
    figure, imshow(img), imagesc, hold on
    himage = imshow(Lrgb);  
    set(himage, 'AlphaData', 0.3);
    print('-dpng','maskimage','-r300')
    
    disp('Press the "Enter" key to continue...')
    pause

%     fid=fopen('mask_paramters.txt','w');
%     fprintf(fid,'Parameters used are:\nCenterX=%d\nCenterY=%d\nRadius=%d\n\n',x(1),x(2),x(3));
%     fclose(fid);
    
    save('last_simulation.mat')
end
% end of calibration 

load('last_simulation.mat')

% with the calibration done, to have results with different thresholds
% remove (%) coment and change the threshold number
% intensity_threshold=30;

intensity_threshold= 33*1.25; % max median pixel value of controls + 25%


% Find all the files that end with ".tif"
file_list= dir('*.TIF');
intensities=zeros(1,5);
counter=1;
file_names={};
for i=1: length(file_list)
    name= file_list(i).name;
    if ~isempty(strfind(name, 'nuclei'))
        continue
    end
    file_names{counter}=name;
    img= imread(name);
    masked_img= img(circlePixels);
    % count how many pixels with intensity higher than the threshold are inside the optimal circle
    intensities(counter,1)= sum(masked_img>=intensity_threshold);
    % count how many pixels are inside the optimal circle
    %intensities(counter,2)= sum(circlePixels(:));
    intensities(counter,2)= length(masked_img);
    % Percentage of pixels inside the optimal circle
    intensities(counter,3)= 100*intensities(counter,1)/intensities(counter,2);
    % Mean intensity inside the optimal circle
    intensities(counter,4)= mean(masked_img);
    % Mode of intensities inside the optimal circle
    intensities(counter,5)= median(masked_img);
     
    name
    disp('Percent')
    100*intensities(counter,1)/intensities(counter,2)
    disp('mean')
    mean(masked_img)
    disp('median')
    median(masked_img)
%     pause
    
    
    %Show the images:
    %figure %create a new figure 
%     imshow(img), imagesc, hold on
%     himage = imshow(Lrgb);
%     set(himage, 'AlphaData', 0.3);
%     title(name);
%     hold off
    
    counter=counter+1;
end

% Create the output file
% Column 1: Name of the file
% Column 2: Pixels above the threshold
% Column 3: Total pixles counts (threshold = 0)
% Column 4: Percentage of pixels above the threshold
fid1=fopen('output_c1.txt','w');
fid3=fopen('output_c3.txt','w');
for i=1:length(file_names)
    if ~isempty(strfind(file_names{i}, 'c1'))
        fprintf(fid1,'%s\t%d\t%d\t%2.4g\n', file_names{i}, intensities(i,1), intensities(i,2), intensities(i,3));
    else
        fprintf(fid3,'%s\t%d\t%d\t%2.4g\n', file_names{i}, intensities(i,1), intensities(i,2), intensities(i,3));
    end
end
fclose(fid1);
fclose(fid3);

% fid=fopen('temp_output.txt','w');
% 
% fprintf(fid,'Parameters used are:\nCenterX=%d\nCenterY=%d\nRadius=%d\n\n',x(1),x(2),x(3));
% for i=1:length(file_names)
%     name = file_names{i};
%     img= imread(name);
%     masked_img= img(circlePixels);
%     fprintf(fid,'%s\t%d\t%d\t%2.2g\t', file_names{i}, intensities(i,1), intensities(i,2), intensities(i,3));
% %     for j = 20:2.5:35
% %         fprintf(fid,'theshold=%d %2.2g\t', j,100*sum(masked_img>=j)/intensities(i,2));
% %     end
%     
%     for j = [30,50,90]
%         fprintf(fid,'theshold=%d %2.2g\t', j,100*sum(masked_img>=j)/intensities(i,2));
%     end
%     
%     
%     
% %     for j = 0:25:255
% %         fprintf(fid,'theshold=%d %2.2g\t', j,100*sum(masked_img<=j)/intensities(i,2));
% %     end
%     
% %     for j = 25:25:100
% %         for k = 255:25:150
% %             fprintf(fid,'theshold=[%d,%] %2.2g\t', j,100*sum(masked_img<=j)/intensities(i,2));
% %         end
% %     end
%     
%     fprintf(fid,'\n');
% end
% 
% fclose(fid);