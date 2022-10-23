
%% clear
clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%% necessary folders
%maximum and processed border images segmented (handcorrection.tif within
%numbered folder)
b_dir =[filedir, '/borders'];
%average projection of cytoskeleton component/intracellular marker
ske_dir = [filedir, '/cytoskeleton_average'];
%average projection of protein in the cell borders
bor_dir = [filedir, '/borders_average'];

%Folder to save information about cells
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];

Intensity_corr = [0,0];

%Reading 16-bit average intensity projection files (border protein)
cd(bor_dir);
files_tif = dir('*.tif');

%%introduce here the number of columns in the summary table
summary = zeros(numel(files_tif),11);
%%write a loop that gets data
for k=1:numel(files_tif)
    
    bd_dir = [b_dir,'/', num2str(k)];
    cd(bd_dir);
    I=imbinarize(imread('handCorrection.tif'));%reads the segmentation mask
    I = I(:,:,1); %you can now use imread
    %eliminate incomplete borders by making the image borders equal to
    %zero
    I(1,:) = 0;
    I(:,1) = 0;
    I(end,:) = 0;
    I(:,end) = 0;
    IC = imcomplement(I);
    IC2 = imclearborder(IC,4);
    
    %show two images together using "imshowpair(I,IC,'montage')"
    %why 4 and not 8, is it because the mask is only 1 pixel thickness?
    %generates complement image and remove incomplete cells
    IC3 = imdilate(IC2, strel('diamond', 1));
    IC4 = IC3 - IC2;
    %generates the complement
    [ic1,ic2] = find(IC4==1);
    B =[ic1,ic2];
    %finds the position of the pixels that are borders
    cd(ske_dir);
    %%use this line instead of the next two to use the original skeleton image
    %%without adjusting (necessary for tubulin antibody)
    % %     SignalSke = double(imread([num2str(k),'.tif']));
    
    SignalSke_original = imread([num2str(k),'.tif']);
    SignalSke=imadjust(SignalSke_original);
    cd(bor_dir);
    SignalBor = imread([num2str(k),'.tif']);
    SignalBor = imadjust(SignalBor);
    
    thresh = 2.5*graythresh(SignalBor);
    if thresh > 1
        thresh = 0.8;
    end
    Puncta = bwareaopen(imbinarize(SignalBor, thresh),10);
    
    Image1 = figure;
    imshowpair(SignalBor,Puncta,'montage');
    cd(im_dir);
    image_filename = [num2str(k),'_Puncta.tif'];
    print(Image1, '-dtiff', '-r150', image_filename);
    
    IC5 = imdilate(IC4, strel('disk', 3,0));
    %generate the cytoplasmic mask
    IC6 = imcomplement(IC3);
    IC7 = IC5 + IC6;
    Signal_Cyto=imcomplement(IC7);
    Signal_P = Puncta .* IC5;
    %Signal_P is the default puncta mask
    Signal_PD = imdilate(Signal_P, strel('disk', 8,0));
    %Signal_PD is the dilated puncta mask
    Signal_nonP = imsubtract(IC5,Signal_PD);
    %Signal_nonP is now borders - (dilated puncta)
    Signal_nonP(Signal_nonP==-1)=0;
    %this corrects -1 from previous operation
    Signal_nonPD = imdilate(Signal_nonP, strel('disk', 8,0));
    Signal_PD = imsubtract(Signal_PD,Signal_P);
    %Signal_PD is the mask for the halos around puncta
    Signal_nonPD = imsubtract(Signal_nonPD,Signal_nonP);
    %Signal_nonPD is the mask for the halos around non-puncta area
    Signal_PD = Signal_PD(:);
    Signal_nonPD = Signal_nonPD(:);
    Signal_line = SignalSke(:);
    Signal_line_ds = SignalBor(:);
    
    Signal_Ske_P = mean(Signal_line(Signal_P == 1));
    Signal_Ske_nonP = mean(Signal_line(Signal_nonP == 1));
    Signal_Ske_halo = mean(Signal_line(Signal_PD == 1));
    Signal_Ske_nonP_halo = mean(Signal_line(Signal_nonPD == 1));
    Signal_Ske_cyto = mean(Signal_line(Signal_Cyto == 1));
    
    Signal_Bor_P = mean(Signal_line_ds(Signal_P == 1));
    Signal_Bor_nonP = mean(Signal_line_ds(Signal_nonP == 1));
    Signal_Bor_halo = mean(Signal_line_ds(Signal_PD == 1));
    Signal_Bor_nonP_halo = mean(Signal_line_ds(Signal_nonPD == 1));
    Signal_Bor_cyto = mean(Signal_line_ds(Signal_Cyto == 1));
    
    summary(k,:) = [k, Signal_Ske_P, Signal_Ske_halo, Signal_Ske_nonP, Signal_Ske_nonP_halo, Signal_Ske_cyto,...
        Signal_Bor_P, Signal_Bor_halo, Signal_Bor_nonP, Signal_Bor_nonP_halo, Signal_Bor_cyto];
    
    CC = bwconncomp(Signal_P,8);
    S = regionprops(CC,'Area');
    IntensitySke = regionprops(CC, SignalSke, 'MeanIntensity');
    IntensityBor = regionprops(CC, SignalBor, 'MeanIntensity');
    
    Intensity_corr = [Intensity_corr; [IntensityBor.MeanIntensity]', [IntensitySke.MeanIntensity]'];
end


cd(filedir);
Intensity_corr(1,:) = [];
[r, p] = corrcoef(Intensity_corr(:,1), Intensity_corr(:,2));

Intensity_corr = [Intensity_corr; 0,0; r(1,2), p(1,2)];
correlation = array2table(Intensity_corr);
correlation.Properties.VariableNames = {'PCP', 'Signal'};
writetable(correlation,'punctaintensities.csv');

all = array2table(summary);
all.Properties.VariableNames = {'Wing', 'Tubulin_puncta', 'Tubulin_puncta_halo', 'Tubulin_nopuncta', 'Tubulin_nonpuncta_halo', 'Tubulin_cyto',...
    'Ds_puncta', 'Ds_puncta_halo', 'Ds_nopuncta', 'Ds_nopuncta_halo', 'Ds_cyto'};

summary(:,2:6) = summary(:,2:6)/mean(summary(:,2));
summary(:,7:11) = summary(:,7:11)/mean(summary(:,7));
all2 = array2table(summary);
all2.Properties.VariableNames = {'Wing', 'Tubulin_puncta', 'Tubulin_puncta_halo', 'Tubulin_nopuncta', 'Tubulin_nonpuncta_halo', 'Tubulin_cyto',...
    'Ds_puncta', 'Ds_puncta_halo', 'Ds_nopuncta', 'Ds_nopuncta_halo', 'Ds_cyto'};
writetable(all,'proximity.csv');
writetable(all2,'proximity_normalised.csv');
cd(currdir);



close all;
clear variables;
clc;