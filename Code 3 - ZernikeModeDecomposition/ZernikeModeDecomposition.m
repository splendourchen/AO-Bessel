%% Zernike mode decomposition Script
% Get AO pattern
% addpath(genpath('H:\Projects\Codes\C2_AO_Bessel\Functions\'));
currentFolder = pwd;
[file,FilePath] = uigetfile([currentFolder, '*.tif']);
fileName=[FilePath file];
AOPhase=imread(fileName,'tif');% 
AOPhase_length=AOPhase/255*0.94; % in unit um
AOPhase_length=imresize(AOPhase_length,[512 512],'cubic');

%% Zernike decomposition
ZernikeModeN=55;
wavefront=AOPhase_length;
centerPos=[256 256];
pupilSize=256*2;
wavefront_fit=AOPhase_length(centerPos(1)-pupilSize/2+1:centerPos(1)+pupilSize/2,centerPos(2)-pupilSize/2+1:centerPos(2)+pupilSize/2);
zernikeCoeff = ZernikeDecomposition(wavefront, centerPos, pupilSize,ZernikeModeN);


%% Reconstruct wavefront

reconstructWF = zeros(pupilSize);
for i = 1:ZernikeModeN
    [Mode, ~] = zernike_fun(i,pupilSize);                                % Mode: 51x51 wavefront for each Zernike mode
    reconstructWF = reconstructWF + Mode * zernikeCoeff(i);   % Zernike_coef*mode=phase --> Zer_coef=phase/mode 
end

[DefocusMode, ~]=zernike_fun(5,pupilSize); 
DefocusWF=DefocusMode* zernikeCoeff(5);

reconstructWFNoDefocus = zeros(pupilSize);
for i = 1:ZernikeModeN
    if (i~=5&&i~=1&&i~=2&&i~=3)
    [Mode, ~] = zernike_fun(i,pupilSize);                                % Mode: 51x51 wavefront for each Zernike mode
    reconstructWFNoDefocus = reconstructWFNoDefocus + Mode * zernikeCoeff(i);   % Zernike_coef*mode=phase --> Zer_coef=phase/mode 
    end
end

%% Plot
Wavelength=0.94;%um
cRange=[0 0.5];
axis1=figure(1);
bar(zernikeCoeff);
xticks([1:1:10,13]);
xticklabels({'piston','V-tilt','H-tilt','O-Asti','Defocus','V-Asti','V-Trefoil','V-Coma','H-Coma','O-Trefoil','SA'})
xlabel('ANSI standard Zernike modes');ylabel('coefficients');
xlim([4 15]);
ylim([-0.1 0.1]);
set(axis1,'color','w');

axis2=figure(2);
subplot(2,2,1)
imagesc(wavefront_fit./Wavelength);
caxis(cRange);
xlabel('x (pixels)');ylabel('y (pixels)');
h = colorbar;
set(get(h,'title'),'string','wave');

subplot(2,2,2)
imagesc(reconstructWF./Wavelength);
caxis(cRange);
xlabel('x (pixels)');ylabel('y (pixels)');
h = colorbar;
set(get(h,'title'),'string','wave');

subplot(2,2,3)
imagesc(DefocusWF./Wavelength);
caxis(cRange);
xlabel('x (pixels)');ylabel('y (pixels)');
h = colorbar;
set(get(h,'title'),'string','wave');

subplot(2,2,4)
imagesc(reconstructWFNoDefocus./Wavelength);
caxis(cRange);
xlabel('x (pixels)');ylabel('y (pixels)');
h = colorbar;
set(get(h,'title'),'string','wave');

set(axis2,'color','w');
%% Save reconstruct pattern
% wavefrontSave=reconstructWF./0.94.*255;
% imwrite(wavefrontSave, [fileName(1:end-4), '_reconstructWF.tif']);


wavefrontSave=single(reconstructWFNoDefocus./0.94.*255);
t=creatTiffFile([fileName(1:end-4), '_NoDefocus.tif'],wavefrontSave);
t.write(wavefrontSave);
t.close();

wavefrontSave=single(reconstructWF./0.94.*255);
t=creatTiffFile([fileName(1:end-4), '_ReConsWF.tif'],wavefrontSave);
t.write(wavefrontSave);
t.close();

saveas(axis1,[fileName(1:end-4), '_ZernikeCoeff'],'fig');
saveas(axis1,[fileName(1:end-4), '_ZernikeCoeff'],'png');

saveas(axis2,[fileName(1:end-4), '_ZernikeImg'],'fig');
saveas(axis2,[fileName(1:end-4), '_ZernikeImg'],'png');

function t=creatTiffFile(fileName,data)
t = Tiff(fileName,'w');
% Setup tags
% Lots of info here:
% http://www.mathworks.com/help/matlab/ref/tiffclass.html
tagstruct.ImageLength     = size(data,1);
tagstruct.ImageWidth      = size(data,2);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)
end