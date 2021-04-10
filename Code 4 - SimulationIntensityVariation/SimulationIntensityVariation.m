clear all;
% Select simulation mode: Gaussian beam or Bessel beam; 
% Add phase error or not; focal AO or pupil AO
flag_Bessel=1;
flag_Gaussian=0;
flag_addPhaseError=1;
flag_focalAO=1;
flag_pupilAO=0;

Abb_Z=200*1e3;
SlitX1=5107;
SlitX2=5580;
SlitX10=5207;
SlitX20=5480;
currentFolder=pwd;
Resultpath=[currentFolder '\UniformInt2_S'];
load([currentFolder '\field_int.mat']);


%% field at back pupil plane
phase1=imread([currentFolder '\AstiAO_Cal.tif']);
phase1=pi/255*double(phase1)-pi;

phase=ones(size(phase1));
phase=2*pi/255*double(phase)-pi;
[pixel_num,~]=size(phase);
%% SLM parameters2
SLM.pitch=15; %um
SLM.pixelNumber=512;
SLM.size=SLM.pitch*SLM.pixelNumber; %um
upper_bound=255;  
lower_bound=0;
delta_degree=2*pi/(upper_bound-lower_bound);
laser_sig=3.262; 
refind=1.33;
lambda=940*1e-3; 
obj.mag=25;
obj.NA=1.05;
obj.f=180/obj.mag; 
obj.D=2*obj.f*obj.NA; 
mag=2; 
k=2*pi/lambda; 

%coordinate
field_x=(-pixel_num/2:pixel_num/2-1)*SLM.pitch*2;
field_y=(-pixel_num/2:pixel_num/2-1)*SLM.pitch*2; 
[field_xx,field_yy]=meshgrid(field_x,field_y);
field_rr=sqrt(field_xx.^2+field_yy.^2);  


%% Uniform intensity
intensity=zeros(pixel_num);
if(flag_Bessel)
intensity(field_rr>SlitX10&field_rr<SlitX20)=1;
else
intensity(:)=1;
end
field=intensity.*exp(1.j*phase);

%% Propogate field and add phase error other than the pupil
zeroPadding=4;
intensity_new=zeros(4*size(intensity));
phase_new=zeros(4*size(phase));
phase_error=zeros(4*size(phase));

intensity_new((end-length(intensity))/2+1:(end+length(intensity))/2,(end-length(intensity))/2+1:(end+length(intensity))/2)=intensity;
phase_new((end-length(phase))/2+1:(end+length(phase))/2,(end-length(phase))/2+1:(end+length(phase))/2)=phase;
phase_error((end-length(phase1))/2+1:(end+length(phase1))/2,(end-length(phase1))/2+1:(end+length(phase1))/2)=phase1;

N_upSample=4096*2;
phase_error=imresize(phase_error,[N_upSample N_upSample]);
field_upSample_Intensity=imresize(intensity_new,[N_upSample N_upSample]);
field_upSample_Phase=imresize(phase_new,[N_upSample N_upSample]);
field_upSample=field_upSample_Intensity.*exp(1.j*field_upSample_Phase);
field_temp=propIR(field_upSample,field_rr(end/2,end)*2,lambda,Abb_Z);
%%Add phase error
if(flag_addPhaseError)
phase_Z=angle(field_temp)+phase_error;
else
phase_Z=angle(field_temp);
end
intensity_Z=abs(field_temp);
field_temp_Z=intensity_Z.*exp(1.j*phase_Z);

figure(1);
subplot(1,3,1);
imagesc(intensity_Z(length(intensity_Z)*3/8+1:length(intensity_Z)*5/8,length(intensity_Z)*3/8+1:length(intensity_Z)*5/8));
subplot(1,3,2);
imagesc(angle(field_temp(length(intensity_Z)*3/8+1:length(intensity_Z)*5/8,length(intensity_Z)*3/8+1:length(intensity_Z)*5/8)));
subplot(1,3,3);
imagesc(phase_Z(length(intensity_Z)*3/8+1:length(intensity_Z)*5/8,length(intensity_Z)*3/8+1:length(intensity_Z)*5/8));
field_temp2=propIR(field_temp_Z,field_rr(end/2,end)*2,lambda,-Abb_Z);
field_back=field_temp2(length(field_temp2)*3/8+1:length(field_temp2)*5/8,length(field_temp2)*3/8+1:length(field_temp2)*5/8);
figure(2);
subplot(2,2,1);
imagesc(field_x,field_y,abs(field))
subplot(2,2,2);
imagesc(field_x,field_y,angle(field))

subplot(2,2,3);
imagesc(field_x,field_y,abs(field_back))
subplot(2,2,4);
x=-1023:1:1024;
y=-1023:1:1024;
[XX,YY]=meshgrid(x,y);
rr=sqrt(XX.^2+YY.^2);
rr_mask=ones(size(rr));
rr_mask(rr>760-27)=0;%Bessel narrow slit
rr_mask(rr<670+27)=0;
phase_back=angle(field_back);
phase_back(~rr_mask)=0;
imagesc(field_x,field_y,phase_back)

%% Construc the final field
int_tem=imresize(abs(field_back),[512,512]);
int_tem1=zeros([512,512]);
int_tem2=zeros([512,512]);
int_tem3=zeros([512,512]);

int_tem1(field_rr>SlitX10&field_rr<SlitX20)=int_tem(field_rr>SlitX10&field_rr<SlitX20);
int_tem2(field_rr>SlitX10&field_rr<SlitX20)=AB(field_rr>SlitX10&field_rr<SlitX20);
AB1=AB*sqrt(sum(sum(int_tem1^2))/sum(sum(int_tem2^2)));
int_tem3(field_rr>SlitX10&field_rr<SlitX20)=AB1(field_rr>SlitX10&field_rr<SlitX20);
%%
if(flag_Gaussian)
 field=int_tem.*exp(1.j*1);   
elseif(flag_focalAO)
 field=int_tem3.*exp(1.j*1);   
end

field_x2=field_x;
field_y2=field_y;
field2=field;
if(flag_Bessel)
x=-2:0.05:2; %Bessel
y=-2:0.05:2;
z=-50:5:50;
else
x=-1:0.05:1; %Gaussian
y=-1:0.05:1;
z=-2:0.2:2;
end

for a=1:length(z)
    PSF=Calc_Annular_Field_Integrals_V2(x, y, z(a),field2, field_x2*1e-3, field_y2*1e-3, lambda,refind,obj.f);
    PSF=uint16(squeeze(PSF)*1e-9);
    imwrite(PSF,[Resultpath,num2str(z(a)),'.tif']);
end

img=imread([Resultpath,num2str(z(1)),'.tif']);
imwrite(img, [Resultpath,'Stack_0.05xy_1z.tif'], 'tif', 'WriteMode', 'overwrite','Compression', 'none');

for a=2:length(z)
    img=imread([Resultpath,num2str(z(a)),'.tif']);
    imwrite(img, [Resultpath,'Stack_0.05xy_1z.tif'], 'tif', 'WriteMode', 'append','Compression', 'none');
end

















