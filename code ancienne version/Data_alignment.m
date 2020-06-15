clear all; clc

 cd '/home/vicente.enguix/Desktop'
 FolderName='AP05_1';

imgFilesList = dir([FolderName filesep 'img_0*.bin']);
aiFilesList = dir([FolderName filesep 'ai_*.bin']);

header=memmapfile([FolderName filesep imgFilesList(1).name],'Offset',0,'Format',{'int32',5,'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
nx=header.Data.header(2);
ny=header.Data.header(3);

frameFormat = {'uint64',3, 'framej'; 'uint16', [double(nx), double(ny)], 'imgj'};

im_res=[nx, ny];
SizeImage=nx*ny*2+8;

k=1;
M=[];
i=1;
n=1;
%% Read single slice red data
   img = memmapfile([FolderName filesep imgFilesList(i).name],'Offset',5*4,'Format',frameFormat);
    ii=1;
    a=img.Data(ii).imgj;
    a=imrotate(a,90);
    %figure(1),imshow(a,[])
  

%% Skeleton
    B = fibermetric(a,70,'ObjectPolarity','dark');
    BW2 = im2bw(B,0.0005);
    BW3 = bwskel(BW2);
    Area = regionprops(BW3,'Area');
    BW4 = bwareaopen(BW3,1000);
    
    %imshow(BW4)
   
%% Select bregma and lambda
display('select bregma and lambda')
 imshow(imfuse(a,BW4 ,'falsecolor','ColorChannels',[0 1 2]));
 [x,y]=ginput(2);close();
 


 %% Image centering visual check
 centro_x=(size(a,1))/2;
 centro_y=(size(a,2))/2;
 mover_x=[x(1)-centro_x];
 mover_y=[y(1)-centro_y];
 

 a2=imtranslate(a,[-mover_x, -mover_y]);%ahora bregma esta en el centro
 c_contiguo=abs(x(1)-x(2));
 c_opuesto=y(1)-y(2);%si lambda esta mas bajo, valor negativo
 tan_alpha=(c_opuesto/c_contiguo);
 angulo_rad=atan(tan_alpha);
 angulo_grad=(angulo_rad*360)/(6.28);
 img=imrotate(a2,angulo_grad,'bilinear','crop');

 imshow(img,[])
 
 line_x=[centro_x centro_x-300];
 line_y=[centro_y centro_y];
 
 hold on
 plot(centro_x, centro_y,'r*')
 plot(centro_x-300, centro_y,'b*')
 plot(line_x, line_y,'b-')
 legend('image center','point in medial line','medial line')

 %% Realign all the data test
 %i=1;
 for i = 1:size(imgFilesList,1)
 imagen = memmapfile([FolderName filesep imgFilesList(i).name],'Offset',5*4,'Format',frameFormat);
   
 for ii=1:size(imagen.Data,1)
     data=imagen.Data(ii).imgj;
     data=imrotate(data,90);
     data=imtranslate(data,[-mover_x, -mover_y]);
     data=imrotate(data,angulo_grad,'bilinear','crop');
     %imagen.Data(ii).imgj=data;
%      figure(1),imshow(imagen.Data(ii).imgj,[])
figure(1),imshow(data,[])
hold on
 plot(centro_x, centro_y,'r*')
 plot(centro_x-300, centro_y,'b*')
 plot(line_x, line_y,'b-')
 legend('image center','point in medial line','medial line')


    end
 end
