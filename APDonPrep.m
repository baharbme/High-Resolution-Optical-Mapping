% clear;clc;close all;
% load 'APDs.mat';
function APDonPrep (APD,BGFrame, flipoption,titlestr)
f2=figure;

APD = imresize(APD,3);
MinAPD = 2;
MaxAPD = max(max(APD));
% Levels = ceil(MaxAPD./2);
Levels=20;
delta_row=2;
delta_col=2;
color_const = 20;
StepAPD = (MaxAPD-MinAPD)/Levels;
map = colormap(jet(Levels+1));
MAP=map;
MAP(1,:)=[0,0,0];
map=map([Levels+1:-1:1],:);

%  colormap(map);

% steps = floor(color_const/Levels); % color increase steps
% steps = floor(color_const/Levels); % color increase steps

% BGFrame = imresize(RawData(:,:,50),3);
if flipoption == 1
    BGFrame = fliplr(BGFrame);
end
figure(f2);
% subplot(2,1,1);
% imshow(mat2gray(BGFrame));
h = mat2gray(BGFrame);
H(:,:,1)=h;
H(:,:,2)=h;
H(:,:,3)=h;
% scrsz = get(0,'ScreenSize');
% f1=figure('Position',[2*scrsz(3)/3 scrsz(4)/3 10 scrsz(4)/2]);imshow(H);
% f2=figure('Position',[2*scrsz(3)/3 2*scrsz(4)/3 10 scrsz(4)/2],'Visible','off');
imshow(H); hold on;

for k=2:Levels+1
    IMG = (APD > MinAPD+k*StepAPD) ;
    s=size(IMG);
    
    for row = 1:delta_row:s(1)
        for col=1:delta_col:s(2)
            if IMG(row,col)
                break;
            end
        end
        contour = bwtraceboundary(IMG, [row, col],'e',8,50000,'clockwise');
        %         color=steps*k+(color_const-steps*Levels);
        color=k;
        %
        if(~isempty(contour))
            figure(f2);
            hold on;
            q = fill(contour(:,2),contour(:,1),MAP(color,:),'LineWidth',1,'LineSmoothing','off','EdgeColor','none');
            axis image;
            axis off;
        end
    end
end
            axis off;
hold off
figure(f2); title(titlestr)
figure(f2);colorbar('YTick',[linspace(0,1,2)], 'YTickLabel',{' 0 ms',strcat(num2str(MaxAPD),' ms')});
colormap(map);
% title(titlestr)

% subplot(2,1,2);
% imshow(imresize(APD,3),MAP);colorbar
% colorbar('Ticks',[2,max(max(APD))],...
%          'TickLabels',{'0 ms',strcat(num2str(ceil(nanmax(nanmax(APD)))),'ms')})
% figure(f2);colorbar(map);

