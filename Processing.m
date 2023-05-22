% % --- Executes on button press in pushbutton38.
% function pushbutton38_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton38 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % RoiOption = get(handles.popupmenu2,'Value'); %% 1 for freehand - 2 for impoly - 3 for imellipse - 4 for imrect
% clc
% FullPath = get(handles.edit29,'String');
% load(strcat(FullPath,'.mat'));
% info = imfinfo(strcat(FullPath,'.tif'));
% frames = numel(info);
% num_images = numel(info);
% bit_depth = info.BitDepth;
% total_frames = get(handles.edit46,'String');
% total_frames = str2double(total_frames);
% total_time = get(handles.edit47,'String');
% total_time = str2double(total_time);
% ScalingOpt = get(handles.checkbox30,'Value');

% Scaling_progess = waitbar(0,'Please wait....');
% for i=1:size(Data,1)
%     waitbar(i/size(Data,1));
%     for j=1:size(Data,2)
%         ScaledData(i,j,:) = (Data(i,j,:)-nanmin(squeeze(Data(i,j,:))))/(nanmax(squeeze(Data(i,j,:)))-nanmin(squeeze(Data(i,j,:))));
%     end
% end
% close (Scaling_progess);
% Data = mat2gray(ScaledData);
Data = mat2gray(RawData);
% if Invchb == 1
    Dmax = 1;
    Dmin = 0;
    Data_Inv = Data-Dmax;
    Data_Inv = abs(Data_Inv);
    Data_Inv = Data_Inv+Dmin;
    %     Sig = squeeze(nanmean(nanmean(Data)));
    %     Siginv = squeeze(nanmean(nanmean(Data_Inv)));
    %     figure(3);cla;plot(Sig); hold on;plot(Siginv); legend('before inverting', 'after inverting')
    Data = Data_Inv;
    Data = mat2gray(Data);%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
%     end
    
% end

%     Dmax = 1;
%     Dmin = 0;
%     Data_Inv = Data-Dmax;
%     Data_Inv = abs(Data_Inv);
%     Data_Inv = Data_Inv+Dmin;
%     %     Sig = squeeze(nanmean(nanmean(Data)));
%     %     Siginv = squeeze(nanmean(nanmean(Data_Inv)));
%     %     figure(3);cla;plot(Sig); hold on;plot(Siginv); legend('before inverting', 'after inverting')
%     Data = Data_Inv;
%     Data = mat2gray(Data);

% LPFchb = get(handles.checkbox14,'Value');
% MAchb = get(handles.checkbox15,'Value');
% SMchb = get(handles.checkbox17,'Value');
% Invchb = get(handles.checkbox18,'Value');
% PCAchb = get(handles.checkbox19,'Value');
% preproschb = get(handles.checkbox20,'Value');
% bcorrchb = get(handles.checkbox21,'Value');

% if LPFchb == 1
    frame_rate=1000*total_frames/total_time;
%     cutoff_frequency = str2num(get(handles.edit79,'String'));
%     filter_order = str2num(get(handles.edit80,'String'));
    normalized_cutoff_frequency = cutoff_frequency/(frame_rate/2);
    [b,a] = butter(filter_order,normalized_cutoff_frequency,'low');
    FilteredData = filter(b,a,Data,[],3);
    FilteredData(:,:,1:50)=repmat(FilteredData(:,:,51),[1,1,50]);
%     windowWidth = 11; % Whatever you want.
% kernel = ones(windowWidth,1) / windowWidth;
% FilteredData = filter(kernel, 1, FilteredData); %moving average
%     hdl = waitbar(0,['IMAGE STACK TEMPORAL PROCESSING']);
% median = 40
%     for i=1:size(Data,1)
%         parfor j=1:size(Data,2)
%             if ~isnan(Mask(i,j))
%                     PIXEL=double(squeeze(Data(i,j,:)));
% 
%             %LOW PASS
%             [A,B]=butter(6,normalized_cutoff_frequency,'low');% lp filter
% 
%             LP_PIXEL=filtfilt(A,B,PIXEL-mean(PIXEL))+mean(PIXEL);
%             %median filter
%             MED_PIXEL=medfilt1(LP_PIXEL-mean(LP_PIXEL),median)+mean(LP_PIXEL);
%             Data(i,j,:)=MED_PIXEL;
%             end
%             end
%     waitbar(i/size(Data,1));
% end
% close(hdl);

    %     Sig = squeeze(nanmean(nanmean(Data)));
    %     Sigf = squeeze(nanmean(nanmean(FilteredData)));
    %     figure(3);cla; plot(Sig); hold on;plot(Sigf); legend('before LPF', strcat('after',num2str(cutoff_frequency),'Hz LPF'))
    %     anskeep = questdlg('Do you like to keep the changes?','Lowpass filter', 'YES', 'NO', 'YES');
    %     if strcmp(anskeep, 'YES')
    Data = FilteredData;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
%     end
    %     end
% end


Data(:,:,1:100)=repmat(Data(:,:,101),[1,1,100]);
Sig = squeeze(nanmean(nanmean(Data,1),2));
Sig = smooth(squeeze(Sig),0.01);
Sig2 = smooth(squeeze(Sig),0.5);
Sig3 = filter(ones(1,101)/101,1,Sig2);
Sig3(1:101)=Sig3(102);
Sig3 = Sig3+(max(Sig)-max(Sig3));
Sig = Sig-Sig3;
Sig = 100*(Sig-min(Sig))/(max(Sig)-min(Sig));
figure;plot(Sig)
%% correct for bleaching

bleachsignal (1,1,:) = Sig3;
bleachmatrix = repmat(bleachsignal,[size(Data,1),size(Data,2),1]);
CorrectedData = Data - bleachmatrix;
BeforCorrectionData = Data;
Data = CorrectedData;
Data = mat2gray(Data);
%  Sig = squeeze(nanmean(nanmean(Data,1),2));
% Sig = smooth(squeeze(Sig),0.01);
% figure;plot(Sig)
t=1:size(Sig,1);
% t=t*tresolution;
[I,J]=findpeaks(Sig,t,'MinPeakProminence',50,'Annotate','extents')
figure;findpeaks(Sig,t,'MinPeakProminence',50,'Annotate','extents')
hold on; plot(t(J(1)+mean(diff(J))/2-mean(diff(J))/5:J(1)+mean(diff(J))/2+mean(diff(J))/5),...
    Sig(J(1)+mean(diff(J))/2-mean(diff(J))/5:J(1)+mean(diff(J))/2+mean(diff(J))/5),'r')

for cl = 1: length(J)
    bgstart = floor(J(cl)+mean(diff(J))/2-mean(diff(J))/5);
    bgend = floor(J(cl)+mean(diff(J))/2+mean(diff(J))/5);
    if bgend>length(Sig)
        bgend = length(Sig);
    end
   BG(:,:,cl)=nanmean(Data(:,:,bgstart:bgend),3);
    
end
BGFrame = nanmean(BG,3);
% BGFrame = squeeze(BG(:,:,1));

 
    %% subtracting background
%     BGFrame = (mat2gray(nanmean(Data,3))+mat2gray(nanmax(Data,[],3)))./2;
%     BGFrame = mat2gray(nanmax(Data,[],3));
    background = repmat(BGFrame, [1 1 size(Data,3)]);
    Data_fluor = Data-background;
    Data_fluor(Data_fluor<0)=0;
    Data = mat2gray(Data_fluor);
    clear Data_fluor;
%% apply 3D median filter
    MData = medfilt3(Data,'replicate');
    Data =MData;
    GData = imgaussfilt3(Data,'Padding','replicate','FilterDomain','Frequency');
    Data = mat2gray(GData);
    implay(Data)
%            for i = 1:size(Data,3)
%             temp = squeeze(Data(:,:,i));
%             temp = imgaussfilt(temp,2);
%             temp = medfilt2(temp);
%             Data_binned(:,:,i) = temp;
%            end
%            Data= Data_binned;
%         Data = mat2gray(Data);
%     clear Data_binned;
           % binning
%     N=3;
%     [x,y]=meshgrid(-N:N,-N:N);
%     avePattern=exp(-(x.^2/(2*(0.7^2))+y.^2/(2*(0.7^2))));
%     N=size(avePattern,1);
%     Processing_progess = waitbar(0,'Processing....');
%     for k=1:10
%         waitbar(k/10);
%         for i = 1:size(Data,3)
%             temp = squeeze(Data(:,:,i));
%             temp = 1/N/N*conv2(temp,avePattern,'same');
%             Data_binned(:,:,i) = temp;
%         end
%         Data= Data_binned;
%         Data = mat2gray(Data);
%     end
%     clear Data_binned;
% close(Processing_progess)
%     %% spatial smoothing (Moving average)
%     mask = ones(2,2);
%     for i = 1:size(Data,3)
%         C = conv2(Data(:,:,i),mask,'same');
%         Data_MA(:,:,i) = C;
%     end
%     Data = Data_MA;
%     clear Data_MA;
%     Data = mat2gray(Data);
% Data = imgaussfilt(Data,0.5);

% if preproschb == 1
%     %     Correct = get(handles.checkbox9,'Value');
%     startp = get(handles.edit35,'String');
%     startp = str2double(startp);
%     stopp = get(handles.edit36,'String');
%     stopp = str2double(stopp);
%     Data=mat2gray(Data);
%     
%     %% Correct baseline
%     
%     %     if Correct==1
%     %         AVESignal(1,1:frames) = nanmean(nanmean(Data(:,:,1:frames)));
%     %
%     %         %%%%%%Linear trend: wts = [repmat(1/110,100,1)];
%     %         wts = [repmat(1/100,100,1)];
%     %         AVES = conv(AVESignal,wts,'valid');
%     %         BS = str2num(get(handles.edit73,'String'));
%     %         AVES = AVES-BS;
%     %         LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
%     %         AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
%     %         % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
%     %
%     %         for j=1:LAVE
%     %             Data(:,:,j) = Data (:,:,j)-AVES(j);
%     %         end
%     %     end
%     
%     
%     %             figure(4);map = colormap(gray(256));
%     %         for i = 1:size(Data,3)
%     %             temp = Data(:,:,i);
%     %             temp = imresize(temp,2);
%     %             imshow(temp,'Colormap',map);
%     %             title(num2str(i));drawnow
%     %         end
%     %         save(strcat(FullPath,'.mat'), 'Data','RawData');
%     set(handles.popupmenu5,'Value',2)
%     % end
%     %
%     % if PCAchb ==1
%     row=size(Data,1); col=size(Data,2); nFrames = size(Data,3);
%     PCAData = reshape(Data, row*col,nFrames)';
%     PCAData = PCAData - median(PCAData);
%     dec = mdwtdec('c',PCAData,5,'db2');
%     decBIS = chgwdeccfs(dec,'cd',0,1:2);
%     [XD,decDEN,THRESH] = mswden('den',decBIS,'sqtwolog','sln');
%     
%     Xbis = mdwtrec(decBIS);
%     
%     PCAData_reverse = reshape(Xbis', row,col,nFrames);
%     PCAData_reverse2 = reshape(XD', row,col,nFrames);
%     level = 5;
%     wname = 'bior6.8';
%     npc = 'heur';
%     [x_sim, qual, NPC] = wmspca(PCAData',level,wname,npc);
%     close(Processing_progess)
%     PCAData_reverse3 = reshape(x_sim, row,col,nFrames);
%     %     figure(3);cla;subplot(2,2,1);imshow(squeeze(Data(:,:,500))); title('before PCA');
%     %     subplot(2,2,2);imshow(squeeze(PCAData_reverse(:,:,500))); title('after PCA')
%     Sig = squeeze(nanmean(nanmean(Data)));
%     Sigp = squeeze(nanmean(nanmean(PCAData_reverse3)));
%     %     subplot(2,2,3:4);plot(Sig); hold on;plot(Sigp); legend('before PCA', 'after PCA')
%     %     anskeep = questdlg('Do you like to keep the changes?','PCA', 'YES', 'NO', 'YES');
%     %     if strcmp(anskeep, 'YES')
%     Data = PCAData_reverse;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');        
%     end
%     %     end
% end
% 
% if MAchb == 1
%     filter_span = str2num(get(handles.edit82,'String'));
%     FilteredData = filter(ones(1,filter_span)/filter_span,1,Data,[],3);
%     FilteredData(:,:,1:filter_span) = repmat(FilteredData(:,:,filter_span+1),[1,1,filter_span]);
%     %     Sig = squeeze(nanmean(nanmean(Data)));
%     %     Sigf = squeeze(nanmean(nanmean(FilteredData)));
%     %     figure(3);cla;plot(Sig); hold on;plot(Sigf); legend('before MA filter', 'after MA filter')
%     %     anskeep = questdlg('Do you like to keep the changes?','Moving Average filter', 'YES', 'NO', 'YES');
%     %     if strcmp(anskeep, 'YES')
%     Data = FilteredData;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');        
%     end
%     %     end
% end
% if SMchb == 1
%     SMOption = get(handles.popupmenu3,'Value'); %% 1 for smoothing - 2 for Averaging filter
%     % 3 for Gaussian filter - 4 for Median filter
%     % 5 for Adaptive Wiener filter
%     dwindowval = get(handles.popupmenu4,'Value');
%     dwindowlist = [3,5,7,9,11,13,15];
%     dwindow = dwindowlist(dwindowval);
%     % SMOOTH IMAGE:
%     smoothing_progess = waitbar(0,'Smoothing....');
%     for k=1:size(Data,3)
%         waitbar(k/size(Data,3));
%         if SMOption==1
%             pre_smoothed_image = squeeze(Data(:,:,k));
%             for i=1:size(pre_smoothed_image,1)
%                 for j=1:size(pre_smoothed_image,2)
%                     first_included_i_pixel = i-floor(dwindow/2);
%                     if first_included_i_pixel<1
%                         first_included_i_pixel = 1;
%                     end
%                     last_included_i_pixel = i+floor(dwindow/2);
%                     if last_included_i_pixel>size(pre_smoothed_image,1)
%                         last_included_i_pixel = size(pre_smoothed_image,1);
%                     end
%                     
%                     first_included_j_pixel = j-floor(dwindow/2);
%                     if first_included_j_pixel<1
%                         first_included_j_pixel = 1;
%                     end
%                     
%                     last_included_j_pixel = j+floor(dwindow/2);
%                     if last_included_j_pixel>size(pre_smoothed_image,2)
%                         last_included_j_pixel = size(pre_smoothed_image,2);
%                     end
%                     
%                     pixels_to_average = pre_smoothed_image(first_included_i_pixel:last_included_i_pixel,...
%                         first_included_j_pixel:last_included_j_pixel);
%                     
%                     pixels_to_average = pixels_to_average(pixels_to_average>=0);
%                     if length(pixels_to_average)>0
%                         SmoothedData(i,j,k) = nanmean(pixels_to_average);
%                     else
%                         SmoothedData(i,j,k) = Data (i,j,k);
%                     end
%                 end
%             end
%         elseif SMOption==2
%             h = fspecial('average',[dwindow dwindow]);
%             SmoothedData(:,:,k) = imfilter(Data(:,:,k),h,'replicate','same');
%         elseif SMOption==3
%             h = fspecial('gaussian',[dwindow dwindow],2);
%             SmoothedData(:,:,k) = imfilter(Data(:,:,k),h,'replicate','same');
%         elseif SMOption==4
%             SmoothedData(:,:,k) = medfilt2(Data(:,:,k),[dwindow dwindow],'replicate');
%         elseif SMOption==5
%             SmoothedData(:,:,k) = wiener2(Data(:,:,k),[dwindow dwindow]);
%         end
%     end
%     close(smoothing_progess)
%     %     figure(3);cla;subplot(2,2,1);imshow(squeeze(Data(:,:,500))); title('before smoothing');
%     %                   subplot(2,2,2);imshow(squeeze(SmoothedData(:,:,500))); title('after smoothing')
%     %     Sig = squeeze(nanmean(nanmean(Data)));
%     %     Sigs = squeeze(nanmean(nanmean(SmoothedData)));
%     %     subplot(2,2,3:4);plot(Sig); hold on;plot(Sigs); legend('before smoothing', 'after smoothing')
%     
%     %     anskeep = questdlg('Do you like to keep the changes?','Smoothing', 'YES', 'NO', 'YES');
%     %     if strcmp(anskeep, 'YES')
%     Data = SmoothedData;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
%     end
%     %     end
% end
% 
% 
% %         MaskedData = Data.*repmat(Mask,[1,1,frames]);
% %         Data = MaskedData;
% % % if bcorrchb ==1
%     Mask = double(Mask);
%     Mask( Mask == 0 ) = NaN;
%     Data = Data .* Mask;
%         h = mat2gray(squeeze(RawData(:,:,50)));
%     h= imresize(h,2);
% %     flipoption = get(handles.checkbox29,'Value');
% 
%     if flipoption == 1
%         h = fliplr(h);
%     end
% %     [x_coord,y_coord,intensity_val] = impixel(h);
% %     x_coord = floor(x_coord./2); y_coord = floor(y_coord./2);
% %     AVESignal= squeeze(Data (y_coord,x_coord,:));
%     AVESignal = squeeze(nanmean(nanmean(Data,1),2));
% % %     AVESignal
%     figure;plot(AVESignal);
%     h = drawpolyline(gca);
%     P=h.Position
%     X=P(:,1);
%     Y=P(:,2);
%     p = polyfit(X,Y,1);
% X2=1:size(Data,3);
% Y2=polyval(p,X2);
% bleachsignal = Y2;
% % %     bcorrmethod = get(handles.popupmenu6,'Value');
% %     if bcorrmethod==1
% %         %option1
% %         t = (1:length(AVESignal));
% %         opol = 6;
% %         [p,s,mu] = polyfit(t,AVESignal,opol);
% %         bleachsignal = polyval(p,t,[],mu);
% %         AVECorrected = AVESignal - bleachsignal;
% %         
% %         %         AVECorrected = detrend(AVESignal);
% %         
% %     elseif bcorrmethod ==2        %option2:
% %         frame_rate=1000*total_frames/total_time;
% %         cutoff_frequency = 2;
% %         filter_order = 1;
% %         normalized_cutoff_frequency = cutoff_frequency/(frame_rate/2);
% %         [b,a] = butter(filter_order,normalized_cutoff_frequency,'high');
% %         AVECorrected = filter(b,a,AVESignal);
% %         AVECorrected (1:50) = AVECorrected (51);
% %         bleachsignal = AVESignal - AVECorrected;
% %          CorrectedData = filter(b,a,Data,[],3);
% %     elseif bcorrmethod ==3 %option1:
% %         wts = [repmat(1/100,100,1)];
% %         AVES = conv(AVESignal,wts,'valid');
% %         % AVES = AVES-(AVESignal(1)-AVECorrected(1));
% %         LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
% %         AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
% %         % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
% %         Deltay = (AVESignal(1)-AVECorrected(1));
% %         AVES = AVES-(AVESignal(1)-AVECorrected(1));
% %         AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
% %         AVECorrected(900:1000)=AVECorrected(899);
% %         bleachsignal = AVES+(AVESignal(1)-Deltay);
% %         bleachsignal(900:1000)=bleachsignal(899);
% %     end
% %     %     figure;plot(AVESignal,'k');hold on; plot(AVECorrected,'g')
% %     %     legend('Bleached Signal','Corrected Signal')
%     Correcting_progess = waitbar(0,'Correct for Bleaching....');
%     
%     for i=1:size(Data,3)
%         waitbar (i/size(Data,3));
%         CorrectedData(:,:,i) = Data(:,:,i) - bleachsignal(i);
%     end
%     close (Correcting_progess)
%             BeforCorrectionData = Data;
%     Data = CorrectedData;
% %     if size(Data,3)>2000
% %         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
% %     else
% %         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
% %     end
% % end

% %         MaskedData = Data.*repmat(Mask,[1,1,frames]);
% %         Data = MaskedData;
% 
% min_manual_value = 0;
% max_manual_value = 1;
% 
% figure(4);
% if preproschb == 1
% else
% SelectedMap = get(handles.popupmenu5,'Value'); %% 1 for hot - 2 for gray - 3 for jet - 4 for summer
% if SelectedMap == 1
%     map = colormap(hot);
% elseif SelectedMap == 2
%     map = colormap(gray);
% elseif SelectedMap == 3
%     map = colormap(jet);
% elseif SelectedMap == 4
%     map = colormap(summer);
% end
% DataforPlot = mat2gray(Data);
% vraw = VideoWriter([fullfile(FullPath),'_raw','.avi'])
% vraw.FrameRate = 30;
% figure;
% open(vraw)
% for k=1:size(Data,3)
%     Frame = Data(:,:,k);
%     Frame = imresize(Frame,2);
%     figure(4);imshow(Frame,[],'Colormap',map); title (num2str(k));
%     %     caxis([min_manual_value max_manual_value]);
%     freezeColors_Bahar
%     title(num2str(k))
%     drawnow
%     M = getframe(gcf);
%     writeVideo(vraw,M)
% end
% close(vraw)
% Data = BeforCorrectionData;
% % save (strcat(FullPath,'.mat'),'RawData','Data','Mask');
% 
% % end
% Sig = squeeze(nanmean(nanmean(Data)));
% figure;plot(Sig)
% fprintf('Processing Complete!\n');
% 
