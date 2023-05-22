%% Processing and analysis of optical mapping data - Rose Lab - University of Calgary
%% Developed by: Motahareh Moghtadaei (motahareh.moghtadaei@ucalgary.ca) September 2016
%% Includes :
% %               -processing
% %               -colourmaps
% %               -conduction velocity
% %               -optical action potential morphology
%% Version2: Updated on January 30, 2019
%% Updates include:
% %               -The inverted signal is now supported
% %               -The new version is compatible with Matlab 2018a
% %               -The new version is compatible with Mac
%% Version3: Updated on February 24, 2019
%% Updates include:
% %               -Reading the .tif stacks
% %               -Baseline drift correction for bleaching
% %               -Quick representation of the raw data for quick
% %                 assessment of the quality of raw data
%% Version3.1: Updated on March 7, 2019
% %               -Bug fixes
%% Version3.2: Updated on March 19, 2019
% %               -The possibility to normalise the data based on a segment
% %                 with no artifact
% %               -The program calculates the CV even if one of the
% %                 neighboring pixels is unavailable
%% Version4: Updated on May 15, 2019
% %               -PCA processing, LPF, MA, and smoothing filters
% %                 implemented
% %               -ROI selection option
%% Version5: Updated on June 19, 2019
% %               -The option to choose between making the activation map either from a single beat or from average of all beats.
% %               -The conduction velocity values for each pixel in a user defined region of interest (polygon) outputted on an excel sheet.
% %               -Generating APD50, APD70 and APD90 maps and providing values in an excel sheet.
% %               -Generating the representative action potential form either a single pixel (defined by a click) or in a region of interest defined by a polygon.
%% Version6: Updated on November 14, 2019
% %               -The option to choose between calculating the conduction velocity either from a single beat or from average of all beats.
% %               -Calculation of conduction velocity for a region of interest defined by a click at the center and the number of surrounding pixels -instead of the fixed 7 by 7 pixels- (demonstrated by a box).
% %               -The option to choose between generating the APD maps and APD values either from a single beat or from average of all beats.
% %               -The option to generate representative OAPs either from a single beat or from average of all beats (both from a single pixel and user-defined ROI).
% %               -Improved conduction velocity quantification
% %               -Conduction Velocity vector field plot
%% Version6.1: Updated on November 17, 2019
% %               -The option to determine the acceptable range for conduction velocity by the user 
% %               -The option to calculate either Mean or median of the conduction velocity values in the ROI
% %               -Providing the histogram of the conduction velocity values in the ROI
%% Version7: Updated on January 16, 2020
% %               -Modified to load the data from MiCAM
% %               -Modified excel write command for mac
%% Version8: Updated on February 07, 2020
% %               -Bug fixes and MAC compatibility

%% Vesrion11:   Updated on September 01, 2020
% %               -Processing revised to be effective for higer frame rates too
% %               -CV measurment updated 
%% Vesrion12:   Updated on Nov 10, 2020

function varargout = opticalmapping(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @opticalmapping_OpeningFcn, ...
    'gui_OutputFcn',  @opticalmapping_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before opticalmapping is made visible.
% Initializing GUI parameters
function opticalmapping_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to opticalmapping (see VARARGIN)

% Choose default command line output for opticalmapping
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes opticalmapping wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = opticalmapping_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% h = imgmap('HandleVisibility','off');
% set(handles.pushbutton6,'Enable','on')
% Get default command line output from handles structure
% varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%********************** CREATE COLOR MAP  **********************%
function pushbutton5_Callback(hObject, eventdata, handles)
% close(figure(1));close(figure(2));close(figure(3));
clc
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
    FName = get(handles.edit28,'String');
end
if endp>size(Data,3)
    endp=size(Data,3)
end
AVERAGING = get(handles.checkbox2,'Value');

if AVERAGING
    Data = AveragedData;
else
    Data = Data(:,:,startp:endp);
end
% if max(size(RawData,1),size(RawData,2))>=200
%     Data =imresize(Data,0.5);
% end
Data = mat2gray(Data);
% figure; imshow(Data(:,:,50))

% %            binning
    N=3;
    [x,y]=meshgrid(-N:N,-N:N);
    avePattern=exp(-(x.^2/(2*(0.7^2))+y.^2/(2*(0.7^2))));
    N=size(avePattern,1);
    Processing_progess = waitbar(0,'Processing....');
    for k=1:10
        waitbar(k/10);
        for i = 1:size(Data,3)
            temp = squeeze(Data(:,:,i));
            temp = 1/N/N*conv2(temp,avePattern,'same');
            Data_binned(:,:,i) = temp;
        end
        Data= Data_binned;
        Data = mat2gray(Data);
    end
    clear Data_binned;
close(Processing_progess)
    %% spatial smoothing (Moving average)
    mask = ones(2,2);
    for i = 1:size(Data,3)
        C = conv2(Data(:,:,i),mask,'same');
        Data_MA(:,:,i) = C;
    end
    size(Data)
    Mask = squeeze(Mask(:,:,1));
    size(Mask)

Data = mat2gray(Data_MA);
    Data = Data.*repmat(Mask,[1,1,size(Data,3)]);

Data (isnan(Data))=0;


frames = size(Data,3);
pix_val_cutoff = get(handles.edit43,'String');
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;

startf = get(handles.edit44,'String');
startf = str2double(startf)
endf = get(handles.edit45,'String');
endf = str2double(endf)
startp =  1;
endp = frames;
Data_map = Data (:,:,startf-startp+1:endf-startp+1);
SMOption = get(handles.popupmenu9,'Value'); %% 1 for smoothing - 2 for Averaging filter
% 3 for Gaussian filter - 4 for Median filter
% 5 for Adaptive Wiener filter
dwindowval = get(handles.popupmenu10,'Value');
dwindowlist = [3,5,7,9,11,13,15];
dwindow = dwindowlist(dwindowval);
% SMOOTH IMAGE:
smoothing_progess = waitbar(0,'Smoothing....');
for k=1:size(Data_map,3)
    waitbar(k/size(Data_map,3));
    if SMOption==1
        pre_smoothed_image = squeeze(Data_map(:,:,k));
        for i=1:size(pre_smoothed_image,1)
            for j=1:size(pre_smoothed_image,2)
                first_included_i_pixel = i-floor(dwindow/2);
                if first_included_i_pixel<1
                    first_included_i_pixel = 1;
                end
                last_included_i_pixel = i+floor(dwindow/2);
                if last_included_i_pixel>size(pre_smoothed_image,1)
                    last_included_i_pixel = size(pre_smoothed_image,1);
                end
                
                first_included_j_pixel = j-floor(dwindow/2);
                if first_included_j_pixel<1
                    first_included_j_pixel = 1;
                end
                
                last_included_j_pixel = j+floor(dwindow/2);
                if last_included_j_pixel>size(pre_smoothed_image,2)
                    last_included_j_pixel = size(pre_smoothed_image,2);
                end
                
                pixels_to_average = pre_smoothed_image(first_included_i_pixel:last_included_i_pixel,...
                    first_included_j_pixel:last_included_j_pixel);
                
                pixels_to_average = pixels_to_average(pixels_to_average>=0);
                if length(pixels_to_average)>0
                    SmoothedData(i,j,k) = nanmean(pixels_to_average);
                else
                    SmoothedData(i,j,k) = Data_map (i,j,k);
                end
            end
        end
    elseif SMOption==2
        h = fspecial('average',[dwindow dwindow]);
        SmoothedData(:,:,k) = imfilter(Data_map(:,:,k),h,'replicate','same');
    elseif SMOption==3
        h = fspecial('gaussian',[dwindow dwindow],2);
        SmoothedData(:,:,k) = imfilter(Data_map(:,:,k),h,'replicate','same');
    elseif SMOption==4
        SmoothedData(:,:,k) = medfilt2(Data_map(:,:,k),[dwindow dwindow],'replicate');
    elseif SMOption==5
        SmoothedData(:,:,k) = wiener2(Data_map(:,:,k),[dwindow dwindow]);
    end
end
close(smoothing_progess)


%     anskeep = questdlg('Do you like to keep the changes?','Smoothing', 'YES', 'NO', 'YES');
%     if strcmp(anskeep, 'YES')
Data_map = SmoothedData;
Data_map = mat2gray(Data_map);
Data_map = mat2gray(Data_map,[mean(mean(mean(Data_map))),1]);

k=10; % the grey backgroud image

Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);

% xlim = [1,Nframes];

if series==1
    i0=findtime(frames,k)
    % A = strcat(PathName, i0, num2str(k), '.tif');
    A = fullfile(PName, strcat(FName, i0, num2str(k), '.tif'));
    h = imread(A);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
elseif stack==1
    h = squeeze(RawData(:,:, k));
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
end
% if max(size(RawData,1),size(RawData,2))>=200
%     h =imresize(h,0.5);
% end
h = mat2gray(h);
H(:,:,1)=h;
H(:,:,2)=h;
H(:,:,3)=h;
scrsz = get(0,'ScreenSize');
% f1=figure('Position',[2*scrsz(3)/3 scrsz(4)/3 10 scrsz(4)/2]);imshow(H);
f2=figure('Position',[2*scrsz(3)/3 2*scrsz(4)/3 10 scrsz(4)/2],'Visible','off');imshow(imresize(H,3)); hold on;
% f3=figure('Position',[50 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);
% global colorMAP
% colorMAP = figure;

pix_val_min = 0;
pix_val_max = 0;
% max_nframes = 27;
% PathName = get(handles.edit6,'String');
% PathName = 'F:\reentry data\2014-7-22\2-analysis\isomanpstim-3-filter-good\corrected\map1\corrected\numbered\';
%
% A = strcat(PathName, '*.tif');
% list_of_frames = dir(A);
% frames = numel(list_of_frames);
frames = size(Data_map,3);
COLOR=[];
s1 = 'Loading ';
s2 = num2str(frames);
s3 = ' frames from ';
s4 = FullPath;
str = [s1 s2 s3 s4];
disp(str);
str = [ ];
disp(str);
pt=frames;
max_nframes = frames;
% max_nframes = 51;
%Set trace sensitivity value
% pix_val_cutoff = get(handles.edit1,'String');

pix_val_cutoff = str2num(pix_val_cutoff);

%Set frame step size
% frame_step_size = get(handles.edit2,'String');
frame_step_size = '1';
frame_step_size = str2num(frame_step_size);

%Set color spectrum
% color = get(handles.edit3,'String');
color = '256';
color = str2num(color);

%Set delta values for countour boundry detection
% delta_row = get(handles.edit4,'String');
delta_row = '2';
delta_row = str2num(delta_row);

% delta_col = get(handles.edit5,'String');
delta_col = '2';
delta_col = str2num(delta_col);

color_const = color;
% color_const = frames;

% mapVal = get(handles.popupmenu1,'Value');
mapVal = 3;
if mapVal == 1
    map = colormap(gray(color_const));
end
if mapVal == 2
    map = colormap(hsv(color_const));
end
if mapVal == 3
    map = colormap(jet(color_const));
end
if mapVal == 4
    map = colormap(cool(color_const));
end
if mapVal == 5
    map = colormap(winter(color_const));
end

MAP=map;
map=map([color_const:-1:1],:);
% map=MAP(steps*frames+(color_const-steps*frames):-1:steps+(color_const-steps*frames),:);
colormap(map);
% color_lastlayer = frames;
% color_const
% steps = floor(tresolution*(color_const/30))
steps = floor(color_const/max_nframes); % color increase steps
% steps=15;

for k = 1:frames
    I = squeeze(Data_map(:,:,pt));
    I = imresize(I , 3);
    %     I=256*I(:,:,1)/nanmax(nanmax(I(:,:,1)));
    img_param = size(I);
    img_area = img_param(1)*img_param(2);
    
    %Remove saturated data
    if I(1) > 0.7
        for i = 1:img_area
            if I(i) > 0.7
                I(i) = 0;
            end
        end
    end
    %pixel intensity cutoff
    for i = 1:img_area
        if I(i) < pix_val_cutoff
            I(i) = 0;
        end
    end
    %     end
    
    I_sub = I;
    %     B = strcat('C:\Users\John\Desktop\sub\', num2str(pt), '.tif');
    x = size(I_sub);
    wid_I_sub = x(1);
    
    pix_area = x(1) * x(2);
    
    pix_val = I_sub(wid_I_sub);
    
    for j=1:wid_I_sub
        if pix_val < pix_val_min
            pix_val_min = pix_val;
        end
        
        if pix_val > pix_val_max
            pix_val_max = pix_val;
        end
    end
    
    Final_img = I_sub;
    if pt>=1
        s=size(Final_img);
        for row = 1:delta_row:s(1)
            for col=1:delta_col:s(2)
                if Final_img(row,col)
                    break;
                end
            end
            figure(f2);
            hold on;
            contour = bwtraceboundary(Final_img, [row, col],'e',8,50000,'counterclockwise');
            %             figure(f3);
            %             hold on;
            %             contour = bwtraceboundary(Final_img, [row, col], 'e', 8, 500,...
            %                 'counterclockwise');
            %Uncommenting the next line will cause the entire color map range
            %to be utilized.
            %             color = color - floor(color_const/frames);
            
            %Uncommenting the next line causes a the progression of the colors
            %to increase in a standardized step size. Entire color map range might
            %not be utilized.
            
            %             color = (frames*15+1) - k*15
            % color = k*floor(color_const/frames);
            % color=color_const-(k-1)*5;
            %
            color=steps*k+(color_const-steps*frames);
            %
            if(~isempty(contour))
                figure(f2);
                hold on;
                q = fill(contour(:,2),contour(:,1),MAP(color,:),'LineWidth',1,'LineSmoothing','off','EdgeColor','none');
                %                 axis([0 80 -100 10]);
                axis image;
                axis off;
                %                 figure(f3);
                %                 hold on;
                %                 q = fill(contour(:,2),-1*contour(:,1),MAP(color,:),'LineWidth',1,'LineSmoothing','off','EdgeColor','none');
                %                 %                 axis([0 80 -100 10]);
                %                 axis image;
                %                 axis off;
                %                 E = strcat('C:\Users\John\Desktop\ContourFrames\', num2str(pt), '.tif');
            end
        end
    end
    s1 = 'Processed frame ';
    s2 = num2str(pt);
    s = [s1 s2];
    disp(s);
    pt = pt - frame_step_size;
    COLOR=[COLOR, color]   ;
end
figure(f2); title(['" ** total time = "', num2str(frames*tresolution), '" ms ** "'])
% figure(f3); title(['map from "', FName,'" ** total time = "', num2str(frames*tresolution), '" ms ** "', num2str(frames),'" frames were loaded', ' ** frames: ',num2str(startf),' - ',num2str(endf)])
colormap(map);

figure(f2);colorbar('YTick',[linspace(0,1,2)], 'YTickLabel',{' 0 ms',strcat(num2str(ceil(frames*tresolution*100)/100),' ms')});


%***********************************************************************%

%**************************** EXIT BUTTON ******************************%
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to  pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
close all
closereq
%***********************************************************************%

%********************** [Color Map Panel]->[SAVE] **********************%
function pushbutton8_Callback(hObject, eventdata, handles)
clc
global colorMAP
saveimg_path = get(handles.edit19,'String');
saveas(colorMAP,saveimg_path);
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%***********************************************************************%

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a
%        double
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a
%        double

function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% OAP and APD
%****************** Optical Action Potential and duration ******************%
function pushbutton12_Callback(hObject, eventdata, handles)
clc
Averaging=get(handles.checkbox2,'Value')
manual_mode=1;
% Averaging=1;
APDT=50 ;     % choose 50 for APD50, 70 for APD70, and 90 for APD90
manual_RMP=0;
RMP=6.2;
t_RMP=88.82;


series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');

if series == 1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')]
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack == 1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end



total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)
if endp>size(Data,3)
    endp = size(Data,3); 
end
Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
xlim = [1,Nframes];
NumFrame = get(handles.edit62,'String');
NumFrame = str2double(NumFrame);

% i0=findtime(frames,NumFrame)
% A = strcat([fullfile(PName, FName, 'results2','corrected', strcat( i0, num2str(NumFrame))), '.tif']);
%
% h = imread(A);
if series == 1
    load ([fullfile(PName,'results2','corrected','results',FName),'_improved.mat'])
    h = Data(:,:,NumFrame);
    h = imresize(h,2);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
elseif stack ==1
    h = mat2gray(squeeze(RawData(:,:,50)));
    h = imresize(h,2);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
end

% h = Data(:,:,NumFrame);
%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    Sig=[];
    if manual_mode
        figure
        [x_coord,y_coord,intensity_val] = impixel(h);
        x_coord = floor(x_coord./2); y_coord = floor(y_coord./2);
        x= Data (y_coord,x_coord,startp:endp);
    else
        x= Data (29,25,startp:endp);
    end
    x=squeeze(x);
%         Sig = x;
%     Sig = smooth(squeeze(Sig),0.01);
%     Sig2 = smooth(squeeze(Sig),0.5);
%     Sig3 = filter(ones(1,101)/101,1,Sig2);
%     Sig3(1:101)=Sig3(102);
%     Sig3 = Sig3+(max(Sig)-max(Sig3));
%   x=x-Sig3;
% %     figure;plot(x)
    x= smoothdata(x,'sgolay');
%     hold on;plot(x)
    
    if Averaging
       figure; plot(x,'linewidth',1)
       xlabel('frame number')
    end
%     x= smooth(x,0.008);
    frame_rate=1000*total_frames/total_time;
%     Signal=[];
% % smooth x data
Signal(1,:) = x;
% Signal(1,:) = nanmean(Sig);
S2 = smooth(1:size(x),Signal,0.9,'rloess');
Signal = Signal-S2'; % drift removal
window = 5;
mask = ones(1,window)/window;
maY=conv(Signal,mask,'same'); % temporal smoothing
x=maY(1,:);
x= smooth(x,0.01);

    OAP_signal=x-nanmin(x);
    OAP_signal=OAP_signal';
    % %
    t=1:size(OAP_signal,2);
    t=t*tresolution;
    I=t/tresolution;
    OAP = 100*(OAP_signal-min(OAP_signal))/(max(OAP_signal)-min(OAP_signal));
    tt  = t;
    II  = I;
    size(tt)
    size(OAP)
    [I,J]=findpeaks(OAP,II,'MinPeakProminence',50,'Annotate','extents')
    figure;findpeaks(OAP,II,'MinPeakProminence',50,'Annotate','extents')
        Starting_frame = get(handles.edit35,'String');
%     Starting_frame = str2double(Starting_frame);
%     Ending_frame = get(handles.edit36,'String');
%     Ending_frame = str2double(Ending_frame);
%     if Ending_frame>size(Data,3)
%         Ending_frame =size(Data,3)
%     end
    if length(I)>3 && Averaging
        Y=[];
        k=0;
        CL=diff(J)
        halfCL = floor(mean(diff(J))./2);
        for j=1:length(J)-1
            if (J(j)-halfCL > 0) && (floor(J(j+1)+halfCL))<length(OAP)
                k=k+1;
%                 OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2)),OAP(IMax(l+1)-floor(CL./2):IMax(l+1)+floor(CL./2))
% size(OAP)
% Starting_frame+floor(J(j)-halfCL)
% Starting_frame+floor(J(j)+halfCL)
% Starting_frame+floor(J(j+1)-halfCL)
% Starting_frame+floor(J(j+1)+halfCL)
                Y=[Y;OAP(floor(J(j)-halfCL):floor(J(j)+halfCL)),OAP(floor(J(j+1)-halfCL):floor(J(j+1)+halfCL))];
                %             AData(:,:,:,k) = Data(:,:,(Starting_frame+floor(J(j)-halfCL):Starting_frame+floor(J(j)+halfCL)));
            end
        end
        OAP_new = nanmean(Y)-nanmin(nanmean(Y));
        
        tt_new=tt(1:length(OAP_new));
        II_new=II(1:length(OAP_new));
    else
        OAP_new = OAP;
        tt_new=tt;
        II_new=II;
    end
    
%     [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
%     CL=nanmean(diff(tMax))
%     CL=CL/tresolution;
%     if Npeaks>3 && Averaging
%         Y=[];
%         cycles = 0;
%         CL=ceil(CL);
%         for l=1:length(IMax)-1
%             if (IMax(l)-floor(CL./2) > 0) & (IMax(l+1)+floor(CL./2)<length(OAP))
%                 Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2)),OAP(IMax(l+1)-floor(CL./2):IMax(l+1)+floor(CL./2))];
%                 cycles = cycles +1;
%             end
%         end
%         
% %         OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
%                 OAP_new = nanmean(Y)-nanmin(nanmean(Y));
% 
%         tt_new=tt(1:length(OAP_new));
%         II_new=II(1:length(OAP_new));
%     else
%         OAP_new = OAP;
%         tt_new=tt;
%         II_new=II;
%     end
    figure;plot(Y')
    OAP_new=OAP_new-nanmin(OAP_new);
    
    OAP_new = (OAP_new-nanmin(OAP_new))/(nanmax(OAP_new)-nanmin(OAP_new))*100;
    
%     [pks,locs]=findpeaks(OAP_new,tt_new,'MaxPeakProminence',10)
    
%     if Averaging
%     SMspan = 0.1
% else
%     SMspan = 0.02
%     end
%     figure; plot(OAP_new)
% %     OAP_new = smooth(OAP_new,SMspan)';
%     OAP_new = smoothdata(OAP_new);
%         hold on; plot(OAP_new)

%     figure; plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');
%     hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
%     axislimits=nanmax(OAP_new)-nanmin(OAP_new);
%     axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
[yimax, tjmax]=findpeaks(OAP_new,'MinPeakProminence',max(OAP_new)/2,'Annotate','extents');
D=diff(OAP_new);
% figure;plot(OAP_new);hold on;plot(D)
D=smooth(D,0.03)';
Dpart1 = D(1:tjmax);
Dpart1(Dpart1<0.3)=0;
D(1:tjmax)=Dpart1;
% hold on; plot(D)
S=sign(D);
% hold on; plot(S)
SD=diff(S);
SD2=diff(SD);
% hold on; plot(SD)
% hold on; plot(SD2)
SDS=sign(SD);
SDS2=-sign(SD2);
SDS=[0,SDS,0];
SDS2=[0,SDS2,0,0];
[SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(tt_new, SDS, 'model');
SIMax;
scrsz = get(0,'ScreenSize');figure('Position',[scrsz(3)/3 scrsz(4)/1.8 scrsz(3)/2 scrsz(4)/3]);
plot(tt_new,OAP_new,'k','linewidth',1);xlabel('time (ms)');
[Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt_new, OAP_new, 'data');
IMax;

%%%%%%%%%%%%%
% OAPEach=[]; DataAverage=zeros(size(Data,1),size(Data,2),length(-floor(CL./2):floor(CL./2))); cycles=0;
% for k=1:length(IMax)
% %     IMax(k)-floor(CL./2) > 0
% %     IMax(k)+floor(CL./2)<length(OAP_new)
%     if (IMax(k)-floor(CL./2) > 0) & (IMax(k)+floor(CL./2)<length(OAP_new))
%    OAPEach = [OAPEach ; OAP_new(IMax(k)-floor(CL./2):IMax(k)+floor(CL./2))];
%    DataAverage = DataAverage + Data(:,:,IMax(k)-floor(CL./2):IMax(k)+floor(CL./2));
%    cycles = cycles +1;
%     end
%
% end
% DataAverage = DataAverage ./cycles;
% cycles
%
% size(DataAverage)
% Sig = squeeze(nanmean(nanmean(DataAverage)));
% Sig
% figure;plot(Sig');
% size(OAPEach)
% figure;plot(OAPEach')
% figure;
%%%%%%%%%%%%%
a1=find(SIMin==IMax(1));
b1=find(SIMax<IMax(1));
if size(b1,2)>1
    b1=b1(end-1);
end
b2=find(SIMax>=IMax(1));
b2=b2(1);
% APA   = OAP_new(SIMin(a1));
APA = yimax(1);
% t_APA = tt_new(SIMin(a1));
t_APA = tt_new(tjmax(1));
% I_APA = SIMin(a1);
I_APA = tjmax(1);
if manual_RMP
    I_RMP=ceil(t_APA/tresolution);
else
    RMP   = OAP_new(SIMax(b2));
    t_RMP = tt_new(SIMax(b2));
    I_RMP = SIMax(b2);
end
TP = OAP_new(SIMax(b1));
t_TP = tt_new(SIMax(b1));
I_TP = SIMax(b1);
hold on; plot(t_APA,APA,'rs','markerfacecolor','r');
hold on; plot(t_RMP,RMP,'go','markerfacecolor','g');
hold on; plot(t_TP,TP,'c^','markerfacecolor','c');
hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
axislimits=nanmax(OAP_new)-nanmin(OAP_new);
axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
% legend('Action Potential (AP)','Action Potential Amplitude (APA)','Rest Membrane Potential (RMP)','Take-off Potential (TP)' )
%% APD50
AP50= (((50)/100)*(APA-TP))+TP
tt_new(end)
hold on; plot([0, tt_new(end)],[AP50, AP50] , 'k--')

hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
C=AP50*ones(size(OAP_new));
D=diff(sign(OAP_new-C));
hold on; plot(tt_new,sign(OAP_new-C));
hold on; plot(tt_new(2:end),D,'r');
I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
t2a_APD50=tt_new(I2a_APD50);t2b_APD50=tt_new(I2b_APD50);
OAP2a_APD50=OAP_new(I2a_APD50);OAP2b_APD50=OAP_new(I2b_APD50);
T1=((AP50-OAP1a_APD50(1))/(OAP1b_APD50(1)-OAP1a_APD50(1)))*(t1b_APD50(1)-t1a_APD50(1))+t1a_APD50(1);
T2=((AP50-OAP2a_APD50(1))/(OAP2b_APD50(1)-OAP2a_APD50(1)))*(t2b_APD50(1)-t2a_APD50(1))+t2a_APD50(1);
T3=((AP50-OAP1a_APD50(2))/(OAP1b_APD50(2)-OAP1a_APD50(2)))*(t1b_APD50(2)-t1a_APD50(2))+t1a_APD50(2);
T4=((AP50-OAP2a_APD50(2))/(OAP2b_APD50(2)-OAP2a_APD50(2)))*(t2b_APD50(2)-t2a_APD50(2))+t2a_APD50(2);
hold on; plot(T2, AP50,'d', 'markeredgecolor', [0.25,0.75,0.75], 'Markerfacecolor',[0.25,0.75,0.75])
t_APA
nanmin(OAP_new)-5
nanmax(OAP_new)+20
hold on; plot([t_APA,t_APA],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
hold on; plot([T2,T2],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
APD = T2-t_APA
text(20,85,['APD50',' = ', num2str(APD)], 'color', [0.25,0.75,0.75])

%% APD 70

AP50= (((30)/100)*(APA-TP))+TP;
hold on; plot([0, tt_new(end)],[AP50, AP50] , 'k--')
hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
C=AP50*ones(size(OAP_new));
D=diff(sign(OAP_new-C));
hold on; plot(tt_new,sign(OAP_new-C));
hold on; plot(tt_new(2:end),D,'r');
I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
t2a_APD50=tt_new(I2a_APD50);t2b_APD50=tt_new(I2b_APD50);
OAP2a_APD50=OAP_new(I2a_APD50);OAP2b_APD50=OAP_new(I2b_APD50);
T1=((AP50-OAP1a_APD50(1))/(OAP1b_APD50(1)-OAP1a_APD50(1)))*(t1b_APD50(1)-t1a_APD50(1))+t1a_APD50(1);
T2=((AP50-OAP2a_APD50(1))/(OAP2b_APD50(1)-OAP2a_APD50(1)))*(t2b_APD50(1)-t2a_APD50(1))+t2a_APD50(1);
T3=((AP50-OAP1a_APD50(2))/(OAP1b_APD50(2)-OAP1a_APD50(2)))*(t1b_APD50(2)-t1a_APD50(2))+t1a_APD50(2);
T4=((AP50-OAP2a_APD50(2))/(OAP2b_APD50(2)-OAP2a_APD50(2)))*(t2b_APD50(2)-t2a_APD50(2))+t2a_APD50(2);
hold on; plot(T2, AP50,'d', 'markeredgecolor', [1,0.5,0.3], 'Markerfacecolor',[1,0.5,0.3])
hold on; plot([t_APA,t_APA],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
hold on; plot([T2,T2],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
APD = T2-t_APA
text(20,75,['APD70',' = ', num2str(APD)],'color',[1,0.5,0.3])

%% APD90
AP50= (((10)/100)*(APA-TP))+TP;
hold on; plot([0, tt_new(end)],[AP50, AP50] , 'k--')
hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
C=AP50*ones(size(OAP_new));
D=diff(sign(OAP_new-C));
hold on; plot(tt_new,sign(OAP_new-C));
hold on; plot(tt_new(2:end),D,'r');
I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
t2a_APD50=tt_new(I2a_APD50);t2b_APD50=tt_new(I2b_APD50);
OAP2a_APD50=OAP_new(I2a_APD50);OAP2b_APD50=OAP_new(I2b_APD50);
T1=((AP50-OAP1a_APD50(1))/(OAP1b_APD50(1)-OAP1a_APD50(1)))*(t1b_APD50(1)-t1a_APD50(1))+t1a_APD50(1);
T2=((AP50-OAP2a_APD50(1))/(OAP2b_APD50(1)-OAP2a_APD50(1)))*(t2b_APD50(1)-t2a_APD50(1))+t2a_APD50(1);
T3=((AP50-OAP1a_APD50(2))/(OAP1b_APD50(2)-OAP1a_APD50(2)))*(t1b_APD50(2)-t1a_APD50(2))+t1a_APD50(2);
T4=((AP50-OAP2a_APD50(2))/(OAP2b_APD50(2)-OAP2a_APD50(2)))*(t2b_APD50(2)-t2a_APD50(2))+t2a_APD50(2);
hold on; plot(T2, AP50,'d', 'markeredgecolor', [0.7,0.4,1], 'Markerfacecolor',[0.7,0.4,0.5])
hold on; plot([t_APA,t_APA],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
hold on; plot([T2,T2],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
APD = T2-t_APA
text(20,65,['APD90',' = ', num2str(APD)],'color',[0.7,0.4,1])
% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/3 scrsz(4)/10 scrsz(3)/2 scrsz(4)/3])
% plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');
% hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
% axislimits=nanmax(OAP_new)-nanmin(OAP_new);
% axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double

function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%******* Conduction Velocity **********%
function pushbutton17_Callback(hObject, eventdata, handles)

clc
Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end

% AVERAGING = 1;
AVERAGING=get(handles.checkbox2,'Value')
if AVERAGING
    Data = AveragedData;
else
    Data = Data(:,:,startp:endp);
end

% Data_OAP = Data ;
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
total_pixels = get(handles.edit52,'String');
total_pixels = str2double(total_pixels);
if get(handles.checkbox30,'Value')==1
    total_pixels = total_pixels/2;
end
total_distance = get(handles.edit53,'String');
total_distance = str2double(total_distance);
sresolution=total_distance/total_pixels;



xlim = [1,Nframes];
FirstFrame = get(handles.edit62,'String');
FirstFrame = str2double(FirstFrame);
% [x]= signalprep(PathName,FileName,'manual',[54,55], xlim);


%Upstroke Velocity Vector
global UV;
UV = [];

% OAP_chkbx_status = get(handles.checkbox1,'Value');
OAP_chkbx_status =0;

FAMI = [];
FA_50_I = [];
x_coord_array = [];
y_coord_array = [];

%Record Optical Action Potentials
pt = 1;
% array_area = get(handles.edit14,'String');
% array_area = str2num(array_area);
% array_area = 3;
array_area = floor(str2num(get(handles.edit85,'String'))/2);
% FileName = get(handles.edit17,'String');
% PathName = get(handles.edit16,'String');
ImScaleF = 2;
c1=figure;
if series
    NumFrame = '0073';
    A = strcat(PathName, NumFrame, '.tif');
    h = imread(A);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
    figure(c1);
    [x_coord,y_coord,intensity_val] = impixel(h);
elseif stack
    tempframe = squeeze(RawData(:,:,5));
    tempframe = uint16(double(tempframe).*Mask);
    tempframe = imresize(tempframe,ImScaleF);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        tempframe = fliplr(tempframe);
    end
    figure(c1); [x_coord,y_coord,intensity_val] = impixel(tempframe);
    figure(c1); hold on; rectangle('Position',[x_coord-ImScaleF*array_area y_coord-ImScaleF*array_area ImScaleF*array_area ImScaleF*array_area],'EdgeColor','r')
    x_coord = floor(x_coord / ImScaleF);
    y_coord = floor(y_coord / ImScaleF);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x_coord1=x_coord;y_coord1=y_coord;

% x= Data (y_coord,x_coord,:);
%
% hp = impixelinfo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct the 5x5 pixel array
x_init = x_coord - array_area;
y_init = y_coord - array_area;

x_final = x_coord + array_area;
y_final = y_coord + array_area;
% %
x_curr = x_init;
y_curr = y_init;
% %
% %
% % % wait_h = waitbar(0,'Calculating conduction velocity in selected region...');
% %
idx1 = 1;
idx2 = ((array_area*2)+1)^2;
CLs=[];
SIGNAL=[];
slope=[];

x_coord_array = [];
y_coord_array = [];
% figure;
RT_before=0;
c2=figure;
Data = mat2gray(Data);

for ypos = y_curr:y_final
    for xpos = x_init:x_final
        if ~isnan(Mask(ypos,xpos))
            current_x = xpos;
            current_y = ypos;
            current_pixel = [current_x,current_y];
            %             x=ones(1,size(Data,3));
            %             x(1,:)= Data (current_y,current_x,:);
            OAP = squeeze(Data(current_y,current_x,:))-nanmin(squeeze(Data(current_y,current_x,:)));
%             figure;plot(OAP)
            OAP = smooth(OAP,10,'sgolay');
            OAP = OAP';
            I=1:length(OAP);
            t=I*tresolution;
            figure(c2);
             subplot(array_area*2+1,array_area*2+1,idx1); plot(t,OAP,'k','linewidth',2);
             hold on; plot([0, t(end)],[0, 0] , 'k--')
             axis([t(1) t(end) nanmin(OAP)-abs(nanmin(OAP)/4) nanmax(OAP)+abs(nanmax(OAP)/4)])
             drawnow
             
            AverageOAP = squeeze(mean(mean(Data)))-nanmin(squeeze(mean(mean(Data))));
%             figure;plot(AverageOAP);
            if AVERAGING
            OAPth = min(AverageOAP)+(0.1*(max(AverageOAP)-min(AverageOAP)));
            else 
                OAPth = 0.2*max(AverageOAP);
            end
%             OAPth
            I_AverageOAPth =find(AverageOAP>=OAPth);
%             figure;plot(AverageOAP)
%             size(OAP)
%             I_AverageOAPth(1)
            OAP(1:I_AverageOAPth(1)-20)=0;
            OAPinterp = interp(OAP,10);
            OAPinterp(1:2) = 0;
            Iinterp = interp(I,10);
            I_OAPth =find(OAPinterp>=0.1);
            if ~isempty(I_OAPth)
                if Iinterp(I_OAPth(1))~=1
                    Iinterp_OAPth =  Iinterp(I_OAPth(1));
                elseif length(I_OAPth>2)
                    Iinterp_OAPth =  Iinterp(I_OAPth(2));
                end
                rt=Iinterp_OAPth;
                if rt~=0
                    FA_50_I = [FA_50_I,rt]; % rising times
                    x_coord_array = [x_coord_array,xpos];
                    y_coord_array = [y_coord_array,ypos];
                end
                
                title(rt)

            end
            
            
%             tt  = t;
%             II  = I;
%             
%             [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
%             maximums=tMax;
%             Npeaks
%             CL=nanmean(diff(tMax));
%             CL=CL/tresolution;
%             % if (Npeaks>2) && strcmp(type1, 'Averaging')
%             %             if Npeaks>2 && AVERAGING
%             %                 Y=[];
%             %                 cycles = 0;
%             %                 CL=ceil(CL);
%             %                 for l=1:length(IMax)
%             %                     if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP))
%             %                         Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
%             %                         cycles = cycles +1;
%             %                     end
%             %                 end
%             %
%             %                 OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
%             %                 tt_new=tt(1:length(OAP_new));
%             %                 II_new=II(1:length(OAP_new));
%             %             else
%             %                 OAP_new = OAP;
%             %                 tt_new=tt;
%             %                 II_new=II;
%             %             end
%             OAP_new = OAP;
%             tt_new=tt;
%             II_new=II;
%             OAP_new=OAP_new-nanmin(OAP_new);
%             D=diff(OAP_new);
%             S=sign(D);
%             SD=diff(S);
%             SDS=sign(SD);
%             SDS=[0,SDS,0];
%             [SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(tt_new, SDS, 'model');
%             [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt_new, OAP_new, 'data');
%             A1=[];B1=[];B2=[];
%             for k3=1:Npeaks
%                 a1=find(SIMin==IMax(k3));
%                 b1=find(SIMax<=IMax(k3));
%                 if numel(b1)~=0
%                     b1=b1(end);
%                 end
%                 b2=find(SIMax>=IMax(k3));
%                 if numel(b2~=0)
%                     b2=b2(1);
%                 end
%                 A1=[A1,a1];B1=[B1,b1];B2=[B2,b2];
%             end
%             a1=A1;b1=B1;b2=B2;
%             APA   = OAP_new(SIMin(a1));
%             t_APA = tt_new(SIMin(a1));
%             I_APA = SIMin(a1);
%             RMP   = OAP_new(SIMax(b2));
%             t_RMP = tt_new(SIMax(b2));
%             I_RMP = SIMax(b2);
%             TP = OAP_new(SIMax(b1));
%             t_TP = tt_new(SIMax(b1));
%             I_TP = SIMax(b1);
%             APDT=50;        % choose 50 for APD50, 70 for APD70, and 90 for APD90
%             if numel(TP)==numel(APA)
%                 AP50= ((100-APDT)/100)*(APA-TP);
%                 AP50= nanmean(AP50);
%                 C=AP50*ones(size(OAP_new));
%                 D=diff(sign(OAP_new-C));
%                 
%                 I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
%                 I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
%                 
%                 t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
%                 OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
%                 RT = t_TP;
%                 rt=RT(1);
%                 figure(c2);
%                 size(tt_new,2)
%                 size(OAP_new,2)
%                 if size(tt_new)~=size(OAP_new,2)
%                     subplot(array_area*2+1,array_area*2+1,idx1); plot(tt_new(1:end-1),OAP_new,'k','linewidth',2);
%                     %                 ;xlabel('time (ms)');
%                     hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
%                     axis([tt_new(1) tt_new(end) nanmin(OAP_new)-abs(nanmin(OAP_new)/4) nanmax(OAP_new)+abs(nanmax(OAP_new)/4)])
%                     drawnow
%                 else
%                     subplot(array_area*2+1,array_area*2+1,idx1); plot(tt_new,OAP_new,'k','linewidth',2);
%                     %                 ;xlabel('time (ms)');
%                     hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
%                     axis([tt_new(1) tt_new(end) nanmin(OAP_new)-abs(nanmin(OAP_new)/4) nanmax(OAP_new)+abs(nanmax(OAP_new)/4)])
%                     drawnow
%                 end
%                 title(rt)
%                 FA_50_I = [FA_50_I,rt]; % rising times
%                 CLs=[CLs,CL];
%                 %Calculating the slope of discrete points
%                 slope = AP50/rt;
%                 
%                 UV = [UV; slope];
%                 x_coord_array = [x_coord_array,xpos];
%                 y_coord_array = [y_coord_array,ypos];
%                 RT_before=RT;
%             end
            %         waitbar(idx1/idx2);
            idx1 = idx1 + 1;
        end
    end
end
xlabel('time (ms)');
CLs
% % % close(wait_h)
% % % wait_h2 = waitbar(1,'Conduction velocity calculation COMPLETE.');
% % % hold off
% %
% %
%Eliminate outliers in optical action potential
mean_OAP = nanmean(FA_50_I);
q = 1;
clean_OAP=[];
for w = 1:length(FA_50_I)
    if (FA_50_I(w) < mean_OAP-5)
        clean_OAP(w) = mean_OAP;
    else
        clean_OAP(w) = FA_50_I(w);
    end
end

global VelocityPlot
VelocityPlot = figure;

x = x_coord_array';
y = y_coord_array';
z = clean_OAP';
Xcolv = x; % Make X a column vector
Ycolv = y; % Make Y a column vector
Zcolv = z; % Make Z a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term
size(Xcolv);
size(Ycolv);
size(Zcolv);
Coefficients = [Xcolv Ycolv Const]\Zcolv; %Find the coefficients
XCoeff = Coefficients(1); % X coefficient
YCoeff = Coefficients(2); % X coefficient
CCoeff = Coefficients(3); % constant term
%Using the above variables, z = XCoeff * x + YCoeff * y + CCoeff
L=plot3(x,y,z,'ro'); % Plot the original data points
set(L,'Markerfacecolor','k') % Filling in the markers
hold on
[xx, yy]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting

zz = XCoeff * xx + YCoeff * yy + CCoeff;
if size(zz,1)>1 & size(zz,2)>1 & size(zz,1)*size(zz,2)>24
    map = colormap(hsv(256));
    surf(xx,yy,zz) % Plotting the surface
    hold on
    
    hbar = colorbar;
    %mycmap = get(colorbar,'Colormap');
    %set(colorbar,'Colormap',flipud(mycmap));
    initpos = get(hbar,'Position');
    initfontsize = get(hbar,'FontSize');
    set(hbar,'Position',[initpos(1)*1.13 initpos(2)*3.3 initpos(3)*0.5 initpos(4)*0.4],...
        'FontSize',initfontsize*0.75);
    
    %2D-projection
    x_2d = [x_init, x_init, x_final, x_final, x_init];
    y_2d = [y_init, y_final, y_final, y_init, y_init];
    z_2d = [0 0 0 0 0];
    plot3(x_2d, y_2d, z_2d)
    
    %Plane equation: zz = XCoeff * xx + YCoeff * yy + CCoeff;
    [FX,FY] = gradient(zz);
    contour(xx,yy,zz);
    hold on
    
    [xxm, yym]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting
    
    xxs = size(xx);
    xxs = xxs(1)*xxs(2);
    xx_mid = round(xxs/2);
    zzm = zz;
    
    for i = 1:xxs
        xxm(i) = xx(xx_mid);
        yym(i) = yy(xx_mid);
        zzm(i) = zz(xx_mid);
    end
    
    
    
%     quiver(xx,yy,FX,FY);
    quiver3(xx, yy, zz, FX,FY,zz/100);
    
    %quiver3(xxm, yym, zzm, FX,FY,zzm/10);
    %hold on
    %quiver(xxm,yym,FX,FY);
    
    % Fc = get(handles.edit24,'String');
    % Fc = str2num(Fc);
    % Tt = get(handles.edit26,'String');
    % Tt = str2num(Tt);
    Fc=total_frames;
    Tt=total_time;
    Conversion_factor = Fc/Tt;
    
    a = abs(FX(1)); %width converted to mm (assuming 13pixels in 1mm)
    b = abs(FY(1)); %height in frames
    %b = b; %height converted to ms
%     xx
%     yy
%     zz
    a_sq = (a)^2;
    b_sq = (b)^2;
    %Determine the magnitude of the conduction velocity
    
    Times = abs([zz(1,:),zz(2:end,end)',zz(end,1:end-1),zz(2:end-1,1)']-zz(ceil(size(zz,1)./2),ceil(size(zz,2)./2)));
    Distancex = [xx(1,:),xx(2:end,end)',xx(end,1:end-1),xx(2:end-1,1)']-xx(ceil(size(xx,1)./2),ceil(size(xx,2)./2));
    Distancey = [yy(1,:),yy(2:end,end)',yy(end,1:end-1),yy(2:end-1,1)']-yy(ceil(size(yy,1)./2),ceil(size(yy,2)./2));
    locsMax = find(Times==max(Times));
    maxY = Times(locsMax(1));
    a = Distancex(locsMax(1));
    b = Distancey(locsMax(1));
    x_of_max_plane_slope = sqrt(a^2+b^2);
    % x_of_max_plane_slope = sqrt(a_sq+b_sq);
    % y1_plane = abs(zz(1) - zz(size(zz,2)+1));
    % y2_plane = abs(zz(1) - zz(2));
    % y3_plane = abs(zz(1) - zz(size(zz,2)+2));
    %
    %
    % maxY = max(max(y1_plane,y2_plane),y3_plane);
    % % if y1_plane>=y2_plane
    % %     maxY = y1_plane;
    % % else
    % %     maxY = y2_plane;
    % % end
    % %
    % % if y3_plane>=maxY
    % %     maxY = y3_plane;
    % % end
    % %
    %
    % if maxY == 0
    %     maxY = abs(zz(size(zz,1)*size(zz,2))-zz(1));
    % end
    
    max_plane_slope = abs(maxY)/(x_of_max_plane_slope*sresolution); %sresolution mm/pixel
    Conduction_Velocity = 100*(1/max_plane_slope)
    set(handles.text63,'String',[num2str(Conduction_Velocity),'cm/s']);
    
    mTextBox = uicontrol('style','text');
    s1 = 'Conduction velocity = ';
    s2 = num2str(Conduction_Velocity);
    s3 = ' mm/ms.';
    s = [s1 s2 s3];
    set(mTextBox,'String',s,'Units','characters')
    set(mTextBox,'Position',[1 1 15 5])
    
    mTextBox = uicontrol('style','text');
    s1 = 'Average Upstroke Velocity = ';
    s2 = num2str(nanmean(UV));
    s = [s1 s2];
    set(mTextBox,'String',s,'Units','characters')
    set(mTextBox,'Position',[20 1 15 5])
    
    grid on
    
else
    errordlg('Try another pixel','Error: no value found in range');
end

function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double

function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double

function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double

function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double

function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%******************* [Conduction Velocity]->[SAVE] ***********************%
function pushbutton19_Callback(hObject, eventdata, handles)
clc
global VelocityPlot
saveimg_path = get(handles.edit20,'String');
saveas(VelocityPlot,saveimg_path);
%*************************************************************************%


function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double

function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double

function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double

function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double

function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double

function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

%% pre-processing-ordinary
% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
Inverted=get(handles.checkbox18,'Value');
Correct = get(handles.checkbox9,'Value');
if series==1
    clc
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = fullfile(PName,FName)
    B = strcat(PathName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames)
    %     A = strcat(PathName, num2str(frames), '.tif');
    A = strcat(PathName, num2str(frames), '.tif');
    h = imread(A);
    startp = get(handles.edit35,'String');
    startp = str2double(startp);
    stopp = get(handles.edit36,'String');
    stopp = str2double(stopp);
    %% map specifications
    figure;map = colormap(gray(256));close;
    %% making data matrix
    for k = startp:stopp
        i0=findtime(frames,k);
        %         A = strcat(PathName, i0, num2str(k), '.tif');
        A = strcat(PathName,  i0 , num2str(k), '.tif');
        
        h = imread(A);
        Data(:,:,k-startp+1)=h;
    end
    
elseif stack==1
    clc
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = strcat(PName,'.tif');
    startp = get(handles.edit35,'String');
    startp = str2double(startp);
    stopp = get(handles.edit36,'String');
    stopp = str2double(stopp);
    info = imfinfo(PathName);
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
    
    % bit_depth=info.BitDepth;
    % % bit_depth=12
    % bit_depth = 2^(bit_depth-1);
    % matlab_version = version;
    % windows_version = computer;
    %% map specifications
    figure;map = colormap(gray(256));close;
    %% making data matrix
    %     for k = startp:stopp
    %         Data(:,:,k-startp+1) = imread(PathName, k);
    %     end
end
Data=mat2gray(Data);

pixel_max = nanmax(nanmax(nanmax(Data)));
pixel_min = nanmin(nanmin(nanmin(Data)));

if Inverted == 1
    Data = Data-pixel_max;
    Data = abs(Data);
    Data = Data+pixel_min;
end



PName = get(handles.edit29,'String');
FName = get(handles.edit28,'String');
% PathName = fullfile(PName,FName);

DD=[fullfile(PName,FName),'_raw','.mat']
Data_raw = Data;
% save(fullfile(PathName, strcat('_raw','.mat')),'Data_raw')
save([fullfile(PName,FName),'_raw','.mat'],'Data_raw')

clear Data_raw;
%% subtracting background
% background
Nbgarray = get(handles.edit56,'String');
Nbgarray = str2num(Nbgarray)
for cbg = 1:length(Nbgarray) %count background%
    Nbg= Nbgarray(cbg);
    bg = Data(:,:,Nbg-startp+1);
    background = repmat(bg, [1 1 size(Data,3)]);
    Inverted=get(handles.checkbox18,'Value')
    if Inverted
        Data_fluor = abs(background-Data);
    else
        Data_fluor = abs(Data-background);
    end
    Data=Data_fluor;
end
clear Data_fluor;
Data = mat2gray(Data);



% % for MIP of bg frames::
%     % background
% Nbgarray = get(handles.edit56,'String');
% Nbgarray = str2num(Nbgarray)
% for cbg = 1:length(Nbgarray) %count background%
%     Nbg= Nbgarray(cbg)
%     bg(:,:,cbg) = Data(:,:,Nbg-startp+1);
%     bgMIP = nanmax(bg,[],3);
% end
%     background = repmat(bgMIP, [1 1 size(Data,3)]);
%     Inverted=get(handles.checkbox18,'Value')
%     if Inverted
%         Data_fluor = abs(background-Data);
%     else
%         Data_fluor = abs(Data-background);
%     end
%     Data=Data_fluor;
%     clear Data_fluor;
%     cbg=cbg
% % end

%% Correct baseline

if Correct==1
    AVESignal(1,1:num_images-1) = nanmean(nanmean(Data(:,:,1:end)));
    
    %%%%%%Linear trend: wts = [repmat(1/110,100,1)];
    wts = [repmat(1/100,100,1)];
    AVES = conv(AVESignal,wts,'valid');
    BS = str2num(get(handles.edit73,'String'));
    AVES = AVES-BS;
    LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
    AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
    % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
    
    for j=1:LAVE
        Data(:,:,j) = Data (:,:,j)-AVES(j);
    end
    %     stoppp = LAVE;
end


%% binning
N=3;
% avePattern = ones(N,N);
[x,y]=meshgrid(-N:N,-N:N);
avePattern=exp(-(x.^2/(2*(0.7^2))+y.^2/(2*(0.7^2))));
N=size(avePattern,1);
for k=1:10
    parfor i = 1:size(Data,3)
        temp = Data(:,:,i);
        temp = 1/N/N*conv2(temp,avePattern,'same');
        Data_binned(:,:,i) = temp;
    end
    Data= Data_binned;
    Data = mat2gray(Data);
end
clear Data_binned;
%% spatial smoothing (Moving average)
mask = ones(2,2);
parfor i = 1:size(Data,3)
    C = conv2(Data(:,:,i),mask,'same');
    Data_MA(:,:,i) = C;
end
Data = Data_MA;
clear Data_MA;
Data = mat2gray(Data);
if series ==1
    PathName = fullfile(PName,'results1');
    figure;
    for i = 1:size(Data,3)
        temp = Data(:,:,i);
        %     temp = mat2gray(temp);
        i0=findtime(frames,i+startp-1);
        imshow(temp,'Colormap',map);
        title(['" time = ', strcat( i0, num2str(i+startp-1)) , ' ms "'])
        M(i,1) = getframe(gcf);
        m=M(i,1).cdata;
        A = fullfile(PathName, strcat( i0, num2str(i+startp-1), '.tif'));
        %     A = strcat(PathName, i0, num2str(i+startp-1), '.tif');
        imwrite(temp,A)
    end
    Data_ordinary = Data;
    save([fullfile(PName,FName),'_ordinary','.mat'],'Data_ordinary')
    clear Data_ordinary;
    save([fullfile(PathName,'results',FName),'_ordinary','.mat'],'Data')
elseif stack==1
    figure;
    for i = 1:size(Data,3)
        temp = Data(:,:,i);
        %     temp = mat2gray(temp);
        i0=findtime(frames,i+startp-1);
        imshow(temp,'Colormap',map); drawnow
        title(['" time = ', strcat( i0, num2str(i+startp-1)) , ' ms "'])
        %         M(i,1) = getframe(gcf);
        %         m=M(i,1).cdata;
        %         A = fullfile(PathName, strcat( i0, num2str(i+startp-1), '.tif'));
        %         %     A = strcat(PathName, i0, num2str(i+startp-1), '.tif');
        %         imwrite(temp,A)
    end
    if size(Data,3)>2000
        save(strcat(FullPath,'.mat'), 'Data','RawData','-v7.3');
    else
        save(strcat(FullPath,'.mat'), 'Data','RawData');
    end
end


function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double

% Nbg1 = get(hObject,'String');


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');

if series
    [FileName,PathName,FilterIndex] = uigetfile('*.tif','Select the first file from tif series');
    A = strcat(PathName, FileName);
    set(handles.edit29,'String',PathName);
    set(handles.edit28,'String',FileName(1:end-4));
elseif stack
    [FileName,PathName,FilterIndex] = uigetfile('*.tif','Select the tif stack file');
    set(handles.edit29,'String',fullfile(PathName,FileName(1:end-4)));
    set(handles.edit28,'String',FileName(1:end-4));
    info = imfinfo(fullfile(PathName,FileName));
    frames = numel(info);
    bit_depth=info.BitDepth;
    set(handles.edit74,'String',num2str(bit_depth));
    set(handles.edit46,'String',num2str(frames));
    set(handles.edit36,'String',num2str(frames));
    set(handles.edit40,'String',num2str(frames));
    if frames==2000
            set(handles.edit47,'String',num2str(2099));
    end

    FullPath = fullfile(PathName,FileName);
    startp = get(handles.edit35,'String');
    startp = str2double(startp);
    stopp = get(handles.edit36,'String');
    stopp = str2double(stopp);
    %%
    % For upright signal :
    %%
    Loading_progess = waitbar(0,'Loading data....');
    
    for k = startp:stopp
        waitbar(k/(stopp-startp+1))
        RawData(:,:,k) = imread(strcat(FullPath), k);
    end
    RawData(:,:,1) = RawData(:,:,2);
    close (Loading_progess)
    ScalingOpt = 0;
    if size(RawData,1)>200
        scaling = questdlg(strcat('The image size is : ',...
            num2str(size(RawData,1)),'X',num2str(size(RawData,2)),...
            ' . Would you like to scale image to 1/2?'), ...
            'Scaling', ...
            'Yes','No','Yes');
        switch scaling
            case 'Yes'
                set(handles.checkbox30,'Value',1);
                RawData=imresize(RawData,0.5);
                ScalingOpt = 1;
            case 'No'
                set(handles.checkbox30,'Value',0);
                ScalingOpt = 0;
        end
    end
    
    cropping = questdlg(strcat('The image size is : ',...
            num2str(size(RawData,1)),'X',num2str(size(RawData,2)),...
            ' . Would you like to crop the image?'), ...
            'Crop', ...
            'Yes','No','Yes');
        switch cropping
            case 'Yes'
                figure;
                [J,rect] = imcrop(squeeze(mat2gray(RawData(:,:,5))))
                CroppedData = RawData(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);
                RawData = CroppedData;
                clear CroppedData
                croppingOpt = 1;
            case 'No'
                croppingOpt = 0;
        end
        
        
        Data = RawData;
        Data(:,:,1) = Data(:,:,2);
        Mask = ones(size(Data));
        %     Sig = squeeze(nanmean(nanmean(Data)));
        %     figure;plot(Sig); title('Step#1')
        %% Step#2 NORMALIZE DATA
        Data=mat2gray(Data);
        if size(Data,3)>2000
            save(strcat(FullPath(1:end-4),'.mat'),'Data', 'RawData','Mask','ScalingOpt','-v7.3');
        else 
            save(strcat(FullPath(1:end-4),'.mat'),'Data', 'RawData','Mask','ScalingOpt');
        end
        set(handles.edit40,'String',num2str(frames));
end
ArraySize = size(Data)
fprintf('Loading Complete\n');


function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit30.
function edit30_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% processing-smooth
% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
PName = get(handles.edit29,'String');
FName = get(handles.edit28,'String');
PathName = [fullfile(PName,'results1','results')];
load ([fullfile(PathName,FName),'_ordinary','.mat'])
PathName = [fullfile(PName,'results1')];
% B = strcat(PathName,'*.tif');
% list_of_frames = dir(B);
startp = get(handles.edit35,'String');
startp = str2double(startp);
stopp = get(handles.edit36,'String');
stopp = str2double(stopp);
frames = get(handles.edit40,'String');
frames = str2double(frames);
Data = Data(:,:,1:stopp-startp);
%% map specifications
figure;map = colormap(gray(256));close;
%% making data matrix
% for k = startp:stopp
%     i0=findtime(frames,k);
%     A = strcat(PathName, i0, num2str(k), '.tif');
%     h = imread(A);
%     Data(:,:,k-startp+1)=h;
% end
% Data = mat2gray(Data);
%% temporal drift removal
%% AND
%% temporal smoothing
PathName = [fullfile(PName,'results2')];
wait_h = waitbar(0,'Full processing...');
tstart=tic;%start timer
for i=1:size(Data,1)
    i
    parfor j=1:size(Data,2)
        Signal = squeeze(Data(i,j,:))';
        S2 = smooth(1:size(Data,3),Signal,0.9,'rloess');
        Signal = Signal-S2'; % drift removal
        window = 5;
        mask = ones(1,window)/window;
        maY=conv(Signal,mask,'same'); % temporal smoothing
        Data_smooth(i,j,:)=maY(1,:);
    end
    save([fullfile(PathName,'results',FName),'_smooth','.mat'])
    waitbar(i/size(Data,1))
end
close(wait_h)
elapsedtime=toc(tstart);fprintf(['Elapsed time: ',num2str(elapsedtime),' s\n']);

Selected_data = Data_smooth(:,:,str2double(get(handles.edit76,'String')):str2double(get(handles.edit77,'String')));
Data_smooth = Data_smooth./nanmax(nanmax(nanmax(Selected_data)));

save([fullfile(PName,FName),'_smooth','.mat'],'Data_smooth')
Data = Data_smooth;
clear Data_smooth Selected_data;
Data2 = mat2gray(Data);
% map specifications
figure;map = colormap(jet(256));close;
%
% P='C:\Users\Bahar\Documents\bahar\data\frailty project\frail group4 (CNP)\4-12-2014\2\SAN-CNP';
% name = 'SAN-CNP';
v = VideoWriter([fullfile(PathName,'results',FName),'_smooth','.avi'])
v.FrameRate = 10;
figure;
open(v)
for i = 1:size(Data2,3)
    temp = Data2(:,:,i);
    i0=findtime(frames,i+startp-1);
    imshow(temp,'Colormap',map);
    title(['" time = ', strcat( i0, num2str(i+startp-1)) , ' ms "'])
    M(i,1) = getframe(gcf);
    m=M(i,1).cdata;
    %     A = strcat(PathName, i0, num2str(i+startp-1), '.tif');
    A = fullfile(PathName, strcat( i0, num2str(i+startp-1), '.tif'));
    writeVideo(v,M(i,1))
end
close(v)
save([fullfile(PathName,'results',FName),'_smooth','.mat'])
% movie2avi(M,[PathName,'results\',FName,'_smooth','.avi'], 'compression', 'None','fps',10);

% map specifications
figure; map = colormap(gray(256));close;

Selected_data = Data(:,:,str2double(get(handles.edit76,'String')):str2double(get(handles.edit77,'String')));
Data = Data./nanmax(nanmax(nanmax(Selected_data)));

% Data = Data / nanmax(nanmax(nanmax(Data)));
figure; map = colormap(jet(256));close;
PathName = fullfile(PName,'results2','corrected');
v2 = VideoWriter([fullfile(PathName,'results',FName),'_improved','.avi'])
v2.FrameRate = 10;
figure;
open(v2)
for i = 1:size(Data,3)
    temp = Data(:,:,i);
    i0=findtime(frames,i+startp-1);
    imshow(temp,'Colormap',map);
    title(['" time = ', strcat( i0, num2str(i+startp-1)) , ' ms "'])
    M(i,1) = getframe(gcf);
    m=M(i,1).cdata;
    %     A = strcat(PathName, 'jet', i0, num2str(i+startp-1), '.tif');
    A = fullfile(PathName, strcat('jet',i0, num2str(i+startp-1), '.tif'));
    imwrite(m,map,A)
    A = fullfile(PathName, strcat( i0, num2str(i+startp-1), '.tif'));
    imwrite(temp,A)
    writeVideo(v2,M(i,1))
end
close(v2)
Data_improved = Data;
save([fullfile(PName,FName),'_improved','.mat'],'Data_improved')
clear Data_improved;
save([fullfile(PathName,'results',FName),'_improved','.mat'],'Data')
% movie2avi(M,[PathName,'results\',FName,'_improved','.avi'], 'compression', 'None','fps',10);
load ([fullfile(PName,FName),'_raw','.mat'])
load ([fullfile(PName,FName),'_ordinary','.mat'])
load ([fullfile(PName,FName),'_smooth','.mat'])
load ([fullfile(PName,FName),'_improved','.mat'])
save([fullfile(PName,FName),'_processed','.mat'])

function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [FileName,PathName,FilterIndex] = uigetfile;
% A = strcat(PathName, FileName);
% set(handles.edit39,'String',PathName);


function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit53_Callback(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit53 as text
%        str2double(get(hObject,'String')) returns contents of edit53 as a double


% --- Executes during object creation, after setting all properties.
function edit53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit55_Callback(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit55 as text
%        str2double(get(hObject,'String')) returns contents of edit55 as a double


% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox2,'Value')==1
    set(handles.checkbox23,'Value',0)
elseif get(handles.checkbox2,'Value')==0
    set(handles.checkbox23,'Value',1)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=1000*total_frames/total_time;
tresolution=num2str(tresolution);
set(handles.text53,'String',[tresolution, ' fps'])
total_pixels = get(handles.edit52,'String');
total_pixels = str2double(total_pixels);
if get(handles.checkbox30,'Value')==1
    total_pixels = total_pixels/2;
end
total_distance = get(handles.edit53,'String');
total_distance = str2double(total_distance);
sresolution=total_pixels/total_distance;
sresolution=num2str(sresolution);
set(handles.text54,'String',[sresolution,' ppmm'])


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end

AVERAGING = get(handles.checkbox2,'Value');
Data_OAP = Data ;
Data_maximum=nanmax(Data,[],3);
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
total_pixels = get(handles.edit52,'String');
total_pixels = str2double(total_pixels);
if get(handles.checkbox30,'Value')==1
    total_pixels = total_pixels/2;
end
total_distance = get(handles.edit53,'String');
total_distance = str2double(total_distance);
sresolution=total_distance/total_pixels;

Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)

xlim = [1,size(Data,3)];
% FirstFrame = get(handles.edit57,'String');
% FirstFrame = str2double(FirstFrame);
FirstFrame = 50;
% [x]= signalprep(PathName,FileName,'manual',[54,55], xlim);
measuringrange=5;
% measuringrange=floor(str2num(get(handles.edit85,'String'))/2)
% [x]= signalprep(PathName,FileName,'manual',[54,55], xlim);

%Upstroke Velocity Vector
global UV;
UV = [];

% OAP_chkbx_status = get(handles.checkbox1,'Value');
OAP_chkbx_status =0;

FAMI = [];
FA_50_I = [];
x_coord_array = [];
y_coord_array = [];

%Record Optical Action Potentials
pt = 1;
% array_area = get(handles.edit14,'String');
% array_area = str2num(array_area);
% array_area = 3;
array_area = floor(str2num(get(handles.edit85,'String'))/2);
% FileName = get(handles.edit17,'String');
% PathName = get(handles.edit16,'String');


% NumFrame = '0073';
% A = strcat(PathName, NumFrame, '.tif');
% h = imread(A);
ImScaleF = 2;
% figure
% [x_coord,y_coord,intensity_val] = impixel(h);
tempframe = squeeze(RawData(:,:,FirstFrame));
tempframe = uint16(double(tempframe).*Mask);
tempframe = imresize(tempframe,ImScaleF);
c1=figure;
figure(c1);

[x_coord,y_coord,intensity_val] = impixel(tempframe);
hold on; rectangle('Position',[x_coord-ImScaleF*array_area y_coord-ImScaleF*array_area ImScaleF*array_area ImScaleF*array_area],'EdgeColor','r')
x_coord = floor(x_coord / ImScaleF);
y_coord = floor(y_coord / ImScaleF);
x_coord1=x_coord;y_coord1=y_coord;
x= Data (y_coord,x_coord,:);
Signal(1,:) = 100*x/nanmax(x);
S2 = smooth(1:size(x,3),Signal,0.9,'rloess');
Signal = Signal-S2'; % drift removal
window = 5;
mask = ones(1,window)/window;
maY=conv(Signal,mask,'same'); % temporal smoothing
x=maY(1,:);
OAP_signal=x-nanmin(x);
% %
I=xlim(1):xlim(2);
t=I*tresolution;
[Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, x, 'data');
maximums=tMax
difmax=diff(maximums)
CL=nanmean(diff(tMax));
OAP = OAP_signal(1:endp-startp);
tt  = t(startp:endp)-t(startp)+1;
II  = I(startp:endp)-I(startp)+1;

[Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
maximums=tMax;

CL=nanmean(diff(tMax))
CL=CL/tresolution;
if Npeaks>2 && AVERAGING
    Y=[];
    cycles = 0;
    CL=ceil(CL);
    for l=1:length(IMax)
        if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP))
            Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
            cycles = cycles +1;
        end
    end
    
    OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
    tt_new=tt(1:length(OAP_new));
    II_new=II(1:length(OAP_new));
else
    OAP_new = OAP;
    tt_new=tt;
    II_new=II;
end



% figure; plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_V =zeros(size(Data,1),size(Data,2));
wait_h = waitbar(0,'Calculating conduction velocity...');
for k1=y_coord1-measuringrange:y_coord1+measuringrange
    for k2=x_coord1-measuringrange:x_coord1+measuringrange
        x_coord=k2; y_coord=k1;
        %         x= Data (y_coord,x_coord,:);
        if Data_maximum(y_coord,x_coord)>0.1
            
            %         % [x]= signalprep(PathName,FileName,'auto',[x_coord,y_coord], xlim);
            %         % % % smooth x data
            %         Signal(1,:) = x;
            %         S2 = smooth(1:size(x,3),Signal,0.9,'rloess');
            %         Signal = Signal-S2'; % drift removal
            %         window = 5;
            %         mask = ones(1,window)/window;
            %         maY=conv(Signal,mask,'same'); % temporal smoothing
            %         x=maY(1,:);
            %         % % %
            %         t=xlim(1):xlim(2);
            %         t=t*tresolution;
            %         [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, x, 'data');
            %         % maximums=tMax;
            %         % checkresult_d(t, x, Min, tMin, Max, tMax)
            %         % close all
            %
            %         CL=nanmean(diff(tMax));
            %         CL=CL/tresolution;
            %         N=2;
            %         Y=[];
            %         CL=ceil(CL)
            %         for l=1:N
            %             Y=[Y;x((l-1)*CL+1:size(x,2)-(N-l)*CL)];
            %         end
            %
            %         % Y1=x(1:end-4*CL);
            %         % Y2=x(CL+1:end-3*CL);
            %         % Y3=x(2*CL+1:end-2*CL);
            %         % Y4=x(3*CL+1:end-CL);
            %         % Y5=x(4*CL+1:end);
            %         % Y=[Y1;Y2;Y3;Y4;Y5];
            %         % figure; plot(Y')
            %         OAP_signal = nanmean(Y)-nanmin(nanmean(Y));
            %
            %         % OAP_signal=nanmean(Y)-nanmin(nanmean(Y));
            %         % figure; plot(t(1:length(OAP_signal)),Y(1,:)); hold on; plot(t(1:length(OAP_signal)),OAP_signal,'k','linewidth',2); legend('original signal','averaged signal')
            %         % xlabel('time (ms)'); ylabel(' amplitude' ) ; title('OAP')
            %         % hold on; plot([0 t(length(OAP_signal))],[0 0],'k--')
            
            
            
            %Construct the 5x5 pixel array
            x_init = x_coord - array_area;
            y_init = y_coord - array_area;
            
            x_final = x_coord + array_area;
            y_final = y_coord + array_area;
            % %
            x_curr = x_init;
            y_curr = y_init;
            % %
            % %
            % % % wait_h = waitbar(0,'Calculating conduction velocity in selected region...');
            % %
            idx1 = 1;
            idx2 = ((array_area*2)+1)^2;
            CLs=[];
            SIGNAL=[];
            slope=[];
            
            x_coord_array = [];
            y_coord_array = [];
            FA_50_I=[];
            %             figure;
            for ypos = y_curr:y_final
                for xpos = x_init:x_final
                    if ~isnan(Mask(ypos,xpos))
                        current_x = xpos;
                        current_y = ypos;
                        x_coord_array = [x_coord_array,xpos];
                        y_coord_array = [y_coord_array,ypos];
                        current_pixel = [current_x,current_y];
                        % [x]= signalprep(PathName,FileName,'auto',[current_x,current_y], xlim);
                        x=ones(1,size(Data,3));
                        x(1,:)= Data (current_y,current_x,:);
                        
                        
                        
                        OAP_signal=x-nanmin(x);
                        % %
                        I=xlim(1):xlim(2);
                        t=I*tresolution;
                        [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, x, 'data');
                        maximums=tMax;
                        CL=nanmean(diff(tMax));
                        
                        OAP = OAP_signal(1:endp-startp);
                        tt  = t(startp:endp)-t(startp)+1;
                        II  = I(startp:endp)-I(startp)+1;
                        
                        [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
                        maximums=tMax;
                        
                        CL=nanmean(diff(tMax));
                        CL=CL/tresolution;
                        % if (Npeaks>2) && strcmp(type1, 'Averaging')
                        
                        
                        if Npeaks>2 && (Npeaks<20) && AVERAGING
                            Y=[];
                            cycles = 0;
                            CL=ceil(CL);
                            for l=1:length(IMax)
                                if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP))
                                    Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
                                    cycles = cycles +1;
                                end
                            end
                            
                            OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
                            tt_new=tt(1:length(OAP_new));
                            II_new=II(1:length(OAP_new));
                        else
                            OAP_new = OAP;
                            tt_new=tt;
                            II_new=II;
                        end
                        
                        OAP_new=OAP_new-nanmin(OAP_new);
                        %                     D=diff(OAP_new);
                        %                     S=sign(D);
                        %                     SD=diff(S);
                        %                     SDS=sign(SD);
                        %                     SDS=[0,SDS,0];
                        %                     [SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(tt_new, SDS, 'model');
                        [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt_new, OAP_new, 'data');
                        %                     A1=[];B1=[];B2=[];
                        %                     for k3=1:Npeaks
                        %                         a1=find(SIMin==IMax(k3));
                        %                         b1=find(SIMax<=IMax(k3));
                        %                         if numel(b1)~=0
                        %                             b1=b1(end);
                        %                         end
                        %                         b2=find(SIMax>=IMax(k3));
                        %                         if numel(b2~=0)
                        %                             b2=b2(1);
                        %                         end
                        %                         A1=[A1,a1];B1=[B1,b1];B2=[B2,b2];
                        %                     end
                        %                     a1=A1;b1=B1;b2=B2;
                        %                     APA   = OAP_new(SIMin(a1));
                        %                     t_APA = tt_new(SIMin(a1));
                        %                     I_APA = SIMin(a1);
                        %                     RMP   = OAP_new(SIMax(b2));
                        %                     t_RMP = tt_new(SIMax(b2));
                        %                     I_RMP = SIMax(b2);
                        %                     TP = OAP_new(SIMax(b1));
                        %                     t_TP = tt_new(SIMax(b1));
                        %                     I_TP = SIMax(b1);
                        %                     APDT=50;        % choose 50 for APD50, 70 for APD70, and 90 for APD90
                        %                     if numel(TP)~=0
                        %                     AP50= ((100-APDT)/100)*(APA(1)-TP(1));
                        %                     AP50= nanmean(AP50);
                        %                     C=AP50*ones(size(OAP_new));
                        %                     D=diff(sign(OAP_new-C));
                        %
                        %                     I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
                        %                     I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
                        %
                        %                     t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
                        %                     OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
                        RT = tMax(1);
                        rt=RT(1);
                        
                        
                        
                        FA_50_I = [FA_50_I,rt]; % rising times
                        CLs=[CLs,CL];
                        %         SIGNAL=[SIGNAL; OAP_new];
                        %Calculating the slope of discrete points
                        slope = Max(1)/rt;
                        
                        UV = [UV; slope];
                        
                        %         waitbar(idx1/idx2);
                        idx1 = idx1 + 1;
                        %                     else FA_50_I = 0;
                        %                     end
                    end
                end
                
            end
            if FA_50_I~=0
                % % % close(wait_h)
                % % % wait_h2 = waitbar(1,'Conduction velocity calculation COMPLETE.');
                % % % hold off
                % %
                % %
                %Eliminate outliers in optical action potential
                mean_OAP = nanmean(FA_50_I);
                q = 1;
                clean_OAP=[];
                for w = 1:length(FA_50_I)
                    if (FA_50_I(w) < mean_OAP-5)
                        clean_OAP(w) = mean_OAP;
                    else
                        clean_OAP(w) = FA_50_I(w);
                    end
                end
                
                %                 global VelocityPlot
                %                 VelocityPlot = figure;
                
                x = x_coord_array';
                y = y_coord_array';
                z = clean_OAP';
                
                Xcolv = x; % Make X a column vector
                Ycolv = y; % Make Y a column vector
                Zcolv = z; % Make Z a column vector
                Const = ones(size(Xcolv)); % Vector of ones for constant term
                
                Coefficients = [Xcolv Ycolv Const]\Zcolv; %Find the coefficients
                XCoeff = Coefficients(1); % X coefficient
                YCoeff = Coefficients(2); % X coefficient
                CCoeff = Coefficients(3); % constant term
                %Using the above variables, z = XCoeff * x + YCoeff * y + CCoeff
                
                % L=plot3(x,y,z,'ro'); % Plot the original data points
                % set(L,'Markerfacecolor','k') % Filling in the markers
                % hold on
                [xx, yy]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting
                
                zz = XCoeff * xx + YCoeff * yy + CCoeff;
                if size(zz,1)>1 & size(zz,2)>1 & size(zz,1)*size(zz,2)>24
                    % map = colormap(hsv(256));
                    % surf(xx,yy,zz) % Plotting the surface
                    % hold on
                    
                    % hbar = colorbar;
                    %mycmap = get(colorbar,'Colormap');
                    %set(colorbar,'Colormap',flipud(mycmap));
                    % initpos = get(hbar,'Position');
                    % initfontsize = get(hbar,'FontSize');
                    % set(hbar,'Position',[initpos(1)*1.13 initpos(2)*3.3 initpos(3)*0.5 initpos(4)*0.4],...
                    %     'FontSize',initfontsize*0.75);
                    
                    %2D-projection
                    x_2d = [x_init, x_init, x_final, x_final, x_init];
                    y_2d = [y_init, y_final, y_final, y_init, y_init];
                    z_2d = [0 0 0 0 0];
                    %                 plot3(x_2d, y_2d, z_2d)
                    
                    %Plane equation: zz = XCoeff * xx + YCoeff * yy + CCoeff;
                    [FX,FY] = gradient(zz);
                    % contour(xx,yy,zz);
                    % hold on
                    
                    [xxm, yym]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting
                    
                    xxs = size(xx);
                    xxs = xxs(1)*xxs(2);
                    xx_mid = round(xxs/2);
                    zzm = zz;
                    
                    for i = 1:xxs
                        xxm(i) = xx(xx_mid);
                        yym(i) = yy(xx_mid);
                        zzm(i) = zz(xx_mid);
                    end
                    
                    
                    %
                    % quiver(xx,yy,FX,FY);
                    % quiver3(xx, yy, zz, FX,FY,zz/100);
                    
                    %quiver3(xxm, yym, zzm, FX,FY,zzm/10);
                    %hold on
                    %quiver(xxm,yym,FX,FY);
                    
                    % Fc = get(handles.edit24,'String');
                    % Fc = str2num(Fc);
                    % Tt = get(handles.edit26,'String');
                    % Tt = str2num(Tt);
                    Fc=total_frames;
                    Tt=total_time;
                    Conversion_factor = Fc/Tt;
                    
                    a = abs(FX(1)); %width converted to mm (assuming 13pixels in 1mm)
                    b = abs(FY(1)); %height in frames
                    %b = b; %height converted to ms
                    
                    a_sq = (a)^2;
                    b_sq = (b)^2;
                    
                    %Determine the magnitude of the conduction velocity
                    x_of_max_plane_slope = sqrt(a_sq+b_sq);
                    y1_plane = zz(1) - zz(6);
                    y2_plane = zz(1) - zz(2);
                    y3_plane = zz(1) - zz(7);
                    
                    maxY = [];
                    
                    if y1_plane>=y2_plane
                        maxY = y1_plane;
                    else
                        maxY = y2_plane;
                    end
                    
                    if y3_plane>=maxY
                        maxY = y3_plane;
                    end
                    
                    
                    if maxY == 0
                        maxY = zz(25)-zz(1);
                    end
                    max_plane_slope = abs(maxY)/(x_of_max_plane_slope*sresolution); %sresolution mm/pixel
                    Conduction_Velocity = 1/max_plane_slope;
                    
                    C_V(k1,k2)=Conduction_Velocity;
                end
            end
            %         plot3(k1,k2,C_V(k1,k2),'o','markerfacecolor','k'); hold on; drawnow
        end
    end
    waitbar(k1/(y_coord1+measuringrange))
end
close(wait_h)
% save( [PName,'Processings\AP\',FName,'-CV.mat'] , 'C_V','FirstFrame')
% mTextBox = uicontrol('style','text');
%
% A=(1-(C_V>=0.08)).*C_V;
% [x,y]=find(isnan(A));
% A(x,y)=0;
% [X,Y] = meshgrid(1:1:128,1:1:82);
% figure;
% subplot(2,1,1);contourf(A.*100)
% colormap(hot(8)); colorbar
% subplot(2,1,2); imshow(A.*100,colormap(hot(8))); colorbar
% SAN=get(handles.checkbox3,'Value');
% RAA=get(handles.checkbox4,'Value');
% LAA=get(handles.checkbox5,'Value');
% if SAN
%     minval=0.04;
%     maxval=0.1;
% elseif RAA
%     minval=0.2;
%     maxval=0.5;
% elseif LAA
%     minval=0.2;
%     maxval=0.5;
% end
% set(handles.edit59,'String',minval)
% set(handles.edit60,'String',maxval)

minval =get(handles.edit59,'String');
minval = str2double(minval);
maxval =get(handles.edit60,'String');
maxval = str2double(maxval);

A=(1-(C_V>=maxval)).*C_V;
A=(1-(A<minval)).*A;
A=A*100;
[x,y]=find(isnan(A));
A(x,y)=0;
[x,y]=find(A>0);
kcount=0;
for k=1:length(x)
    val(k)=A(x(k),y(k));
    kcount=kcount+1;
end
if kcount==0
    errordlg('Try another pixel','Error: no value found in range');
end
CV1=nanmean(val)
set(handles.text63,'String',[num2str(CV1), ' cm/s'])

function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

set(handles.checkbox4,'Value',0)
set(handles.checkbox5,'Value',0)

minval=0.04;
maxval=0.1;
set(handles.edit59,'String',minval)
set(handles.edit60,'String',maxval)

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox3,'Value',0)
set(handles.checkbox5,'Value',0)
minval=0.15;
maxval=0.5;
set(handles.edit59,'String',minval)
set(handles.edit60,'String',maxval)
% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox3,'Value',0)
set(handles.checkbox4,'Value',0)
minval=0.15;
maxval=0.5;
set(handles.edit59,'String',minval)
set(handles.edit60,'String',maxval)
% Hint: get(hObject,'Value') returns toggle state of checkbox5



function edit59_Callback(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit59 as text
%        str2double(get(hObject,'String')) returns contents of edit59 as a double


% --- Executes during object creation, after setting all properties.
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit60_Callback(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit60 as text
%        str2double(get(hObject,'String')) returns contents of edit60 as a double


% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
set(handles.checkbox2,'Value',1);
set(handles.checkbox23,'Value',0)

AVERAGING=get(handles.checkbox2,'Value')
manual_mode=1;
% Averaging=1;
APDT=50 ;     % choose 50 for APD50, 70 for APD70, and 90 for APD90
manual_RMP=0;
RMP=6.2;
t_RMP=88.82;

if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)
Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
xlim = [1,size(Data,3)];
NumFrame = get(handles.edit62,'String');
NumFrame = str2double(NumFrame);
i0=findtime(frames,NumFrame)
% A = strcat([fullfile(PName, 'results2','corrected', i0), num2str(NumFrame), '.tif']);
% h = imread(A);
if series
    load ([fullfile(PName,'results2','corrected','results',FName),'_improved.mat'])
    h = Data(:,:,NumFrame);
    h= imresize(h,2);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
elseif stack
    h = mat2gray(squeeze(RawData(:,:,50)));
    h= imresize(h,2);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
end


% A = strcat([PName, '\results2\corrected\', NumFrame, '.tif']);
% h = imread(A);
Sig=[];
if manual_mode
    dh = figure;
    [x_coord,y_coord,intensity_val] = impixel(h);
    x_coord = floor(x_coord./2); y_coord = floor(y_coord./2);
    x= Data (y_coord,x_coord,startp:endp);
else
    x= Data (29,25,startp:endp);
end


x=squeeze(x);
frame_rate=1000*total_frames/total_time;
% cutoff_frequency = 0.5;
% filter_order = 2;
% normalized_cutoff_frequency = cutoff_frequency/(frame_rate/2);
% [b,a] = butter(filter_order,normalized_cutoff_frequency,'high');
% xCorrected = filter(b,a,x);
% xCorrected (1:50) = xCorrected (51);
% x = xCorrected;

Signal=[];
% % smooth x data
Signal(1,:) = x;
% Signal(1,:) = nanmean(Sig);
S2 = smooth(1:size(x),Signal,0.9,'rloess');
Signal = Signal-S2'; % drift removal
window = 5;
mask = ones(1,window)/window;
maY=conv(Signal,mask,'same'); % temporal smoothing
x=maY(1,:);
x= smooth(x,0.01)';
OAP_signal=x-nanmin(x);
% %
t=1:size(OAP_signal,2);
t=t*tresolution;
I=t/tresolution;
[Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, x, 'data');
maximums=tMax
difmax= diff(maximums)
% checkresult_d(t, x, Min, tMin, Max, tMax)
CL=nanmean(diff(tMax));

OAP = OAP_signal;
tt  = t;
II  = I;

[Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
maximums=tMax;
% checkresult_d(tt, OAP, Min, tMin, Max, tMax)

CL=nanmean(diff(tMax))
CL=CL/tresolution;
% if Npeaks>2 && AVERAGING
%     Y=[];
%     cycles = 0;
%     CL=ceil(CL);
%     for l=1:length(IMax)
%         if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP))
%             Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
%             cycles = cycles +1;
%         end
%     end
%     
%     OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
%     tt_new=tt(1:length(OAP_new));
%     II_new=II(1:length(OAP_new));
% else
%     OAP_new = OAP;
%     tt_new=tt;
%     II_new=II;
% end
% figure; plot(Y'); hold on; plot(nanmean(Y),'k')

%                 if Npeaks>3 && AVERAGING
%                     Y=[];
%                     cycles = 0;
%                     CL=ceil(CL);
%                     for l=1:length(IMax)-1
%                         if (IMax(l)-floor(CL./2) > 0) & (IMax(l+1)+floor(CL./2)<length(OAP))
%                             Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2)),OAP(IMax(l+1)-floor(CL./2):IMax(l+1)+floor(CL./2))];
%                             cycles = cycles +1;
%                         end
%                     end
% 
%                     %         OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
%                     OAP_new = nanmean(Y)-nanmin(nanmean(Y));
% 
%                     tt_new=tt(1:length(OAP_new));
%                     II_new=II(1:length(OAP_new));
%                 else
%                     OAP_new = OAP;
%                     tt_new=tt;
%                     II_new=II;
%                 end
    OAP = 100*(OAP_signal-min(OAP_signal))/(max(OAP_signal)-min(OAP_signal));

    [I,J]=findpeaks(OAP,II,'MinPeakProminence',50,'Annotate','extents')
    figure;findpeaks(OAP,II,'MinPeakProminence',50,'Annotate','extents')
        Starting_frame = get(handles.edit35,'String');
%     Starting_frame = str2double(Starting_frame);
%     Ending_frame = get(handles.edit36,'String');
%     Ending_frame = str2double(Ending_frame);
%     if Ending_frame>size(Data,3)
%         Ending_frame =size(Data,3)
%     end
    if length(I)>3 && AVERAGING
        Y=[];
        k=0;
        CL=diff(J)
        halfCL = floor(mean(diff(J))./2);
        for j=1:length(J)-1
            if (J(j)-halfCL > 0) && (floor(J(j+1)+halfCL))<length(OAP)
                k=k+1;
%                 OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2)),OAP(IMax(l+1)-floor(CL./2):IMax(l+1)+floor(CL./2))
% size(OAP)
% Starting_frame+floor(J(j)-halfCL)
% Starting_frame+floor(J(j)+halfCL)
% Starting_frame+floor(J(j+1)-halfCL)
% Starting_frame+floor(J(j+1)+halfCL)
floor(J(j)-halfCL)
floor(J(j)+halfCL)
floor(J(j+1)-halfCL)
floor(J(j+1)+halfCL)
size([OAP(floor(J(j))-halfCL:floor(J(j))+halfCL),OAP(floor(J(j+1))-halfCL:floor(J(j+1))+halfCL)])
size(Y)
                Y=[Y;OAP(floor(J(j))-halfCL:floor(J(j))+halfCL),OAP(floor(J(j+1))-halfCL:floor(J(j+1))+halfCL)];
                %             AData(:,:,:,k) = Data(:,:,(Starting_frame+floor(J(j)-halfCL):Starting_frame+floor(J(j)+halfCL)));
            end
        end
        OAP_new = nanmean(Y)-nanmin(nanmean(Y));
        
        tt_new=tt(1:length(OAP_new));
        II_new=II(1:length(OAP_new));
    else
        OAP_new = OAP;
        tt_new=tt;
        II_new=II;
    end
    
    
figure;plot(Y');title('individual beats')
                Signal=[];
              x=OAP_new;
              % % smooth x data
              Signal(1,:) = x;
              % Signal(1,:) = nanmean(Sig);
              S2 = smooth(1:length(x),Signal,0.9,'rloess');
              Signal = Signal-S2'; % drift removal
              window = 5;
              mask = ones(1,window)/window;
              maY=conv(Signal,mask,'same'); % temporal smoothing
              x=maY(1,:);
              x= smooth(x,0.01)';
              OAP_new=x;
OAP_new=OAP_new-nanmin(OAP_new);
OAP_new = OAP_new/nanmax(OAP_new)*100;
% D=diff(OAP_new);
% S=sign(D);
% SD=diff(S);
% SD2=diff(SD);
% SDS=sign(SD);
% SDS2=-sign(SD2);
% SDS=[0,SDS,0];
% SDS2=[0,SDS2,0,0];
% [SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(tt_new, SDS, 'model'); close;
% figure; plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');
% [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt_new, OAP_new, 'data');
% a1=find(SIMin==IMax(1));
% b1=find(SIMax<IMax(1));
% b1=b1(end-1);
% b12=find(SIMax<IMax(2));
% b12=b12(end-1);
% b2=find(SIMax>=IMax(1));
% b2=b2(1);
% APA   = OAP_new(SIMin(a1));
% t_APA = tt_new(SIMin(a1));
% I_APA = SIMin(a1);
% if manual_RMP
%     I_RMP=ceil(t_APA/tresolution);
% else
%     RMP   = OAP_new(SIMax(b2));
%     t_RMP = tt_new(SIMax(b2));
%     I_RMP = SIMax(b2);
% end
% TP = OAP_new(SIMax(b1));
% t_TP = tt_new(SIMax(b1));
% I_TP = SIMax(b1);
% TP2 = OAP_new(SIMax(b12));
% t_TP2 = tt_new(SIMax(b12));
% I_TP2 = SIMax(b12);
% RMPsearch = OAP_new(I_RMP:I_RMP+20);
% [RMP,IRMP] = min(RMPsearch);
% t_RMP = tt(I_RMP+IRMP-1);
% % hold on; plot(t_APA,APA,'rs','markerfacecolor','r');
% % hold on; plot(t_RMP,RMP,'go','markerfacecolor','g');
% % hold on; plot(t_TP,TP,'c^','markerfacecolor','c');
% % hold on; plot(t_TP2,TP2,'m^','markerfacecolor','m');
% hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
% axislimits=nanmax(OAP_new)-nanmin(OAP_new);
% axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
% % legend('Action Potential (AP)','Action Potential Amplitude (APA)','Rest Membrane Potential (RMP)','Take-off Potential (TP)' )
%% DD slope
% DDslope=(TP2-RMP)/(t_TP2-t_RMP)
% hold on; plot([t_TP2,t_RMP],[TP2,RMP],'m--');
close(dh)
% tt_new
% OAP_new
figure; plot(tt_new,OAP_new,'k','linewidth',1);xlabel('time (ms)');
[x,y] = ginput(2);
hold on;plot(x,y)
DDslope=(diff(y))/(diff(x))

text(20,65,['DD slope',' = ', num2str(DDslope)])


function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8



function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit65_Callback(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit65 as text
%        str2double(get(hObject,'String')) returns contents of edit65 as a double


% --- Executes during object creation, after setting all properties.
function edit65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
PName = get(handles.edit29,'String');
FName = get(handles.edit28,'String');
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
% if stack == 1
PathName = strcat(PName,'.tif')
info = imfinfo(PathName);
num_images = numel(info);
bit_depth=info.BitDepth;

user_bit_depth=str2num(get(handles.edit74,'String'));
if user_bit_depth ~= bit_depth
    bit_depth=user_bit_depth
end
% bit_depth=12
bit_depth = 2^(bit_depth-1)
matlab_version = version;
windows_version = computer;
Inverted=get(handles.checkbox18,'Value')
Correct = get(handles.checkbox9,'Value')

AllImage = get(handles.checkbox13,'Value')
SelectRegion = get(handles.checkbox12,'Value')

for k = 2:num_images
    A(:,:,k-1) = imread(PathName, k);
end

startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)

if (windows_version(end)=='N')
    One=1
    A = single(A)/bit_depth;
else
    Two=1
    A = double(A)/bit_depth;
end

pixel_max = nanmax(nanmax(nanmax(A)));
pixel_min = nanmin(nanmin(nanmin(A)));

if Inverted
    A=A-pixel_max;
    A=abs(A);
    A=A+pixel_min;
end

AVESignal(1,1:num_images-1) = nanmean(nanmean(A(:,:,1:end)));
A=(A-nanmin(AVESignal))./nanmax(nanmax(nanmax(A-nanmin(AVESignal))));

% A=(A-nanmin(nanmin(nanmin(A))))./nanmax(nanmax(nanmax(A-nanmin(nanmin(nanmin(A))))));
AVESignal=[];
if AllImage==1
    AVESignal(1,1:num_images-1) = nanmean(nanmean(A(:,:,1:end)));
elseif SelectRegion==1
    AShow=imread(PathName, 10);
    h = mat2gray(AShow);
    H(:,:,1)=h;
    H(:,:,2)=h;
    H(:,:,3)=h;
    figure;imshow(H);title('Select Your Region of Interest')
    [x,y,button] = ginput(2);
    x=round(x); y=round(y);
    rectangle('Position',[x(1),y(1),x(2)-x(1),y(2)-y(1)],'EdgeColor','r',...
        'LineWidth',3)
    AVESignal(1,1:num_images-1) = nanmean(nanmean(A(y(1):y(2),x(1):x(2),1:end)));
    close
end
% figure;plot(AVESignal)
%%%%%%%%%% EDIT WHEN TRYING TO CORRECT THE baseline drift%%%%%%%%%%%%%%%%%%
if Correct==1
    wts = [repmat(1/100,100,1)];
    AVES = conv(AVESignal,wts,'valid');
    % AVES = AVES-(AVESignal(1)-AVECorrected(1));
    LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
    AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
    % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
    
    Deltay = (AVESignal(1)-AVECorrected(1));
    AVES = AVES-(AVESignal(1)-AVECorrected(1));
    AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
    % figure;plot(AVESignal-Deltay,'b');hold on;plot(AVES+(AVESignal(1)-Deltay),'r');plot(AVECorrected-Deltay,'g')
    for j=1:LAVE
        A(:,:,j) = A (:,:,j)-AVES(j);
    end
    
    % %%%%%%Linear trend: wts = [repmat(1/110,100,1)];
    % wts = [repmat(1/10,100,1)];
    % AVES = conv(AVESignal,wts,'valid');
    BS = str2num(get(handles.edit73,'String'));
    % AVES = AVES-BS;
    % LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
    % AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
    % % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
    %
    % for j=1:LAVE
    %     A(:,:,j) = A (:,:,j)-AVES(j);
    % end
    %
else
    LAVE=length(AVESignal); %% Length of the drif signal and corrected signal
    AVECorrected = AVESignal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbg = get(handles.edit56,'String');
Nbg = str2double(Nbg);
bg = A(:,:,Nbg-startp+1);
B = A - bg;
B = B-nanmin(nanmin(nanmin(B)))./nanmax(nanmax(nanmax(B-nanmin(nanmin(nanmin(B))))));

screensize = get( groot, 'Screensize' )
f=figure('visible','on','Position',[screensize(3)/4 screensize(4)/10 screensize(3)/3 4*screensize(4)/5]); %[left bottom width height]
ax = axes('Units','pixels');


for k = 1:LAVE
    C=imgaussfilt(B(:,:,k),1);
    subplot(2,1,1);imshow(C) %could be the mask
    if SelectRegion
        rectangle('Position',[x(1),y(1),x(2)-x(1),y(2)-y(1)],'EdgeColor','r',...
            'LineWidth',1); hold on
    end
    if Correct==1
        subplot(2,1,2);
        %             plot([1:k],AVESignal(1:k)-Deltay,'-k');
        %             hold on;plot(AVES(1:k)+(AVESignal(1)-Deltay),'r');
        plot(AVECorrected(1:k)-Deltay,'k')
        axis([1, ...
            LAVE,...
            nanmin(AVECorrected-Deltay)-0.1*(nanmax(AVECorrected-Deltay)-nanmin(AVECorrected-Deltay)),...
            nanmax(AVECorrected-Deltay)+0.1*(nanmax(AVECorrected-Deltay)-nanmin(AVECorrected-Deltay))]);
        %            plot(AVESignal-Deltay,'b');hold on;plot(AVES+(AVESignal(1)-Deltay),'r');plot(AVECorrected-Deltay,'g')
        %             legend('rawsignal','baseline drift','corrected signal')
    else
        subplot(2,1,2);plot([1:k],AVESignal(1:k),'-k'); hold on
        axis([1 LAVE nanmin(AVESignal)-0.1*(nanmax(AVESignal)-nanmin(AVESignal)) nanmax(AVESignal)+0.1*(nanmax(AVESignal)-nanmin(AVESignal))]);
        
    end
    drawnow
end
% elseif series == 1
%
%     PathName = strcat(fullfile(PName,FName),'.tif');
%
% end
function edit69_Callback(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit69 as text
%        str2double(get(hObject,'String')) returns contents of edit69 as a double


% --- Executes during object creation, after setting all properties.
function edit69_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
PName = get(handles.edit29,'String');
FName = get(handles.edit28,'String');
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
if stack == 1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
    user_bit_depth=str2num(get(handles.edit74,'String'));
    if user_bit_depth ~= bit_depth
        bit_depth=user_bit_depth
    end
    % bit_depth=12
    bit_depth = 2^(bit_depth-1)
    matlab_version = version;
    windows_version = computer;
    Inverted=get(handles.checkbox18,'Value')
    Correct = get(handles.checkbox9,'Value')
    
    AllImage = get(handles.checkbox13,'Value')
    SelectRegion = get(handles.checkbox12,'Value')
    
    %     for k = 2:num_images
    %         A(:,:,k-1) = imread(PathName, k);
    %     end
    A=Data (:,:,1:end);
    
    startp = get(handles.edit35,'String');
    startp = str2double(startp);
    endp = get(handles.edit36,'String');
    endp = str2double(endp)
    
    if (windows_version(end)=='N')
        One=1
        A = single(A)/bit_depth;
    else
        Two=1
        A = double(A)/bit_depth;
    end
    
    pixel_max = nanmax(nanmax(nanmax(A)));
    pixel_min = nanmin(nanmin(nanmin(A)));
    
    if Inverted
        A=A-pixel_max;
        A=abs(A);
        A=A+pixel_min;
    end
    
    AVESignal(1,1:num_images) = nanmean(nanmean(A(:,:,1:end)));
    A=(A-nanmin(AVESignal))./nanmax(nanmax(nanmax(A-nanmin(AVESignal))));
    
    % A=(A-nanmin(nanmin(nanmin(A))))./nanmax(nanmax(nanmax(A-nanmin(nanmin(nanmin(A))))));
    AVESignal=[];
    if AllImage==1
        AVESignal(1,1:num_images) = nanmean(nanmean(A(:,:,1:end)));
    elseif SelectRegion==1
        %         AShow=imread(PathName, 10);
        h = mat2gray(RawData(:,:,50));
        h = imresize(h,2);
        H(:,:,1)=h;
        H(:,:,2)=h;
        H(:,:,3)=h;
        figure;imshow(H);title('Select Your Region of Interest')
        [x,y,button] = ginput(2);
        rectangle('Position',[x(1),y(1),x(2)-x(1),y(2)-y(1)],'EdgeColor','r',...
            'LineWidth',3)
        x=floor(x./2); y=floor(y./2);
        AVESignal(1,1:num_images) = nanmean(nanmean(A(y(1):y(2),x(1):x(2),1:end)));
    end
    % figure;plot(AVESignal)
    %%%%%%%%%% EDIT WHEN TRYING TO CORRECT THE baseline drift%%%%%%%%%%%%%%%%%%
    if Correct==1
        %option1
        t = (1:length(AVESignal));
        opol = 6;
        [p,s,mu] = polyfit(t,AVESignal,opol);
        f_y = polyval(p,t,[],mu);
        
        AVECorrected = AVESignal - f_y;
        
        %         AVECorrected = detrend(AVESignal);
        figure;plot(AVESignal,'k');hold on; plot(AVECorrected,'g')
        legend('rawsignal','corrected signal')
        %option2:
        %         frame_rate = 1000;
        %         %     frame_rate=frames*total_frames/total_time;
        %         %     cutoff_frequency = str2num(get(handles.edit79,'String'));
        %         cutoff_frequency = 2;
        %         %     filter_order = str2num(get(handles.edit80,'String'));
        %         filter_order = 2;
        %         normalized_cutoff_frequency = cutoff_frequency/(frame_rate/2);
        %         [b,a] = butter(filter_order,normalized_cutoff_frequency,'high');
        %         FilteredData = filter(b,a,AVESignal);
        % %         FilteredData(:,:,1:50)=repmat(FilteredData(:,:,51),[1,1,50]);
        %         size(FilteredData)
        % %         AVECorrected = squeeze(nanmean(nanmean(FilteredData)));
        %         size(FilteredData)
        %         figure;plot(AVESignal,'k');hold on; plot(FilteredData,'g')
        %         legend('rawsignal','corrected signal')
        %option1:
        %         wts = [repmat(1/100,100,1)];
        %         AVES = conv(AVESignal,wts,'valid');
        %         % AVES = AVES-(AVESignal(1)-AVECorrected(1));
        %         LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
        %         AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
        %         % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
        %
        %         Deltay = (AVESignal(1)-AVECorrected(1));
        %         AVES = AVES-(AVESignal(1)-AVECorrected(1));
        %         AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
        %         figure;plot(AVESignal-Deltay,'k');hold on;plot(AVES+(AVESignal(1)-Deltay),'r');plot(AVECorrected-Deltay,'g')
        %         % figure;plot(AVESignal,'k');hold on;plot(AVES,'r');plot(AVECorrected,'g')
        %         legend('rawsignal','baseline drift','corrected signal')
    else
        figure;plot(AVESignal,'k')
        legend('rawsignal')
    end
elseif series == 1
    PathName = strcat(fullfile(PName,FName),'.tif');
end
% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9



function edit72_Callback(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit72 as text
%        str2double(get(hObject,'String')) returns contents of edit72 as a double


% --- Executes during object creation, after setting all properties.
function edit72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit73_Callback(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit73 as text
%        str2double(get(hObject,'String')) returns contents of edit73 as a double


% --- Executes during object creation, after setting all properties.
function edit73_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox10,'Value')==1
    set(handles.checkbox11,'Value',0)
    set(handles.text68,'Visible','off')
    set(handles.edit62,'Visible','off')
elseif get(handles.checkbox10,'Value')==0
    set(handles.checkbox11,'Value',1)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox11,'Value')==1
    set(handles.checkbox10,'Value',0)
    set(handles.text68,'Visible','on')
    set(handles.edit62,'Visible','on')
elseif get(handles.checkbox11,'Value')==0
    set(handles.checkbox10,'Value',1)
    set(handles.text68,'Visible','off')
    set(handles.edit62,'Visible','off')
end
% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox12,'Value')==1
    set(handles.checkbox13,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox13,'Value')==1
    set(handles.checkbox12,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox13



function edit74_Callback(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit74 as text
%        str2double(get(hObject,'String')) returns contents of edit74 as a double


% --- Executes during object creation, after setting all properties.
function edit74_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
PName = get(handles.edit29,'String');
FName = get(handles.edit28,'String');
PathName = [fullfile(PName,'results1','results')];
load ([fullfile(PathName,FName),'_ordinary','.mat'])
Selected_data = Data(:,:,str2double(get(handles.edit76,'String')):str2double(get(handles.edit77,'String')));
Data = Data./nanmax(nanmax(nanmax(Selected_data)));
Data_ordinary = Data;
%% (check if it is necessary)
save([fullfile(PName,FName),'_ordinary','.mat'],'Data_ordinary')
clear Data_ordinary;
save([fullfile(PathName,FName),'_ordinary','.mat'],'Data')
%% %% %% %% %% %% %% %% %% %% %% %%
v = VideoWriter([fullfile(PathName,FName),'_Selected','.avi'])
v.FrameRate = 10;
figure;
open(v)
screensize = get( groot, 'Screensize' )
f=figure('visible','on','Position',[screensize(3)/4 screensize(4)/10 screensize(3)/3 2*screensize(4)/5]); %[left bottom width height]
ax = axes('Units','pixels');
map = colormap(gray(256));close;
for i = 1:size(Data,3)
    temp = Data(:,:,i);
    %     i0=findtime(frames,i+startp-1);
    imshow(temp,'Colormap',map);
    %     title(['" time = ', strcat( i0, num2str(i+startp-1)) , ' ms "'])
    M(i,1) = getframe(gcf);
    %     m=M(i,1).cdata;
    % %     A = strcat(PathName, i0, num2str(i+startp-1), '.tif');
    %      A = fullfile(PathName, strcat( i0, num2str(i+startp-1), '.tif'));
    writeVideo(v,M(i,1))
end
close(v)


function edit76_Callback(hObject, eventdata, handles)
% hObject    handle to edit76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit76 as text
%        str2double(get(hObject,'String')) returns contents of edit76 as a double


% --- Executes during object creation, after setting all properties.
function edit76_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit77_Callback(hObject, eventdata, handles)
% hObject    handle to edit77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit77 as text
%        str2double(get(hObject,'String')) returns contents of edit77 as a double


% --- Executes during object creation, after setting all properties.
function edit77_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RoiOption = get(handles.popupmenu2,'Value'); %% 1 for freehand - 2 for impoly - 3 for imellipse - 4 for imrect
if RoiOption == 1
    set(handles.slider1,'Visible','on')
    set(handles.edit84,'Visible','on')
    set(handles.text102,'Visible','on')
    
else
    set(handles.slider1,'Visible','off')
    set(handles.edit84,'Visible','off')
    set(handles.text102,'Visible','off')
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% RoiOption = get(handles.popupmenu2,'Value'); %% 1 for freehand - 2 for impoly - 3 for imellipse - 4 for imrect
clc
FullPath = get(handles.edit29,'String');
load(strcat(FullPath,'.mat'));
info = imfinfo(strcat(FullPath,'.tif'));
frames = numel(info);
num_images = numel(info);
bit_depth = info.BitDepth;
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
ScalingOpt = get(handles.checkbox30,'Value');

% %
% % %% Define a ROI
% % % if RoiOption == 1
% % %     BGFrame = mat2gray(nanmax(RawData,[],3));
% % %     th = get(handles.slider1, 'Value');
% % %     Mask = BGFrame>th;
% % %     Mask = double(Mask);
% % %     Mask( Mask == 0 ) = NaN;
% % %     MaskedData = Data.*repmat(Mask,[1,1,frames]);
% % %     Data = MaskedData;
% % %
% % % % BGFrame = mat2gray(nanmax(RawData,[],3));
% % % %     ansSatisfied='NO';
% % % %     for n=1:frames
% % % %         if strcmp(ansSatisfied,'YES')
% % % %             Data = MaskedData;
% % % %             break
% % % %         elseif strcmp(ansSatisfied,'Cancel masking. I changed my mind!')
% % % %             break
% % % %         elseif strcmp(ansSatisfied,'NO')
% % % %             mask = BGFrame>0.08;
% % % %             figure;imshow(mask); title('MASK');
% % % %             MaskedData = Data.*repmat(mask,[1,1,frames]);;
% % % %             ansSatisfied = questdlg('Do you like the mask?','MASK', 'YES', 'NO', 'Cancel masking. I changed my mind!', 'NO');
% % % %         end
% % % %     end
% % % else
% % %     ansSatisfied='NO';
% % %     for n=1:frames
% % %         if strcmp(ansSatisfied,'YES')
% % %             Data = MaskedData;
% % %             break
% % %         elseif strcmp(ansSatisfied,'Cancel masking. I changed my mind!')
% % %             break
% % %         elseif strcmp(ansSatisfied,'NO')
% % %             figure(1);imshow(imresize(Data(:,:,5),3));
% % %
% % %             if RoiOption==2 %freehand
% % %                 h = imfreehand(gca);
% % %             elseif RoiOption==3 %impoly
% % %                 h = impoly(gca);
% % %             elseif RoiOption==4 %imellipse
% % %                 h = imellipse (gca);
% % %             elseif RoiOption==5 %imrect
% % %                 h = imrect (gca);
% % %             end
% % %             Mask =createMask(h);
% % %             Mask = imresize(Mask,1/3);
% % %                 Mask = double(Mask);
% % %                 Mask( Mask == 0 ) = NaN;
% % %             MaskedData = Data .*repmat(Mask,[1,1,frames]);
% % %         end
% % %         ansSatisfied = questdlg('Do you like the mask?','MASK', 'YES', 'NO', 'Cancel masking. I changed my mind!', 'NO');
% % %     end
% % %     figure(2);imshow(imresize(Data(:,:,50),3));
% % %     save (strcat(FullPath,'.mat'),'RawData','Data','Mask');
% % % end
% % %
% Scaling_progess = waitbar(0,'Please wait....');
% for i=1:size(Data,1)
%     waitbar(i/size(Data,1));
%     for j=1:size(Data,2)
%         ScaledData(i,j,:) = (Data(i,j,:)-nanmin(squeeze(Data(i,j,:))))/(nanmax(squeeze(Data(i,j,:)))-nanmin(squeeze(Data(i,j,:))));
%     end
% end
% close (Scaling_progess);
% Data = mat2gray(ScaledData);
Data = mat2gray(Data);

LPFchb = get(handles.checkbox14,'Value');
MAchb = get(handles.checkbox15,'Value');
SMchb = get(handles.checkbox17,'Value');
Invchb = get(handles.checkbox18,'Value');
% PCAchb = get(handles.checkbox19,'Value');
preproschb = get(handles.checkbox20,'Value');
fastpreproschb = get(handles.checkbox31,'Value');

bcorrchb = get(handles.checkbox21,'Value');
if Invchb == 1
    Dmax = 1;
    Dmin = 0;
    Data_Inv = Data-Dmax;
    Data_Inv = abs(Data_Inv);
    Data_Inv = Data_Inv+Dmin;
    %     Sig = squeeze(nanmean(nanmean(Data)));
    %     Siginv = squeeze(nanmean(nanmean(Data_Inv)));
    %     figure(3);cla;plot(Sig); hold on;plot(Siginv); legend('before inverting', 'after inverting')
    Data = Data_Inv;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
%     end
    
end

if LPFchb == 1
    frame_rate=1000*total_frames/total_time;
    
    cutoff_frequency = str2num(get(handles.edit79,'String'));
    filter_order = str2num(get(handles.edit80,'String'));
    normalized_cutoff_frequency = cutoff_frequency/(frame_rate/2);
    [b,a] = butter(filter_order,normalized_cutoff_frequency,'low');
    FilteredData = filter(b,a,Data,[],3);
    %     FilteredData(:,:,1:50)=repmat(FilteredData(:,:,51),[1,1,50]);
    %     anskeep = questdlg('Do you like to keep the changes?','Lowpass filter', 'YES', 'NO', 'YES');
    %     if strcmp(anskeep, 'YES')
    FilteredData(:,:,1:cutoff_frequency)=repmat(FilteredData(:,:,cutoff_frequency+1),[1,1,cutoff_frequency]);
    Sig = squeeze(nanmean(nanmean(Data)));
    Sigf = squeeze(nanmean(nanmean(FilteredData)));
    figure(3);cla; plot(Sig); hold on;plot(Sigf); legend('before LPF', strcat('after',num2str(cutoff_frequency),'Hz LPF'));
    xlabel('frame');
    Data = FilteredData;
%     Data(:,:,1:cutoff_frequency*1.2)=repmat(Data(:,:,cutoff_frequency*1.2+1),1,1,cutoff_frequency*1.2);
    %     if size(Data,3)>2000
    %         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
    %     else
    %         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
    %     end
    %     end
end
% figure;imshow(squeeze(Data(:,:,50)));
if fastpreproschb == 1
    Sig = squeeze(nanmean(nanmean(Data,1),2));
    Sig = smooth(squeeze(Sig),0.01);
    Sig3 = smooth(1:length(Sig),squeeze(Sig),0.9,'loess');

%     Sig2 = smooth(squeeze(Sig),0.5);
%     Sig3 = filter(ones(1,101)/101,1,Sig2);
%     Sig3(1:101)=Sig3(102);
%     Sig3 = Sig3+(max(Sig)-max(Sig3));
    Sig = Sig-Sig3;
    Sig = 100*(Sig-min(Sig))/(max(Sig)-min(Sig));
     %% correct for bleaching
    
    bleachsignal (1,1,:) = Sig3;
    bleachmatrix = repmat(bleachsignal,[size(Data,1),size(Data,2),1]);
    CorrectedData = Data - bleachmatrix;
    BeforCorrectionData = Data;
    Data = mat2gray(CorrectedData);
    %%
    t=1:size(Sig,1);
%     time_ms = (1:length(Sig))./(total_frames/total_time);
%     t=t*tresolution;
    [I,J]=findpeaks(Sig,t,'MinPeakProminence',30,'Annotate','extents');
    figure;findpeaks(Sig,t,'MinPeakProminence',30,'Annotate','extents');
    hold on; plot(t(J(1)+mean(diff(J))/2-mean(diff(J))/7:J(1)+mean(diff(J))/2+mean(diff(J))/7),...
        Sig(J(1)+mean(diff(J))/2-mean(diff(J))/7:J(1)+mean(diff(J))/2+mean(diff(J))/7),'r');
    xlabel('frame');
    background2 = [];
    for cl = 1: length(J)
        bgstart = floor(J(cl)+mean(diff(J))/2-mean(diff(J))/7);
        bgend = floor(J(cl)+mean(diff(J))/2+mean(diff(J))/7);
        if bgend>length(Sig)
            bgend = length(Sig);
        end
        BG(:,:,cl)=nanmean(Data(:,:,bgstart:bgend),3);
        if cl == 1
            bgstart = 1;
        end
        if (cl < length(J))&&((floor(J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/7-1)<size(Data,3)))
%         floor(J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/5-bgstart)
%         bgstart
%         J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/5-1
        background2 (:,:,bgstart:floor(J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/7-1)) = repmat(squeeze(BG(:,:,cl)),[1 1 floor(J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/7-bgstart)]);
        elseif (cl == length(J))||(floor(J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/7-1)>size(Data,3))
%                 floor(J(cl+1)+mean(diff(J))/2-mean(diff(J+1))/5-bgstart)
%         bgstart
%         size(Data,3)
        background2 (:,:,bgstart:size(Data,3)) = repmat(squeeze(BG(:,:,cl)),[1 1 floor(size(Data,3)-bgstart+1)]);
        end

    end
%     size(background2)
%     size(Data)
    BGFrame = nanmean(BG,3);

% BGFrame = squeeze(BG(:,:,1));
%    figure;imshow(squeeze(BGFrame));
%     %% correct for bleaching
%     
%     bleachsignal (1,1,:) = Sig3;
%     bleachmatrix = repmat(bleachsignal,[size(Data,1),size(Data,2),1]);
%     CorrectedData = Data - bleachmatrix;
%     BeforCorrectionData = Data;
%     Data = mat2gray(CorrectedData);
    
    %% subtracting background
    %     BGFrame = (mat2gray(nanmean(Data,3))+mat2gray(nanmax(Data,[],3)))./2;
    %     BGFrame = mat2gray(nanmax(Data,[],3));
    background = repmat(BGFrame, [1 1 size(Data,3)]);
%     BG = nanmean(Data(:,:,2048:2118),3);
%     BGFrame = nanmean(BG,3);
%     background2 = repmat(BGFrame, [1 1 size(Data,3)]);
    Data_fluor = Data-background2;
%         Data_fluor(Data_fluor<0)=0;

    Data = mat2gray(Data_fluor);
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
%     Sig = squeeze(nanmean(nanmean(Data,1),2));
% %     Sig = smooth(squeeze(Sig),0.01);
%     Sig3 = smooth(1:length(Sig),squeeze(Sig),0.9,'loess');
%     Sig = Sig-Sig3;
%     Sig = 100*(Sig-min(Sig))/(max(Sig)-min(Sig));
%      %% correct for bleaching
%     bleachsignal (1,1,:) = Sig3;
%     bleachmatrix = repmat(bleachsignal,[size(Data,1),size(Data,2),1]);
%     CorrectedData = Data - bleachmatrix;
%     BeforCorrectionData = Data;
%     Data = mat2gray(CorrectedData);
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%

clear Data_fluor;
clear CorrectedData;
Processing_progess = waitbar(0,'Processing....');

for k=1:10
    waitbar(k/10);
    Data= imgaussfilt(Data,0.67);
    Data = mat2gray(Data);
end
close(Processing_progess)
FilteredData = filter(ones(1,2)/2,1,Data,[],3);
FilteredData(:,:,1:2) = repmat(FilteredData(:,:,2+1),[1,1,2]);
Data=FilteredData;
Data = mat2gray(Data);
clear FilteredData;

    %% apply 3D median filter
    MData = medfilt3(Data,'replicate');
    Data = mat2gray(MData);
clear MData;

end
%     figure;imshow(squeeze(Data(:,:,50)));

if preproschb == 1
    %     Correct = get(handles.checkbox9,'Value');
    startp = get(handles.edit35,'String');
    startp = str2double(startp);
    stopp = get(handles.edit36,'String');
    stopp = str2double(stopp);
    Data=mat2gray(Data);
    %% subtracting background
    BGFrame = (mat2gray(nanmean(Data,3))+mat2gray(nanmax(Data,[],3)))./2;
    background = repmat(BGFrame, [1 1 size(Data,3)]);
    Data_fluor = Data-background;
    Data = mat2gray(Data_fluor);
    clear Data_fluor;
    
    %% Correct baseline
    
    %     if Correct==1
    %         AVESignal(1,1:frames) = nanmean(nanmean(Data(:,:,1:frames)));
    %
    %         %%%%%%Linear trend: wts = [repmat(1/110,100,1)];
    %         wts = [repmat(1/100,100,1)];
    %         AVES = conv(AVESignal,wts,'valid');
    %         BS = str2num(get(handles.edit73,'String'));
    %         AVES = AVES-BS;
    %         LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
    %         AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
    %         % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
    %
    %         for j=1:LAVE
    %             Data(:,:,j) = Data (:,:,j)-AVES(j);
    %         end
    %     end
    
    
    %% binning
    N=3;
    [x,y]=meshgrid(-N:N,-N:N);
    avePattern=exp(-(x.^2/(2*(0.7^2))+y.^2/(2*(0.7^2))));
    N=size(avePattern,1);
    Processing_progess = waitbar(0,'Processing....');
    for k=1:10
        waitbar(k/10);
        for i = 1:size(Data,3)
            temp = Data(:,:,i);
            temp = 1/N/N*conv2(temp,avePattern,'same');
            Data_binned(:,:,i) = temp;
        end
        Data= Data_binned;
        Data = mat2gray(Data);
    end
    clear Data_binned;
    %% spatial smoothing (Moving average)
    mask = ones(2,2);
    for i = 1:size(Data,3)
        C = conv2(Data(:,:,i),mask,'same');
        Data_MA(:,:,i) = C;
    end
    Data = Data_MA;
    clear Data_MA;
    Data = mat2gray(Data);
    %             figure(4);map = colormap(gray(256));
    %         for i = 1:size(Data,3)
    %             temp = Data(:,:,i);
    %             temp = imresize(temp,2);
    %             imshow(temp,'Colormap',map);
    %             title(num2str(i));drawnow
    %         end
    %         save(strcat(FullPath,'.mat'), 'Data','RawData');
%     set(handles.popupmenu5,'Value',2);
    % end
    %
    % if PCAchb ==1
    row=size(Data,1); col=size(Data,2); nFrames = size(Data,3);
    PCAData = reshape(Data, row*col,nFrames)';
    PCAData = PCAData - median(PCAData);
    dec = mdwtdec('c',PCAData,5,'db2');
    decBIS = chgwdeccfs(dec,'cd',0,1:2);
    [XD,decDEN,THRESH] = mswden('den',decBIS,'sqtwolog','sln');
    
    Xbis = mdwtrec(decBIS);
    
    PCAData_reverse = reshape(Xbis', row,col,nFrames);
    PCAData_reverse2 = reshape(XD', row,col,nFrames);
    level = 5;
    wname = 'bior6.8';
    npc = 'heur';
    [x_sim, qual, NPC] = wmspca(PCAData',level,wname,npc);
    close(Processing_progess);
    PCAData_reverse3 = reshape(x_sim, row,col,nFrames);
    %     figure(3);cla;subplot(2,2,1);imshow(squeeze(Data(:,:,500))); title('before PCA');
    %     subplot(2,2,2);imshow(squeeze(PCAData_reverse(:,:,500))); title('after PCA')
    Sig = squeeze(nanmean(nanmean(Data)));
    Sigp = squeeze(nanmean(nanmean(PCAData_reverse3)));
    %     subplot(2,2,3:4);plot(Sig); hold on;plot(Sigp); legend('before PCA', 'after PCA')
    %     anskeep = questdlg('Do you like to keep the changes?','PCA', 'YES', 'NO', 'YES');
    %     if strcmp(anskeep, 'YES')
    Data = PCAData_reverse;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');        
%     end
    %     end
end

if MAchb == 1
    filter_span = str2num(get(handles.edit82,'String'));
    FilteredData = filter(ones(1,filter_span)/filter_span,1,Data,[],3);
    FilteredData(:,:,1:filter_span) = repmat(FilteredData(:,:,filter_span+1),[1,1,filter_span]);
    %     Sig = squeeze(nanmean(nanmean(Data)));
    %     Sigf = squeeze(nanmean(nanmean(FilteredData)));
    %     figure(3);cla;plot(Sig); hold on;plot(Sigf); legend('before MA filter', 'after MA filter')
    %     anskeep = questdlg('Do you like to keep the changes?','Moving Average filter', 'YES', 'NO', 'YES');
    %     if strcmp(anskeep, 'YES')
    Data = FilteredData;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');        
%     end
    %     end
end
if SMchb == 1
    SMOption = get(handles.popupmenu3,'Value'); %% 1 for smoothing - 2 for Averaging filter
    % 3 for Gaussian filter - 4 for Median filter
    % 5 for Adaptive Wiener filter
    dwindowval = get(handles.popupmenu4,'Value');
    dwindowlist = [3,5,7,9,11,13,15];
    dwindow = dwindowlist(dwindowval);
    % SMOOTH IMAGE:
    smoothing_progess = waitbar(0,'Smoothing....');
    for k=1:size(Data,3)
        waitbar(k/size(Data,3));
        if SMOption==1
            pre_smoothed_image = squeeze(Data(:,:,k));
            for i=1:size(pre_smoothed_image,1)
                for j=1:size(pre_smoothed_image,2)
                    first_included_i_pixel = i-floor(dwindow/2);
                    if first_included_i_pixel<1
                        first_included_i_pixel = 1;
                    end
                    last_included_i_pixel = i+floor(dwindow/2);
                    if last_included_i_pixel>size(pre_smoothed_image,1)
                        last_included_i_pixel = size(pre_smoothed_image,1);
                    end
                    
                    first_included_j_pixel = j-floor(dwindow/2);
                    if first_included_j_pixel<1
                        first_included_j_pixel = 1;
                    end
                    
                    last_included_j_pixel = j+floor(dwindow/2);
                    if last_included_j_pixel>size(pre_smoothed_image,2)
                        last_included_j_pixel = size(pre_smoothed_image,2);
                    end
                    
                    pixels_to_average = pre_smoothed_image(first_included_i_pixel:last_included_i_pixel,...
                        first_included_j_pixel:last_included_j_pixel);
                    
                    pixels_to_average = pixels_to_average(pixels_to_average>=0);
                    if length(pixels_to_average)>0
                        SmoothedData(i,j,k) = nanmean(pixels_to_average);
                    else
                        SmoothedData(i,j,k) = Data (i,j,k);
                    end
                end
            end
        elseif SMOption==2
            h = fspecial('average',[dwindow dwindow]);
            SmoothedData(:,:,k) = imfilter(Data(:,:,k),h,'replicate','same');
        elseif SMOption==3
            h = fspecial('gaussian',[dwindow dwindow],2);
            SmoothedData(:,:,k) = imfilter(Data(:,:,k),h,'replicate','same');
        elseif SMOption==4
            SmoothedData(:,:,k) = medfilt2(squeeze(Data(:,:,k)),[dwindow dwindow]);
        elseif SMOption==5
            SmoothedData(:,:,k) = wiener2(Data(:,:,k),[dwindow dwindow]);
        end
    end
    close(smoothing_progess);
    %     figure(3);cla;subplot(2,2,1);imshow(squeeze(Data(:,:,500))); title('before smoothing');
    %                   subplot(2,2,2);imshow(squeeze(SmoothedData(:,:,500))); title('after smoothing')
    %     Sig = squeeze(nanmean(nanmean(Data)));
    %     Sigs = squeeze(nanmean(nanmean(SmoothedData)));
    %     subplot(2,2,3:4);plot(Sig); hold on;plot(Sigs); legend('before smoothing', 'after smoothing')
    
    %     anskeep = questdlg('Do you like to keep the changes?','Smoothing', 'YES', 'NO', 'YES');
    %     if strcmp(anskeep, 'YES')
    Data = SmoothedData;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
%     end
    %     end
end


%         MaskedData = Data.*repmat(Mask,[1,1,frames]);
%         Data = MaskedData;
% if Invchb == 1
%     Dmax = 1;
%     Dmin = 0;
%     Data_Inv = Data-Dmax;
%     Data_Inv = abs(Data_Inv);
%     Data_Inv = Data_Inv+Dmin;
%     %     Sig = squeeze(nanmean(nanmean(Data)));
%     %     Siginv = squeeze(nanmean(nanmean(Data_Inv)));
%     %     figure(3);cla;plot(Sig); hold on;plot(Siginv); legend('before inverting', 'after inverting')
%     Data = Data_Inv;
% %     if size(Data,3)>2000
% %         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
% %     else
% %         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
% %     end
%     
% end
if bcorrchb ==1
%     Mask = double(Mask);
%     Mask( Mask == 0 ) = NaN;
%     Data = Data .* Mask;
    AVESignal(1,1:size(Data,3)) = nanmean(nanmean(Data(:,:,1:size(Data,3))));
    bcorrmethod = get(handles.popupmenu6,'Value');
    if bcorrmethod==1
        %option1
%         t = (1:length(AVESignal));
%         opol = 6;
%         [p,s,mu] = polyfit(t,AVESignal,opol);
%         bleachsignal = polyval(p,t,[],mu);
%         AVECorrected = AVESignal - bleachsignal;
%         
%         %         AVECorrected = detrend(AVESignal);
  Sig = smooth(squeeze(AVESignal),0.01);
  Sig3 = smooth(1:length(Sig),squeeze(Sig),0.9,'loess');
    
%     Sig2 = smooth(squeeze(Sig),0.5);
%     Sig3 = filter(ones(1,101)/101,1,Sig2);
%     Sig3(1:101)=Sig3(102);
    bleachsignal = Sig3+(max(Sig)-max(Sig3));
  AVECorrected = AVESignal - bleachsignal;
        
    elseif bcorrmethod ==2        %option2:
        frame_rate=1000*total_frames/total_time;
        cutoff_frequency = 0.5;
        filter_order = 2;
        normalized_cutoff_frequency = cutoff_frequency/(frame_rate/2);
        [b,a] = butter(filter_order,normalized_cutoff_frequency,'high');
        AVECorrected = filter(b,a,AVESignal);
        AVECorrected (1:50) = AVECorrected (51);
        bleachsignal = AVESignal - AVECorrected;
    elseif bcorrmethod ==3 %option1:
        wts = [repmat(1/100,100,1)];
        AVES = conv(AVESignal,wts,'valid');
        % AVES = AVES-(AVESignal(1)-AVECorrected(1));
        LAVE=nanmin(length(AVESignal),length(AVES)); %% Length of the drif signal and corrected signal
        AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
        % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
        Deltay = (AVESignal(1)-AVECorrected(1));
        AVES = AVES-(AVESignal(1)-AVECorrected(1));
        AVECorrected = AVESignal(1:LAVE)-AVES(1:LAVE);
        AVECorrected(900:1000)=AVECorrected(899);
        bleachsignal = AVES+(AVESignal(1)-Deltay);
        bleachsignal(900:1000)=bleachsignal(899);
    end
%         figure;plot(AVESignal,'k');hold on; plot(AVECorrected,'g')
%         legend('Bleached Signal','Corrected Signal')
%     Correcting_progess = waitbar(0,'Correct for Bleaching....');
%       bleachsignal (1,1,:) = bleachsignal; 
%     for i=1:size(Data,3)
%         waitbar (i/size(Data,3));
%         CorrectedData(:,:,i) = Data(:,:,i) - bleachsignal(i);
%     end
%     close (Correcting_progess)
    
       bleachsignal1 (1,1,:) = bleachsignal;
    bleachmatrix = repmat(bleachsignal1,[size(Data,1),size(Data,2),1]);
    CorrectedData = Data - bleachmatrix;

    %         BeforCorrectionData = Data;
    Data = CorrectedData;
%     if size(Data,3)>2000
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
%     else
%         save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
%     end
end
Data = Data -min(min(min(Data)));
Data = mat2gray(Data,[mean(mean(mean(Data)))./2,1]);
if size(Data,3)>2000
    save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
else
    save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
end


%         MaskedData = Data.*repmat(Mask,[1,1,frames]);
%         Data = MaskedData;

% min_manual_value = 0;
% max_manual_value = 1;
% 
% figure(4);
% % if preproschb == 1
% % else
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
% % DataforPlot = mat2gray(Data);
% vraw = VideoWriter([fullfile(FullPath),'_raw','.avi'])
% vraw.FrameRate = 30;
% % figure;
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
% save (strcat(FullPath,'.mat'),'RawData','Data','Mask');
implay(mat2gray(Data),60)
% end
Sig = squeeze(nanmean(nanmean(Data)));
figure;plot(Sig);xlabel('frame');
fprintf('Processing Complete!\n');

% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15



function edit79_Callback(hObject, eventdata, handles)
% hObject    handle to edit79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit79 as text
%        str2double(get(hObject,'String')) returns contents of edit79 as a double


% --- Executes during object creation, after setting all properties.
function edit79_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit80_Callback(hObject, eventdata, handles)
% hObject    handle to edit80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit80 as text
%        str2double(get(hObject,'String')) returns contents of edit80 as a double


% --- Executes during object creation, after setting all properties.
function edit80_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit82_Callback(hObject, eventdata, handles)
% hObject    handle to edit82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit82 as text
%        str2double(get(hObject,'String')) returns contents of edit82 as a double


% --- Executes during object creation, after setting all properties.
function edit82_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit83_Callback(hObject, eventdata, handles)
% hObject    handle to edit83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit83 as text
%        str2double(get(hObject,'String')) returns contents of edit83 as a double


% --- Executes during object creation, after setting all properties.
function edit83_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function checkbox11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
th = get(handles.slider1, 'Value');
set(handles.edit84,'String',num2str(th));
FullPath = get(handles.edit29,'String');
load(strcat(FullPath,'.mat'));
BGFrame = mat2gray(nanmax(RawData,[],3));

Mask = BGFrame>th;
Mask = double(Mask);
% Mask( Mask == 0 ) = NaN;
flipoption = get(handles.checkbox29,'Value');
if flipoption == 1
    Mask = fliplr(Mask);
    Data = fliplr(Data);
end
MaskedData = Data.*repmat(Mask,[1,1,size(Data,3)]);
figure(10);imshow(Mask); title('MASK');
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit84_Callback(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
th = get(handles.edit84, 'String');
set(handles.slider1,'Value',str2num(th));
FullPath = get(handles.edit29,'String');
load(strcat(FullPath,'.mat'));
BGFrame = mat2gray(nanmax(RawData,[],3));

th = get(handles.slider1, 'Value');
Mask = BGFrame>th;
Mask = double(Mask);
% Mask( Mask == 0 ) = NaN;
flipoption = get(handles.checkbox29,'Value');
if flipoption == 1
    Mask = fliplr(Mask);
    Data = fliplr(Data);
end
MaskedData = Data.*repmat(Mask,[1,1,size(Data,3)]);
figure(10);imshow(Mask); title('MASK');
% Hints: get(hObject,'String') returns contents of edit84 as text
%        str2double(get(hObject,'String')) returns contents of edit84 as a double


% --- Executes during object creation, after setting all properties.
function edit84_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
RoiOption = get(handles.popupmenu2,'Value'); %% 1 for freehand - 2 for impoly - 3 for imellipse - 4 for imrect
FullPath = get(handles.edit29,'String');
load(strcat(FullPath,'.mat'));
info = imfinfo(strcat(FullPath,'.tif'));
frames = numel(info);
ScalingOpt = get(handles.checkbox30,'Value');
%% Define a ROI
if RoiOption == 1
    %     BGFrame = mat2gray(nanmax(RawData,[],3));
    %     th = get(handles.slider1, 'Value');
    %     Mask = BGFrame>th;
    %     Mask = double(Mask);
    %     Mask( Mask == 0 ) = NaN;
    %     MaskedData = Data.*repmat(Mask,[1,1,frames]);
    %     Data = MaskedData;
    BGFrame = mat2gray(nanmax(RawData,[],3));
    ansSatisfied='NO';
    for n=1:1000
        if strcmp(ansSatisfied,'YES')
            flipoption = get(handles.checkbox29,'Value');
            if flipoption == 1
                Data = fliplr(Data);
            end
            MaskedData = Data.*repmat(Mask,[1,1,size(Data,3)]);
            Data = MaskedData;
            if size(Data,3)>2000
                save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
            else
                save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
            end
            break
        elseif strcmp(ansSatisfied,'Cancel masking. I changed my mind!')
            break
        elseif strcmp(ansSatisfied,'NO')
            th = get(handles.slider1, 'Value');
            Mask = BGFrame>th;
            Mask = double(Mask);
            %             Mask( Mask == 0 ) = NaN;
            flipoption = get(handles.checkbox29,'Value');
            if flipoption == 1
                Mask = fliplr(Mask);
                %                 Data = fliplr(Data);
            end
            figure;imshow(Mask); title('MASK');
            %             MaskedData = Data.*repmat(Mask,[1,1,frames]);
            ansSatisfied = questdlg('Do you like the mask?','MASK', 'YES', 'Cancel masking. I changed my mind!', 'YES');
        end
    end
else
    ansSatisfied='NO';
    for n=1:1000
        if strcmp(ansSatisfied,'YES')
            if flipoption == 1
                Data = fliplr(Data);
            end
            MaskedData = Data .*repmat(Mask,[1,1,size(Data,3)]);
            Data = MaskedData;
            if size(Data,3)>2000
                save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt','-v7.3');
            else
                save (strcat(FullPath,'.mat'),'RawData','Data','Mask','ScalingOpt');
            end
            break
        elseif strcmp(ansSatisfied,'Cancel masking. I changed my mind!')
            break
        elseif strcmp(ansSatisfied,'NO')
            tempframe = mat2gray(imresize(RawData(:,:,50),3));
            
            flipoption = get(handles.checkbox29,'Value');
            if flipoption == 1
                tempframe = fliplr(tempframe);
                %                 Data = fliplr(Data);
            end
            
            figure(1);imshow(tempframe);
            
            if RoiOption==2 %freehand
                h = imfreehand(gca);
            elseif RoiOption==3 %impoly
                h = impoly(gca);
            elseif RoiOption==4 %imellipse
                h = imellipse (gca);
            elseif RoiOption==5 %imrect
                h = imrect (gca);
            end
            Mask =createMask(h);
            Mask = imresize(Mask,1/3);
            Mask = double(Mask);
            %             Mask( Mask == 0 ) = NaN;
            figure;imshow(Mask); title('MASK');
            
            %             flipoption = get(handles.checkbox29,'Value');
            %             if flipoption == 1
            %                     Mask = fliplr(Mask);
            %                 Data = fliplr(Data);
            %             end
            %             MaskedData = Data .*repmat(Mask,[1,1,frames]);
        end
        ansSatisfied = questdlg('Do you like the mask?','MASK', 'YES', 'NO', 'Cancel masking. I changed my mind!', 'NO');
    end
    %     figure(10);imshow(imresize(Data(:,:,50),3));
end

fprintf('Masking Complete!\n');

% --- Executes on button press in pushbutton47
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
FullPath = get(handles.edit29,'String');
load(strcat(FullPath,'.mat'));
AVERAGING = get(handles.checkbox2,'Value');
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp);

if AVERAGING
    AveragedData(find(isnan(AveragedData)))=0;
    AveragedData = mat2gray(AveragedData);
%     AveragedData = AveragedData .*Mask;
    if size(AveragedData,1)>200
        AveragedData = imresize(AveragedData,0.5);
    end
    AveragedData = mat2gray(AveragedData);
    implay(imresize(AveragedData,3))
else
    Data(find(isnan(Data)))=0;
    Data = mat2gray(Data);
%     Data = Data .*Mask;
    if size(Data,1)>200
        Data = imresize(Data,0.5);
    end
    Data = Data(:,:, startp:endp);
    implay(imresize(Data,3))
end

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
% set(handles.checkbox2,'Value',1);
% set(handles.checkbox23,'Value',0);
% Averaging=get(handles.checkbox2,'Value');
Averaging=1
manual_mode=1;
% Averaging=1;
APDT=50 ;     % choose 50 for APD50, 70 for APD70, and 90 for APD90
manual_RMP=0;
RMP=6.2;
t_RMP=88.82;


series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');

if series == 1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack == 1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end



total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp);
if endp>size(Data,3)
    endp = size(Data,3); 
end
Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
xlim = [1,size(Data,3)];
NumFrame = get(handles.edit62,'String');
NumFrame = str2double(NumFrame);
ScalingOpt = get(handles.checkbox30,'Value');

frame_rate=1000*total_frames/total_time;
t=1:length(startp:endp);
t=t*tresolution;
I=t/tresolution;
filter_span=5;
FilteredData = filter(ones(1,filter_span)/filter_span,1,Data,[],3);
FilteredData(:,:,1:filter_span) = repmat(FilteredData(:,:,filter_span+1),[1,1,filter_span]);
Data=FilteredData;
Sig = squeeze(nanmean(nanmean(Data)))';
Sig = 100*(Sig-min(Sig))/(max(Sig)-min(Sig));
[Min, tMin, IMin, Max, tMax, IMax, NPEAKS]=fpeaks([1:size(Sig,2)],Sig, 'data');
Sig=[]; APD50=[];APD70=[];APD90=[];
tic
h_1 = waitbar (0 , 'Please wait...');
for x_coord = 1:size(Data,2)
    waitbar(x_coord/size(Data,2));
    for y_coord = 1:size(Data,1)
        % for x_coord = 14:14
        %     for y_coord = 22:22
        if isnan(Mask(y_coord,x_coord))==0
            x = squeeze(Data (y_coord,x_coord,startp:endp))';
            x=(x-nanmin(x))/(nanmax(x)-nanmin(x));
            %             if isempty(find(isnan(x)))
            
            %     xCorrected = filter(b,a,x);
            %     xCorrected (1:50) = xCorrected (51);
            %     x = xCorrected;
            %     Signal=[];
            %     % % smooth x data
            %     Signal(1,:) = x;
            %     % Signal(1,:) = nanmean(Sig);
            %     S2 = smooth(1:size(x),Signal,0.9,'rloess');
            %     Signal = Signal-S2'; % drift removal
            %     window = 5;
            %     mask = ones(1,window)/window;
            %     maY=conv(Signal,mask,'same'); % temporal smoothing
            %     x=maY(1,:);
            OAP_signal=x-nanmin(x);
            % %
%             size(t)
%             size(x)
%             figure;plot(t,x)
            [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, x, 'data');
            %                 IMax
            Npeaks;
            NPEAKS;
            %                 if Npeaks == NPEAKS
            if Npeaks>=3
                % maximums=tMax
                % difmax= diff(maximums)
                %     checkresult_d(t, x, Min, tMin, Max, tMax)
                CL=nanmean(diff(tMax));
                
                OAP = OAP_signal;
                tt  = t;
                II  = I;
                
                [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
                %                 IMax
                %     close
                maximums=tMax;
                % checkresult_d(tt, OAP, Min, tMin, Max, tMax)
                
                CL=nanmean(diff(tMax));
                %     CLstr=num2str(CL);
                %     set(handles.text55,'String',[CLstr, ' ms'])
                %
                CL=CL/tresolution;
                OAP_new = OAP;
                tt_new=tt;
                II_new=II;
                OAP_new=OAP_new-nanmin(OAP_new);
                OAP_new = OAP_new/nanmax(OAP_new)*100;
%                     D=diff(OAP_new);
%                     % figure;plot(OAP_new);hold on;plot(D)
%                     D=smooth(D,0.03)';
%                     Dpart1 = D(1:tjmax);
%                     Dpart1(Dpart1<0.3)=0;
%                     D(1:tjmax)=Dpart1;
%                     S=sign(D);
%                     SD=diff(S);
%                     SD2=diff(SD);
%                     SDS=sign(SD);
%                     SDS2=-sign(SD2);
%                     SDS=[0,SDS,0];
%                     SDS2=[0,SDS2,0,0];
%                     [SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(tt_new, SDS, 'model');
%                     % SIMax
%                     % scrsz = get(0,'ScreenSize');figure('Position',[scrsz(3)/3 scrsz(4)/1.8 scrsz(3)/2 scrsz(4)/3])
                    % plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');
                    [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt_new, OAP_new, 'data');
                    %                 IMax
                    %                 if Npeaks>2
                    if Npeaks>2 && Averaging
                        Y=[];
                        cycles = 0;
                        CL=ceil(CL);
                        for l=1:length(IMax)
                            if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP_new))
                                Y=[Y;OAP_new(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
                                cycles = cycles +1;
                            end
                        end
                        
                        OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
                        tt_new=tt(1:length(OAP_new));
                        II_new=II(1:length(OAP_new));
                    end
                    if cycles<2 || Npeaks<2 || ~Averaging
                        OAP_new = OAP;
                        tt_new=tt;
                        II_new=II;
                    end
                    OAP_new=OAP_new-nanmin(OAP_new);
                    OAP_new = OAP_new/nanmax(OAP_new)*100;
                                    [yimax, tjmax]=findpeaks(OAP_new,'MinPeakProminence',max(OAP_new)/2,'Annotate','extents');
                if ~isempty(yimax)

                    D=diff(OAP_new);
                    % figure;plot(OAP_new);hold on;plot(D)
                    D=smooth(D,0.03)';
                    Dpart1 = D(1:tjmax);
                    Dpart1(Dpart1<0.3)=0;
                    D(1:tjmax)=Dpart1;
                    S=sign(D);
                    SD=diff(S);
                    SD2=diff(SD);
                    SDS=sign(SD);
                    SDS2=-sign(SD2);
                    SDS=[0,SDS,0];
                    SDS2=[0,SDS2,0,0];
                    [SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(tt_new, SDS, 'model');
                    %     SIMax
                    
                    %                 cla(figure(2));plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');
                    [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt_new, OAP_new, 'data');
                    % IMax
                    a1=find(SIMin==IMax(1));
                    b1=find(SIMax<IMax(1));
                    if size(b1,2)>1
                        b1=b1(end-1);
                    end
                    b2=find(SIMax>=IMax(1));
                    if size(b2>1)
                    b2=b2(1);
                    % APA   = OAP_new(SIMin(a1));
                    APA = yimax(1);
                    % t_APA = tt_new(SIMin(a1));
%                     tjmax
                    t_APA = tt_new(tjmax(1));
                    % I_APA = SIMin(a1);
                    I_APA = tjmax(1);
                    if manual_RMP
                        I_RMP=ceil(t_APA/tresolution);
                    else
                        RMP   = OAP_new(SIMax(b2));
                        t_RMP = tt_new(SIMax(b2));
                        I_RMP = SIMax(b2);
                    end
                    TP = OAP_new(SIMax(b1));
                    t_TP = tt_new(SIMax(b1));
                    I_TP = SIMax(b1);
                    end
                    %                 hold on; plot(t_APA,APA,'rs','markerfacecolor','r');
                    %                 hold on; plot(t_RMP,RMP,'go','markerfacecolor','g');
                    %                 hold on; plot(t_TP,TP,'c^','markerfacecolor','c');
                    %                 hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
                    axislimits=nanmax(OAP_new)-nanmin(OAP_new);
                    %                 axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
                    %                 legend('Action Potential (AP)','Action Potential Amplitude (APA)','Rest Membrane Potential (RMP)','Take-off Potential (TP)' )
                    %% APD50
                    if ~isempty(TP)
                        AP50= (((50)/100)*(APA-TP))+TP;
                        % hold on; plot([0, tt_new(end)],[AP50, AP50] , 'k--')
                        % hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
                        C=AP50*ones(size(OAP_new));
                        D=diff(sign(OAP_new-C));
                        %                                 Npeaks
                        %                                 sum(find(D==2))
                        %                                 sum(find(D==-2))
                        if Npeaks>1&&sum(find(D==2))>1&&sum(find(D==-2))>1
                            % hold on; plot(tt_new,sign(OAP_new-C));
                            % hold on; plot(tt_new(2:end),D,'r');
                            I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
                            I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
                            t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
                            OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
                            t2a_APD50=tt_new(I2a_APD50);t2b_APD50=tt_new(I2b_APD50);
                            OAP2a_APD50=OAP_new(I2a_APD50);OAP2b_APD50=OAP_new(I2b_APD50);
                            T1=((AP50-OAP1a_APD50(1))/(OAP1b_APD50(1)-OAP1a_APD50(1)))*(t1b_APD50(1)-t1a_APD50(1))+t1a_APD50(1);
                            T2=((AP50-OAP2a_APD50(1))/(OAP2b_APD50(1)-OAP2a_APD50(1)))*(t2b_APD50(1)-t2a_APD50(1))+t2a_APD50(1);
                            %     T3=((AP50-OAP1a_APD50(2))/(OAP1b_APD50(2)-OAP1a_APD50(2)))*(t1b_APD50(2)-t1a_APD50(2))+t1a_APD50(2);
                            %     T4=((AP50-OAP2a_APD50(2))/(OAP2b_APD50(2)-OAP2a_APD50(2)))*(t2b_APD50(2)-t2a_APD50(2))+t2a_APD50(2);
                            %                     hold on; plot(T2, AP50,'d', 'markeredgecolor', [0.25,0.75,0.75], 'Markerfacecolor',[0.25,0.75,0.75])
                            %                     hold on; plot([t_APA,t_APA],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
                            %                     hold on; plot([T2,T2],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
                            APD50(y_coord,x_coord) = T2-t_APA;
                            %                     text(20,85,['APD50',' = ', num2str(APD50(y_coord,x_coord))], 'color', [0.25,0.75,0.75])
                        end
                        %% APD 70
                        
                        AP50= (((30)/100)*(APA-TP))+TP;
                        % hold on; plot([0, tt_new(end)],[AP50, AP50] , 'k--')
                        % hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
                        C=AP50*ones(size(OAP_new));
                        D=diff(sign(OAP_new-C));
                        if Npeaks>1&&sum(find(D==2))>1&&sum(find(D==-2))>1
                            
                            % hold on; plot(tt_new,sign(OAP_new-C));
                            % hold on; plot(tt_new(2:end),D,'r');
                            I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
                            I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
                            t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
                            OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
                            t2a_APD50=tt_new(I2a_APD50);t2b_APD50=tt_new(I2b_APD50);
                            OAP2a_APD50=OAP_new(I2a_APD50);OAP2b_APD50=OAP_new(I2b_APD50);
                            T1=((AP50-OAP1a_APD50(1))/(OAP1b_APD50(1)-OAP1a_APD50(1)))*(t1b_APD50(1)-t1a_APD50(1))+t1a_APD50(1);
                            T2=((AP50-OAP2a_APD50(1))/(OAP2b_APD50(1)-OAP2a_APD50(1)))*(t2b_APD50(1)-t2a_APD50(1))+t2a_APD50(1);
                            %     T3=((AP50-OAP1a_APD50(2))/(OAP1b_APD50(2)-OAP1a_APD50(2)))*(t1b_APD50(2)-t1a_APD50(2))+t1a_APD50(2);
                            %     T4=((AP50-OAP2a_APD50(2))/(OAP2b_APD50(2)-OAP2a_APD50(2)))*(t2b_APD50(2)-t2a_APD50(2))+t2a_APD50(2);
                            %                     hold on; plot(T2, AP50,'d', 'markeredgecolor', [1,0.5,0.3], 'Markerfacecolor',[1,0.5,0.3])
                            %                     hold on; plot([t_APA,t_APA],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
                            %                     hold on; plot([T2,T2],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
                            APD70(y_coord,x_coord) = T2-t_APA;
                            %                     text(20,75,['APD70',' = ', num2str(APD70(y_coord,x_coord))],'color',[1,0.5,0.3])
                        end
                        %% APD90
                        AP50= (((10)/100)*(APA-TP))+TP;
                        % hold on; plot([0, tt_new(end)],[AP50, AP50] , 'k--')
                        % hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
                        C=AP50*ones(size(OAP_new));
                        D=diff(sign(OAP_new-C));
                        if Npeaks>1&&sum(find(D==2))>1&&sum(find(D==-2))>1
                            
                            % hold on; plot(tt_new,sign(OAP_new-C));
                            % hold on; plot(tt_new(2:end),D,'r');
                            I1b_APD50=find(D== 2);I1a_APD50=I1b_APD50+1;
                            I2b_APD50=find(D==-2);I2a_APD50=I2b_APD50+1;
                            t1a_APD50=tt_new(I1a_APD50);t1b_APD50=tt_new(I1b_APD50);
                            OAP1a_APD50=OAP_new(I1a_APD50);OAP1b_APD50=OAP_new(I1b_APD50);
                            t2a_APD50=tt_new(I2a_APD50);t2b_APD50=tt_new(I2b_APD50);
                            OAP2a_APD50=OAP_new(I2a_APD50);OAP2b_APD50=OAP_new(I2b_APD50);
                            T1=((AP50-OAP1a_APD50(1))/(OAP1b_APD50(1)-OAP1a_APD50(1)))*(t1b_APD50(1)-t1a_APD50(1))+t1a_APD50(1);
                            T2=((AP50-OAP2a_APD50(1))/(OAP2b_APD50(1)-OAP2a_APD50(1)))*(t2b_APD50(1)-t2a_APD50(1))+t2a_APD50(1);
                            %     T3=((AP50-OAP1a_APD50(2))/(OAP1b_APD50(2)-OAP1a_APD50(2)))*(t1b_APD50(2)-t1a_APD50(2))+t1a_APD50(2);
                            %     T4=((AP50-OAP2a_APD50(2))/(OAP2b_APD50(2)-OAP2a_APD50(2)))*(t2b_APD50(2)-t2a_APD50(2))+t2a_APD50(2);
                            %                     hold on; plot(T2, AP50,'d', 'markeredgecolor', [0.7,0.4,1], 'Markerfacecolor',[0.7,0.4,0.5])
                            %                     hold on; plot([t_APA,t_APA],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
                            %                     hold on; plot([T2,T2],[nanmin(OAP_new)-5 nanmax(OAP_new)+20],'k--')
                            APD90(y_coord,x_coord) = T2-t_APA;
                            %                     text(20,65,['APD90',' = ', num2str(APD90(y_coord,x_coord))],'color',[0.7,0.4,1])
                            %     scrsz = get(0,'ScreenSize');
                            %     figure('Position',[scrsz(3)/3 scrsz(4)/10 scrsz(3)/2 scrsz(4)/3])
                            %     plot(tt_new,OAP_new,'k','linewidth',2);xlabel('time (ms)');
                            %     hold on; plot([0, tt_new(end)],[0, 0] , 'k--')
                            %     axislimits=nanmax(OAP_new)-nanmin(OAP_new);
                            %     axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
                        end
                    end
                end
            end
        end
    end
end
close (h_1)
toc(tic)
BGFrame = imresize(RawData(:,:,50),3);

APD50(APD50<1) = 0;
APD70(APD70<1) = 0;
APD90(APD90<1) = 0;
APD50p = imgaussfilt(APD50,3);
APD70p = imgaussfilt(APD70,3);
APD90p = imgaussfilt(APD90,3);
map50 = jet(ceil(max(max(APD50p))));
map50 (1,:,:) = [0 0 0];
map70 = jet(ceil(max(max(APD70p))));
map70 (1,:,:) = [0 0 0];
map90 = jet(ceil(max(max(APD90p))));
map90 (1,:,:) = [0 0 0];
figure;imshow(imresize(APD50p,3),map50);colorbar;title('APD50 (ms)')
figure;imshow(imresize(APD70p,3),map70);colorbar;title('APD70 (ms)')
figure;imshow(imresize(APD90p,3),map90);colorbar;title('APD90 (ms)')
flipoption = get(handles.checkbox29,'Value');
onPrepchb = get(handles.checkbox32, 'Value');
if onPrepchb
    APDonPrep (APD50p,BGFrame, flipoption,'APD50 (ms)')
    APDonPrep (APD70p,BGFrame, flipoption,'APD70 (ms)')
    APDonPrep (APD90p,BGFrame, flipoption,'APD90 (ms)')
end
% % % x=linspace(1,size(APD50p,2),size(APD50p,2));
% % % y=linspace(1,size(APD50p,1),size(APD50p,1));
% % % [X,Y]=meshgrid(x,y);
% % % figure;imshow(BGFrame);hold on; contour(X,Y,APD50p,map50);colorbar;title('APD50 (ms)')
if size(Data,3)>2000
    save(strcat(FullPath,'APDs.mat'),'APD50p','APD70p','APD90p','RawData','Data','ScalingOpt','-v7.3');
else
    save(strcat(FullPath,'APDs.mat'),'APD50p','APD70p','APD90p','RawData','Data','ScalingOpt')
end
% figure;imshow(mat2gray(BGFrame));hold on; contour(imresize(APD50p,3));colorbar;title('APD50 (ms)')
APD50(APD50==0)=nan;
APD70(APD70==0)=nan;
APD90(APD90==0)=nan;
exportchb = get(handles.checkbox24, 'Value');
if exportchb
    javaaddpath('C:/Yair/Utils/JExcelAPI/jxl.jar')
    javaaddpath('C:/Yair/Utils/JExcelAPI/MXL.jar')
    javaaddpath('poi_library/poi-3.8-20120326.jar');
    javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
    javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
    javaaddpath('poi_library/xmlbeans-2.3.0.jar');
    javaaddpath('poi_library/dom4j-1.6.1.jar');
    javaaddpath('poi_library/stax-api-1.0.1.jar');

    xlwrite(strcat(FullPath, 'APD50.xls'),APD50,1)
    xlwrite(strcat(FullPath, 'APD70.xls'),APD70,1)
    xlwrite(strcat(FullPath, 'APD90.xls'),APD90,1)
    %%%% OR:
    filename = strcat(FullPath, 'APDs.xls');
    xlswrite(filename,APD50,1)
    xlswrite(filename,APD70,2)
    xlswrite(filename,APD90,3)
%     
%     ActivexServer = actxserver('Excel.Application'); % # open Activex server
%     ExcelFile = ActivexServer.Workbooks.Open(filename); % # open file
%     ExcelFile.Worksheets.Item(1).Name = 'APD50'; % # rename 1st sheet
%     ExcelFile.Worksheets.Item(2).Name = 'APD70';
%     ExcelFile.Worksheets.Item(3).Name = 'APD90';
%     ExcelFile.Save % # save to the same file
%     ExcelFile.Close(false)
%     ActivexServer.Quit
end


% figure;contourf (APD50)
% figure;contourf (APD70)
% figure;contourf (APD90)

% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox23,'Value')==1
    set(handles.checkbox2,'Value',0)
elseif get(handles.checkbox23,'Value')==0
    set(handles.checkbox2,'Value',1)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox23



% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end
ScalingOpt = get(handles.checkbox30,'Value');
% AVERAGING = 1;
AVERAGING=get(handles.checkbox2,'Value');
if AVERAGING
    Data = AveragedData;
end
DataSize= size(Data)

% size(Data);
%**********
SMOption = get(handles.popupmenu3,'Value'); %% 1 for smoothing - 2 for Averaging filter
    % 3 for Gaussian filter - 4 for Median filter
    % 5 for Adaptive Wiener filter
    dwindowval = get(handles.popupmenu4,'Value');
    dwindowlist = [3,5,7,9,11,13,15];
    dwindow = dwindowlist(dwindowval);
    % SMOOTH IMAGE:
    smoothing_progess = waitbar(0,'Please wait....');
    for k=1:size(Data,3)
        waitbar(k/size(Data,3));
        if SMOption==1
            pre_smoothed_image = squeeze(Data(:,:,k));
            for i=1:size(pre_smoothed_image,1)
                for j=1:size(pre_smoothed_image,2)
                    first_included_i_pixel = i-floor(dwindow/2);
                    if first_included_i_pixel<1
                        first_included_i_pixel = 1;
                    end
                    last_included_i_pixel = i+floor(dwindow/2);
                    if last_included_i_pixel>size(pre_smoothed_image,1)
                        last_included_i_pixel = size(pre_smoothed_image,1);
                    end
                    
                    first_included_j_pixel = j-floor(dwindow/2);
                    if first_included_j_pixel<1
                        first_included_j_pixel = 1;
                    end
                    
                    last_included_j_pixel = j+floor(dwindow/2);
                    if last_included_j_pixel>size(pre_smoothed_image,2)
                        last_included_j_pixel = size(pre_smoothed_image,2);
                    end
                    
                    pixels_to_average = pre_smoothed_image(first_included_i_pixel:last_included_i_pixel,...
                        first_included_j_pixel:last_included_j_pixel);
                    
                    pixels_to_average = pixels_to_average(pixels_to_average>=0);
                    if length(pixels_to_average)>0
                        SmoothedData(i,j,k) = nanmean(pixels_to_average);
                    else
                        SmoothedData(i,j,k) = Data (i,j,k);
                    end
                end
            end
        elseif SMOption==2
            h = fspecial('average',[dwindow dwindow]);
            SmoothedData(:,:,k) = imfilter(Data(:,:,k),h,'replicate','same');
        elseif SMOption==3
            h = fspecial('gaussian',[dwindow dwindow],2);
            SmoothedData(:,:,k) = imfilter(Data(:,:,k),h,'replicate','same');
        elseif SMOption==4
            SmoothedData(:,:,k) = medfilt2(Data(:,:,k),[dwindow dwindow],'replicate');
        elseif SMOption==5
            SmoothedData(:,:,k) = wiener2(Data(:,:,k),[dwindow dwindow]);
        end
    end
    close(smoothing_progess)
    %     figure(3);cla;subplot(2,2,1);imshow(squeeze(Data(:,:,500))); title('before smoothing');
    %                   subplot(2,2,2);imshow(squeeze(SmoothedData(:,:,500))); title('after smoothing')
    %     Sig = squeeze(nanmean(nanmean(Data)));
    %     Sigs = squeeze(nanmean(nanmean(SmoothedData)));
    %     subplot(2,2,3:4);plot(Sig); hold on;plot(Sigs); legend('before smoothing', 'after smoothing')
    
    %     anskeep = questdlg('Do you like to keep the changes?','Smoothing', 'YES', 'NO', 'YES');
    %     if strcmp(anskeep, 'YES')
    Data = SmoothedData;
DataSize= size(Data)

%**********



%%binning%%
% binsize = 2;
% CVData = imresize(AveragedData,1/binsize);
% RawData = imresize(RawData,1/binsize);
% Data = imresize(Data,1/binsize);
% Mask = imresize(Mask,1/binsize);
%%%%%%%%

% Data_OAP = Data ;
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
total_pixels = get(handles.edit52,'String');
total_pixels = str2double(total_pixels);
if get(handles.checkbox30,'Value')==1
    total_pixels = total_pixels/2;
end
total_distance = get(handles.edit53,'String');
total_distance = str2double(total_distance);
sresolution=total_distance/total_pixels;

Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp);
binsize = get(handles.edit87,'String');
binsize = str2double (binsize);
vectorfield = get(handles.checkbox28,'Value');
xlim = [1,size(Data,3)];
MinCV = get(handles.edit88 ,'String');
MinCV = str2double(MinCV);
MaxCV = get(handles.edit89 ,'String');
MaxCV = str2double(MaxCV);
% FirstFrame = get(handles.edit62,'String');
% FirstFrame = str2double(FirstFrame);
% OAP_chkbx_status =0;
% array = [];
%Record Optical Action Potentials
% pt = 1;
% array_area = 2;

Sig = squeeze(nanmean(nanmean(Data)))';
% size(Sig);
[Min, tMin, IMin, Max, tMax, IMax, NPEAKS]=fpeaks([1:size(Sig,2)],Sig, 'data');

array_area = floor(str2num(get(handles.edit85,'String'))/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RoiOption = get(handles.popupmenu7,'Value'); %% 1 for freehand - 2 for impoly - 3 for imellipse - 4 for imrect
% size(Mask)
%% Define a ROI
% d1 = figure;
if RoiOption == 1
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    MASK = Mask;
%     figure(d1); imshow(mat2gray(double(imresize(RawData(:,:,5),3))).*imresize(Mask,3));
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
else
    d1 = figure;
    ansSatisfied='NO';
    for n=1:1000
        if strcmp(ansSatisfied,'YES')
            pos = getPosition (h);
            tempframe = mat2gray(double(imresize(RawData(:,:,5),3)));
            flipoption = get(handles.checkbox29,'Value');
            if flipoption == 1
                tempframe = fliplr(tempframe);
            end
            figure(d1); imshow(tempframe.*imresize(Mask,3));
            hold on; plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],'g','LineWidth',2)
            break
        elseif strcmp(ansSatisfied,'Cancel!')
            break
        elseif strcmp(ansSatisfied,'NO')
            tempframe = mat2gray(double(imresize(RawData(:,:,5),3)));
            flipoption = get(handles.checkbox29,'Value');
            if flipoption == 1
                tempframe = fliplr(tempframe);
            end
            figure(d1); imshow(tempframe.*imresize(Mask,3));
            h = impoly(gca);
            MASK = createMask(h);
            MASK = imresize(MASK,1/3);
            MASK = double(MASK);
            MASK( MASK == 0 ) = NaN;
        end
        ansSatisfied = questdlg('Calulate CV in this ROI?','MASK', 'YES', 'NO', 'Cancel!', 'NO');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CV = zeros(size(MASK));
wait_h = waitbar(0,'Calculating conduction velocity in the ROI...');


%%**********************
DataSize= size(Data)

CVData = mat2gray(Data);

wait1=waitbar(0,'Please wait...');
for x = 1:size(CVData,1)
    for y = 1:size(CVData,2)
        if ~isnan(MASK(x,y))
%             OAP = (squeeze(CVData(x,y,:))-nanmin(squeeze(CVData(x,y,:))))/(nanmax(squeeze(CVData(x,y,:)))-nanmin(squeeze(CVData(x,y,:))));
            OAP = squeeze(CVData(x,y,:))-nanmin(squeeze(CVData(x,y,:)));
            OAP = smooth(OAP,10,'sgolay');
            OAP = OAP';
            I = 1 : length(OAP);
            t = I * tresolution;
%             AverageOAP = squeeze(mean(mean(CVData)))-nanmin(squeeze(mean(mean(CVData))));
            AverageOAP = (squeeze(nanmean(nanmean(CVData)))-nanmin(squeeze(nanmean(nanmean(CVData)))))/(nanmax(squeeze(nanmean(nanmean(CVData))))-nanmin(squeeze(nanmean(nanmean(CVData)))));
            MASK(x,y);
            OAPth = min(AverageOAP)+(0.1*(max(AverageOAP)-min(AverageOAP)));
            I_AverageOAPth =find(AverageOAP>=OAPth);
            OAP(1:I_AverageOAPth(1)-20)=0;
            OAPinterp = interp(OAP,10);
            OAPinterp(1:2) = 0;
            Iinterp = interp(I,10);
            I_OAPth =find(OAPinterp>=0.1);
            if ~isempty(I_OAPth)
                if Iinterp(I_OAPth(1))~=1
                    Iinterp_OAPth =  Iinterp(I_OAPth(1));
                elseif length(I_OAPth>2)
                    Iinterp_OAPth =  Iinterp(I_OAPth(2));
                end
                rrto_vectorplot(x,y)=Iinterp_OAPth;
            end
        end
    end
end
close(wait1)

Data = CVData;
DataSize= size(Data)

rtto = zeros(size(Data,1),size(Data,2));
rtto (1:size(rrto_vectorplot,1),1:size(rrto_vectorplot,2)) = rrto_vectorplot;
rtto = imgaussfilt(rtto,2);
rttosize = size(rtto)
rrto_vectorplotsize = size(rrto_vectorplot)
% RTTO = zeros(size(Data,1),size(Data,2));
for y_coord = array_area+1:size(Data,1)-array_area
    waitbar(y_coord/(size(Data,1)-2*array_area))
    for x_coord = array_area+1:size(Data,2)-array_area
        if ~isnan(MASK(y_coord,x_coord))
            %Construct the 5x5 pixel array
            x_init = x_coord - array_area;
            y_init = y_coord - array_area;
            
            x_final = x_coord + array_area;
            y_final = y_coord + array_area;
            % %
            idx1 = 1;
            idx2 = ((array_area*2)+1)^2;
            CLs=[];
            SIGNAL=[];
            slope=[];
            FAMI = [];
            FA_50_I = [];
            x_coord_array = [];
            y_coord_array = [];
            UV = [];
            RT_before=0;
            Times2=zeros(y_final-y_init+1,y_final-y_init+1);
            for ypos = y_init:y_final
                for xpos = x_init:x_final
                    current_x = xpos;
                    current_y = ypos;
                    current_pixel = [current_x,current_y];
                    rt = rtto(current_y,current_x);
                    %                     rt = rrto_vectorplot(current_y,current_x);
                    if rt~=0
                        FA_50_I = [FA_50_I,rt]; % rising times
                        x_coord_array = [x_coord_array,xpos];
                        y_coord_array = [y_coord_array,ypos];
                    end
                end
            end
            if numel(FA_50_I)~=0
                %Eliminate outliers in optical action potential
                size(FA_50_I);
                mean_OAP = nanmean(FA_50_I);
                clean_OAP=[];
                for w = 1:length(FA_50_I)
                    w;
                    FA_50_I(w);
                    if (FA_50_I(w) < mean_OAP-5)||(FA_50_I(w) > mean_OAP+5)
                        clean_OAP(w) = mean_OAP;
                    else
                        clean_OAP(w) = FA_50_I(w);
                    end
                end
                x = x_coord_array';
                y = y_coord_array';
                z = clean_OAP';
                Xcolv = x; % Make X a column vector
                Ycolv = y; % Make Y a column vector
                Zcolv = z; % Make Z a column vector
                if size(Zcolv,1)>=6
                    OPTIONS = fitoptions('poly22');
                    [fitresult, gof] = fit([Xcolv,Ycolv],Zcolv,'poly22',OPTIONS );
                    %              xx=reshape(Xcolv,[array_area*2+1,array_area*2+1])';
                    %              yy=reshape(Ycolv,[array_area*2+1,array_area*2+1])';
                    [xx, yy]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting
                    zz=fitresult(xx,yy);
                    
                    %             mean_OAP
                    %             FA_50_I
                    %             Const = ones(size(Xcolv)); % Vector of ones for constant term
                    %             Coefficients = [Xcolv Ycolv Const]\Zcolv; %Find the coefficients
                    %             XCoeff = Coefficients(1); % X coefficient
                    %             YCoeff = Coefficients(2); % X coefficient
                    %             CCoeff = Coefficients(3); % constant term
                    %             [xx, yy]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting
                    %
                    %             zz = XCoeff * xx + YCoeff * yy + CCoeff;
                    %             zz
                    Times2;
                    %             FA_50_I
                    %             size(zz)
                    if size(zz,1)>1 & size(zz,2)>1
                        [FX,FY] = gradient(zz);
                        [xxm, yym]=meshgrid(nanmin(x):1:nanmax(x),nanmin(y):1:nanmax(y)); % Generating a regular grid for plotting
                        
                        xxs = size(xx);
                        xxs = xxs(1)*xxs(2);
                        xx_mid = round(xxs/2);
                        zzm = zz;
                        
                        for i = 1:xxs
                            xxm(i) = xx(xx_mid);
                            yym(i) = yy(xx_mid);
                            zzm(i) = zz(xx_mid);
                        end
                        Fc=total_frames;
                        Tt=total_time;
                        Conversion_factor = Tt/Fc;
                        
                        a = abs(FX(1)); %width converted to mm (assuming 13pixels in 1mm)
                        b = abs(FY(1)); %height in frames
                        %b = b; %height converted to ms
                        
                        a_sq = (a)^2;
                        b_sq = (b)^2;
                        
                        
                        
                        %Determine the magnitude of the conduction velocity
                        
                        Times = ([zz(1,:),zz(2:end,end)',zz(end,1:end-1),zz(2:end-1,1)']-zz(ceil(size(zz,1)./2),ceil(size(zz,2)./2)));
                        Distancex = [xx(1,:),xx(2:end,end)',xx(end,1:end-1),xx(2:end-1,1)']-xx(ceil(size(xx,1)./2),ceil(size(xx,2)./2));
                        Distancey = [yy(1,:),yy(2:end,end)',yy(end,1:end-1),yy(2:end-1,1)']-yy(ceil(size(yy,1)./2),ceil(size(yy,2)./2));
                        
                        %                 Times     = [zz(array_area,[array_area:array_area+2]),zz(array_area+1,array_area+1),zz(array_area-1,array_area-1),zz(array_area+2,[array_area:array_area+2])]-zz(array_area+1,array_area+1);
                        %                 Distancex = [xx(array_area,[array_area:array_area+2]),xx(array_area+1,array_area+1),xx(array_area-1,array_area-1),xx(array_area+2,[array_area:array_area+2])]-xx(array_area+1,array_area+1);
                        %                 Distancey = [yy(array_area,[array_area:array_area+2]),yy(array_area+1,array_area+1),yy(array_area-1,array_area-1),yy(array_area+2,[array_area:array_area+2])]-yy(array_area+1,array_area+1);
                        %                     MaxTimes = max(Times);
                        locsMax = find(Times==max(Times));
                        
                        DistMaxTime = sqrt(Distancex(locsMax).^2 + Distancey(locsMax).^2);
                        locsDistMin = find(DistMaxTime == min(DistMaxTime));
                        locsMaxTimeMinDist = locsMax(locsDistMin);
                        if numel(locsMaxTimeMinDist) ~= 0
                            
                            maxY = Times(locsMaxTimeMinDist(1));
                            %                     a = Distancex(locsMaxTimeMinDist(1));
                            %                     b = Distancey(locsMaxTimeMinDist(1));
                            %                     x_of_max_plane_slope = sqrt(a^2+b^2);
                            x_of_max_plane_slope = DistMaxTime(locsDistMin(1));
                            u = Distancex(locsMaxTimeMinDist(1))/(maxY.*tresolution);
                            v = Distancey(locsMaxTimeMinDist(1))/(maxY.*tresolution);
                        else
                            maxY = abs(zz(1) - zz(size(zz,2)*size(zz,1)));
                            x_of_max_plane_slope =sqrt(size(zz,1)^2+size(zz,2)^2);
                            u = size(zz,1)/(maxY.*tresolution);
                            v = size(zz,2)/(maxY.*tresolution);
                        end
                        %                     maxY
                        %                     max_plane_slope = abs(maxY)/(x_of_max_plane_slope*sresolution); %sresolution mm/pixel
                        max_plane_slope = abs(maxY.*tresolution)/(x_of_max_plane_slope*sresolution); %sresolution mm/pixel
                        
                        Conduction_Velocity = 1/max_plane_slope;
                        CV (y_coord,x_coord) = Conduction_Velocity;
                        U(y_coord,x_coord) = u; %x velocity
                        V(y_coord,x_coord) = v; %y velocity
                        X(y_coord,x_coord) = x_coord; %x location
                        Y(y_coord,x_coord) = y_coord; %y location
                    end
                end
            end
        end
    end
end
cvsize = size(CV)
save (strcat(FullPath,'CV.mat'),'rtto','CV','ScalingOpt')
% toc(timeval)
close(wait_h)
CV = CV*100;
CV(CV==0) = NaN;
CV(CV>200) = NaN;
cvfig = figure;
histogram(CV,100)
title(strcat('General Histogram of Conduction Velocity values in the ROI','\newline','withought threshold(CV<200 cm/s)'))
xlabel('Conduction Velocity (cm/s)');
ylabel('Count (pixels)')
% CV(CV==0) = NaN;
% CV(CV>100) = 100;
% meanCV = nanmean(nanmean(CV));
filename = strcat(FullPath, 'CV.xlsx');
% cvsize=size(CV)
xlswrite(filename,CV,1)
    javaaddpath('C:/Yair/Utils/JExcelAPI/jxl.jar')
    javaaddpath('C:/Yair/Utils/JExcelAPI/MXL.jar')
    javaaddpath('poi_library/poi-3.8-20120326.jar');
    javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
    javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
    javaaddpath('poi_library/xmlbeans-2.3.0.jar');
    javaaddpath('poi_library/dom4j-1.6.1.jar');
    javaaddpath('poi_library/stax-api-1.0.1.jar');

    xlwrite(strcat(FullPath, 'ConductionVelocity.xls'),CV,1)

% CV(find(isnan(CV))) = 0;
% CV = imgaussfilt(CV,4);
% figure;imshow(imresize(CV,3),jet(ceil(max(max(CV)))));colorbar;title('Conduction Velocity (cm/s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Data,3)>2000
    save (strcat(FullPath,'CVallData.mat'),'Data','RawData','MASK','Mask','rrto_vectorplot','sresolution','tresolution','binsize',...
        'rtto','CV','ScalingOpt','-v7.3')
else
    save (strcat(FullPath,'CVallData.mat'),'Data','RawData','MASK','Mask','rrto_vectorplot','sresolution','tresolution','binsize',...
        'rtto','CV','ScalingOpt')
end
% MASK = Mask;
% sresolution=10/235;
% tresolution = 2099/2000;
dwindow = 5;

%%%%%%%%%%%%%% Test here: %%%%%%%%%%%%%%%
% CV(CV>100) = 100;
h = fspecial('gaussian',[dwindow dwindow],2);
CV = imfilter(CV,h,'replicate','same');
CVsize = size(CV)
if vectorfield
    %% Vector field plot
%     MASKsize = size(MASK)
%     if size(Mask,1)>200
%         MASK = imresize(MASK,0.5);
%         CV = imresize(CV,0.5);
%         rrto_vectorplot = imresize(rrto_vectorplot,0.5);
%     end
%     CVsize = size(CV)

    %% CV measurement
    CV = Scale2D(CV,binsize);
    rrto_vectorplot = Scale2D(rrto_vectorplot,binsize);
    MASK = Scale2D(MASK,binsize);
    FrameYLim=[1,size(MASK,1)];
    FrameXLim=[1,size(MASK,2)];
    SIGNAL = MASK;
    ImgScale = 5;
    CV(CV>40) = 40;

    plotvectors=1;
    normvectorscale=30;
    vectorscale=0.1;
    plotsumvector=1;
    sumvectorscale=10;
    
    
    %%%
    %%%
    %vector size and linewidth
    vectorwidth=1;%[pts] vector line width [2.5]
    sumvectorwidth=4;
    unitvectorwidth=4;
    normsumvector=0;%normalizes sum vector so that all sumvectors have the same length
    
    %COLORS
    framecolor=[0 0 0];%color of frame
    frameweight=3;
    vectorcolor=[1 1 0];
    sumvectorcolor=[1 0 0];
    unitvectorcolor=[0.5 0.5 0.5];
    
   
    
    %%%
    %%%
    [VFRAMEi,VFRAMEj,RMSEFRAME,XFRAMEINTERP_0,YFRAMEINTERP_0,SIGNALFRAMEINTERP_0]=fitpixels3BM(rrto_vectorplot,FrameXLim,FrameYLim,50,MASK);
    %VFRAMEi: conduction velocity in i-direction (line) [pixel/frame]
    %VFRAMEj: conduction velocity in j-direction (column) [pixel/frame]
    %XFRAMEINTERP_0: X coordinates of frame used for CV calculation
    %YFRAMEINTERP_0: Y coordinates of frame
    
%     VFRAMEi=imresize(VFRAMEi,0.5);
%     VFRAMEj=imresize(VFRAMEj,0.5);
%     RMSEFRAME=imresize(RMSEFRAME,0.5);
%     XFRAMEINTERP_0=imresize(XFRAMEINTERP_0,0.5);
%     YFRAMEINTERP_0=imresize(YFRAMEINTERP_0,0.5);
%     SIGNALFRAMEINTERP_0=imresize(SIGNALFRAMEINTERP_0,0.5);
    %% Litte fit statistics
    %framepixels: total number of active pixels in frame
    framepixels=0;
    rmsepixels=0;
    maxrmse = 0.75;
    RMSELIST=[];%to determine smallest and largest RMSE
    for i=1:size(SIGNALFRAMEINTERP_0,1)
        for j=1:size(SIGNALFRAMEINTERP_0,2)
            if SIGNALFRAMEINTERP_0(i,j)==1
                framepixels=framepixels+1;
                if 0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse
                    rmsepixels=rmsepixels+1;
                    RMSELIST=[RMSELIST;RMSEFRAME(i,j)];
                end
            end
        end
    end
    SIGNALFRAME=SIGNAL(FrameYLim(1):FrameYLim(end),FrameXLim(1):FrameXLim(end));%original SIGNALS
    framearea=floor(bwarea(SIGNALFRAME))*sresolution*sresolution;%estimate of tissue area in frame
    area=floor(bwarea(SIGNAL))*sresolution*sresolution;%estimate total excitable tissue area
    fprintf(['number of active pixels available for fit: ',num2str(framepixels),'\n']);
    fprintf(['active pixels with RMSE<',num2str(maxrmse),': ',num2str(rmsepixels),'\n']);
    fprintf(['max RMSE: ',num2str(max(RMSELIST)),'\n']);
    fprintf(['min RMSE: ',num2str(min(RMSELIST)),'\n']);
    %% determine ROI polygon coordinates in interpolated frame space
    %make sure that POLYROI from last run is wiped out
    
    %% calibrate velocities
    %DEFINITIONS: VELVECTOR=[VFRAMEj(i,j);VFRAMEi(i,j)];
    % CALVFRAMEi=VFRAMEi*sresolution*10^3/tresolution;%[mm/s]
    % CALVFRAMEj=VFRAMEj*sresolution*10^3/tresolution;%[mm/s]
    %% calibrate velocities
    %DEFINITIONS: VELVECTOR=[VFRAMEj(i,j);VFRAMEi(i,j)];
    CALVFRAMEi=VFRAMEi*(sresolution*0.1)/(tresolution*1.0e-3);%[cm/s]
    CALVFRAMEj=VFRAMEj*(sresolution*0.1)/(tresolution*1.0e-3);%[cm/s]
    %% filter vectors
    %determine mean and standard deviation of all velocities with RMSE<maxrmse
    VLIST=[];
    for i=1:size(SIGNALFRAMEINTERP_0,1)
        for j=1:size(SIGNALFRAMEINTERP_0,2)
            if SIGNALFRAMEINTERP_0(i,j)==1 && (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
                VLIST=[VLIST;norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2)];
            end
        end
    end
    meanvelall=mean(VLIST);
    stdvelall=std(VLIST);
    %filter vectors
    VELINTERVAL=[meanvelall-5*stdvelall,meanvelall+5*stdvelall];
    VPLOTBIN=zeros(size(SIGNALFRAMEINTERP_0));%binary matrix of filtered pixels
    VPLOTLIST=[];%list of filtered velocities
    CVplot=[];
    for i=1:size(SIGNALFRAMEINTERP_0,1)
        for j=1:size(SIGNALFRAMEINTERP_0,2)
            if SIGNALFRAMEINTERP_0(i,j)==1
%                 if (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
                    v=norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2);
                    %if VELINTERVAL(1)<v && v<VELINTERVAL(2)
                    if v>0
                        VPLOTBIN(i,j)=1;
                        VPLOTLIST=[VPLOTLIST;v];
                        CVplot(i,j) = v;
                    end
%                 end
            end
        end
    end
    
    
    % %show numerical results
    % meanvelocity=mean(VPLOTLIST);
    % stdvelocity=std(VPLOTLIST);
    %
    % fprintf(['mean velocity of selected vectors: (',num2str(meanvelocity),' +- ',num2str(stdvelocity/sqrt(length(VPLOTLIST))),') mm/s \n']);
    % fprintf(['number of selected vectors: ',num2str(length(VPLOTLIST)),' out of total number: ',num2str(length(VLIST)),'\n']);
    % % fprintf(['window size (w x h) :  ',num2str(frameumwidth),' x ',num2str(frameumheight),' um \n']);
    % %% plot velocity vectors
    
    CALVFRAMEi2=CALVFRAMEi;
    CALVFRAMEi2(CALVFRAMEi==0)=NaN;
    CALVFRAMEi2=(CALVFRAMEi2-nanmin(nanmin(CALVFRAMEi2))./nanmax(nanmax(CALVFRAMEi-nanmin(nanmin(CALVFRAMEi2)))));
    CALVFRAMEj2=CALVFRAMEj;
    CALVFRAMEj2(CALVFRAMEj==0)=NaN;
    CALVFRAMEj2=(CALVFRAMEj2-nanmin(nanmin(CALVFRAMEj2))./nanmax(nanmax(CALVFRAMEj-nanmin(nanmin(CALVFRAMEj2)))));
    % CALVFRAMEi(abs(CALVFRAMEi)>100)=100*sign(CALVFRAMEi(abs(CALVFRAMEi)>100));
    % CALVFRAMEj(abs(CALVFRAMEj)>100)=100*sign(CALVFRAMEj(abs(CALVFRAMEj)>100));
    % V = sqrt(CALVFRAMEi.^2+CALVFRAMEj.^2);
    % mean(mean(V))
    % max(max(V))
    % min(min(V))
    % median(median(V))
    % figure;quiver(YFRAMEINTERP_0,XFRAMEINTERP_0,CALVFRAMEi2,CALVFRAMEj2)
%     save ('test.mat','CALVFRAMEi','CALVFRAMEj','VPLOTLIST')
%     size(CV)
%     CV = interp2(CV);
%     size(CV)
    wait1=waitbar(0,'Generating the map...');
    figure(100);  cla;
    tempframe = mat2gray(double(imresize(RawData(:,:,5),ImgScale)));
            flipoption = get(handles.checkbox29,'Value');
% flipoption=1;
            if flipoption == 1
                tempframe = fliplr(tempframe);
            end
%             if size(RawData,1)>200
%                 tempframe = imresize(tempframe,2);
%             end
    figure(100); imshow(tempframe);
    hold on
    for i=1:1:size(VPLOTBIN,1)
        waitbar(i/size(VPLOTBIN,1))
        for j=1:1:size(VPLOTBIN,2)
            if VPLOTBIN(i,j)==1
%                 x=XFRAMEINTERP_0(1,j)-floor(binsize/2);
                x=XFRAMEINTERP_0(1,j)-ceil(binsize/2);
                y=YFRAMEINTERP_0(i,1)-binsize;
                vx=CALVFRAMEj(i,j);
                vy=CALVFRAMEi(i,j);
%                             vxn=vx/norm([vx,vy],2);
%                             vyn=vy/norm([vx,vy],2);
%                             vxn=vx/max(max(abs(V)));
%                             vyn=vy/max(max(abs(V)));
%                 vxn=vx/100;
%                             vyn=vy/100;
                vxn=(vx/norm([vx,vy],2))*(CV(i,j)/40);
                vyn=(vy/norm([vx,vy],2))*(CV(i,j)/40);
%                 figure(100); hold on;
                if normvectorscale==0
                    %plot vector [vx,vy] at position [x,y] NOT NORMALIZED
                    quiver(x*ImgScale*binsize,y*ImgScale*binsize,vx,vy,vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',unitvectorwidth);
                else
                    %plot vector [vx,vy] at position [x,y] NORMALIZED
                    quiver(x*ImgScale*binsize,y*ImgScale*binsize,vxn,vyn,normvectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',unitvectorwidth);
                end
            end
        end
    end
    hold off
    close(wait1)
    
    % [I,J]=find(MASK==1);
    % FRAMECENTER = [(J(1)+J(end))/2,(I(1)+I(end))/2];
    % %% plot mean vector
    % if plotsumvector==1;%plot sum velocity vector
    %     %calculate mean velocity vector
    %     VECTORLISTxy=[];%to calculate mean vector
    %     for i=1:size(VPLOTBIN,1)
    %         for j=1:size(VPLOTBIN,2)
    %             if VPLOTBIN(i,j)==1
    %                 VECTORLISTxy=[VECTORLISTxy;[CALVFRAMEj2(i,j),CALVFRAMEi2(i,j)]];
    %             end
    %         end
    %     end
    %     %     MEANVECTOR=[sum(VECTORLISTxy(:,1)),sum(VECTORLISTxy(:,2))]./length(VECTORLISTxy);
    %     MEANVECTOR=nanmean(VECTORLISTxy);
    %     MEANVECTOR=MEANVECTOR*ImgScale*binsize./(norm(MEANVECTOR))
    %     %plot sum vector into main isochronal map
    %     figure(1);hold on
    %     R=FRAMECENTER;
    %     quiver(R(1)*ImgScale*binsize,R(2)*ImgScale*binsize,MEANVECTOR(1),MEANVECTOR(2),sumvectorscale,'LineWidth',sumvectorwidth,'Color',sumvectorcolor,'MaxHeadSize',1);
    %
    % end
    
    
end


% CV(CV>=MaxCV) = 0;
% CV(CV<=MinCV)  = 0;
% CV(CV==0)   = NaN;
% MEAN_CV_bounded = nanmean(nanmean(CV))
% % CV(CV>=80) = 0;
% % CV(CV<=10) = 0;
% % CV(CV==0)  = NaN;
% % MEAN_CV_10_80  = nanmean(nanmean(CV))
% set(handles.text109,'String',[num2str(MEAN_CV_bounded), 'cm/s']);
% % set(handles.text55,'String',[CLstr, ' ms']);
% 
% % save ('cvtestdata.mat','CV')
% % save ('cvtestdata.mat','CV')
% 
% % save ('CVRTTO.mat','rtto','CV')

    fprintf('Conduction Velocity Measurement Complete!\n');
% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat']);
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end
% size(Data)
total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
ScalingOpt = get(handles.checkbox30,'Value');
% AVERAGING = 1;
AVERAGING=get(handles.checkbox2,'Value');
if AVERAGING 
    Starting_frame = get(handles.edit35,'String');
    Starting_frame = str2double(Starting_frame);
    Ending_frame = get(handles.edit36,'String');
    Ending_frame = str2double(Ending_frame);
    if Ending_frame>size(Data,3)
        Ending_frame =size(Data,3);
    end
    
    Sig = squeeze(nanmean(nanmean(Data(:,:,Starting_frame:Ending_frame),1),2));
    Sig = smooth(squeeze(Sig),0.01);
    Sig = 100*(Sig-min(Sig))/(max(Sig)-min(Sig));
    time=1:size(Sig,1);
    
    %     Sig = squeeze(nanmean(nanmean(Data(:,:,Starting_frame:Ending_frame))));
    %     Sig=Sig-nanmin(Sig);
    %     Sig = Sig'/nanmax(Sig)*100;
    %         size(Sig)
    %     time=1:size(Sig,2);
    [I,J]=findpeaks(Sig,time,'MinPeakProminence',20,'Annotate','extents');
    figure;findpeaks(Sig,time,'MinPeakProminence',20,'Annotate','extents');
    AData=[];
    k=0;
    CL=diff(J);
    Cycle_length = mean(CL)*tresolution
    halfCL = floor(mean(diff(J))./2);
    for j=1:length(J)
        if (J(j)-halfCL > 0) && (J(j)+halfCL<length(Sig))
            k=k+1;
            AData(:,:,:,k) = Data(:,:,(Starting_frame+floor(J(j)-halfCL):Starting_frame+floor(J(j)+halfCL)));
        end
    end
    AveragedData = nanmean(AData,4);
    AveragedData = mat2gray(AveragedData);
    figure;plot( squeeze(nanmean(nanmean(AData,1),2)));
    hold on; plot(squeeze(nanmean(nanmean(AveragedData,1),2)),'k','Linewidth',2); legend;
    
    %
    %     [Min, tMin, IMin, Max, tMax, IMax, Npeaks] = fpeaks(time, Sig, 'data');
    %     %     figure;plot(time,Sig)
    %     CL=nanmean(diff(tMax));
    %     %                 CL=CL/tresolution;
    %     Y=[];
    %     cycles = 0;
    %     CL=ceil(CL)
    %     wait1=waitbar(0,'Please wait...');
    %     for l=1:length(IMax)
    %         if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(Sig))
    %             Y(:,:,:,l)=Data(:,:,IMax(l)-floor(CL./2):IMax(l)+floor(CL./2));
    %             cycles = cycles +1;
    %         end
    %     end
    %
    %     AveragedData = nanmean(Y,4);
    %     size(AveragedData)
    %     AveragedData(find(isnan(AveragedData)))=0;
    %     AveragedData =mat2gray(AveragedData);
    %     AveragedData = AveragedData .*Mask;
    %     AveragedData =mat2gray(AveragedData);
    if size(Data,3)>2000
        save (strcat(FullPath,'.mat'),'RawData','Data','Mask','AveragedData','ScalingOpt','-v7.3');
    else
        save (strcat(FullPath,'.mat'),'RawData','Data','Mask','AveragedData','ScalingOpt');
    end
    %     close(wait1)
    set(handles.edit93,'String',num2str(size(AveragedData,3)));
    set(handles.edit92,'String',num2str(1));
    fprintf('Averaging Complete!\n');
else
       figure; plot(squeeze(nanmean(nanmean(Data,1),2)),'k','Linewidth',2);

end

function edit85_Callback(hObject, eventdata, handles)
% hObject    handle to edit85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit85 as text
%        str2double(get(hObject,'String')) returns contents of edit85 as a double


% --- Executes during object creation, after setting all properties.
function edit85_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
set(handles.checkbox2,'Value',1);
set(handles.checkbox23,'Value',0);
RoiOption = get(handles.popupmenu8,'Value'); %% 1 for signel pixel - 2 for impoly
AVERAGING=get(handles.checkbox2,'Value')
manual_mode=1;
% Averaging=1;
% APDT=50 ;     % choose 50 for APD50, 70 for APD70, and 90 for APD90
% manual_RMP=0;
% RMP=6.2;
% t_RMP=88.82;

total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp)
Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
% xlim = [1,Nframes];
NumFrame = get(handles.edit62,'String');
NumFrame = str2double(NumFrame);


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
    load ([fullfile(PName,'results2','corrected','results',FName),'_improved.mat'])
    h = Data(:,:,NumFrame);
    h= imresize(h,2);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
    h = mat2gray(squeeze(RawData(:,:,50)));
    h= imresize(h,2);
    flipoption = get(handles.checkbox29,'Value');
    if flipoption == 1
        h = fliplr(h);
    end
end
% i0=findtime(frames,NumFrame)
% A = strcat([fullfile(PName, 'results2','corrected', i0), num2str(NumFrame), '.tif']);
% h = imread(A);
% if series
%     load ([fullfile(PName,'results2','corrected','results',FName),'_improved.mat'])
%     h = Data(:,:,NumFrame);
%     h= imresize(h,2);
%     flipoption = get(handles.checkbox29,'Value');
%     if flipoption == 1
%         h = fliplr(h);
%     end
% elseif stack
%     h = mat2gray(squeeze(RawData(:,:,50)));
%     h= imresize(h,2);
%     flipoption = get(handles.checkbox29,'Value');
%     if flipoption == 1
%         h = fliplr(h);
%     end
% end

if AVERAGING
    Data = AveragedData;
else
    Data = Data(:,:,startp:endp);
end
Data = mat2gray(Data);
%% Define a ROI
d2 = figure;
if RoiOption == 1
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    Sig=[];
    if manual_mode
        [x_coord,y_coord,intensity_val] = impixel(h);
        x_coord = floor(x_coord./2); y_coord = floor(y_coord./2);
        x= Data (y_coord,x_coord,:);
    else
        x= Data (29,25,startp:endp);
    end
    close(d2)
    x=squeeze(x);
% %     x= smooth(x,0.01);
%     frame_rate=1000*total_frames/total_time;
%     Signal=[];
%     % % smooth x data
%     Signal(1,:) = x;
%     % Signal(1,:) = nanmean(Sig);
%     S2 = smooth(1:size(x),Signal,0.9,'rloess');
%     Signal = Signal-S2'; % drift removal
%     window = 5;
%     mask = ones(1,window)/window;
%     maY=conv(Signal,mask,'same'); % temporal smoothing
%     x=maY(1,:);
    x= smooth(x,0.01);

    OAP_signal=x-nanmin(x);
    % %
        t=1:size(OAP_signal,1);

%     t=1:size(OAP_signal,2);
    t=t*tresolution;
    I=t/tresolution;
    OAP = OAP_signal;
    tt  = t;
    II  = I;
    
%     [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
%     CL=nanmean(diff(tMax))
%     CL=CL/tresolution;
%     
%     if Npeaks>2 && AVERAGING
%         Y=[];
%         cycles = 0;
%         CL=ceil(CL);
%         for l=1:length(IMax)
%             if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP))
%                 Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
%                 cycles = cycles +1;
%             end
%         end
%         
%         OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
%         tt_new=tt(1:length(OAP_new));
%         II_new=II(1:length(OAP_new));
%     else
%         OAP_new = OAP;
%         tt_new=tt;
%         II_new=II;
%     end
        OAP_new = OAP;
        tt_new=tt;
        II_new=II;
        OAP_new=OAP_new-nanmin(OAP_new);
%     OAP_new = OAP_new/nanmax(OAP_new)*100;
%     OAP_new = smooth(OAP_new,0.08);
% size(x)
% figure;plot(x-nanmin(x));
    figure; plot(tt_new,OAP_new,'k','linewidth',1);xlabel('time (ms)');
    hold on; plot([0, tt_new(end)],[0.1, 0.1] , 'k--')
%     axislimits=nanmax(OAP_new)-nanmin(OAP_new);
%     axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
    
    
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
else
    ansSatisfied='NO';
    for n=1:1000
        if strcmp(ansSatisfied,'YES')
            pos = getPosition (hpoly);
            h = mat2gray(squeeze(RawData(:,:,50)));
            h= imresize(h,3);
            flipoption = get(handles.checkbox29,'Value')
            if flipoption == 1
                h = fliplr(h);
            end
            figure(d2); imshow(mat2gray(h.*imresize(Mask,3)));
            hold on; plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],'g','LineWidth',2)
            break
        elseif strcmp(ansSatisfied,'Cancel!')
            break
        elseif strcmp(ansSatisfied,'NO')
            h = mat2gray(nanmax(RawData,[],3));
            h= imresize(h,3);
            flipoption = get(handles.checkbox29,'Value')
            if flipoption == 1
                h = fliplr(h);
            end
            figure(d2); imshow(mat2gray(h.*imresize(Mask,3)));
            hpoly = impoly(gca);
            MASK = createMask(hpoly);
            MASK = imresize(MASK,1/3);
            MASK = double(MASK);
            MASK( MASK == 0 ) = NaN;
        end
        ansSatisfied = questdlg('Generate OAP from this ROI?','MASK', 'YES', 'NO', 'Cancel!', 'NO');
    end
    close(d2)
    DataOAP = Data .* MASK;
    DataOAP = Data .* MASK;
    x = squeeze(nanmean(nanmean(DataOAP)));
    frame_rate=1000*total_frames/total_time;
%     Signal=[];
%     % % smooth x data
%     Signal(1,:) = x
%     % Signal(1,:) = nanmean(Sig);
%     S2 = smooth(1:size(x),Signal,0.9,'rloess');
%     Signal = Signal-S2'; % drift removal
%     window = 5;
%     mask = ones(1,window)/window;
%     maY=conv(Signal,mask,'same'); % temporal smoothing
%     x=maY(1,:);
    x= smooth(x,0.01);
    OAP_signal=x-nanmin(x);
    % %
    t=1:size(OAP_signal,1);
    t=t*tresolution;
    I=t/tresolution;
    OAP = OAP_signal;
    tt  = t;
    II  = I;
    
%     [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(tt, OAP, 'data');
%     CL=nanmean(diff(tMax))
%     CL=CL/tresolution;
%     if Npeaks>2 && AVERAGING
%         Y=[];
%         cycles = 0;
%         CL=ceil(CL);
%         for l=1:length(IMax)
%             if (IMax(l)-floor(CL./2) > 0) & (IMax(l)+floor(CL./2)<length(OAP))
%                 Y=[Y;OAP(IMax(l)-floor(CL./2):IMax(l)+floor(CL./2))];
%                 cycles = cycles +1;
%             end
%         end
%         
%         OAP_new = repmat(nanmean(Y)-nanmin(nanmean(Y)),1,2);
%         tt_new=tt(1:length(OAP_new));
%         II_new=II(1:length(OAP_new));
%     else
%         OAP_new = OAP;
%         tt_new=tt;
%         II_new=II;
%     end
        OAP_new = OAP;
        tt_new=tt;
        II_new=II;
    OAP_new=OAP_new-nanmin(OAP_new);
%     OAP_new = OAP_new/nanmax(OAP_new)*100;
%     OAP_new = smooth(OAP_new,0.08)
    figure; plot(tt_new,OAP_new,'k','linewidth',1);xlabel('time (ms)');
    hold on; plot([0, tt_new(end)],[0.1, 0.1] , 'k--')
%     axislimits=nanmax(OAP_new)-nanmin(OAP_new);
%     axis([tt_new(1) tt_new(end) nanmin(OAP_new)-(5*axislimits)/100 nanmax(OAP_new)+(5*axislimits)/100])
    
end




% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox27


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end

% AVERAGING = 1;
AVERAGING=get(handles.checkbox2,'Value');
if AVERAGING
    Data = AveragedData;
end


total_frames = get(handles.edit46,'String');
total_frames = str2double(total_frames);
total_time = get(handles.edit47,'String');
delay = get(handles.edit91,'String');
total_time = str2double(total_time)-str2double(delay)+1;
tresolution=total_time/total_frames;
total_pixels = get(handles.edit52,'String');
total_pixels = str2double(total_pixels);
if get(handles.checkbox30,'Value')==1
    total_pixels = total_pixels/2;
end
total_distance = get(handles.edit53,'String');
total_distance = str2double(total_distance);
sresolution=total_distance/total_pixels;

Nframes = get(handles.edit40,'String');
Nframes = str2double(Nframes);
startp = get(handles.edit35,'String');
startp = str2double(startp);
endp = get(handles.edit36,'String');
endp = str2double(endp);

xlim = [1,size(Data,3)];
% FirstFrame = get(handles.edit62,'String');
% FirstFrame = str2double(FirstFrame);
% OAP_chkbx_status =0;
% array = [];
%Record Optical Action Potentials
% pt = 1;
% array_area = 2;
Sig = squeeze(nanmean(nanmean(Data)))';
[Min, tMin, IMin, Max, tMax, IMax, NPEAKS]=fpeaks([1:size(Sig,2)],Sig, 'data');

array_area = floor(str2num(get(handles.edit85,'String'))/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RoiOption = get(handles.popupmenu7,'Value'); %% 1 for freehand - 2 for impoly - 3 for imellipse - 4 for imrect
%
% %% Define a ROI
% d1 = figure;
% if RoiOption == 1
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
%     MASK = Mask;
%     figure(d1); imshow(mat2gray(double(imresize(RawData(:,:,5),3))).*imresize(Mask,3));
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
% else
%     ansSatisfied='NO';
%     for n=1:1000
%         if strcmp(ansSatisfied,'YES')
%             pos = getPosition (h);
%             figure(d1); imshow(mat2gray(double(imresize(RawData(:,:,5),3))).*imresize(Mask,3));
%             hold on; plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],'g','LineWidth',2)
%             break
%         elseif strcmp(ansSatisfied,'Cancel!')
%             break
%         elseif strcmp(ansSatisfied,'NO')
%             figure(d1); imshow(mat2gray(double(imresize(RawData(:,:,5),3))).*imresize(Mask,3));
%             h = impoly(gca);
%             MASK = createMask(h);
%             MASK = imresize(MASK,1/3);
%             MASK = double(MASK);
%             MASK( MASK == 0 ) = NaN;
%         end
%         ansSatisfied = questdlg('Calulate CV in this ROI?','MASK', 'YES', 'NO', 'Cancel!', 'NO');
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binsize = 1;
CVData = imresize(mat2gray(Data),1/binsize);
MASK = imresize(Mask,1/binsize);
MASK = MASK;
% if RoiOption == 1
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
%     MASK = MASK;
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%
% else
% MASK = imresize(MASK,1/binsize);
% end
wait1=waitbar(0,'Please wait...');
% figure(1); imshow(mat2gray(squeeze(RawData(:,:,50))));
for x = 1:size(CVData,1)
    for y = 1:size(CVData,2)
        if ~isnan(MASK(x,y))
            OAP = squeeze(CVData(x,y,:))-nanmin(squeeze(CVData(x,y,:)));
            OAP = OAP';
            %             OAP = smooth(OAP,5,'sgolay',1);
            I = 1 : length(OAP);
            t = I * tresolution;
            AverageOAP = squeeze(mean(mean(CVData)))-nanmin(squeeze(mean(mean(CVData))));
            I_AverageOAPth =find(AverageOAP>=0.1);
            OAP(1:I_AverageOAPth(1)-20)=0;
            OAPinterp = interp(OAP,10);
            OAPinterp(1:2) = 0;
            Iinterp = interp(I,10);
            I_OAPth =find(OAPinterp>=0.1);
            if ~isempty(I_OAPth)
                if Iinterp(I_OAPth(1))~=1
                    Iinterp_OAPth =  Iinterp(I_OAPth(1));
                elseif length(I_OAPth>2)
                    Iinterp_OAPth =  Iinterp(I_OAPth(2));
                end
                rrto_vectorplot(x,y)=Iinterp_OAPth;
            end
            %             figure(1) ; hold on ; plot(binsize*y , binsize*x , '*');
            %             figure(2) ; cla;
            %             figure(2) ; plot(t , OAP)
            %             [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, OAP, 'data');
            %             if numel(Min)~=0
            %                 D=diff(OAP);
            %                 S=sign(D);
            %                 SD=diff(S);
            %                 SDS=sign(SD);
            %                 SDS=[0,SDS,0];
            %                 [SMin, StMin, SIMin, SMax, StMax, SIMax, SNpeaks]=fpeaks(t, SDS, 'model');
            %                 [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks(t, OAP, 'data');
            %                 A1=[];B1=[];B2=[];
            %                 for k3=1:Npeaks
            %                     a1=find(SIMin==IMax(k3));
            %                     b1=find(SIMax<=IMax(k3));
            %                     if numel(b1)~=0
            %                         b1=b1(end);
            %                     end
            %                     b2=find(SIMax>=IMax(k3));
            %                     if numel(b2~=0)
            %                         b2=b2(1);
            %                     end
            %                     A1=[A1,a1];B1=[B1,b1];B2=[B2,b2];
            %                 end
            %                 a1=A1;b1=B1;b2=B2;
            %                 APA   = OAP(SIMin(a1));
            %                 t_APA = t(SIMin(a1));
            %                 I_APA = SIMin(a1);
            %                 APAval = max(APA);
            %                 APAtime = t_APA(APA == APAval);
            %                 figure(2);hold on; plot(APAtime,APAval,'or')
            %                 RMP   = OAP(SIMax(b2));
            %                 t_RMP = t(SIMax(b2));
            %                 I_RMP = SIMax(b2);
            %                 figure(2);hold on; plot(t_RMP,RMP,'og')
            %                 TP = OAP(SIMax(b1));
            %                 t_TP = t(SIMax(b1));
            %                 I_TP = SIMax(b1)
            %                 RTind = find(t_TP<APAtime(1));
            %                 if ~isempty(RTind)
            %                     RT = t_TP(RTind(end));
            %                     RTval = TP((RTind(end)));
            %                     figure(2);hold on; plot(RT,RTval,'om')
            %                     if numel(TP)==numel(APA)
            %                         %                     rt(x,y)=RT(1);
            %                         %                       rrto_vectorplot(x,y)=I_APA(1);
            %                         rrto_vectorplot(x,y)=RTind(end);
            %                     end
            %                 end
            %             end
            %                     pause(1)
        end
    end
end
close(wait1)
% save ('CVrt.mat','rrto_vectorplot')

%% Vector field plot

%% CV measurement
FrameYLim=[1,size(MASK,1)];
FrameXLim=[1,size(MASK,2)];
SIGNAL = MASK;
ImgScale = 10;

plotvectors=1;
normvectorscale=30;
vectorscale=0.1;
plotsumvector=1;
sumvectorscale=10;


%%%
%%%
%vector size and linewidth
vectorwidth=1;%[pts] vector line width [2.5]
sumvectorwidth=4;
unitvectorwidth=4;
normsumvector=0;%normalizes sum vector so that all sumvectors have the same length

%COLORS
framecolor=[0 0 0];%color of frame
frameweight=3;
vectorcolor=[1 1 0];
sumvectorcolor=[1 0 0];
unitvectorcolor=[0.5 0.5 0.5];

%%%
%%%
[VFRAMEi,VFRAMEj,RMSEFRAME,XFRAMEINTERP_0,YFRAMEINTERP_0,SIGNALFRAMEINTERP_0]=fitpixels3BM(rrto_vectorplot,FrameXLim,FrameYLim,50,MASK);
%VFRAMEi: conduction velocity in i-direction (line) [pixel/frame]
%VFRAMEj: conduction velocity in j-direction (column) [pixel/frame]
%XFRAMEINTERP_0: X coordinates of frame used for CV calculation
%YFRAMEINTERP_0: Y coordinates of frame
%% Litte fit statistics
%framepixels: total number of active pixels in frame
framepixels=0;
rmsepixels=0;
maxrmse = 0.75;
RMSELIST=[];%to determine smallest and largest RMSE
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1
            framepixels=framepixels+1;
            if 0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse
                rmsepixels=rmsepixels+1;
                RMSELIST=[RMSELIST;RMSEFRAME(i,j)];
            end
        end
    end
end
SIGNALFRAME=SIGNAL(FrameYLim(1):FrameYLim(end),FrameXLim(1):FrameXLim(end));%original SIGNALS
framearea=floor(bwarea(SIGNALFRAME))*sresolution*sresolution;%estimate of tissue area in frame
area=floor(bwarea(SIGNAL))*sresolution*sresolution;%estimate total excitable tissue area
fprintf(['number of active pixels available for fit: ',num2str(framepixels),'\n']);
fprintf(['active pixels with RMSE<',num2str(maxrmse),': ',num2str(rmsepixels),'\n']);
fprintf(['max RMSE: ',num2str(max(RMSELIST)),'\n']);
fprintf(['min RMSE: ',num2str(min(RMSELIST)),'\n']);
%% determine ROI polygon coordinates in interpolated frame space
%make sure that POLYROI from last run is wiped out

%% calibrate velocities
%DEFINITIONS: VELVECTOR=[VFRAMEj(i,j);VFRAMEi(i,j)];
% CALVFRAMEi=VFRAMEi*sresolution*10^3/tresolution;%[mm/s]
% CALVFRAMEj=VFRAMEj*sresolution*10^3/tresolution;%[mm/s]
%% calibrate velocities
%DEFINITIONS: VELVECTOR=[VFRAMEj(i,j);VFRAMEi(i,j)];
CALVFRAMEi=VFRAMEi*(sresolution*0.1)/(tresolution*1.0e-3);%[cm/s]
CALVFRAMEj=VFRAMEj*(sresolution*0.1)/(tresolution*1.0e-3);%[cm/s]
%% filter vectors
%determine mean and standard deviation of all velocities with RMSE<maxrmse
VLIST=[];
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1 && (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
            VLIST=[VLIST;norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2)];
        end
    end
end
meanvelall=mean(VLIST);
stdvelall=std(VLIST);
%filter vectors
VELINTERVAL=[meanvelall-5*stdvelall,meanvelall+5*stdvelall];
VPLOTBIN=zeros(size(SIGNALFRAMEINTERP_0));%binary matrix of filtered pixels
VPLOTLIST=[];%list of filtered velocities
CVplot=[];
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1
            if (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
                v=norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2);
                %if VELINTERVAL(1)<v && v<VELINTERVAL(2)
                if v>0
                    VPLOTBIN(i,j)=1;
                    VPLOTLIST=[VPLOTLIST;v];
                    CVplot(i,j) = v;
                end
            end
        end
    end
end


% %show numerical results
% meanvelocity=mean(VPLOTLIST);
% stdvelocity=std(VPLOTLIST);
%
% fprintf(['mean velocity of selected vectors: (',num2str(meanvelocity),' +- ',num2str(stdvelocity/sqrt(length(VPLOTLIST))),') mm/s \n']);
% fprintf(['number of selected vectors: ',num2str(length(VPLOTLIST)),' out of total number: ',num2str(length(VLIST)),'\n']);
% % fprintf(['window size (w x h) :  ',num2str(frameumwidth),' x ',num2str(frameumheight),' um \n']);
% %% plot velocity vectors

CALVFRAMEi2=CALVFRAMEi;
CALVFRAMEi2(CALVFRAMEi==0)=NaN;
CALVFRAMEi2=(CALVFRAMEi2-nanmin(nanmin(CALVFRAMEi2))./nanmax(nanmax(CALVFRAMEi-nanmin(nanmin(CALVFRAMEi2)))));
CALVFRAMEj2=CALVFRAMEj;
CALVFRAMEj2(CALVFRAMEj==0)=NaN;
CALVFRAMEj2=(CALVFRAMEj2-nanmin(nanmin(CALVFRAMEj2))./nanmax(nanmax(CALVFRAMEj-nanmin(nanmin(CALVFRAMEj2)))));
CALVFRAMEi(abs(CALVFRAMEi)>100)=100*sign(CALVFRAMEi(abs(CALVFRAMEi)>100));
CALVFRAMEj(abs(CALVFRAMEj)>100)=100*sign(CALVFRAMEj(abs(CALVFRAMEj)>100));
V = sqrt(CALVFRAMEi.^2+CALVFRAMEj.^2);
% mean(mean(V))
% max(max(V))
% min(min(V))
% median(median(V))
% figure;quiver(YFRAMEINTERP_0,XFRAMEINTERP_0,CALVFRAMEi2,CALVFRAMEj2)
% save ('test.mat','CALVFRAMEi','CALVFRAMEj','VPLOTLIST')
wait1=waitbar(0,'Generating the map...');
figure(100);  cla;
figure(100); imshow(mat2gray(double(imresize(RawData(:,:,5),ImgScale))));
hold on
for i=1:1:size(VPLOTBIN,1)
    waitbar(i/size(VPLOTBIN,1))
    for j=1:1:size(VPLOTBIN,2)
        if VPLOTBIN(i,j)==1
            x=XFRAMEINTERP_0(1,j);
            y=YFRAMEINTERP_0(i,1);
            vx=CALVFRAMEj(i,j);
            vy=CALVFRAMEi(i,j);
            %             vxn=vx/norm([vx,vy],2);
            %             vyn=vy/norm([vx,vy],2);
            %             vxn=vx/max(max(abs(V)));
            %             vyn=vy/max(max(abs(V)));
            vxn=vx/100;
            vyn=vy/100;
            figure(100); hold on;
            if normvectorscale==0
                %plot vector [vx,vy] at position [x,y] NOT NORMALIZED
                quiver(x*ImgScale*binsize,y*ImgScale*binsize,vx,vy,vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',unitvectorwidth);
            else
                %plot vector [vx,vy] at position [x,y] NORMALIZED
                quiver(x*ImgScale*binsize,y*ImgScale*binsize,vxn,vyn,normvectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',unitvectorwidth);
            end
        end
    end
end
hold off
close(wait1)

% [I,J]=find(MASK==1);
% FRAMECENTER = [(J(1)+J(end))/2,(I(1)+I(end))/2];
% %% plot mean vector
% if plotsumvector==1;%plot sum velocity vector
%     %calculate mean velocity vector
%     VECTORLISTxy=[];%to calculate mean vector
%     for i=1:size(VPLOTBIN,1)
%         for j=1:size(VPLOTBIN,2)
%             if VPLOTBIN(i,j)==1
%                 VECTORLISTxy=[VECTORLISTxy;[CALVFRAMEj2(i,j),CALVFRAMEi2(i,j)]];
%             end
%         end
%     end
%     %     MEANVECTOR=[sum(VECTORLISTxy(:,1)),sum(VECTORLISTxy(:,2))]./length(VECTORLISTxy);
%     MEANVECTOR=nanmean(VECTORLISTxy);
%     MEANVECTOR=MEANVECTOR*ImgScale*binsize./(norm(MEANVECTOR))
%     %plot sum vector into main isochronal map
%     figure(1);hold on
%     R=FRAMECENTER;
%     quiver(R(1)*ImgScale*binsize,R(2)*ImgScale*binsize,MEANVECTOR(1),MEANVECTOR(2),sumvectorscale,'LineWidth',sumvectorwidth,'Color',sumvectorcolor,'MaxHeadSize',1);
%
% end




function edit87_Callback(hObject, eventdata, handles)
% hObject    handle to edit87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit87 as text
%        str2double(get(hObject,'String')) returns contents of edit87 as a double


% --- Executes during object creation, after setting all properties.
function edit87_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28



function edit88_Callback(hObject, eventdata, handles)
% hObject    handle to edit88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit88 as text
%        str2double(get(hObject,'String')) returns contents of edit88 as a double


% --- Executes during object creation, after setting all properties.
function edit88_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit89_Callback(hObject, eventdata, handles)
% hObject    handle to edit89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit89 as text
%        str2double(get(hObject,'String')) returns contents of edit89 as a double


% --- Executes during object creation, after setting all properties.
function edit89_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
FullPath = get(handles.edit29,'String');
load (strcat(FullPath,'CV.mat'));
CV=CV.*100;
MinCV = get(handles.edit88 ,'String');
MinCV = str2double(MinCV);
MaxCV = get(handles.edit89 ,'String');
MaxCV = str2double(MaxCV);
CV(CV>=MaxCV) = 0;
CV(CV<=MinCV)  = 0;
CV(CV==0)   = NaN;
figure(101); histogram(CV,50);
title('Histogram of Conduction Velocity values in the ROI')
xlabel('Conduction Velocity (cm/s)');
ylabel('Count (pixels)')
cvOption = get(handles.popupmenu12,'Value'); %% 1 for smoothing - 2 for Averaging filter
if cvOption == 1
    Est_CV_bounded = nanmean(nanmean(CV));
elseif cvOption == 2
    Est_CV_bounded = nanmedian(nanmedian(CV));
end

set(handles.text109,'String',[num2str(Est_CV_bounded), ' cm/s']);


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
th = get(handles.slider1, 'Value');
set(handles.edit84,'String',num2str(th));
FullPath = get(handles.edit29,'String');
load(strcat(FullPath,'.mat'));
BGFrame = mat2gray(nanmax(RawData,[],3));
% figure; imshow(squeeze(Data(:,:,50)));

Mask = BGFrame>th;
Mask = double(Mask);
% Mask( Mask == 0 ) = NaN;
flipoption = get(handles.checkbox29,'Value');
            if flipoption == 1
                Mask = fliplr(Mask);
                Data = fliplr(Data);
            end
MaskedData = Data.*repmat(Mask,[1,1,size(Data,3)]);
figure(10);imshow(Mask); title('MASK');
% figure; imshow(squeeze(Data(:,:,50)));
% Hint: get(hObject,'Value') returns toggle state of checkbox29


% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30


% --- Executes on button press in pushbutton52.
function pushbutton52_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');
if series
    [FileName,PathName,FilterIndex] = uigetfile('*.tif','Select the first file from tif series');
    A = strcat(PathName, FileName);
    set(handles.edit29,'String',PathName);
    set(handles.edit28,'String',FileName(1:end-4));
elseif stack
    [FileName,PathName,FilterIndex] = uigetfile('*.tif','Select the tif stack file');
    set(handles.edit29,'String',fullfile(PathName,FileName(1:end-4)));
    set(handles.edit28,'String',FileName(1:end-4));
    info = imfinfo(fullfile(PathName,FileName));
    frames = numel(info);
    bit_depth=info.BitDepth;
    set(handles.edit74,'String',num2str(bit_depth));
end
% if ScalingOpt == 1
%     set(handles.checkbox30,'Value',1);
% end


% --- Executes on button press in checkbox31.
function checkbox31_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox31


% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox32



function edit91_Callback(hObject, eventdata, handles)
% hObject    handle to edit91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit91 as text
%        str2double(get(hObject,'String')) returns contents of edit91 as a double


% --- Executes during object creation, after setting all properties.
function edit91_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
series = get(handles.checkbox11,'Value');
stack  = get(handles.checkbox10,'Value');


if series ==1
    PName = get(handles.edit29,'String');
    FName = get(handles.edit28,'String');
    PathName = [fullfile(PName,'results2','corrected','results')];
    load ([fullfile(PathName,FName),'_improved','.mat'])
    B = strcat(PName,'*.tif');
    list_of_frames = dir(B);
    frames = numel(list_of_frames);
elseif stack ==1
    FullPath = get(handles.edit29,'String');
    load(strcat(FullPath,'.mat'));
    info = imfinfo(strcat(FullPath,'.tif'));
    frames = numel(info);
    num_images = numel(info);
    bit_depth=info.BitDepth;
end
PName = get(handles.edit29,'String');
fst = str2num(get(handles.edit92,'String'));
fen = str2num(get(handles.edit93,'String'));

ScaleOpt = get(handles.popupmenu13,'Value');

v = VideoWriter([PName,'_activation','.avi']);
v.FrameRate = 5;
% v.Quality = 95;
figure(9); cla;
open(v)
BFrame=0.5*((mat2gray(RawData(:,:,1))-mean(mean(mat2gray(RawData(:,:,1)).*double(~Mask)))));
% Scale = 0.5/max(mean(mean(AveragedData-mean(mean(AveragedData(:,:,1))),1)));
% Scale = 1/max(mean(mean(AveragedData,1)));
% AveragedData = (AveragedData-min(mean(mean(AveragedData))))./max(mean(mean(AveragedData-min(mean(mean(AveragedData))))));
% figure; plot(squeeze(mean(mean(AveragedData))));
AVERAGING=get(handles.checkbox2,'Value');
if AVERAGING
    AveragedData = AveragedData(:,:,fst:fen);
else
    AveragedData = Data(:,:,fst:fen);
end

% Scale=min(1/max(mean(mean(AveragedData,1))),2)
% Scale = 1;
Scale = ScaleOpt
for i = 1:size(AveragedData,3)
    %     temp = imresize(mat2gray(RawData(:,:,1))+AveragedData(:,:,i),4);
    frame = AveragedData(:,:,i)-mean(mean(AveragedData(:,:,1)));
    %     frame = AveragedData(:,:,i);
    
    %     Scale = 1/max(mean(mean(AveragedData,1)));
    temp = imresize((BFrame+Scale.*frame).*double(Mask),4);
    figure(9);imshow(temp,'Colormap',hot(256));
    M = getframe(figure(9));
    writeVideo(v,M)
end
close(v)
fprintf('Video Saved!\n');




function edit92_Callback(hObject, eventdata, handles)
% hObject    handle to edit92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit92 as text
%        str2double(get(hObject,'String')) returns contents of edit92 as a double


% --- Executes during object creation, after setting all properties.
function edit92_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit93_Callback(hObject, eventdata, handles)
% hObject    handle to edit93 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit93 as text
%        str2double(get(hObject,'String')) returns contents of edit93 as a double


% --- Executes during object creation, after setting all properties.
function edit93_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit93 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
