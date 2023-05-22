% [x_coord,y_coord,intensity_val] = impixel(Data(:,:,55));
% OAP_signal= Data (y_coord,x_coord,:);
OAP_signal=OAP_signal';
wts = [repmat(1/100,100,1)];
        AVES = conv(OAP_signal,wts,'valid');
        % AVES = AVES-(AVESignal(1)-AVECorrected(1));
        LAVE=min(length(OAP_signal),length(AVES)); %% Length of the drif signal and corrected signal
        OAPCorrected = OAP_signal(1:LAVE)-AVES(1:LAVE);
        % figure;plot(AVESignal,'b');hold on;plot(AVES,'r');plot(AVECorrected,'g')
        figure;plot(OAP_signal)
        hold on;plot(OAPCorrected,'r')
%                 [Min, tMin, IMin, Max, tMax, IMax, Npeaks]=fpeaks([1:size(OAPCorrected)], OAPCorrected, 'data');
                type='data';
                t=[1:size(OAPCorrected)];
                x=OAPCorrected;
if strcmp(type, 'data')||strcmp(type, 'video'); % that is the data
    % % Changed these two lines
%     [PKS1,t1] = findpeaks(x); % maximums of the data
%     [PKS1_2,t1_2] = findpeaks([-10000 PKS1 -100000]);
[PKS1,t1] = findpeaks(x); % maximums of the data
[PKS1_2,t1_2] = findpeaks([-10000 PKS1 -100000]);
    if strcmp(type, 'data')
        a=find(PKS1_2>=max(PKS1)/2);
    elseif strcmp(type, 'video')
        videot = 1;
        a=find(PKS1_2>=videot);
    end
    PKS1_3 = PKS1_2(a);
    t1_3 = t(t1(t1_2(a)-1));
    t1_3_diff = diff(t1_3);
    dMax = PKS1_3  ;        % maximums of the data
    tdMax = t1_3    ;       % times of maximums of the data
    IdMax = t1(t1_2(a)-1);  % indices of maximums of the data
    Ndpeaks=length(dMax);   % number of action potentials in the data
    % find the peaks that are too close to each other and choose the bigger
    % value az the main peak
    [r1,c1]=find(abs(t1_3_diff)<=30);
    if numel(c1)~=0
        l=1;
        PKS1_4=[];
        I1_4=[];
        t1_4=[];
        c1_b = c1(l);
        c1_a = c1(l);
        for k=1:Ndpeaks
            % k=2;
            if (k~=c1_b)&&(k~=c1_b+1)&&(k~=c1_b-1)&&(k~=c1_a)
                PKS1_4 = [PKS1_4, PKS1_3(k)];
                I1_4 = [I1_4, IdMax(k)];
                t1_4 = [t1_4, t1_3(k)];
                KK1=k;
                C1 = c1_b;
            elseif k==c1_a
                [pks1,loc1] = max([PKS1_3(k),PKS1_3(k+1)]);
                k1=k;
                
                PKS1_4 = [PKS1_4, pks1];
                I1_4 = [I1_4, IdMax(k+loc1-1)];
                t1_4 = [t1_4, t1_3(k+loc1-1)];
                if l<length(c1)
                    l=l+1;
                    c1_a = c1(l);
                    c1_b = c1(l-1);
                elseif l==length(c1)
                    c1_b = c1(l);
                end
            end
        end
        %
        dMax = PKS1_4;          % maximums of the data
        tdMax = t1_4 ;          % times of maximums of the data
        IdMax = I1_4;           % indices of maximums of the data
        Ndpeaks=length(dMax);   % number of action potentials in the data
    end
    % find minimums of the signal
    if strcmp(type, 'data')
        [dMin,IdMin]=min(x(1 : IdMax(1)));
        tdMin=t(IdMin);
        for j=1:Ndpeaks-1
            [dmin,I] = min(x(IdMax(j) : IdMax(j+1)));
            IdMin = [IdMin, I+IdMax(j)-1];
            dMin=[dMin, dmin];
            tdMin=[tdMin, t(I+IdMax(j)-1)];
        end
        [dmin,I] = min(x(IdMax(j+1) : end));
        if numel(I)~=0 && numel(IdMax)~=0
        IdMin = [IdMin, I+IdMax(j+1)-1] ;       % indices of minimums of the data
        dMin=[dMin, dmin];                      % minimums of the data
        tdMin=[tdMin, t(I+IdMax(j+1)-1)];       % times of minimums of the data
        end
    elseif strcmp(type, 'video')
        x=-x;
        % % modified:
%         [dMin,IdMin] = findpeaks(x(1 : IdMax(1)));
[IdMin,dMin] = findpeaks(x(1 : IdMax(1)));
% %
        tdMin=t(IdMin(end));
        dMin=-dMin(end);
        IdMin=IdMin(end);
        for j=1:Ndpeaks-1
            % Modified:
%             [dmin,I] = findpeaks(x(IdMax(j) : IdMax(j+1)));
            [I,dmin] = findpeaks(x(IdMax(j) : IdMax(j+1)));
            IdMin = [IdMin, I(1)+IdMax(j)-1];
            IdMin = [IdMin, I(end)+IdMax(j)-1];
            dMin=[dMin, -dmin(1)];
            dMin=[dMin, -dmin(end)];
            tdMin=[tdMin, t(I(1)+IdMax(j)-1)];
            tdMin=[tdMin, t(I(end)+IdMax(j)-1)];
        end
        % Modified:
%         [dmin,I] = findpeaks(x(IdMax(j+1) : end));
        [I,dmin] = findpeaks(x(IdMax(j+1) : end));
        IdMin = [IdMin, I(1)+IdMax(j+1)-1];        % indices of minimums of the data
        dMin=[dMin, -dmin(1)];                      % minimums of the data
        tdMin=[tdMin, t(I(1)+IdMax(j+1)-1)];       % times of minimums of the data
        
    end
    Min    =    dMin;
    tMin   =    tdMin;
    IMin   =    IdMin;
    Max    =    dMax;
    tMax   =    tdMax;
    IMax   =    IdMax;
    Npeaks =    Ndpeaks;
    
    
elseif strcmp(type , 'model')

    % % modified these two lines
%     [mMax,ImMax] = findpeaks(x);  % maximums of the piece of model and their indices
%     [mMin,ImMin] = findpeaks(-x); % minimums of the piece of model and their indices
    [ImMax,mMax] = findpeaks(x);  % maximums of the piece of model and their indices
    [ImMin,mMin] = findpeaks(-x); % minimums of the piece of model and their indices
    % %
    mMin    = -mMin;
    tmMin   = t(ImMin);    % times of minimums of the model
    tmMax   = t(ImMax);    % times of maximums of the model
    Nmpeaks = length(mMax); % number of action potentials in the model
    Min=mMin;
    tMin= tmMin;
    IMin=ImMin;
    Max=mMax;
    tMax=tmMax;
    IMax=ImMax;
    Npeaks=Nmpeaks;
end