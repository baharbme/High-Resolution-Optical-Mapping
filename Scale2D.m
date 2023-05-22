function Y = Scale2D(X,binsize) 
if length(size(X)) == 2
    X = imgaussfilt(X,0.5);
    Y = nan(floor(size(X,1)/binsize),floor(size(X,2)/binsize));
    for i = 1 : size(Y,1)
        for j = 1 : size(Y,2)
            Y(i,j) = nanmean(nanmean(X((i-1)*binsize+1:i*binsize , (j-1)*binsize+1:j*binsize)));
        end
    end
elseif length(size(X)) == 3
    Xall = X;
    Yall = nan(floor(size(X,1)/binsize),floor(size(X,2)/binsize),size(X,3));
    for k=1:size(X,3)
        X=squeeze(Xall(:,:,k));
        X = imgaussfilt(X,0.5);
        Y = nan(floor(size(X,1)/binsize),floor(size(X,2)/binsize));
        for i = 1 : size(Y,1)
            for j = 1 : size(Y,2)
                Y(i,j) = nanmean(nanmean(X((i-1)*binsize+1:i*binsize , (j-1)*binsize+1:j*binsize)));
            end
        end
        Yall(:,:,k)=Y;
    end
    Y=Yall;
end
