function i0=findtime(frames,k)
    if frames > 999
        if k < 10
            i0 = '000';
        elseif k >= 10 && k < 100
            i0 = '00';
        elseif k >= 100 && k < 1000
            i0 = '0';
        else
            i0 = '';
        end
    end
    %>99 frames
    if frames > 99 && frames<=999
        if k < 10
            i0 = '00';
        elseif k >= 10 && k < 100
            i0 = '0';
        else
            i0='';
        end
    end
    
    if frames < 100
        if k < 10
            i0 = '0';
        else
            i0='';
        end
    end
    