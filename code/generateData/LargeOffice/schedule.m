function sp = schedule(ctrlsp, TOD, DOW)

if strcmp('ClgSP', ctrlsp)
    
    %     "Weekdays" : [6,26.7,
    % 				  22,24.0,
    % 				  24,26.7],
    %     "Saturday" : [6,26.7,
    % 				  18,24.0,
    % 				  24,26.7],
    %     "AllOtherDays" : [24,26.7]

    if strcmp('Weekdays', DOW)
        
        if TOD <= 6*3600, sp = 26.7;
        elseif 6*3600 < TOD && TOD <= 22*3600, sp = 24;
        elseif 22*3600 < TOD && TOD <= 24*3600, sp = 26.7;
        end;
        
    elseif strcmp('Saturday', DOW)
        
        if TOD <= 6*3600, sp = 26.7;
        elseif 6*3600 < TOD && TOD <= 18*3600, sp = 24;
        elseif 18*3600 < TOD && TOD <= 24*3600, sp = 26.7;
        end;
        
    elseif strcmp('AllOtherDays', DOW)
        
        if TOD <= 24*3600, sp = 26.7;
        end;
        
    end
    
elseif strcmp('LgtSP', ctrlsp)
   
    %         "Weekdays" : [5,0.05,
    % 				  7,0.1,
    % 				  8,0.3,
    % 				  17,0.9,
    % 				  18,0.7,
    % 				  20,0.5,
    % 				  22,0.3,
    % 				  23,0.1,
    % 				  24,0.05],
    %     "Saturday" : [6,0.5,
    % 				  8,0.1,
    % 				  14,0.5,
    % 				  17,0.15,
    % 				  24,0.05],
    %     "AllOtherDays" : [24,0.05]

    if strcmp('Weekdays', DOW)
        
        if TOD <= 5*3600, sp = 0.05;
        elseif 5*3600 < TOD && TOD <= 7*3600, sp = 0.1;
        elseif 7*3600 < TOD && TOD <= 8*3600, sp = 0.3;
        elseif 8*3600 < TOD && TOD <= 17*3600, sp = 0.9;
        elseif 17*3600 < TOD && TOD <= 18*3600, sp = 0.7;
        elseif 18*3600 < TOD && TOD <= 20*3600, sp = 0.5;
        elseif 20*3600 < TOD && TOD <= 22*3600, sp = 0.3;
        elseif 22*3600 < TOD && TOD <= 23*3600, sp = 0.1;
        elseif 23*3600 < TOD && TOD <= 24*3600, sp = 0.05;
        end;
        
    elseif strcmp('Saturday', DOW)
        
        if TOD <= 6*3600, sp = 0.05;
        elseif 6*3600 < TOD && TOD <= 8*3600, sp = 0.1;
        elseif 8*3600 < TOD && TOD <= 14*3600, sp = 0.5;
        elseif 14*3600 < TOD && TOD <= 17*3600, sp = 0.15;
        elseif 17*3600 < TOD && TOD <= 24*3600, sp = 0.05;
        end;
        
    elseif strcmp('AllOtherDays', DOW)
        
        if TOD <= 24*3600, sp = 0.05;
        end;
        
    end
    
else
    
    error('illegal setpoint');
    
end