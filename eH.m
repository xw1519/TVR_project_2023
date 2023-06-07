% historgram of a signal

function [pxx2,f2] = eH(data,fs,smth, isplot)

    for nchan=1:size(data,1)
        % power Filtering
        [pxx2,f2] = pwelch(detrend((data(nchan,:))),fs*smth,1,fs,fs);
        if isplot
            plot( f2, 10*log10( pxx2 ) ,'linewidth',2);
            hold on
            xlim( [ 0 fs/2 ] )
            ylim([-40,30])
        end
        
        %xlabel('Frequency (Hz)')
        %ylabel('PSD (dB/Hz)')
    end
end
