
function [jitter_local,jitter_local_absolute,jitter_RAP,jitter_PPQ5,shimmer_rel,shimmer_local_absolute,shimmer_APQ3,shimmer_APQ5,shimmer_dB]=JitterShimmer(y,fs)
% y = audioread('1.wav',[1 20000]);
% fs = 44000;
windowLength=0.020*fs;
step= 0.010*fs;
y = y / max(abs(y));
y = filter(1, [1 -0.92], y);

curPos = 1;
L = length(y);
numOfFrames = floor((L-windowLength)/step) + 1; %L=20000 WL= 1100 S=440
H = hamming(windowLength);
dt=1/fs;
jitter_local = [];
jitter_local_absolute = [];
jitter_RAP = [];
jitter_PPQ5 = [];
shimmer_rel = [];
shimmer_local_absolute = [];
shimmer_APQ3 = [];
shimmer_APQ5 = [];
shimmer_dB = [];
tstart = 1;
for i=1:numOfFrames
    window = H.*(y(curPos:curPos+windowLength-1)); % 880x1 which is a single frame
    time = (1/fs)*length (window); % 1 time constant
    t=linspace(tstart,tstart + time,length(window)); % 880 values
    %y=sin(2*pi*100*t);
    %find peaks in signal
    [maxs,~]=peakdet(window,dt);
    %r=peakdet(y,0.6,t);
    %plot(t,y);
    %hold on;
    
    %plot only maximums and minimums onto plot
    % plot(maxs(:,1),maxs(:,2),'g*');
    % plot(mins(:,1),mins(:,2),'r*');
    % hold off;
    
    %maximum peak values time in column 1
    peaktime=maxs(:,1);
    %peak maximum values in column 2
    peakval=maxs(:,2);                  % V(i)
    
    %% Pitch period
    pitchperiods=abs(diff(peaktime));   %P(i) in formula
    pchdiff=zeros(length(pitchperiods)-1,1);   % (N-1) in formula
    
    %difference b/w successive pitch periods
    for k=1:length(pitchperiods)-1             %(N-1)
        pchdiff(k)=abs(pitchperiods(k)-pitchperiods(k+1));  % |P(i)- P(i+1)|
    end
    
    avgpchdiff=mean(pchdiff);           % [Summation |P(i)-P(i+1)|]/(N-1)
    avgpch=mean(pitchperiods);          % [Summation P(i)] / N
    
    jitt_local=(avgpchdiff/avgpch);           %{([Summation |P(i)-P(i+1)|]/(N-1)) / [[Summation P(i)] / N ])}
    jitter_local=[jitter_local;jitt_local];
    
    jitt_local_absolute= mean(pchdiff);
    jitter_local_absolute=[jitter_local_absolute; jitt_local_absolute];
    
    pchdiff2 = zeros(length(pitchperiods)-1,1);
    sumlocpitchper=[];
    for k=1:length(pitchperiods)-1             %i varies from (1:(N-1)) in formula
        if k<2
            sumlocpitchper= mean([pitchperiods(1);pitchperiods(2)]);
            
        else
            localpchperiods=zeros(3,1);
            for i=(k-1):(k+1)                       %j varies from((i-1)to (i+1))
                localpchperiods =(pitchperiods(i));
            end
            sumlocpitchper =mean(localpchperiods);
        end
        
        pchdiff2(k)=abs(pitchperiods(k)-sumlocpitchper);  % |P(i)- (sum(P(j)/3))|
    end
    
    avgpchdiff2=mean(pchdiff2);                             %sum from i=1 to (N-1) (|P(i)- (sum(P(j)/3))|)
    jitt_RAP = (avgpchdiff2 / avgpch);
    jitter_RAP = [jitter_RAP; jitt_RAP];
    
    pchdiff3 = zeros(length(pitchperiods)-1,1);
    sumlocpitchper = [];
    for k=2:length(pitchperiods)-2                  %i varies from (2:(N-2)) in formula
        if k<3
            sumlocpitchper= mean([pitchperiods(1);pitchperiods(2);pitchperiods(3);pitchperiods(4)]);
        else
            localpchperiods=zeros(5,1);
            for i =(k-2):(k+2)                       %j varies from((i-2)to (i+2))
                localpchperiods =(pitchperiods(i));
            end
            sumlocpitchper = mean(localpchperiods);
        end
        
        pchdiff3(k) = abs(pitchperiods(k)-sumlocpitchper);        % |P(i)- (sum(P(j)/5))|
    end
    
    avgpchdiff3 = (sum(pchdiff3))/(length(pitchperiods)-1);
    
    jitt_PPQ5 = (avgpchdiff3 / avgpch);
    jitter_PPQ5 = [jitter_PPQ5; jitt_PPQ5];
    
    pkdiff = abs(diff(peakval));          %|V(I)-V(i+1)|
    avgpkdiff = mean(pkdiff);             %[Summation |V(I)-V(i+1)|]/(N-1)
    avgpks = mean(peakval);               %[Summation V(i)] / N
    
    shim_rel = avgpkdiff/avgpks;              %{([Summation |V(i)-V(i+1)|]/(N-1)) / [[Summation V(i)] / N ])}
    shimmer_rel = [shimmer_rel; shim_rel];
    
    shimm_local_abso = mean(pkdiff);
    shimmer_local_absolute = [shimmer_local_absolute; shimm_local_abso];
    
    pkdiff2 = zeros(length(peakval)-1,1);
    sumlocpekper = [];
    for k=1:length(peakval)-1             %i varies from (1:(N-1)) in formula
        if k<2
            sumlocpekper= mean([peakval(1);peakval(2)]);
            
        else
            localpekperiods=zeros(3,1);
            for i=(k-1):(k+1)                       %j varies from((i-1)to (i+1))
                localpekperiods =(peakval(i));
            end
            sumlocpekper =mean(localpekperiods);
        end
        pkdiff2(k)=abs(peakval(k)-sumlocpekper);      % |V(i)- (sum(V(j))/3))|
    end
    avgpekdiff2=mean( pkdiff2);                     %sum from i=1 to (N-1) (|V(i)- (sum(V(j))/3))|)
    
    shim_APQ3 = (avgpekdiff2 / avgpks);
    shimmer_APQ3 = [shimmer_APQ3; shim_APQ3];
    
    pkdiff3 = zeros(length(peakval)-1,1);
    sumlocpekper=[];
    for k=2:length(peakval)-2               %i varies from (2:(N-2)) in formula
        if k<3
            sumlocpekper= mean([peakval(1);peakval(2);peakval(3);peakval(4)]);
        else
            localpekperiods=zeros(5,1);
            for i=(k-2):(k+2)
                localpekperiods =(peakval(i));
            end
            sumlocpekper =mean(localpekperiods);
        end
        pkdiff3(k)=abs(peakval(k)-sumlocpekper);     % |V(i)- (sum(V(j))/5))|
    end
    avgpekdiff3=mean( pkdiff3);                     %sum from i=1 to (N-1) (|V(i)- (sum(V(j))/5))|)
    
    shim_APQ5 = (avgpekdiff3 / avgpks);
    shimmer_APQ5= [shimmer_APQ5; shim_APQ5];
    
    peklog = zeros(length(peakval)-1,1);
    for k=1:length(peakval)-1
        peklog(k)=[20 * log(peakval(k+1)/peakval(k))];
    end
    peklog = abs(peklog);
   
    shim_dB = mean(peklog);
    shimmer_dB=[shimmer_dB; shim_dB];
    
    curPos = curPos + step;
end
