function formants = formant(x,fs)

   ns = length(x);

% error checking on the signal level, remove the DC bias

x = x - mean(x);    x = x / std(x);
x = filter(1, [1 -0.92], x);
% use a 30msec segment, choose a segment every 20msec
% that means the overlap between segments is 10msec
wt = 11.3;
ov = 10;  

win_sam  = fix( wt * fs * 10^-3 );
overlap = fix(win_sam * ov / 100 );
sam_shift  = win_sam - overlap;
nframe = fix(ns / sam_shift) -1;

% get the formants from each segmented frame

frame =[];
avgF0 =[];
start = 1;
stop = win_sam;
formants=zeros(nframe,4);
i = 1;
window= hamming(win_sam);
while stop < ns
    seg = x(start: stop, 1);
    win_seg = seg .* window;
    % get Linear prediction filter
    ncoeff=2+fs/1000;      % rule of thumb for formant estimation
    a = lpc(win_seg, ncoeff);
    % find frequencies by root-solving
  
    r=roots(a);                  % find roots of polynomial a
    r=r(imag(r)>0.01);           % only look for roots >0Hz up to fs/2
    ffreq=sort(atan2(imag(r),real(r))*fs/(2*pi));  % convert to Hz and sort
    formants(i, :) = ffreq(1:4)';
    start = start + sam_shift ;
    stop  = start + win_sam - 1 ;
    i = i + 1;
end

end

