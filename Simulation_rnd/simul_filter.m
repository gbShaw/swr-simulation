function X = simul_filter(X,fs,fl,fh,order)

wp=[fl/(fs/2) fh/(fs/2)];
b=fir1(order,wp,blackman(order+1));
X = filtfilt(b,1,X);
