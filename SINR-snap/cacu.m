function snr=cacu(woo)
ti=1;
sinr(ti)=sinr(ti)+10*log10((woo'*Rs*woo)/(woo'*Rjn*woo));
 snr=sinr(ti);
end