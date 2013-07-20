clc;
close all;
[s1 fs1]=wavread('s1.wav'); %the mono signal.change here to change the signal
[m n]=size(s1)
figure,plot(s1);
title('original signal');
j=1;
overlap=fs1*0.01
framesize=fs1*.025

%dividing into frames
for i=1:overlap:m
    if(i+framesize<=m)
          frame(j,1:framesize)=s1(i:i+framesize-1)';
          j=j+1;
    else
        frame(j,1:framesize)=[s1(i:m)' zeros(1,i+framesize-m-1)];
        if(mod(j,2)==0)
            break;
        else
            j=j+1;
        end
    end
end

number_of_samples=j   ;
[frsi frsi2]=size(frame(j,:));
k=1;
ham_window=hamming(framesize); %
%now we compute the fourier transform of the signal and compute the
%periodogram
for i=1:number_of_samples
    x=fft(frame(i)*ham_window,(fs1/2));
    fourier(i,1:(fs1/2))=x(1:(fs1/2));
    periodogram(i,1:(fs1/2))=(abs(fourier(i,1:(fs1/2))).^2)/framesize;
end
for i=1:(fs1/2)
    y(i)=i;
end
figure,plot(y,periodogram);
title('periodogram');
%now we compute the mel filters
freq_l=300;
freq_u=fs1/2;
mel_l=1125*log(1+(freq_l/700)); %lower mel freq
mel_u=1125*log(1+(freq_u/700)); %upper mel freq
%we will  create a 26 filterbank
mel_increment=(mel_u-mel_l)/27;
for i=1:28
    mel_f(i)=mel_l+(i-1)*mel_increment;
    freq(i)=700*(exp(mel_f(i)/1125)-1);
end

filter=zeros(26,1:(fs1/2));%initialize the mel filterbank
%  now the mel filterbanks(i.e triangular filterbanks) will created
j=1;
for i=2:27
    for k=1:(fs1/2)
        if(k<freq(i-1))
            fliter(j,k)=0;
        elseif((freq(i-1)<=k) && (k<=freq(i)))
              
               filter(j,k)=(k-freq(i-1))/(freq(i)-freq(i-1));
        elseif((freq(i)<=k) && (k<=freq(i+1)))
                filter(j,k)=(freq(i+1)-k)/(freq(i+1)-freq(i));
        elseif(k>freq(i+1))
                filter(j,k)=0;
        end
    end
    j=j+1;
end
  for i=1:fs1/2
      k(i)=i;
  end
 
  figure,plot(k,filter(1,:),k,filter(2,:),k,filter(3,:),k,filter(4,:),k,filter(5,:),k,filter(6,:),k,filter(7,:),k,filter(8,:),k,filter(9,:),k,filter(10,:),k,filter(11,:),k,filter(12,:),k,filter(13,:),k,filter(14,:),k,filter(15,:),k,filter(16,:),k,filter(17,:),k,filter(18,:),k,filter(19,:),k,filter(20,:),k,filter(21,:),k,filter(22,:),k,filter(23,:),k,filter(24,:),k,filter(25,:),k,filter(26,:));
  title('mel filterbank 1');
  figure,plot(k,filter(j-1,:));
  title('mel filterbank 20');

%now we calculate the filterbank energies 
figure,plot(periodogram);
title('periodogram');
n1=size(filter,2)
temp_frame=zeros(1,n);
n2=size(periodogram,2)
for i=1:number_of_samples
    for j=1:26
        energy(i,j)=periodogram(i,:)*filter(j,:)';
    end
end
[m n]=size(energy);
figure,plot(energy)
title('energy');
energy2=log(energy);
figure,plot(energy2);
title('log of energy');
%calculating the cepstral coefficients i.e taking the DCT of the log of the
%energy coefficients.This decorrelates the matrix.
for i=1:number_of_samples
    cep1(i,1:26)=dct(energy2(i,:));
end
    % taking in account the low 12-13 coefficients for MFCC'S
for i=1:number_of_samples
    mfcc(i,1:12)=cep1(i,15:26);
end
disp('final MFCC calculated')
[m n]=size(mfcc);
figure,plot(mfcc(1,:));
title('MFCC');
    
