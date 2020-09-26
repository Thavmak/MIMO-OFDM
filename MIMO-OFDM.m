%---SIMULATION OF AN OFDM AWGN/LTI CHANNEL---
pkg load communications;
clear all;
close all;
clc;

%-----------PARAMETERS-----------------------
N        = 10^4;       %Bits for transmission 
OS       = 1;          %Oversampling
Rb       = 10^4;       %Bit Rate (bps)
T        = 1;          %How frequent h changes (sec)
Eb       = 1;          %Energy of bit
%SNR     = 0;
SNR      = [0:5:15];   %SNR in dB
%M       = 2;
M        = [2 4 8 16]; %Modulation (BPSK, QPSK, 8PSK, 16PSK)
QAM      = 1;          %apply 16QAM Modulation instead 16PSK (true=1, false=0)
c        = 0;          %Channel estimation inaccuracy (set c=0 if channel estimation is correct)
tx_ant   = 1;          %number of Tx antennas
rx_ant   = 1;          %number of Rx antennas
channelt = 1;          %select if OFDM AWGN (channelt=0), OFDM LTI (channelt=1), LTI (channelt=2)
sub      = 64;         %number of subcarriers
cp       = 16;         %cyclic prefix
%---------------------------------------------


%-------CONSTRUCT h VECTOR (channel)----------

%how many bits am I able to transmit per antenna between changes of h
tbits=Rb*T;
%total transmissions needed to sent all the bits across the channel
tnum=round(N/tbits);
%how many transmissions per antenna
tant=round(tnum/tx_ant);
%how many different h will affect every transmission
hsize=tnum*rx_ant;
%Variance(ó^2) of h
varh=1/sqrt(2);
%Standard deviation(ó) of h
stdh=sqrt(varh/2);
%create h vector
if (channelt==2) %LTI
  h=stdh*(randn(1,hsize)+1i*randn(1,hsize));
else %OFDM LTI
  h=[0.9+0.9i 0.6+0.6i 0.3+0.3i];
end

%---------ENCODED DATA-----------------

%create a vector with random 0 & 1 bits
encodedData=randi([0 1],1,N);
%create data streams for every transmission
enc_split=vec2mat(encodedData,tnum);


for i=1:length(M) %runs for every Modulation scheme
  
  %--------CONSTRUCT THE CONSTELLATION DIAGRAM-------
  
  n=log2(M(i)); %n = bits/symbol
  complvalv=[]; %initialize constellation points table
  
  if (M(i)==16 && QAM==1) %use 16QAM Modulation
    
    x2=2/sqrt(10);
    x3=sqrt(4-x2^2);
    table=[x2 x3 -x2 -x3];
    for xval=1:n
      for yval=1:n
        complval=table(xval)+1i*table(yval);
        %complvalv contains all the constellation points
        complvalv=[complvalv; complval];
      end
    end
    
  else %use M-PSK Modulation
    
    x1=0:M(i)-1;
    %divide 2pi(or 360 degrees) in M equal parts
    %complvalv contains all the constellation points
    complvalv=cosd((360/M(i))*x1)+1i*sind((360/M(i))*x1);
    complvalv=complvalv.';
    
  end
  
  Ic=real(complvalv); 
  Qc=imag(complvalv);
  
  %------------SYMBOL ENCODER-----------------
  
  for l=1:length(SNR) %run for every value of SNR
    
    wrong_symbolsv=[];
    goodpackets=0;
    cnt=1; %counter of transmissions
    
    for q=1:tant %transmissions begin
      
      xv=[]; %initialize signal x
      hm=zeros(rx_ant,tx_ant); %initialize channel h vector/matrix
      sentv=[]; %initialize sentv
      
      for q2=1:tx_ant %run for # signals that will be transmitted from every Tx
        
        %construct signal x
        %take a data stream for this transmission
        enc_h=enc_split(:,cnt);
        %make a matrix with n bits per row
        %every row corresponds to a symbol
        encMat=vec2mat(enc_h,n);
        %convert every row/symbol to a decimal number 
        x4=1:(tbits/n);
        valMat=bi2de(encMat(x4,:),'left-msb'); 
        %match every decimal number to a complex value
        if (M(i)==16 && QAM==1) %16QAM
          sent=complvalv(valMat+1); 
        else %MPSK
          sent=cosd((360/M(i))*valMat)+1i*sind((360/M(i))*valMat);
        end
        sentv=[sentv, sent];
        
        %I & Q signal
        Is=real(sent);
        Qs=imag(sent);
        
        %-----------TX FILTER---------------------
        
        %----OVERSAMPLE------
        %length here is the size, signal will have after the OS
        Is1=zeros(1,length(Is)*OS); 
        Qs1=zeros(1,length(Qs)*OS); 
        u=1;
        
        for w=1:length(Is)
          %add symbols at specified positions
          Is1(u)=Is1(u)+Is(w); 
          Qs1(u)=Qs1(u)+Qs(w); 
          u=u+OS;
        end
        
        %Amplitude of bit after Oversampling
        AmplOS=sqrt(Eb/OS); 
        %impulse response of TX filter
        ir=ones(1,OS)*AmplOS; 
        %signal x passing through TX filter
        IsF=conv(Is1,ir);
        QsF=conv(Qs1,ir);
        x=IsF+1i*QsF;
        %add x to a matrix of tx rows
        xv=[xv; x];
        
        if (channelt~=2) %OFDM AWGN/LTI 
          
          %--------SERIAL/PARALLEL----------
          
          %divide signal into OFDM symbols
          xv=vec2mat(xv,sub);          
          %delete last row
          xv(size(xv,1),:)=[];
          
          %---IDFT (Inverse Discrete Fourier Transform)----
          x_ffts=fftshift(xv,2);
          x_idft=ifft(x_ffts,sub,2);
          
          %-----ADD CP (cyclic prefix)----- 
          xcp2=[];
          xi2=[];
          for cpn=1:size(xv,1) %for every OFDM symbol
            %take every OFDM symbol
            xi=x_idft(cpn,:);
            %copy the last cp symbols
            xcp=xi(sub-cp+1:sub);
            %add them to the front of the OFDM symbol 
            xi2=[xcp xi];
            %new x with cp added
            xcp2=[xcp2; xi2];
          end
         
          xv=xcp2;
          
          %---PSD PLOT----
          %{
          Fs=10000;
          xv_psd=xv.';
          xv_psd=reshape(xv_psd,[1,numel(xv_psd)]);
          xv_psd=xcorr(xv_psd,'biased');
          xv_psd=abs(xv_psd);
          freq = -Fs/2:Fs/length(xv_psd):Fs/2-(Fs/length(xv_psd));
          figure(234,'name','Power Spectral Density plot','NumberTitle','off')
          plot(freq, 10*log10(xv_psd));
          title("PSD of OFDM Signal (Symbol Rate 10KBps)");
          ylabel('PSD(dB/Hz)');
          xlabel('Frequency (Hz)');
          %}
          %---------------
          
        end
        
        %-------CHANNEL------------
        if (channelt==2) %LTI
          %prepare h
          hv=[]; %empty h
          %create h vector/matrix for this transmission
          for qt=((cnt-1)*rx_ant)+1:cnt*rx_ant
            hv=[hv h(qt)];
            hm(:, q2)= hv;
          end
        end
        
        cnt=cnt+1; %count transmission
        
      end %signal/s that will be transmitted is/are ready
      
      %prepare noise
      noisev=[]; %initialize noise
      for rxn=1:rx_ant %construct noise
        if (channelt==2)
          %number of noise variables
          L=length(x);
          %SNRdB to linear scale
          SNRlin=10^(SNR(l)/10);
          %Calculate symbol energy
          Es=sum(abs(x).^2)/(L);
          %Find the noise spectral density/variance(ó^2) 
          N0=Es/SNRlin;
          %Standard deviation(ó) for AWGN
          noiseSigma=sqrt(N0/2);
          %compute noise
          noise=noiseSigma*(randn(1,L)+1i*randn(1,L));
          %add noise to a matrix of tx rows
          noisev=[noisev; noise];         
        else
          %number of noise variables
          if (channelt==1)
            L=numel(xcp2)+size(xcp2,1)*2;
          else
            L=numel(xcp2);
          end
          %SNRdB to linear scale
          SNRlin=10^(SNR(l)/10);
          %Calculate symbol energy
          %Es=1/sqrt(sub);
          Es=1/sub;
          %Find the noise spectral density/variance(ó^2) 
          N0=Es/SNRlin;
          %Standard deviation(ó) for AWGN
          noiseSigma=sqrt(N0/2);
          %compute noise
          noise=noiseSigma*(randn(1,L)+1i*randn(1,L));
          if (channelt==1)
            noisev=vec2mat(noise,sub+cp+2);
          else
            noisev=vec2mat(noise,sub+cp);
          end
        end
      end %end of noise construction
      
      %add h and noise to the signal
      %y contains signals received by all the antennas at the receiver 
      if (channelt==2) %LTI
        y=hm*xv+noisev;
      elseif (channelt==0) %OFDM AWGN
        y=xv+noisev;
      else %OFDM LTI
        y=[];
        for ofdms=1:size(xv,1)
          %take each OFDM symbol
          xvt=xv(ofdms,:);
          %convolute it with h of channel
          xvt2=conv(xvt,h);
          %take noise
          nst=noisev(ofdms,:);
          %add noise to signal
          yt=xvt2+nst;
          y=[y; yt];
        end 
      end
      
      if (channelt~=2) %OFDM AWGN/LTI Channel
        
        %---REMOVE CYCLIC PREFIX----
        ycpm=[];
        for cpr=1:size(y,1)
          yscp=y(cpr,:);
          %remove the cp (first cp symbols) and the last 2 from the convolution
          ycp=yscp(cp+1:sub+cp);
          ycpm=[ycpm; ycp];
        end
        
        %---DFT (Discrete Fourier Transform)----
        y_fft=fft(ycpm,sub,2);
        y_iffts=fftshift(y_fft,2);
        
        %---PARALLEL/SERIAL--------
        ytemp=y_iffts;
        y_iffts=y_iffts.';
        y=reshape(y_iffts,[1,numel(y_iffts)]);
        
      end
      
      Iy=real(y);
      Qy=imag(y);
      
      %---RX FILTER---------------------
      
      yrv=[]; %initialize received signal
      for rxf=1:rx_ant
        
        %convolute the received signal
        IyF1=conv(Iy(rxf,:),ir); 
        QyF1=conv(Qy(rxf,:),ir);
        %y received signal
        yr=IyF1+1i*QyF1;
        yrv=[yrv; yr];
        
      end
      
      %trying to mitigate the c error
      %mhc=mean(h+c);
      %ece=real(mhc);
      ece=0; %no error correction
      
      if (channelt==2) %LTI channel
        
        %---EQUALIZATION---------------
        
        %Maximal Ratio Combining (MRC)(1x1 & 1x2)
        if (tx_ant==1)
          Hn=0; %numerator of h^-1
          Hd=0; %denominator of h^-1
          for lp=1:rx_ant
            Hn=Hn+(conj(hm(lp,:)+(c-ece))*yrv(lp,:));
            Hd=Hd+(abs(hm(lp,:)+(c-ece)).^2);
          end
          yeq=Hn/Hd;
          
        %Least Squares channel inversion (LS)(2x2)
        else
          W=((hm'*hm)^-1)*hm';
          yeq=W*yrv;
        end %end of equalization
        
        %uncomment to disable EQUALIZATION
        %yeq=yrv;
        
      elseif (channelt==0) %OFDM AWGN
        
        yeq=yrv;
        
      else %OFDM LTI
        
        Hn1=[];
        Hn = fftshift(fft(h,64,2));
        for cnt1=1:size(ytemp,1);
          Hn1= [Hn1; Hn];
        end
        yeq=ytemp./Hn1;
        yeq1=yeq;
        yeq=yeq.';
        yeq=reshape(yeq,[1,numel(yeq)]);
        
      end
      
      IyF=real(yeq);
      QyF=imag(yeq);
      
      %---SAMPLING------------------------
      
      for sc=1:tx_ant %sample & decide for every signal in Rx
        
        if (channelt==2) %LTI
          
          %sample each signal separately
          IyF1=IyF(sc,:);
          QyF1=QyF(sc,:);
          IyS=IyF1(OS:OS:OS*tbits/n); 
          QyS=QyF1(OS:OS:OS*tbits/n);
          
        else %OFDM AWGN/LTI
          
          IyS=IyF;
          QyS=QyF;
          
        end
        
        %---DECISION------------------------
        
        distance=[]; %initialize distance matrix
        minDist=[]; %initialize minimum distance vector
        %for every symbol received
        for k=1:length(IyS) 
          
          %from every constellation point
          for j=1:M(i) 
            
            %find the distance
            distance(k,j)=sqrt((IyS(k)-Ic(j))^2+(QyS(k)-Qc(j))^2);
            
          end
          
        end
        
        %now find the smallest distance for every symbol
        minDist = min(distance,[],2);
        
        closest=[]; %initialize closest to constellation point
        decision=[]; %initialize decision vector
        %for every symbol
        for k=1:length(IyS)
          
          %out of which constellation point was the smallest distance
          closest(k)=find(distance(k,:)==minDist(k)); 
          decision(k)=complvalv(closest(k));
          
        end
        
        %---BER CALCULATIONS------------------
        
        %find how many symbols were wrong in this transmission
        wrong_symbols=0;
        sent1=sentv(:,sc);
        sent1=sent1(1:numel(IyS));
        sent1=sent1.';
        %compare sent symbols with decision
        result=(decision~=sent1);
        %count all wrong symbols
        wrong_symbols=sum(result);
        %keep wrong symbols from every transmission
        wrong_symbolsv=[wrong_symbolsv; wrong_symbols];
        
        %if there are no errors count it as goodpacket
        if (wrong_symbols==0)
          goodpackets=goodpackets+1;
        end
        
        %{
        %SCATTER PLOT for SNR vector
        if (M(i)==2)
          movegui(figure(i,'name',['BPSK (c=',num2str(c),')'],'NumberTitle','off'),'northwest');
        elseif (M(i)==4)
          movegui(figure(i,'name',['QPSK (c=',num2str(c),')'],'NumberTitle','off'),'north');
        elseif (M(i)==8)
          movegui(figure(i,'name',[num2str(M(i)),'PSK (c=',num2str(c),')'],'NumberTitle','off'),'northeast');
        else
          movegui(figure(i,'name',[num2str(M(i)),'PSK (c=',num2str(c),')'],'NumberTitle','off'),'southwest');
        end
        subplot(2,2,l)
        scatter(IyS,QyS,1,"b",'filled');
        hold on
        title(['SNR=',num2str(SNR(l)),'dB']);
        %}
        
      end
      
    end %end of transmission loop
    %hold off
    
    %find how many symbols were wrong in total  
    total_wrong=0;
    total_wrong=sum(wrong_symbolsv);
    
    %total symbols sent
    %total_sent=round(tbits/n)*tnum; 
    total_sent=numel(result);
    %Symbol Error Rate (SER)
    SER=total_wrong/total_sent;
    %Calculate BER values for every SNR
    BER(l)=SER/n;
    
    if (M(i)==2||M(i)==4)
      %theoretical BER for B/QPSK, AWGN channel
      BQPSK(l)=0.5*erfc(sqrt(10.^(SNR(l)/10)));
    else
      %theoretical BER for MPSK, AWGN channel
      MPSK(l)=1/n*erfc(sqrt((10.^(SNR(l)/10))*n)*sin(pi/M(i)));
    end
    
    %---THROUGHPUT CALCULATIONS-----------
    
    throughput=Rb*tx_ant*n; %in bps
    %how many symbols in a goodpacket? Rb*n*T
    total_good_bits= goodpackets*Rb*T;
    goodput=(total_good_bits/N)*throughput;
    good_perc=(total_good_bits/N)*100;
    
  end %end of SNR loop
  
  %BER
  %{-
  %BER GRAPH for every modulation 
  movegui(figure(11,'name',['BER values of a ',num2str(tx_ant),'x',num2str(rx_ant),' system'],'NumberTitle','off'),'south');
  warning ("off", "Octave:negative-data-log-axis");
  if (M(i)==2)
    semilogy(SNR,BER,'-ob');
    %hold on
    %semilogy(SNR,BQPSK,'--or');
  elseif (M(i)==4)   
    semilogy(SNR,BER,'-or');
    %semilogy(SNR,BQPSK,'--or');
  elseif (M(i)==8) 
    semilogy(SNR,BER,'-og');
    %semilogy(SNR,MPSK,'--og');
  else 
    semilogy(SNR,BER,'-oy');
    %semilogy(SNR,MPSK,'--oy');
  end
  hold on
  %}
  
end %end of Modulation loop

%{-
%uncomment for BER GRAPH to work 
xlabel({'SNR=2Eb/N0','(in dB)'});
ylabel('BER');
%grid on
%legend('BPSK','QPSK','Diff. enc. B/QPSK','8PSK','Diff. enc. 8PSK','16PSK','Diff. enc. 16PSK','Location','southwest')
if (QAM==1) %use 16QAM Modulation
  %legend('16QAM')
  legend('BPSK','QPSK','8PSK','16QAM','Location','southwest')
else
  %legend('BPSK')
  legend('BPSK','QPSK','8PSK','16PSK','Location','southwest')
end
hold off
%}