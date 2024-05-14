Rb = 100;
fc = 10001;       #(1000 + 1): plus one so cos(x) != 1
num_bits = 10;

ts = 1/ Rb;
t = 0:ts: (num_bits - ts);

bit_stream = randi([0, 1], 1, num_bits);
bit_stream = repelem(bit_stream, 100);

oscillator_phase = [0,30,60,90];
s = 1;

fs = 1/ts;
df = 0.1;
if (rem(length(t),2)==0)
f=- (0.5 * fs) : df : (0.5*fs - df) ;
else
f=- (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df) ;
end


figure(1);
plot (t,bit_stream);
title('bit stream');
ylim([-1.5,1.5]);
grid on;


###################################################################### Transmitter  ###########################################################################
############################################################################################################################################################


################## ASK modulation ##################
####################################################

c   = cos(2*pi*fc*t);
ask = bit_stream .* c  ;


figure(2);
plot (t,ask);
title('ASK Modulated Signal');
xlabel('t');
ylabel('Amplitude');
ylim([-1.5, 1.5]);
grid on;


#################### ASK spectrum ####################
######################################################
ASK = fftshift(fft(ask))*ts;
figure(3);
stem (f,abs(ASK));
title('ASK spectrum');
xlabel('t');
ylabel('Amplitude');
grid on;

################## FSK modulation ##################
####################################################

freq1 = 33*Rb+1;       #plus 1 so cos(x) != 1
freq2 = 11*Rb+7;       #plus 7 so cos(x) != 1

c1 =   cos(2*pi*freq1 * t);
c2 =   cos(2*pi*freq2 * t);


fsk = zeros(1, length(t)^2);
counter =  1:length(t):length(t)^2 + 1;
for i = 1:length(bit_stream)
    if bit_stream(i) == 0
        fsk(   counter(i) : counter(i)+length(t)-1    ) = c1;
    else
        fsk(   counter(i) : counter(i)+length(t)-1    ) = c2;
    end
    counter(i+1) = counter(i+1) - 1;
end


figure(4);
plot (fsk);
title('FSK Modulated Signal');
xlabel('t');
ylabel('Amplitude');
ylim([-1.5, 1.5]);
grid on;

#################### FSK spectrum ####################
######################################################
FSK = fftshift(fft(fsk))*ts;

figure(5);
stem (repelem(f, 1000),abs(FSK));
title('FSK spectrum');
xlabel('t');
ylabel('Amplitude');
grid on;

########### Polar NRZ Encoding ###########
#########################################


polar_nrz = 2 * bit_stream - 1;

figure(6);
plot (t,polar_nrz);
title('polar NRZ');

ylim([-1.5, 1.5]);
grid on;

################## PSK modulation ##################
####################################################


c3 = cos(2*pi*fc*t);
c4 = cos(2*pi*fc*t+pi);


psk = zeros(1, length(t)^2);
counter1 =  1:length(t):length(t)^2 + 1;
for i = 1:length(polar_nrz)
    if polar_nrz(i) == -1
        psk(   counter1(i) : counter1(i)+length(t)-1    ) = c3;
    else
        psk(   counter1(i) : counter1(i)+length(t)-1    ) = c4;
    end
    counter1(i+1) = counter1(i+1) - 1;
end


figure(7);
plot (psk);
title('PSK Modulated Signal');
xlabel('t');
ylabel('Amplitude');
ylim([-1.5, 1.5]);
grid on;


#################### PSK spectrum ####################
######################################################

PSK = fftshift(fft(psk))*ts;

figure(8);
stem (repelem(f, 1000),abs(PSK));
title('PSK spectrum');
xlabel('t');
ylabel('Amplitude');
grid on;

###################################################################### Receiver  ###########################################################################
############################################################################################################################################################


#################### ASK  ############################
######################################################
c_1    = cos(2*pi*fc*t + oscillator_phase(s));
g_ask = ask .* c_1 ;

G_ASK = fftshift(fft(g_ask))*ts;
H1= abs(f)< fc;

ask_received = real(ifft(ifftshift(H1.*G_ASK))/ts);



figure(9);
plot (t,ask_received);
title('ASK Received Signal');
xlabel('t');
ylabel('Amplitude');
ylim([-1.5,1.5]);
grid on;




#################### FSK  ############################
######################################################

c_2 = -cos(2*pi*freq1*t + oscillator_phase(s)) -cos(2*pi*freq2*t + oscillator_phase(s));





fsk_received  = zeros(1,length(fsk)) ;

H2 = abs(f)< (2*freq2)-100;


for i = 1:length(bit_stream)

      segment = fsk  (     (i-1)*length(t)+1 : i*length(t) );
      fsk_received   (     (i-1)*length(t)+1 : i*length(t) ) = real(ifft(ifftshift(H2.* (fftshift(fft(segment .* c_2 )*ts)))/ts));

end




figure(10);
plot (t,fsk_received(1:length(bit_stream):end));
title('FSK Received Signal');
xlabel('t');
ylabel('Amplitude');
ylim([-1.5,1.5]);

grid on;






#################### PSK  ############################
######################################################
c_3 =   cos(2*pi*fc*t + pi  + oscillator_phase(s) );
psk_received  = zeros(1,length(psk)) ;



H3 = abs(f)< fc;


for i = 1:length(polar_nrz)

    segment = psk  (     (i-1)*length(t)+1 : i*length(t) );
    psk_received   (     (i-1)*length(t)+1 : i*length(t) ) = real(ifft(ifftshift(H3.* (fftshift(fft(segment .* c_3)*ts)))     /ts));

end




figure(11);
plot (t,psk_received (1:length(bit_stream):end)) ;
title('PSK Received Signal');
xlabel('t');
ylabel('Amplitude');
ylim([-1.5,1.5]);

grid on;
















