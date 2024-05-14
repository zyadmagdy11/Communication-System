########################### Q8 ###########################
##########################################################

###################### Redefine x(t1)################
fs = 100;       # Given
df = 0.01;      # Given

ts1 = 0.04/fs;
t1 = -2:ts1:2-ts1;

if (rem(2,2)==0)
f= - (0.5 * fs) : df : (0.5*fs - df) ;
else
f= - (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df);
end

b1 = find (t1 ==-1);
b2 = find (t1 == 1);

x1 = t1(1:b1)+2;
x2 = ones(1,b2-b1-1);
x3 = 2-t1(b2:end);

x = [x1 , x2 , x3];


###################### Redefine m(t2)################

fm = 1;

T = 100/fm;
ts2 = 1/fs;
N = ceil(T/ts2);

if (rem(N,2)==0)
f= - (0.5 * fs) : df : (0.5*fs - df) ;
else
f= - (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df);
end

t2 = 0:ts2:((N-1)*ts2);

m = cos(2*pi*1*t2);
m(t2>6) = 0;


############################################################################

###################### x(t): LPF ######################
H =   abs(f)<3;
X =  (fftshift(fft(x))*ts1) .* H ;
x =  real(ifft(ifftshift(X))/ts1);


################### x(t) modulation #################

c1 = cos(2*pi*20*t2);               #The form of t inside carrier must be t = 0:ts:((N-1)*ts); not t1 = -2:ts1:2-ts1;
s1 = x .* c1;                       #This is Why We have t1 and t2. pass t1 to c1 and see the diffrance in fig(1)

S1 = fftshift(fft(s1))*ts1;
figure(1);
stem(f,abs(S1));   #t1: The Signal isn't shifted 20Hz
grid on;


###################### m(t): LPF ######################
H = abs(f)<3;
M =  (fftshift(fft(m))*ts2) .* H ;
m =  real(ifft(ifftshift(M))/ts2);




#################### Q9: LSB

###################### m(t) modulation ######################

###################### Fc of c2 ######################

BW_m =  abs(f(find(abs(M) > 0,1,"first")));

figure(2);

stem(f,abs(X));
hold on;
stem(f+BW_m+2.5,abs(M));              #Zoom to see Empty band
legend("X(f),M(f)");
grid on;

fc2 = BW_m+2.5;     #5.49 Hz

###################### Q10 ######################

c2 = sin(2*pi*fc2*t2);
s2 = m .* c2;
S2 = fftshift(fft(s2))*ts2;

figure(3);
stem(f,abs(S2));
grid on;


###################### Band Pass Filter to LSB of m(t) LSB ######################
H = zeros(size(f));
H(f>(fc2-BW_m) & f<(fc2))=1;
H(f<-(fc2-BW_m) & f>-(fc2))=1;

S2 = S2 .* H;

figure(4);
stem(f,abs(S2));
grid on;


###################### Q11 ######################

s = s1 + s2;
t = -10:0.002:10-0.002;

figure(5);
plot(t,s);
grid on;


S = fftshift(fft(m));

figure(6);
stem(f,abs(S));                      #stem(f,abs(S1+S2));
grid on;

###################### Q12 ######################


######### m(t) #########
g2 = s .* c2;

H = abs(f)<(BW_m +2.5);


g2 =  real(ifft(ifftshift(fftshift(fft(g2)) .* H)));
figure(7);
plot(t2,g2);
hold on;
plot(t2,m,'r');
hold off;



######### x(t) #########

g1 = s .* c1;
H = abs(f)<10;
g1 =  real(ifft(ifftshift(fftshift(fft(g1)) .* H)));
figure(8);
plot(t1,g1);
hold on;
plot(t1,x,'r');
hold off;
clc;


















