fs = 100;
df = 0.01;

ts = 0.04/fs;
t = -2:ts:2-ts;

if (rem(length(t),2)==0)
f= - (0.5 * fs) : df : (0.5*fs - df) ;
else
f= - (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df);
end

b1 = find (t ==-1);
b2 = find (t == 1);

x1 = t(1:b1)+2;
x2 = ones(1,b2-b1-1);
x3 = 2-t(b2:end);

x = [x1 , x2 , x3];


######### Q1 #########
figure(1);
plot(t,x);
ylim([0, 1.5]);
xlabel("t");ylabel("x(t)");
title('Time Domain Representation of x(t)');
grid on;

######### Q2 #########
figure(2);

X_analytical = (-1./(2*(pi^2)*(f.^2))) .* (cos(4*pi*f) - cos(2*pi*f)) ;

stem(f,abs(X_analytical),'r')
xlabel("f"); ylabel("X(f)");

legend("X(f) Analytical Solution");
title('Analytical Fourier Transform of x(t)');
grid on;

######### Q3 #########
X_numerical = fftshift(fft(x))*ts;

figure(3);

stem(f,abs(X_numerical));
hold on;
stem(f+100,abs(X_analytical));

legend("X(f) Numerical Solution" , "X(f) Analytical Solution");
xlabel("Numerical & Analytical are not exactly identical",'FontWeight', 'bold');
title('Comparison of Numerical and Analytical Fourier Transforms');
grid on;


######### Q4 #########

P = abs(X_numerical).^2;

P_max = max(P);
P_max_index = find(P == P_max)
threshold = 0.05 * P_max;

BW = abs(f(P_max_index + find( P(P_max_index:end) <= threshold,1) -1))


######### Q5 #########

H= abs(f)< 1;
X_filtered = X_numerical .* H;
figure(4);
stem(f,abs(X_filtered))
xlabel("f"); ylabel("X_filtered(f)");
grid on

x_filtered = real(ifft(ifftshift(X_filtered))/ts);
figure(5);
plot(t,x_filtered)
xlabel("t"); ylabel("x_filtered(t)");
grid on


######### Q6 #########

H= abs(f)< 0.3;
X_filtered = X_numerical .* H;
figure(6);
stem(f,abs(X_filtered))
xlabel("f"); ylabel("X_filtered(f)");
grid on

x_filtered = real(ifft(ifftshift(X_filtered))/ts);
figure(7);
plot(t,x_filtered)
xlabel("t"); ylabel("x_filtered(t)");
grid on


close all;

########################### Q7 ###########################
##########################################################

fs = 100;
fm = 1;
df = 0.01;

T = 100/fm;
ts = 1/fs;
N = ceil(T/ts);

if (rem(N,2)==0)
f= - (0.5 * fs) : df : (0.5*fs - df) ;
else
f= - (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df);
end

t = 0:ts:((N-1)*ts);

############## Q1 #############

m = cos(2*pi*1*t);
m(t>6) = 0;

figure(1);
plot(t,m);
title("====>> m(t) <<====")
grid on;

############## Q2 #############

M_Analytical = (3* (exp(3*j*(2*pi-2*pi*f)) .* sinc(6-6*f)) + 3*(exp(-3*j*(2*pi+2*pi*f)) .* sinc(6+6*f)));

figure (2);
stem(f,abs(M_Analytical));
title("====>> M_Analytical(f) <<====");
grid on;

############## Q3 #############

M_numerical = fftshift(fft(m))*ts;
figure(3);
stem(f,abs(M_numerical));
title("========>> M_numerical(f) <<========")
grid on;

figure(4);
stem(f,abs(M_Analytical));
hold on;
stem(f+100,abs(M_numerical));
legend("M_Analytical FT" , "M_numerical FT");
grid on;

############## Q4 #############

P = abs(M_Analytical).^2;

P_max = max(P);
P_max_index = 4901;
threshold = 0.05 * P_max;

BW = abs(f(P_max_index + find( P(P_max_index:end) <= threshold,1) -1))




































































































