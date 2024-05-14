clear
clc

num_bits = 64;
bit_stream = randi([0, 1], 1, num_bits);


ts = 0.005;

t = 0:ts:num_bits-ts;



########### AMI_rz line ###########
###################################





ami_rz_code = zeros(1,num_bits);
last_voltage = 1; % Initial voltage state
for i = 1:num_bits
    if bit_stream(i) == 1
        ami_rz_code(i) = last_voltage;
        last_voltage = -last_voltage;
    end
end

ami_signal = repelem(ami_rz_code, 200);


for i = 100:100:length(ami_signal)
    if ami_signal(i) == 1 || ami_signal(i) == -1
        ami_signal(i:i+100) = 0;
    end
end





########### Polar NRZ Encoding ###########
####################################


polar_nrz_code = 2 * bit_stream - 1;


polar_nrz_signal = repelem(polar_nrz_code, 200);


########### Tb Vertical line ###########
########################################
Tb = t(201);
line = linspace(-1.5,1.5,20);
t_line = zeros(1,length(line));


########### AMI Plot ###########
################################
figure(1);
subplot(2, 1, 1);
plot(t,ami_signal, 'b', 'LineWidth', 1);

title('AMI Line Coding (Time Domain)');
xlabel('Time');
ylabel('Voltage');
ylim([-1.5, 1.5]);
hold on

for i = 1:num_bits
    plot(t_line, line ,'|','color',"k",'LineWidth', 0.1);
    t_line = t_line + Tb;    #Return to zero
end
t_line = zeros(1,length(line));
grid on;


########### Polar NRZ Plot ###########
######################################

subplot(2, 1, 2);
plot(t, polar_nrz_signal, 'r', 'LineWidth', 1);
title('Polar NRZ Line Coding (Time Domain)');
xlabel('Time');
ylabel('Voltage');
ylim([-1.5, 1.5]);
hold on;

for i = 1:num_bits
    plot(t_line, line ,'|','color',"k",'LineWidth', 0.1);
    t_line = t_line + Tb;         #non-return to zero
end
grid on;








########### Spectral domain plots ###########
#############################################

ami_ft       = fftshift(fft(ami_signal))*ts;
polar_nrz_ft = fftshift(fft(polar_nrz_signal))*ts;

fs = 64;
df = 0.005;
if (rem(length(t),2)==0)
f= - (0.5 * fs) : df : (0.5*fs - df) ;
else
f= - (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df);
end


figure(2);
subplot(2, 1, 1);
stem(f, abs(ami_ft), 'b');
title('AMI Line Coding (Spectral Domain)');
xlabel('Frequency');
ylabel('Magnitude');
grid on;

subplot(2, 1, 2);
stem(f, abs(polar_nrz_ft), 'r');
title('Polar NRZ Line Coding (Spectral Domain)');
xlabel('Frequency');
ylabel('Magnitude');
grid on;









