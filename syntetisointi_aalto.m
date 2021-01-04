
%
% Syntetisoi sanan aalto syntetisoimalla sanan eri kirjainten foneemit erikseen.
% 
%

% Mikko Seppi
% 10.10.2020

clear all

%

[y,Fs] = audioread('aalto.wav');
f = 200;
%soundsc(y,Fs)

%Foneemien pituuksien maaritys
Ta=0.5;
Tl=0.15;
Tt=0.005;
To =0.25;

%perustaajuus lahde signaalille
F= 200;
%naytteenottotaajuuden maaritys
fs = Fs;

%% A foneemin syntetisointi

%Formananttin keskitaajuudet
fc1=720;
fc2 = 1240;
fc3 = 2455;


G=1;

%kaistanleveydet
Bw1=130;
Bw2=70;
Bw3=160;

%maaritete‰‰n R arvot
R1=exp(-pi*Bw1/fs);
R2=exp(-pi*Bw2/fs);
R3=exp(-pi*Bw3/fs);

%maaritete‰‰n Wc arvot
wc1= 2*pi*fc1/fs;
wc2=2*pi*fc2/fs;
wc3=2*pi*fc3/fs;

%suodattimen b(kaikilla sama)
b= [G 0 0];

%suodattimien a komponentit
a1= [1 -2*R1*cos(wc1) R1^2];
a2 = [1 -2*R2*cos(wc2) R2^2];
a3 = [1 -2*R3*cos(wc3) R3^2];


%Soudattimien spektrien kuvaajat
w = logspace(-1,1);

%%fc1
figure(1)
subplot(2,2,1);
hh = freqs(b,a1,w);
mag=abs(hh);
loglog(w,mag)

xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('a suodatin taajuusvaste fc1')

%fc2
subplot(2,2,2);
hh = freqs(b,a2,w);
mag=abs(hh);
loglog(w,mag)


xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('a suodatin taajuusvaste fc2')

%fc3
subplot(2,2,3);
hh = freqs(b,a3,w);
mag=abs(hh);
loglog(w,mag)


xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('a suodatin taajuusvaste fc3')


%suodatettavan signaalin luominen

nsamps = fs*Ta; %%naytteiden maara

 
w0T = 2*pi*F/fs; 

nharm = floor((fs/2)/F); 

s = zeros(1,nsamps);

n = 0:(nsamps-1);

%luodaan itse signaali
for i=1:nharm
    s = s + cos(i*w0T*n);
end
s = s/max(s);

%kuvaajat glottipulssin aika ja taajustasosta
figure(2)
subplot(2,1,1)
yg = fft(s);
fg = (0:nsamps-1)*(fs/nsamps);     
power = abs(yg).^2/nsamps; 
plot(fg,power)
title('glottipulssi taajuustasossa')
subplot(2,1,2)
plot(0:nsamps-1,s)
title('glottipulssi aikatasossa')


%filtteroidaan signaali, jotta siita tulee foneemi a:lle
est_a = filter(b,a1,filter(b,a2,filter(b,a3,s)));




%% o foneemin syntetisointi

G=1;

%Formananttin keskitaajuudet
fc1=905;
fc2 =515;
fc3=2430;

%kaistanleveydet
Bw1=130;
Bw2=70;
Bw3=160;

%maaritete‰‰n R arvot
R1=exp(-pi*Bw1/fs);
R2=exp(-pi*Bw2/fs);
R3=exp(-pi*Bw3/fs);

%maaritete‰‰n Wc arvot
wc1= 2*pi*fc1/fs;
wc2= 2*pi*fc2/fs;
wc3= 2*pi*fc3/fs;

%suodattimen b(kaikilla sama)
b= [G 0 0];

%suodattimien a komponentit
a1= [1 -2*R1*cos(wc1) R1^2];
a2 = [1 -2*R2*cos(wc2) R2^2];
a3 = [1 -2*R3*cos(wc3) R3^2];

%Soudattimien spektrien kuvaajat
w = logspace(-1,1);

%%fc1
figure(3)
subplot(2,2,1)
hh = freqs(b,a1,w);
mag=abs(hh);

loglog(w,mag)

xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('o suodatin spektri fc1')

%fc2
subplot(2,2,2)
hh = freqs(b,a2,w);
mag=abs(hh);
loglog(w,mag)


xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('o suodatin spektri fc2')

%fc3
subplot(2,2,3)
hh = freqs(b,a3,w);
mag=abs(hh);

loglog(w,mag)

xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('o suodatin spektri fc3')



%suodatettavan signaalin luominen

nsamps = fs*To; %%naytteiden maara

w0T = 2*pi*F/fs; 

nharm = floor((fs/2)/F); 

s = zeros(1,nsamps);

n = 0:(nsamps-1);

%luodaan itse signaali
for i=1:nharm
    s = s + cos(i*w0T*n);
end

s = s/max(s);

%filtteroidaan signaali, jotta siita tulee foneemi o:lle
est_o = filter(b,a1,filter(b,a2,filter(b,a3,s)));



% luetaan l_long.waw, jota kaytetaan avuksi l j t syntetsointiin.
[kon,FS] = audioread('l_long.wav');


%% l syntetisointi

%signaalin luonti
nsamps = fs*Tl; %%naytteiden maara

w0T = 2*pi*F/fs; 

nharm = floor((fs/2)/F); 

s = zeros(1,nsamps);

n = 0:(nsamps-1);


for i=1:nharm
    s = s + cos(i*w0T*n);
end

s = s/max(s);





%luodaan lpc functiota ja l_long.waw tiedostoa kayttaen suodatin l foneemin
%luomiseen
[a,g] = lpc(kon,3);

%kuvaaja suodattimen spektrista
figure(4)
hh = freqs(1,a,w);
mag=abs(hh);

loglog(w,mag)

xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('l suodatin spektri')

%filtteroidaan signaali, jotta siita tulee foneemi l:lle
est_l =  filter(1,a,s);


%% t syntetisointi

%signaalin luonti
nsamps = fs*Tt; %%naytteiden maara

w0T = 2*pi*F/fs; 

nharm = floor((fs/2)/F); 

s = zeros(1,nsamps);

n = 0:(nsamps-1);

%luodaan itse signaali
for i=1:nharm
    s = s + cos(i*w0T*n);
end

s = s/max(s);

%kaytetaan samaa suodatinta kuin l foneemien tekemisessa
est_t = filter(1,a,s);

%luodaan kohinaa
ran =3* randn(1,length(est_t));

%lisataan kohina suodatettuun signaaliin, jotta saadaan t:een foneemi
est_t = ran + est_t;

%viela lisataan nollia eteen loppu tuloksen parantamiseksi
est_t= horzcat(zeros(1,1000),est_t);





%% Lopputulos

%Kestiarvoistetaan viela foneemit
est_a =  est_a./sum(est_a);
est_o =  est_o./sum(est_o);
est_t =  est_t./sum(est_t);
est_l =  est_l./sum(est_l);

%yhdistetaan luodut foneemit toisiinsa
aani= horzcat((est_a),horzcat(est_l,horzcat(est_t,est_o)));
soundsc(aani,fs)




figure(5)
subplot(2,1,1)
plot(1:length(y),y)
title('Malli')
subplot(2,1,2)
plot(1:length(aani),aani)
title('Syntetisoimalla luotu')


%aani kirjoitetaan tiedostoon uusiaalto.wav
%audiowrite('uusiaalto.wav',aani,fs);

