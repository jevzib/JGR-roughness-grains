% debye decomposition after Nordsiek and Weller, 2008
clear

%% read simulation results
rfile = 'sigma_abc_20_4_0,8e-6_MW.txt';
sfile = 'surf_abc_20_4_0,8e-6.txt';
vfile = 'vol_abc_20_4_0,8e-6.txt';

a = readtable(rfile);
% result file contains two columns: column 1 stores angular frequencies
% (w) and column 2 holds complex conductivity (Y)
w = a.Var1;
Y = a.Var2;

LFcut = 1e5;
HFcut = 1e7;
n_peaks = 3;

% read volume file
a  = readmatrix(vfile);
vp = min(a(4,:)); % particle volume
vd = max(a(4,:)); % domain volume

% read surface area file
a  = readmatrix(sfile);
sa = a(1,1);

%% Adjust volumetric particle content of simulation
if ~isnan(vd)
    nusim = vp/vd;
else
    L = 5e-5;

    nusim  = vp/(2*pi*L.^3);    
end

nuscale  = 0.7; % targetted volumetric content

% Convert to correct volumetric particle content
fG0 = 1/nusim*(Y-1)./(Y+2);
Y = (1+2*nuscale*fG0)./(1-nuscale*fG0);


% for denormalization
sigma0 = 0.0096485; % fluid conductivity
eps    = 80 * 8.854e-12; % fluid permittivity
% Y_f = sigma0 + 1i.*eps.*w;
Y_f = sigma0;
Y = Y .* Y_f;

%% Fit
% convert to complex resistivity
Z = Y.^-1;

% normalize
R0 = sqrt(real(Z(1))^2 + imag(Z(1))^2);

Rnorm = (R0 - Z)./R0;
RnormLF = Rnorm(w<LFcut);
RnormHF = Rnorm(w<HFcut);

RnLF = [real(RnormLF)' imag(RnormLF)'];
RnHF = [real(RnormHF)' imag(RnormHF)'];

% discrete relaxation times corresponding to bandwith of target spectrum
% with an additional decade at each end
% taus = logspace(-log10(max(w))-1, -log10(min(w))+1, 700);
taus = logspace(-log10(max(w))-1, -log10(min(w))+1, 700);

% system matrices - Ar: real part | Ai: imaginary part
ArLF = (w(w<LFcut).*taus).^2./(1+(w(w<LFcut).*taus).^2);
AiLF = (w(w<LFcut).*taus)./(1+(w(w<LFcut).*taus).^2);

ArHF = (w(w<HFcut).*taus).^2./(1+(w(w<HFcut).*taus).^2);
AiHF = (w(w<HFcut).*taus)./(1+(w(w<HFcut).*taus).^2);

% solve for m (partial chargabilities) with least squares method
[mLF,resnormLF,residualLF] = lsqnonneg([ArLF; AiLF], RnLF(:));
[mHF,resnormHF,residualHF] = lsqnonneg([ArHF; AiHF], RnHF(:));

% impedance with obtained parameters
ZddLF = R0.*(1 - sum(mLF.*(1-1./(1+1i.*w'.*taus'))));
ZddHF = R0.*(1 - sum(mHF.*(1-1./(1+1i.*w'.*taus'))));

% average relaxation time
% tau = exp( sum(mLF.*log(taus'))./sum(mLF) );

%% interpolate time constants and chargeabilities of neighbours
[mrowLF, ~, mvalLF] = find(mLF);
[mrowHF, ~, mvalHF] = find(mHF);

j       = 1;
tau_LF = zeros(1,floor(length(mrowLF)/2));
m_LF   = zeros(1,floor(length(mrowLF)/2));

for i = 1:length(mrowLF)   
    % check for direct neighbour
    if i+1 <= length(mrowLF) && mrowLF(i+1) == mrowLF(i)+1
        % form average of two neighbouring taus and ms
        tau_LF(j) = exp((mvalLF(i)*log(taus(mrowLF(i)))+mvalLF(i+1)*log(taus(mrowLF(i+1)))) ...
                 /(mvalLF(i)+mvalLF(i+1)));
        m_LF(j) = mvalLF(i) + mvalLF(i+1);
        j = j+1;
    elseif i > 1 && mrowLF(i-1) == mrowLF(i)-1
        % do nothing, already averaged
    else
        % no direct neighbour, keep tau and m
        tau_LF(j) = taus(mrowLF(i));
        m_LF(j) = mvalLF(i);
        j = j+1;
    end
end

j       = 1;
tau_HF = zeros(1,floor(length(mrowHF)/2));
m_HF   = zeros(1,floor(length(mrowHF)/2));

for i = 1:length(mrowHF)   
    % check for direct neighbour
    if i+1 <= length(mrowHF) && mrowHF(i+1) == mrowHF(i)+1
        % form average of two neighbouring taus and ms
        tau_HF(j) = exp((mvalHF(i)*log(taus(mrowHF(i)))+mvalHF(i+1)*log(taus(mrowHF(i+1)))) ...
                 /(mvalHF(i)+mvalHF(i+1)));
        m_HF(j) = mvalHF(i) + mvalHF(i+1);
        j = j+1;
    elseif i > 1 && mrowHF(i-1) == mrowHF(i)-1
        % do nothing, already averaged
    else
        % no direct neighbour, keep tau and m
        tau_HF(j) = taus(mrowHF(i));
        m_HF(j) = mvalHF(i);
        j = j+1;
    end
end


tau_int = [tau_HF(tau_HF<LFcut^-1) tau_LF(tau_LF>LFcut^-1)];
m_int   = [m_HF(tau_HF<LFcut^-1) m_LF(tau_LF>LFcut^-1)];

Zdd_ = R0.*(1 - sum(m_int'.*(1-1./(1+1i.*w'.*tau_int'))));

%% mean relaxation time and total chargeability without MW polarization
thershold_MW = 1e-5;

m_t = sum(m_int(tau_int>thershold_MW));
t_m = exp( sum( m_int(tau_int>thershold_MW).*log(tau_int(tau_int>thershold_MW)) ) ./ m_t );

% first polarization peak relaxation time
% tau_1 = tau_int(m_int==max(m_int(tau_int>1e-6)));

tau_n = zeros(1,n_peaks);
for i = 1:n_peaks
    if i == 1
        tau_n(i) = tau_int(m_int==max(m_int(tau_int>thershold_MW)));
    else
        ti = tau_int>thershold_MW & tau_int<10^(log10(tau_n(i-1))-1);
        tau_n(i)= tau_int(m_int==max(m_int(i)));
    end
end

% normalized chargeability 
m_n = m_t * sqrt(real(Y(1))^2 + imag(Y(1))^2 );

%% radii
r_vref = (3/4/pi.*vp).^(1/3); % calculated radius from particle volume
r_sref = sqrt(1/4/pi.*sa); % calculated radius from perticle surface area
r_diff = r_vref - r_sref; % abs. difference
r_reld = abs(r_diff./r_sref); % relative difference

%% specific surface area
s_tot = sa./vp;    % spec. surface area per unit volume
s_por = s_tot./(1-nuscale);  % spec. surface area per unit pore vol.  
s_pormu = s_por*1e-6;               % scaled S_por to 1/mum


%% saving fitresults
name = split(rfile,'.');
name = name{1,:};

fit.file    = name;
fit.w       = w;
fit.s_sim   = Y;
fit.s_fit   = Zdd_.^-1;
fit.m_all   = mvalLF;
fit.m_int   = m_int;
fit.m_t     = m_t;
fit.m_n     = m_n;
fit.t_all   = taus(mrowLF);
fit.t_int   = tau_int;
fit.t_m     = t_m;
fit.tau_1   = tau_n(1);
fit.tau_n   = tau_n;
fit.nusim   = nusim;
fit.nuscale = nuscale;
fit.vp      = vp;
fit.vd      = vd;
fit.surfa   = sa;
fit.w1Hz    = w(77);
fit.sig1Hz  = imag(Y(77));
fit.sigma0  = sqrt(real(Y(1))^2 + imag(Y(1))^2 );
fit.r_vref  = r_vref;
fit.r_sref  = r_sref;
fit.r_diff  = r_diff;
fit.r_reld  = r_reld;
fit.s_tot   = s_tot;
fit.s_por   = s_por;
fit.s_pormu = s_pormu;

save([name '_LFHFfit_iew.mat'], '-struct','fit');

%% plots
figure('Position', [200 200 800 800])

subplot(5,1,1:3)
yyaxis left
semilogx(w,real(Y),'.')
hold on
% semilogx(w,real(ZddLF.^-1),'b:','LineWidth', 1)
% semilogx(w,real(ZddHF.^-1),'r:','LineWidth', 1)
semilogx(w,real(Zdd_.^-1),'LineWidth', 2)
xlabel('\omega')
ylabel('\sigma''')
set(gca,'XScale','log')
set(gca,'Xlim',[min(w) max(w)])
set(gca,'FontWeight', 'bold','FontSize',10)
hold off

yyaxis right
loglog(w,imag(Y),'.')
hold on
% loglog(w,imag(ZddLF.^-1),'b:','LineWidth', 1)
% loglog(w,imag(ZddHF.^-1),'r:','LineWidth', 1)
loglog(w,imag(Zdd_.^-1),'LineWidth', 2)
ylabel('\sigma''''')
set(gca,'YScale','log')
hold off
lstr = {'\sigma'' data','\sigma'' fit','\sigma'''' data', '\sigma'''' fit'};
legend(lstr, 'Location', 'northwest')
title(['\' replace(name,'_',' ')])
set(gca,'FontWeight', 'bold','FontSize',10)
grid on

subplot(5,1,4)
stem(w, abs(imag(Zdd_.^-1)+imag(Y')) ./ imag(Zdd_.^-1) ,'.')
xlabel('\omega')
ylabel('rel. Error \sigma''''')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'Xlim',[min(w) max(w)])
set(gca, 'Ylim', [1e-5 1])
set(gca,'FontWeight', 'bold','FontSize',10)
grid on

subplot(5,1,5)
stem(tau_int,m_int)
hold on
% plot(tau_1,10^(log10(m_int(tau_int==tau_1))+1),'rv','MarkerFaceColor','r')
plot(tau_n,fliplr(10.^(log10(m_int(ismember(tau_int,tau_n)))+1)),'rv','MarkerFaceColor','r')
xlabel('\tau')
ylabel('m')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'Ylim',[min(m_int) 1])
set(gca,'Xlim',[max(w)^-1 min(w)^-1])
set (gca,'Xdir','reverse')
set(gca,'FontWeight', 'bold','FontSize',10)
grid on