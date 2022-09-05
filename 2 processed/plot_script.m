files = dir('*.mat');

for i = 1:length(files)
    load(files(i).name);
end

%% plot options
lw = 2; % LineWidth
fs = 16; % FontSize

r_model = logspace(-7,4,10);
tau_modelsl = r_model.^2/2/1.3e-10;
tau_modeldl = r_model.^2/2/1.3e-9;

%% m_n-S_por relation
spormn = logspace(-1,2,4);
mnlin1 = 7.*logspace(-3,0,4);
mnlin2 = 2.*logspace(-2,1,4); 
mnlin3 = 3.*logspace(-1,2,4); 

%% sigma''@1Hz-S_por relation
sigma2 = logspace(-4,0,5);
spors2 = logspace(-2,2,5);


%% plots
% tau - r plot
loglog(1e6.*[spheres.r_sref],[spheres.t_m],'o','MarkerSize',8)
hold on
loglog(1e6.*[fractals.r_vref],[fractals.tau_1],'+','MarkerSize',8)
loglog(1e6.*[Nvars.r_vref],[Nvars.tau_1],'bp','MarkerSize',8)
loglog(1e6.*[bvars.r_vref],[bvars.tau_1],'g*','MarkerSize',8)
loglog(1e6.*r_model,tau_modelsl,'k--')
loglog(1e6.*r_model,tau_modeldl,'k-')
hold off

xlabel('grain radius /\mum')
ylabel('\tau /s')
set(gca, 'XLim',[1e-1 1e2])
set(gca, 'YLim',[1e-5 1e1])
grid on

set(gca,'FontSize',fs)
set(gcf,'Position',[400 150 900 700])

lstr1 = {'smooth grains','fractal grains', ...
         'N_{var} rough grain', ...
         'b_{var} rough grain', ...
         '\tau_{grain,Stern}','\tau_{grain,diffuse}'};
legend(lstr1,'location','northwest')

% S_por - m_n plot
figure
loglog([spheres.s_pormu],1e3*[spheres.m_n],'o','MarkerSize',8)
hold on
loglog([fractals.s_pormu],1e3*[fractals.m_n],'+','MarkerSize',8)
loglog([Nvars.s_pormu],1e3*[Nvars.m_n],'p','MarkerSize',8)
loglog([bvars.s_pormu],1e3*[bvars.m_n],'*','MarkerSize',8)
plot(spormn,mnlin1,'k-',spormn,mnlin2,'k--')
hold off

xlabel('S_{por} /1/\mum')
ylabel('m_n /mS/m')
set(gca, 'XLim',[1e-1 1e2])
set(gca, 'YLim',[1e-3 1e2])
grid on

set(gca,'FontSize',fs)
set(gcf,'Position',[400 150 900 700])

lstr2 = {'smooth grains','fractal grains', ...
    'N_{var} rough grains','b_{var} rough grains'};
legend(lstr2,'location','northwest')


% S_por - s'' plot
figure
loglog([spheres.s_pormu],1e3*[spheres.sig1Hz],'o','MarkerSize',8)
hold on
loglog([fractals.s_pormu],1e3*[fractals.sig1Hz],'+','MarkerSize',8)
loglog([Nvars.s_pormu],1e3*[Nvars.sig1Hz],'p','MarkerSize',8)
loglog([bvars.s_pormu],1e3*[bvars.sig1Hz],'*','MarkerSize',8)
plot(spors2,sigma2,'k')
hold off

xlabel('S_{por} /1/\mum')
ylabel('\sigma'''' /mS/m')
set(gca, 'XLim',[1e-1 1e2])
set(gca, 'YLim',[1e-5 1e0])
grid on

set(gca,'FontSize',fs)
set(gcf,'Position',[400 150 900 700])

lstr3 = {'smooth grains','fractal grains', ...
    'N_{var} rough grains','b_{var} rough grains'};
legend(lstr3,'location','northwest')