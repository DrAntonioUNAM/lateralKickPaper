%%% Biomechanical and Signal-Based Characterization of Karate Lateral 
%%% Kicks Using Videogrammetry Analysis – Kinematic Processing Pipeline
%
% Author: Dr. Luis Antonio Aguilar
%
% Description:
%   This script processes 2D kinematic trajectories exported from Kinovea
%   to compute velocity, acceleration, jerk, biomechanical metrics, and
%   time–frequency representations using wavelet analysis.
%
% Input:
%   - Exactly 3 XLS/XLSX files per condition (Kinovea exports)
% Output:
%   - Structured variables with aligned kinematics and summary metrics
%
% Notes:
%   - Raw data are never modified
%   - Processing is fully automatic after folder definition
%
%
%%
clear all
clc
%% Parameters

% --- Temporal normalization
N_interp = 500;     % uniform samples

% ---Sampling of video (processed in (Kinovea)
fs       = 30;      % Hz 

% --- Onset detection
onset_win  = 35; % sliding window length
onset_prct = 99; %percentile for noise-based threshold

% --- Savitzky-Golay
sgolayOrder = 3;  
sgolayFrame = 21; 

% --- Butterworth filter
butter_order = 4;
butter_fc    = 1;   % Hz

%--- Onset detection parameters (trial-dependent)
%P1-Corta
%
k0_vec  = [10 10 8];     % number of initial samples used for noise estimation
win_vec = [30 30 25];    % persistence window length (samples)
offset_hold_vec = [8 8 1];   % uno por patada
%}
%P3-media
%{
k0_vec  = [80 200 200];     % number of initial samples used for noise estimation
win_vec = [5 5 5];    % persistence window length (samples)
offset_hold_vec = [8 8 1];   % uno por patada
%}
%p1larga
%{
k0_vec  = [10 2 10];     % number of initial samples used for noise estimation
win_vec = [30 2 35];    % persistence window length (samples)
offset_hold_vec = [8 8 1];   % uno por patada
%}

%p3larga
%{
k0_vec  = [80 22 100];     % number of initial samples used for noise estimation
win_vec = [5 14 5];    % persistence window length (samples)
offset_hold_vec = [8 8 1];   % uno por patada
%}

%% 1.- File ingestion and standardization
rootFolder   = "D:\programacionKarate";
sourceFolder = fullfile(rootFolder,"Data","KINOVEA","FRONT","LARGA","P3");

workingFolder = "PARTICIPANTE_3";
tempFolder    = fullfile("temp","posiciones","Larga");

savePath     = fullfile(rootFolder, workingFolder);
saveTempPath = fullfile(rootFolder, workingFolder, tempFolder);

baseFigPath = "D:\programacionKarate\karate-lateral-kick-dataset\imagenes";

folderXXX = "LARGA";   % <-- editable
folderYYY = "P3";   % <-- editable

figSavePath = fullfile(baseFigPath, folderXXX, folderYYY);

% --- Create working directory
if ~exist(savePath,'dir')
    [status, message, ~] = mkdir(savePath);
    if status == 0
        error(message)
    end
end

% --- Create temporary directory
if ~exist(saveTempPath,'dir')
    [status, message, ~] = mkdir(saveTempPath);
    if status == 0
        error(message)
    end
end

% Crear carpeta si no existe
if ~exist(figSavePath,'dir')
    mkdir(figSavePath);
end

% --- Search Kinovea files
filesXLS  = dir(fullfile(sourceFolder,"*.xls"));
filesXLSX = dir(fullfile(sourceFolder,"*.xlsx"));
files     = [filesXLS; filesXLSX];
% --- Remove Excel temporary files
files = files(~startsWith({files.name},"~$"));
% --- Safety check
assert(numel(files) == 3, ...
    "Expected exactly 3 files in %s", sourceFolder);
% --- Sort files alphabetically for reproducibility
[~, idx] = sort({files.name});
files = files(idx);
% --- Copy and rename files to standardized format
for i = 1:3
    newName = "P3_intento_" + i + ".xlsx";
    src  = fullfile(files(i).folder, files(i).name);
    dest = fullfile(saveTempPath, newName);
    copyfile(src, dest);
end
% --- Create datastore for batch processing
ds = datastore(saveTempPath);
disp("Standardized files:")
disp(ds.Files)


%% 2.- Time normalization and interpolation

% --- Initialize structure for raw data
dataRaw = struct([]);
% --- Read each standardized file
for i = 1:numel(ds.Files)

    % Detect import options and preserve original variable names
    opts = detectImportOptions(ds.Files{i});
    opts.VariableNamingRule = 'preserve';
    opts.SelectedVariableNames = {'min','seg','milseg','x','y'};

    % Read table
    T = readtable(ds.Files{i}, opts);

    % Store raw data
    dataRaw(i).participant = "P3";
    dataRaw(i).trial       = i;
    dataRaw(i).table       = T;
end
% --- Initialize structure for temporally normalized data
DataTemporal = struct([]);
for i = 1:numel(dataRaw)

    % Original table
    T = dataRaw(i).table;

    % --- Time vector (seconds)
    t = T.min * 60 + T.seg + T.milseg / 100;

    % Re-reference time to t = 0
    t = t - t(1);

    % --- Position signals
    x = T.x;
    y = T.y;

    % --- Uniform time base
    t_u = linspace(0, t(end), N_interp)';

    % --- Interpolation (shape-preserving)
    x_u = interp1(t, x, t_u, 'pchip');
    y_u = interp1(t, y, t_u, 'pchip');

    % --- Temporal consistency diagnostics
    dt = diff(t);
    fprintf('Trial %d: mean dt = %.4f s | std = %.4f s\n', ...
            i, mean(dt), std(dt));

    % --- Store normalized data
    DataTemporal(i).t = t_u;
    DataTemporal(i).x = x_u;
    DataTemporal(i).y = y_u;
end
% --- Quick visualization (sanity check)
%figure01
figure
hold on
for i = 1:numel(DataTemporal)
    plot(DataTemporal(i).t, DataTemporal(i).x, 'LineWidth', 1.2)
end
xlabel('Time (s)')
ylabel('Position X (pixels)')
legend('Trial 1','Trial 2','Trial 3','Location','best')
grid on
hold off

figName = sprintf("figure01");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);



%% 3.- Automatic onset detection


for i = 1:3

    t = DataTemporal(i).t;
    x = DataTemporal(i).x;
    y = DataTemporal(i).y;

    % ---- velocity
    dt = mean(diff(t));
    vx = gradient(x,dt);
    vy = gradient(y,dt);
    v  = sqrt(vx.^2 + vy.^2);

    % ---- displacement (internal exploration)
    dx = diff(x);
    dy = diff(y);
    dpos = sqrt(dx.^2 + dy.^2);

    % ---- trial-specific parameters
    k0  = k0_vec(i);
    win = win_vec(i);

    % ---- noise estimation
    thV = prctile(v(1:k0), 99);

    % ---- persistence-based detector
    cond = movmean(v > thV, win) > 0.8;
    idx0 = find(cond,1,'first');

    % ---- store results
    DataTemporal(i).idx_onset = idx0;
    DataTemporal(i).t_onset  = t(idx0);
    DataTemporal(i).v        = v;
    DataTemporal(i).dpos     = dpos;
    DataTemporal(i).thV      = thV;
    DataTemporal(i).k0       = k0;
    DataTemporal(i).win      = win;
end

% --- Time realignment relative to onset
for i = 1:3
    idx0 = DataTemporal(i).idx_onset;
    DataTemporal(i).t_evt = ...
        DataTemporal(i).t - DataTemporal(i).t(idx0);
end
%figure02
figure
hold on
for i = 1:3
    plot(DataTemporal(i).t_evt, DataTemporal(i).x, 'LineWidth',1.5)
end
xline(0,'k--','Onset')
xlabel('Time relative to onset (s)')
ylabel('X position')
legend('Trial 1','Trial 2','Trial 3','Location','best')
grid on
hold off
figName = sprintf("figure02");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);


%figure03
figure
hold on
for i = 1:3
    plot(DataTemporal(i).t_evt, DataTemporal(i).v, 'LineWidth',1.5)
end
xline(0,'k--','Start')
xlabel('Time relative to onset (s)')
ylabel('Velocity')
legend('Trial  1','Trial 2','Trial 3','Location','best')
grid on
hold off
figName = sprintf("figure03");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

%% 4.- Kinematic computation (SG and Butterworth)

% 4.- Kinematic computation (Savitzky–Golay + Butterworth)

for i = 1:3

    % --- Time relative to onset
    t  = DataTemporal(i).t_evt;
    dt = mean(diff(t));

    % --- Relative positions (recommended reference frame)
    idx0  = DataTemporal(i).idx_onset;
    x_rel = DataTemporal(i).x - DataTemporal(i).x(idx0);
    y_rel = DataTemporal(i).y - DataTemporal(i).y(idx0);

    % --- Displacement magnitude (relative)
    disp_rel = sqrt(x_rel.^2 + y_rel.^2);

    % --- Savitzky–Golay smoothing (geometric smoothing)
    x_s = sgolayfilt(x_rel, sgolayOrder, sgolayFrame);
    y_s = sgolayfilt(y_rel, sgolayOrder, sgolayFrame);

    % --- Safety: remove possible gaps before filtering
    x_s = fillmissing(x_s,'linear');
    y_s = fillmissing(y_s,'linear');

    % --- Butterworth LOW-PASS (spectral control)
    fs_eff = fs;  % Kinovea sampling frequency
    [b,a]  = butter(butter_order, butter_fc/(fs_eff/2), 'low');

    % Zero-phase filtering of positions
    x_f = filtfilt(b,a,x_s);
    y_f = filtfilt(b,a,y_s);

    % --- Velocity (px/s)
    vx = gradient(x_f, dt);
    vy = gradient(y_f, dt);
    v  = sqrt(vx.^2 + vy.^2);

    % --- Acceleration (px/s^2)
    ax = gradient(vx, dt);
    ay = gradient(vy, dt);
    a  = sqrt(ax.^2 + ay.^2);

    % --- Jerk (px/s^3)
    jx = gradient(ax, dt);
    jy = gradient(ay, dt);
    j  = sqrt(jx.^2 + jy.^2);

    % --- Store results (traceable pipeline)
    DataTemporal(i).x_rel = x_rel;
    DataTemporal(i).y_rel = y_rel;
    DataTemporal(i).disp  = disp_rel;

    DataTemporal(i).x_s = x_s;
    DataTemporal(i).y_s = y_s;

    DataTemporal(i).x_f = x_f;
    DataTemporal(i).y_f = y_f;

    DataTemporal(i).vx = vx; DataTemporal(i).vy = vy; DataTemporal(i).v = v;
    DataTemporal(i).ax = ax; DataTemporal(i).ay = ay; DataTemporal(i).a = a;
    DataTemporal(i).jx = jx; DataTemporal(i).jy = jy; DataTemporal(i).j = j;
end

% Trajectory (filtered positions)
%figure04
figure
hold on
for i = 1:3
    plot(DataTemporal(i).x_f, DataTemporal(i).y_f, 'LineWidth',1.5)
end
axis equal
xlabel('Relative X (pixels)')
ylabel('Relative Y (pixels)')
legend('Trial 1','Trial 2','Trial 3','Location','best')
title('Filtered trajectory (SG + Butterworth)')
grid on
hold off
figName = sprintf("figure04");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

% Velocity
%figure05
figure
hold on
for i = 1:3
    plot(DataTemporal(i).t_evt, DataTemporal(i).v, 'LineWidth',1.5)
end
xline(0,'k--','Onset')
xlabel('Time relative to onset (s)')
ylabel('Velocity (px/s)')
legend('Trial 1','Trial 2','Trial 3','Location','best')
title('Resultant velocity (SG + Butterworth)')
grid on
hold off
figName = sprintf("figure05");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

% Acceleration
%figure06
figure
hold on
for i = 1:3
    plot(DataTemporal(i).t_evt, DataTemporal(i).a, 'LineWidth',1.5)
end
xline(0,'k--','Onset')
xlabel('Time relative to onset (s)')
ylabel('Acceleration (px/s^2)')
legend('Trial 1','Trial 2','Trial 3','Location','best')
title('Resultant acceleration (filtered)')
grid on
hold off
figName = sprintf("figure06");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

% Jerk
%figure07
figure
hold on
for i = 1:3
    plot(DataTemporal(i).t_evt, DataTemporal(i).j, 'LineWidth',1.5)
end
xline(0,'k--','Onset')
xlabel('Time relative to onset (s)')
ylabel('Jerk (px/s^3)')
legend('Trial 1','Trial 2','Trial 3','Location','best')
title('Resultant jerk (filtered)')
grid on
hold off
figName = sprintf("figure07");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);


%% 5.- Gesture segmentation and metrics

for i = 1:3

    v   = DataTemporal(i).v;
    t   = DataTemporal(i).t_evt;
    thV = DataTemporal(i).thV;

    win  = DataTemporal(i).win;          
    idx0 = DataTemporal(i).idx_onset;
    offset_hold = offset_hold_vec(i);    

    % --- actividad (misma lógica)
    active = movmean(v > thV, win) > 0.8;

    % --- buscar secuencia estable de inactividad
    inactive = movmean(active == 0, offset_hold) > 0.9;

    % --- buscar fin DESPUÉS del onset
    idxEndRel = find(inactive(idx0:end), 1, 'first');

    if isempty(idxEndRel)
        idxEnd = numel(t);
    else
        idxEnd = idx0 + idxEndRel - 1;
    end

    % --- blindaje
    if idxEnd <= idx0
        idxEnd = numel(t);
    end

    % --- guardar
    DataTemporal(i).idx_end = idxEnd;
    DataTemporal(i).t_end  = t(idxEnd);
    DataTemporal(i).offset_hold = offset_hold;
end


for i = 1:3
    DataTemporal(i).win_active = ...
        DataTemporal(i).idx_onset : DataTemporal(i).idx_end;
end

for i = 1:3
    w = DataTemporal(i).win_active;
    t = DataTemporal(i).t_evt(w);

    DataTemporal(i).duration = t(end) - t(1);
end

for i = 1:3

    t = DataTemporal(i).t_evt;
    v = DataTemporal(i).v;

    idx0  = DataTemporal(i).idx_onset;
    idxEnd = DataTemporal(i).idx_end;
    w = DataTemporal(i).win_active;
    
    %figure08, 09, 10
    figure
    hold on

    % Señal completa
    plot(t, v, 'Color',[0.7 0.7 0.7], 'LineWidth',1)

    % Ventana activa
    plot(t(w), v(w), 'b', 'LineWidth',2)

    % Inicio y fin
    xline(t(idx0), 'k--','Start','LineWidth',1.5)
    xline(t(idxEnd), 'r--','End','LineWidth',1.5)

    xlabel('Time relative to onset (s)')
    ylabel('Velocity')
    title(sprintf('Automatic gesture segmentation – Trial %d', i))
    grid on
    hold off
    figName = sprintf("figure", i+7);
    exportgraphics(gcf, ...
        fullfile(figSavePath, figName + ".png"), ...
        'Resolution',300);

end


%% 6. Biomechanical metrics extraction (active window only)

for i = 1:3

    w = DataTemporal(i).win_active;

    % --- señales dentro del gesto
    t = DataTemporal(i).t_evt(w);
    v = DataTemporal(i).v(w);
    a = DataTemporal(i).a(w);
    j = DataTemporal(i).j(w);
    d = DataTemporal(i).disp(w);

    % --- duración efectiva
    DataTemporal(i).dur = t(end) - t(1);

    % --- picos y tiempos
    [DataTemporal(i).v_max, iv] = max(v);
    DataTemporal(i).t_v_max = t(iv);

    [DataTemporal(i).a_max, ia] = max(a);
    DataTemporal(i).t_a_max = t(ia);

    [DataTemporal(i).j_max, ij] = max(j);
    DataTemporal(i).t_j_max = t(ij);

    [DataTemporal(i).disp_max, id] = max(d);
    DataTemporal(i).t_disp_max = t(id);

    % --- RMS (control global del gesto)
    DataTemporal(i).v_rms = rms(v);
    DataTemporal(i).a_rms = rms(a);
    DataTemporal(i).j_rms = rms(j);
end

T_metrics = table( ...
    (1:3)', ...
    [DataTemporal.duration]', ...
    [DataTemporal.disp_max]', ...
    [DataTemporal.t_disp_max]', ...
    [DataTemporal.v_max]', ...
    [DataTemporal.t_v_max]', ...
    [DataTemporal.a_max]', ...
    [DataTemporal.t_a_max]', ...
    [DataTemporal.j_max]', ...
    [DataTemporal.t_j_max]', ...
    [DataTemporal.v_rms]', ...
    [DataTemporal.a_rms]', ...
    [DataTemporal.j_rms]', ...
    'VariableNames',{ ...
        'Kick','Duration_s', ...
        'DispMax_px','TimeDispMax_s', ...
        'Vmax_px_s','TimeVmax_s', ...
        'Amax_px_s2','TimeAmax_s', ...
        'Jmax_px_s3','TimeJmax_s', ...
        'Vrms_px_s','Arms_px_s2','Jrms_px_s3'} );
%figure11
figure
hold on
for i = 1:3
    w = DataTemporal(i).win_active;
    plot(DataTemporal(i).t_evt(w), DataTemporal(i).a(w), 'LineWidth',1.8)
end
xline(0,'k--','Onset')
xlabel('Time relative to onset (s)')
ylabel('Acceleration (px/s^2)')
title('Acceleration – active window')
legend('Trial 1','Trial 2','Trial 3','Location','best')
grid on
hold off
figName = sprintf("figure11");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

%figure12
figure
hold on
for i = 1:3
    w = DataTemporal(i).win_active;
    plot(DataTemporal(i).t_evt(w), DataTemporal(i).j(w), 'LineWidth',1.8)
end
xline(0,'k--','Onset')
xlabel('Time relative to onset (s)')
ylabel('Jerk (px/s^3)')
title('Jerk – active window')
legend('Trial 1','Trial 2','Trial 3','Location','best')
grid on
hold off
figName = sprintf("figure12");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

%% 7. Time–frequency analysis (CWT) – active window only

% 7.1 Continuous Wavelet Transform (Morlet)
for i = 1:3

    % --- señal y tiempo SOLO del gesto
    w = DataTemporal(i).win_active;
    t = DataTemporal(i).t_evt(w);
    a = DataTemporal(i).a(w);
    a = a - mean(a);      % eliminar offset DC


    % --- frecuencia de muestreo efectiva
    fs_eff = fs;   % 30 Hz (Kinovea)

    % --- CWT Morlet
    [cfs, frq] = cwt(a, fs_eff, 'amor');

    % --- guardar resultados
    DataTemporal(i).cwt.cfs = cfs;
    DataTemporal(i).cwt.frq = frq;
    DataTemporal(i).cwt.t   = t;

    % --- Visualización
    %figure13, 14, 15
    figure
    surface(t, frq, abs(cfs))
    shading interp
    axis tight
    axis square
    set(gca,'yscale','log')
    xlabel('Time relative to onset (s)')
    ylabel('Frequency (Hz)')
    title(sprintf('Morlet CWT – Acceleration (Trial %d)', i))
    colorbar
    figName = sprintf("figure%d", i+12);
    exportgraphics(gcf, ...
        fullfile(figSavePath, figName + ".png"), ...
        'Resolution',300);
end
%% 7.3 Wavelet Mexican Hat (control morfológico)

for i = 1:3

    % --- ventana activa y señal filtrada
    w = DataTemporal(i).win_active;
    a = DataTemporal(i).a(w);
    a = a - mean(a);   % eliminar DC

    fs_eff = fs;
    Ts = 1/fs_eff;

    % --- rango de pseudo-frecuencia
    fmin = 0.5;
    fmax = 10;

    fc = centfrq('mexh');

    s_min = fc/(fmax*Ts);
    s_max = fc/(fmin*Ts);
    scales = logspace(log10(s_min), log10(s_max), 64);
    
    fprintf('Trial %d | length(a) = %d | rms(a) = %.3f\n', ...
        i, numel(a), rms(a))

    wt = cwtft({a, Ts}, 'wavelet','mexh', 'scales', scales);
    freq = fc./(scales*Ts);

    % --- guardar
    DataTemporal(i).cwt_mexh.cfs  = wt.cfs;
    DataTemporal(i).cwt_mexh.freq = freq;
    DataTemporal(i).cwt_mexh.t    = DataTemporal(i).t_evt(w);

    % --- visualización
    %figure16,17,18
    figure
    imagesc(DataTemporal(i).t_evt(w), freq, abs(wt.cfs))
    axis square
    axis xy
    set(gca,'yscale','log')
    xlabel('Time relative to onset (s)')
    ylabel('Pseudo-frequency (Hz)')
    title(sprintf('Mexican Hat CWT – Acceleration (Trial %d)', i))
    colorbar
    figName = sprintf("figure%d", i+15);
    exportgraphics(gcf, ...
        fullfile(figSavePath, figName + ".png"), ...
        'Resolution',300);

end

%%

% 7.4 Energía tiempo–frecuencia (descriptor cuantitativo)
for i = 1:3
    E_tf = abs(DataTemporal(i).cwt.cfs).^2;
    DataTemporal(i).E_tf = E_tf;

    % Energía integrada en el tiempo
    DataTemporal(i).E_t = sum(E_tf, 1);
end


% 7.5 Entropía de Shannon (complejidad dinámica)
%figure15
figure
hold on
for i = 1:3
    plot(DataTemporal(i).cwt.t, DataTemporal(i).E_t, 'LineWidth',1.8)
end
xlabel('Time relative to onset (s)')
ylabel('Wavelet energy')
title('Time–frequency energy (Morlet)')
legend('Trial 1','Trial 2','Trial 3','Location','best')
grid on
hold off
figName = sprintf("figure19");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);

%%

%7.6 Frecuencia dominante (trayectoria espectral)
for i = 1:3

    E_tf = DataTemporal(i).E_tf;

    % normalización por tiempo
    P = E_tf ./ sum(E_tf,1);

    % evitar log(0)
    P(P==0) = eps;

    % entropía de Shannon en el tiempo
    H = -sum(P .* log(P), 1);

    DataTemporal(i).H_t = H;
end

%figure16
figure
hold on
for i = 1:3
    plot(DataTemporal(i).cwt.t, DataTemporal(i).H_t, 'LineWidth',1.8)
end
xlabel('Time relative to onset (s)')
ylabel('Shannon entropy')
title('Dynamic complexity of the movement')
legend('Trial 1','Trial 2','Trial 3','Location','best')
grid on
hold off
figName = sprintf("figure20");
exportgraphics(gcf, fullfile(figSavePath, figName + ".png"), ...
               'Resolution', 300);


