clear
close all
clc

%% ===================== settings =====================
stimulus2plot = 'looming';    % looming, sweeping, dimm, white loom, cricket, wntone2
species       = 0;            % 0 = mouse, 1 = Rhabdomys
datapath      = 'N:\Katja';

ZTlist  = [2 10 14 22];
NDlist  = [0 1 3];
NDlight = [0.8 0.66 0.38];              % light level during the test window for ND0, ND1, ND3
NDcol   = [0 0 1; 0 1 1; 0 0 0];        % plot colors for ND0, ND1, ND3

main4        = {'ND3LD14','ND0LD14','ND0LD02','ND3LD02'};
main4_labels = {'ZT14 dim','ZT14 bright','ZT2 bright','ZT2 dim'};

post = 240+30 : 240+150;                % 30-150 min after light onset (used for the model distributions)

% fixed model parameters
mel_weight     = 0.2;    % weight of melanopsin in ipRGC activation
rodcone_weight = 2;      % weight of rods/cones in ipRGC activation
mfilter        = exp(-(0:0.2:1.2));   % melanopsin onset transient

if species == 0
    spc = 'BEHV'; clockshift = 0;   ipnorm = 3.5; threshold = 30;
else
    spc = 'Rhab'; clockshift = 0.6; ipnorm = 3;   threshold = 50;
end

%% ===================== behavioural data =====================
load(fullfile(datapath, ['all', spc, '2_', stimulus2plot]))   % alldata, conditions, plotrangestarts, animalsex, ...

isND        = strncmp(cellstr(conditions), 'ND', 2);
Y           = alldata(isND,:);
cond_data   = cellstr(conditions(isND,:));
starts_data = plotrangestarts(isND);

% sex of each ND trial (1 = female, 2 = male)
if numel(animalsex) == size(alldata,1)
    sex_data = animalsex(isND);
elseif numel(animalsex) == size(Y,1)
    sex_data = animalsex;
else
    error('animalsex does not match the number of trials in alldata');
end
sex_data = sex_data(:);

tmp = find(Y > 120); Y(tmp) = Y(tmp-1);   % (as in original)

% escape latency: first time after stimulus with 5 consecutive samples above threshold
% onsets in the first 10 frames get 10/60 only if there is a clear acceleration;
% otherwise the first later escape (new run starting at frame 11 or later) is used
% W = pre-escape speed: median speed 60 to 10 frames before the detected onset
L = nan(size(Y,1),1);
W = nan(size(Y,1),1);
for y = 1:size(Y,1)
    above = Y(y,:) > threshold;
    early_rejected = false;
    for f = 1:size(Y,2)-243
        if all(above(239+f:243+f))
            if f < 11
                if early_rejected, continue, end
                pre_mean  = mean(Y(y,180:239));   % 1 s before stimulus
                post_peak = max(Y(y,240:359));    % 0-2 s after stimulus
                if post_peak > 1.5*pre_mean && post_peak > threshold
                    L(y) = 10/60;
                    W(y) = median(Y(y,239+f-60:239+f-10));
                    break
                end
                early_rejected = true;            % no clear acceleration: keep searching
            else
                % after a rejected early escape, skip frames that still belong to that run
                if early_rejected && above(238+f), continue, end
                L(y) = f/60;
                W(y) = median(Y(y,239+f-60:239+f-10));
                break
            end
        end
    end
end

%% ===================== fixed circadian traces =====================
x      = -pi : 2*pi/(1440-1) : pi;
SCN    = -sin(x);              SCN    = [SCN(1321:end) SCN SCN];      % starts at ZT22
clockg = -sin(x + clockshift); clockg = [clockg(1321:end) clockg clockg];
SCN    = (SCN + 1)/2;
clockg = (clockg + 1)/2;

%% ===================== model =====================
nZT = numel(ZTlist); nND = numel(NDlist);
model_latency = cell(nZT, nND);
cond_names    = cell(nZT, nND);

figFull  = figure('Name','Model traces, whole duration');
figRange = figure('Name','Model traces, plotrange');
for iz = 1:nZT
    ZT = ZTlist(iz);
    for n = 1:nND
        ND = NDlist(n);
        cond_names{iz,n} = sprintf('ND%dLD%02d', ND, ZT);

        % --- light ---
        light = ones(1,3000)*0.1;
        light(121:840)   = 0.9;
        light(1561:2280) = 0.8;
        start = (ZT+2)*60;
        light(start:start+240) = NDlight(n);
        plotrange = start-239 : start+240;

        % --- melanopsin: contrast-scaled onset transient, then 30% plateau while light > 0.75 ---
        melanopsin = zeros(1,3000);
        id = 1;
        while id < 3000
            idc = find(light(id:end) > 0.75, 1) + id - 1;
            if isempty(idc)
                break
            end
            contrast = (light(idc) - light(idc-1)) / light(idc-1);
            idc2 = find(light(idc:end) < 0.75, 1) + idc - 1;
            if isempty(idc2)
                idc2 = 3000; 
            end
            melanopsin(idc:idc+6)  = contrast * mfilter;
            melanopsin(idc+7:idc2) = contrast * 0.3;
            id = idc2 + 1;
        end

        % --- rods/cones: step-wise tracking of light changes ---
        rodscones = ones(1,3000)*light(1);
        I_now = light(1);
        id = 1;
        while id < 3000
            idc = find(diff(light(id:end)) ~= 0, 1) + id - 1;
            if isempty(idc)
                break
            end
            I_now = (I_now + light(idc+1) - light(idc-1)) * 0.9;
            rodscones(idc:end) = I_now;
            id = idc + 2;
        end

        % --- ipRGC activation and light history (fast rise, slow decay) ---
        ipRGC = mel_weight*melanopsin + rodcone_weight*rodscones;
        mel_del = zeros(1,3000);
        mel_del(1) = ipRGC(1);
        for i = 2:3000
            if ipRGC(i) > mel_del(i-1)
                if ipRGC(i) > ipRGC(i-1)
                    a = 0.7;
                else
                    a = 1/5; 
                end
            else
                a = 1/25;
            end
            mel_del(i) = mel_del(i-1) + a*(ipRGC(i) - mel_del(i-1));
        end

        % --- latency ---
        mel_to_SCN = mel_del(plotrange) / ipnorm;
        if species == 0
            clockgene_current = exp(13 * (1 - 0.22*clockg(plotrange)) .* (1 - clockg(plotrange)) ...
                                    .* (1 - exp(-max(5*mel_to_SCN - 1, 0))));
            base_latency = 0.3; expo = 0.243;
        else
            clockgene_current = exp(-5.5 * clockg(plotrange) .* (1 - exp(-max(7.5*mel_to_SCN - 1, 0))) ...
                                    + 8 * (clockg(plotrange) - SCN(plotrange)));
            base_latency = 0.75; expo = 0.275;
        end
        SCN_current = (SCN(plotrange) .* clockgene_current) ./ ...
                      (SCN(plotrange) .* clockgene_current + (1 - SCN(plotrange)));
        latency = base_latency .* (SCN_current ./ (1 - SCN_current)).^expo;

        model_latency{iz,n} = latency;

         % --- trace plots for the 4 main conditions ---
        k = find(strcmp(main4, cond_names{iz,n}));
        if ~isempty(k)
            st = starts_data(strcmp(cond_data, cond_names{iz,n}));
            leg = {'SCN','clockg','mel\_to\_SCN','trial starts','clockgene\_current','SCN\_current'};

            % whole duration
            figure(figFull), subplot(2,2,k)
            semilogy(SCN), hold on
            semilogy(clockg)
            semilogy(mel_del/ipnorm)
            semilogy(plotrange(st), SCN(plotrange(st)), '*k')
            semilogy(plotrange, clockgene_current, '--m')
            semilogy(plotrange, SCN_current, '-k')
            title(main4_labels{k})
            if k == 1, legend(leg, 'Location','best'), end

            % from 2 h before experiment start (ZT-2) to the end of plotrange, x-axis in ZT
            figure(figRange), subplot(2,2,k)
            sel = 121:numel(plotrange);
            ztx = ZT + (plotrange(sel) - start)/60;
            st2 = st(st >= sel(1));
            semilogy(ztx, SCN(plotrange(sel))), hold on
            semilogy(ztx, clockg(plotrange(sel)))
            semilogy(ztx, mel_to_SCN(sel))
            semilogy(ZT + (plotrange(st2) - start)/60, SCN(plotrange(st2)), '*k')
            semilogy(ztx, clockgene_current(sel), '--m')
            semilogy(ztx, SCN_current(sel), '-k')
            xt = ZT-2 : ZT+4;
            set(gca, 'XTick', xt, 'XTickLabel', string(mod(xt,24)))
            xlim([ZT-2 ZT+4])
            xlabel('ZT (h)')
            title(main4_labels{k})
            if k == 1, legend(leg, 'Location','best'), end
        end
    end
end

%% ===================== data vs model: split violins (4 main conditions) =====================
figure('Name','Data vs model'), hold on
halfwidth = 0.4;
sides = [-1 1];                              % left = data, right = model
cols  = {[0 0.45 0.74], [0.5 0.5 0.5]};
for k = 1:4
    [iz, n] = find(strcmp(cond_names, main4{k}));
    idx = find(strcmp(cond_data, cond_names{iz,n}) & ~isnan(L));
if isempty(idx)
    vals = {[], model_latency{iz,n}(post)'};
else
    vals = {L(idx), model_latency{iz,n}(starts_data(idx))'};
end
    % idx  = find(strcmp(cond_data, main4{k}) & ~isnan(L));
    % vals = {L(idx), model_latency{iz,n}(starts_data(idx))'};   % model sampled at each trial's time
    for s = 1:2
        v = vals{s}(:);
        if isempty(v)
            continue
        end
        if numel(v) > 1 && max(v) > min(v)
            yi = linspace(min(v), max(v), 100);
            f  = ksdensity(v, yi);
            f  = f/max(f) * halfwidth;
            fill([k + sides(s)*f, k*ones(1,100)], [yi, fliplr(yi)], cols{s}, ...
                 'FaceAlpha', 0.25, 'EdgeColor', cols{s})
            q = quantile(v, [0.25 0.5 0.75]);
            wq = interp1(yi, f, q);
            plot([k k + sides(s)*wq(2)], [q(2) q(2)], '-',  'Color', cols{s}, 'LineWidth', 2)
            plot([k k + sides(s)*wq(1)], [q(1) q(1)], '--', 'Color', cols{s})
            plot([k k + sides(s)*wq(3)], [q(3) q(3)], '--', 'Color', cols{s})
        end
        scatter(k + sides(s)*(0.03 + 0.25*rand(size(v))), v, 12, cols{s}, 'filled', 'MarkerFaceAlpha', 0.6)
    end
end
h1 = fill(nan, nan, cols{1}, 'FaceAlpha', 0.25, 'EdgeColor', cols{1});
h2 = fill(nan, nan, cols{2}, 'FaceAlpha', 0.25, 'EdgeColor', cols{2});
legend([h1 h2], {'data','model'}, 'Location', 'best')
set(gca, 'XTick', 1:4, 'XTickLabel', main4_labels, 'TickDir', 'out')
ylabel('latency')
xlim([0.4 4.6])
box off


%% ===================== model only: bar graph, all conditions =====================
med = zeros(nZT, nND); q25 = med; q75 = med;
for iz = 1:nZT
    for n = 1:nND
        q = quantile(model_latency{iz,n}(post), [0.25 0.5 0.75]);
        q25(iz,n) = q(1); med(iz,n) = q(2); q75(iz,n) = q(3);
    end
end

figure('Name','Model latency, all conditions'), hold on
ZTcol = [0.55 0.8 1; 1 0.6 0.1; 0 0.3 0.8; 0.5 0.2 0.6];   % ZT2 light blue, ZT10 orange, ZT14 blue, ZT22 purple
w  = 0.18;                                   % bar width
hb = gobjects(1, nZT);
for n = 1:nND
    for iz = 1:nZT
        xc = n + (iz - (nZT+1)/2) * (w + 0.03);
        hb(iz) = fill(xc + [-w -w w w]/2, [q25(iz,n) q75(iz,n) q75(iz,n) q25(iz,n)], ZTcol(iz,:), ...
                      'FaceAlpha', 0.35, 'EdgeColor', ZTcol(iz,:));
        plot(xc + [-w w]/2, [med(iz,n) med(iz,n)], '-', 'Color', ZTcol(iz,:), 'LineWidth', 2)
    end
end
set(gca, 'XTick', 1:nND, 'XTickLabel', strcat('ND', string(NDlist)), 'TickDir', 'out')
legend(hb, strcat('ZT', string(ZTlist)), 'Location', 'best')
ylabel('model latency (median, 25–75%)')
xlim([0.5 nND + 0.5])
box off

%% ===================== data vs model: bar graph, all conditions =====================
figure('Name','Data vs model latency, all conditions','position',[402.6000 236.2000 413.6000 420])
hold on
wb = 0.08;                                   % width of each bar
for n = 1:nND
    for iz = 1:nZT
        xc = n + (iz - (nZT+1)/2) * 0.22;    % centre of this condition
        % data (left bar), only if this condition exists in the data
        v = L(strcmp(cond_data, cond_names{iz,n}) & ~isnan(L));
        if ~isempty(v)
            q = quantile(v, [0.25 0.5 0.75]);
            fill(xc - wb/2 - 0.005 + [-wb -wb wb wb]/2, [q(1) q(3) q(3) q(1)], ZTcol(iz,:), ...
                 'FaceAlpha', 0.6, 'EdgeColor', ZTcol(iz,:))
            plot(xc - wb/2 - 0.005 + [-wb wb]/2, [q(2) q(2)], '-k', 'LineWidth', 2)
        end
        % model (right bar)
        fill(xc + wb/2 + 0.005 + [-wb -wb wb wb]/2, [q25(iz,n) q75(iz,n) q75(iz,n) q25(iz,n)], ZTcol(iz,:), ...
             'FaceAlpha', 0.15, 'EdgeColor', ZTcol(iz,:), 'LineStyle', '--')
        plot(xc + wb/2 + 0.005 + [-wb wb]/2, [med(iz,n) med(iz,n)], '-', 'Color', ZTcol(iz,:), 'LineWidth', 2)
    end
end
hz = gobjects(1, nZT);
for iz = 1:nZT, hz(iz) = fill(nan, nan, ZTcol(iz,:), 'FaceAlpha', 0.6, 'EdgeColor', ZTcol(iz,:)); end
hd = fill(nan, nan, [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', [0.5 0.5 0.5]);
hm = fill(nan, nan, [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');
legend([hz hd hm], [cellstr(strcat('ZT', string(ZTlist))), {'data','model'}], 'Location', 'best')
set(gca, 'XTick', 1:nND, 'XTickLabel', strcat('ND', string(NDlist)), 'TickDir', 'out')
ylabel('latency (median, 25–75%)')
xlim([0.5 nND + 0.5])
box off

%% ===================== data vs model: split violins, all conditions =====================
NDorder = [3 0 1];                           % plotting order of ND groups
gap     = 1.4;                               % spacing between ND groups
hw      = 0.09;                              % max half-width of each violin
ybreak  = 3.5;                               % y-axis break
datacol = [0.6 0.6 0.6];                     % data violins: gray

ytop = max(L);
for c = 1:numel(model_latency), ytop = max(ytop, max(model_latency{c}(post))); end
ytop = max(ytop*1.05, ybreak + 0.5);

fig   = figure('Name','Data vs model latency, all conditions (violins)');
axBot = axes(fig, 'Position', [0.1 0.1 0.85 0.64]); hold(axBot, 'on')
axTop = axes(fig, 'Position', [0.1 0.79 0.85 0.15]); hold(axTop, 'on')
for ax = [axBot axTop]
    for g = 1:numel(NDorder)
        n = find(NDlist == NDorder(g));
        for iz = 1:nZT
            xc   = 1 + (g-1)*gap + (iz - (nZT+1)/2) * 0.22;
            vals = {L(strcmp(cond_data, cond_names{iz,n}) & ~isnan(L)), ...   % left: data
                    model_latency{iz,n}(post)'};                              % right: model
            for s = 1:2
                v = vals{s}(:);
                if isempty(v), continue, end
                if s == 1, fc = datacol;     lc = [0 0 0];     fa = 0.5;
                else,      fc = ZTcol(iz,:); lc = ZTcol(iz,:); fa = 0.35;
                end
                if numel(v) > 1 && max(v) > min(v)
                    yi = linspace(min(v), max(v), 100);
                    f  = ksdensity(v, yi);
                    f  = f/max(f) * hw;
                    fill(ax, [xc + sides(s)*f, xc*ones(1,100)], [yi, fliplr(yi)], fc, ...
                         'FaceAlpha', fa, 'EdgeColor', fc)
                    q  = quantile(v, [0.25 0.5 0.75]);
                    wq = interp1(yi, f, q);
                    plot(ax, [xc xc + sides(s)*wq(2)], [q(2) q(2)], '-',  'Color', lc, 'LineWidth', 2)
                    plot(ax, [xc xc + sides(s)*wq(1)], [q(1) q(1)], '--', 'Color', lc)
                    plot(ax, [xc xc + sides(s)*wq(3)], [q(3) q(3)], '--', 'Color', lc)
                end
                if s == 1
                    scatter(ax, xc - (0.01 + 0.06*rand(size(v))), v, 8, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.6)
                end
            end
        end
    end
end

xlims = [0.5, 1 + (numel(NDorder)-1)*gap + 0.5];
set(axBot, 'YLim', [0 ybreak], 'XLim', xlims, 'XTick', 1 + (0:numel(NDorder)-1)*gap, ...
    'XTickLabel', strcat('ND', string(NDorder)), 'TickDir', 'out', 'Box', 'off')
set(axTop, 'YLim', [ybreak ytop], 'XLim', xlims, 'XTick', [], 'XColor', 'none', 'TickDir', 'out', 'Box', 'off')
ylabel(axBot, 'latency')
annotation(fig, 'line', [0.09 0.11], [0.735 0.755])   % axis-break marks
annotation(fig, 'line', [0.09 0.11], [0.775 0.795])

hz = gobjects(1, nZT);
for iz = 1:nZT, hz(iz) = fill(axBot, nan, nan, ZTcol(iz,:), 'FaceAlpha', 0.35, 'EdgeColor', ZTcol(iz,:)); end
hd = fill(axBot, nan, nan, datacol, 'FaceAlpha', 0.5, 'EdgeColor', datacol);
legend(axBot, [hd hz], [{'data'}, cellstr(strcat('model ZT', string(ZTlist)))], 'Location', 'northwest')

figure('Name','Data vs model latency, all conditions (violins v2)'), hold on
for g = 1:numel(NDorder)
    n = find(NDlist == NDorder(g));
    for iz = 1:nZT
        xc   = 1 + (g-1)*gap + (iz - (nZT+1)/2) * 0.22;
        vals = {L(strcmp(cond_data, cond_names{iz,n}) & ~isnan(L)), ...   % left: data
                model_latency{iz,n}(post)'};                              % right: model
        for s = 1:2
            v = vals{s}(:);
            if numel(v) < 2 || max(v) == min(v), continue, end
            if s == 1, fc = datacol;     lc = [0 0 0];     fa = 0.5;
            else,      fc = ZTcol(iz,:); lc = ZTcol(iz,:); fa = 0.35;
            end
            q = quantile(v, [0.25 0.5 0.75]);
            if s == 1
                yi = linspace(q(1), q(3), 100);      % data: only the 25-75% range
            else
                yi = linspace(min(v), max(v), 100);  % model: full range
            end
            f = ksdensity(v, yi);
            f = f/max(f) * hw;
            fill([xc + sides(s)*f, xc*ones(1,100)], [yi, fliplr(yi)], fc, 'FaceAlpha', fa, 'EdgeColor', fc)
            plot([xc xc + sides(s)*interp1(yi, f, q(2))], [q(2) q(2)], '-', 'Color', lc, 'LineWidth', 2)
        end
    end
end
set(gca, 'XLim', xlims, 'XTick', 1 + (0:numel(NDorder)-1)*gap, ...
    'XTickLabel', strcat('ND', string(NDorder)), 'TickDir', 'out', 'Box', 'off')
ylim([0 inf])
ylabel('latency')

hz = gobjects(1, nZT);
for iz = 1:nZT, hz(iz) = fill(nan, nan, ZTcol(iz,:), 'FaceAlpha', 0.35, 'EdgeColor', ZTcol(iz,:)); end
hd = fill(nan, nan, datacol, 'FaceAlpha', 0.5, 'EdgeColor', datacol);
legend([hd hz], [{'data (25-75%)'}, cellstr(strcat('model ZT', string(ZTlist)))], 'Location', 'northwest')

%% ===================== measured data: latency and pre-escape speed per condition =====================
% same 4 conditions and order as the data vs model split violins
plot_conds  = main4;          % {'ND3LD14','ND0LD14','ND0LD02','ND3LD02'}
plot_labels = main4_labels;   % {'ZT14 dim','ZT14 bright','ZT2 bright','ZT2 dim'}
nC = numel(plot_conds);

meas      = {L, W};
meas_name = {'escape latency (s)', 'pre-escape speed (cm/s)'};
for m = 1:2
    V = meas{m};
    figure('Name', ['Measured ', meas_name{m}]), hold on
    for c = 1:nC
        inC = strcmp(cond_data, plot_conds{c}) & ~isnan(V);
        fm  = inC & sex_data == 1;
        ml  = inC & sex_data == 2;
        if any(fm), swarmchart(c*ones(nnz(fm),1), V(fm), 20, 'm', 'XJitterWidth', 0.6), end
        if any(ml), swarmchart(c*ones(nnz(ml),1), V(ml), 20, 'b', 'XJitterWidth', 0.6), end
        plot([c-0.3 c+0.3], median(V(fm))*[1 1],  '-', 'Color', 'm',         'LineWidth', 2)
        plot([c-0.3 c+0.3], median(V(ml))*[1 1],  '-', 'Color', 'b',         'LineWidth', 2)
        plot([c-0.3 c+0.3], median(V(inC))*[1 1], '-', 'Color', [0 0 0.3],   'LineWidth', 2)
    end
    hf = plot(nan, nan, 'om');
    hm = plot(nan, nan, 'ob');
    ha = plot(nan, nan, '-', 'Color', [0 0 0.3], 'LineWidth', 2);
    legend([hf hm ha], {'female (line = median)', 'male (line = median)', 'median all'}, 'Location', 'best')
    set(gca, 'XTick', 1:nC, 'XTickLabel', plot_labels, 'TickDir', 'out', 'TickLabelInterpreter', 'none')
    xtickangle(45)
    xlim([0.5 nC + 0.5])
    ylabel(['measured ', meas_name{m}])
    box off

    print_stats(meas_name{m}, V, cond_data, sex_data, plot_conds, plot_labels)
end

%% ===================== local functions =====================
function print_stats(name, V, cond_data, sex_data, conds, labels)
% medians, 25-75th percentiles and pairwise rank-sum tests for all, female and male mice
grp_name = {'all', 'female', 'male'};
grp_mask = {true(size(V)), sex_data == 1, sex_data == 2};

fprintf('\n==================== measured %s ====================\n', name);
fprintf('\nMedian [25th - 75th percentile] (n trials)\n');
fprintf('%-12s', 'condition');
for g = 1:3, fprintf('  %-34s', grp_name{g}); end
fprintf('\n');
for c = 1:numel(conds)
    fprintf('%-12s', labels{c});
    for g = 1:3
        v = V(strcmp(cond_data, conds{c}) & grp_mask{g} & ~isnan(V));
        if isempty(v)
            s = 'no data';
        else
            q = quantile(v, [0.25 0.5 0.75]);
            s = sprintf('%.3f [%.3f - %.3f] (n=%d)', q(2), q(1), q(3), numel(v));
        end
        fprintf('  %-34s', s);
    end
    fprintf('\n');
end

fprintf('\nWilcoxon rank-sum test, p-values (uncorrected for multiple comparisons)\n');
fprintf('%-12s    %-12s', 'condition A', 'condition B');
for g = 1:3, fprintf('  %-10s', grp_name{g}); end
fprintf('\n');
for a = 1:numel(conds)-1
    for b = a+1:numel(conds)
        fprintf('%-12s vs %-12s', labels{a}, labels{b});
        for g = 1:3
            va = V(strcmp(cond_data, conds{a}) & grp_mask{g} & ~isnan(V));
            vb = V(strcmp(cond_data, conds{b}) & grp_mask{g} & ~isnan(V));
            if isempty(va) || isempty(vb)
                p = NaN;
            else
                p = ranksum(va, vb);
            end
            fprintf('  %-10.4g', p);
        end
        fprintf('\n');
    end
end
end