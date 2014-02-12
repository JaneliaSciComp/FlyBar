function FlyBar(filePath)
    if nargin < 1 || isempty(filePath) || ~ischar(filePath)
        % Prompt the user to pick the input file.
        [fileName, parentPath] = uigetfile('*.*');
        filePath = fullfile(parentPath, fileName);
    end
    
    %% Parameters
    startSize = 20;
    tunnelSize = 100;
    endSize = 25;
    totalSize = startSize + tunnelSize + endSize;
    heatmapSize = [20 30];
    tunnelWidth = 6;    % max(fly_x) - min(fly_x)
    
    heatmapScale = [heatmapSize(1) / tunnelWidth, heatmapSize(2) / totalSize];
    
    xTickPos = [-startSize / 2,  tunnelSize / 2,  tunnelSize + endSize / 2];
    heatXTickPos = (xTickPos + startSize) * heatmapScale(2) + 0.5;
    
    %% Read in the data.
    fileID = fopen(filePath);
    try
        textscan(fileID, '%s %*[^\n]', 1); % skip over the header row
        while ~feof(fileID)
            fields = textscan(fileID, '%d\t%d\t%s\t', 1);
            run = fields{1};
            tunnel = fields{2};
            param = fields{3}{1};
            data(run, tunnel).(param) = sscanf(fgetl(fileID), '%f'); %#ok<*AGROW>
        end
    catch ME
        fclose(fileID);
        rethrow(ME);
    end
    fclose(fileID);
    
    runs = size(data, 1);
    tunnels = size(data, 2);
    
    %% Generate heatmaps and time spent in each section for each run/tunnel.
    heatmap = zeros(runs, tunnels, heatmapSize(1), heatmapSize(2));
    for r = 1:runs
        for t = 1:tunnels
            pauseTime = [diff(data(r, t).time); 1];
%            heatX = int16(round(data(r, t).fly_x * 2)) + 7;
            heatX = int16(round((data(r, t).fly_x + tunnelWidth / 2) * heatmapScale(1)));
%            heatY = int16(round(data(r, t).fly_y / 5)) + 5;
            heatY = int16(round((data(r, t).fly_y + startSize) * heatmapScale(2)));
            timeSpent(r, t).start = 0; %#ok<*AGROW>
            timeSpent(r, t).middle = 0;
            timeSpent(r, t).end = 0;
            speed(r, t).start = [];
            speed(r, t).middle = [];
            speed(r, t).end = [];
            for i = 1:length(pauseTime)
                heatmap(r, t, heatX(i), heatY(i)) = heatmap(r, t, heatX(i), heatY(i)) + pauseTime(i);
                
                if i > 1
                    dist = sqrt((data(r, t).fly_x(i) - data(r, t).fly_x(i - 1))^2 + (data(r, t).fly_y(i) - data(r, t).fly_y(i - 1))^2);
                end
                
                if data(r, t).fly_y(i) < 0
                    timeSpent(r, t).start = timeSpent(r, t).start + pauseTime(i);
                    if i > 1
                        speed(r, t).start(end + 1) = dist / pauseTime(i);
                    end
                elseif data(r, t).fly_y(i) < 100
                    timeSpent(r, t).middle = timeSpent(r, t).middle + pauseTime(i);
                    if i > 1
                        speed(r, t).middle(end + 1) = dist / pauseTime(i);
                    end
                else
                    timeSpent(r, t).end = timeSpent(r, t).end + pauseTime(i);
                    if i > 1
                        speed(r, t).end(end + 1) = dist / pauseTime(i);
                    end
                end
            end
        end
    end
    
    maxTime = max(cellfun(@(x) max(x), {data.time}));
    
    %% Create the figure
    figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    axesWidth = 1.0 / (tunnels + 1);
    axesHeight = 1.0 / (runs + 1);
    axesPad = 0.025;
    axesHandles = [];
    
    heatColormap = hot(256);
    
    % Allow the top and bottom percentile of non-zero values in the heatmap to saturate.
    heatSat = 0.025;
    heatSort = sort(heatmap(heatmap(:) > 0));
    heatMin = heatSort(uint16(heatSat * length(heatSort) + 0.5));
    heatMax = heatSort(uint16((1.0 - heatSat) * length(heatSort) + 0.5));
    
    %% Plot the individual tracks.
    colormap('Winter');
    for r = 1:runs
        for t = 1:tunnels
            axesHandles(r, t) = axes('Units', 'normalized', ...
                 'Position', [(t - 1) * axesWidth + axesPad, 1.0 - r * axesHeight + axesPad, axesWidth - axesPad * 2, axesHeight - axesPad * 2]); %#ok<LAXES>
            color_line(data(r, t).fly_y', data(r, t).fly_x', data(r, t).time' / maxTime);
            line([-startSize, tunnelSize + endSize, tunnelSize + endSize, -startSize, -startSize], ...
                 [-3 -3 3 3 -3], 'Color', 'k', 'LineWidth', 2);
            line([0 0], [-3 3], 'Color', 'k');
            line([tunnelSize tunnelSize], [-3 3], 'Color', 'k');
            set(gca, 'XLim', [-startSize, tunnelSize + endSize], ...
                     'XTick', xTickPos, ...
                     'XTickLabel', {sprintf('%.0f s\n%.1f mm/s', timeSpent(r, t).start, mean(speed(r, t).start)), ...
                                    sprintf('%.0f s\n%.1f mm/s', timeSpent(r, t).middle, mean(speed(r, t).middle)), ...
                                    sprintf('%.0f s\n%.1f mm/s', timeSpent(r, t).end, mean(speed(r, t).end))}, ...
                     'TickLength', [0 0], ...
                     'YLim', [-3 3], ...
                     'YTick', []);
            my_xticklabels(gca, xTickPos, get(gca, 'XTickLabel'));
            if r == 1
                title(sprintf('Tunnel %d', t));
            end
            if t == 1
                set(get(gca, 'YLabel'), 'String', sprintf('Run %d', r), ...
                                        'FontSize', get(get(gca, 'Title'), 'FontSize'));
            end
        end
    end
    
    %% Draw the average heat map for each run in the rightmost column.
    for r = 1:runs
        axesHandles(r, tunnels + 1) = axes('Units', 'normalized', ...
             'Position', [1.0 - axesWidth + axesPad, 1.0 - r * axesHeight + axesPad, axesWidth - axesPad * 2, axesHeight - axesPad * 2]); %#ok<LAXES>
        
        % Calculate the mean heatmap of the run then scale and saturate the result to fit the colormap.
        heatImg = flipud(squeeze(mean(heatmap(r, :, :, :), 2)));
        heatImg = uint16((heatImg - heatMin) / (heatMax - heatMin) * length(heatColormap) + 0.5);
        image(ind2rgb(heatImg, heatColormap));
        
        set(gca, 'XTick', heatXTickPos, ...
                 'XTickLabel', {sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(r, :).start]) / tunnels, mean([speed(r, :).start])), ...
                                sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(r, :).middle]) / tunnels, mean([speed(r, :).middle])), ...
                                sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(r, :).end]) / tunnels, mean([speed(r, :).end]))}, ...
                 'TickLength', [0 0], ...
                 'YTick', [], ...
                 'XAxisLocation', 'bottom');
        axis tight;
        set(gca, 'YDir', 'normal');
        my_xticklabels(gca, heatXTickPos, get(gca, 'XTickLabel'));
        set(gca, 'YDir', 'reverse');
        if r == 1
            title('Run Averages');
        end
    end
    
    %% Draw the average heat map for each tunnel in the bottom row.
    for t = 1:tunnels
        axesHandles(runs + 1, t) = axes('Units', 'normalized', ...
             'Position', [(t - 1) * axesWidth + axesPad, 0.0 + axesPad, axesWidth - axesPad * 2, axesHeight - axesPad * 2]); %#ok<LAXES>
        
        % Calculate the mean heatmap of the tunnel then scale and saturate the result to fit the colormap.
        heatImg = flipud(squeeze(mean(heatmap(:, t, :, :), 1)));
        heatImg = uint16((heatImg - heatMin) / (heatMax - heatMin)  * length(heatColormap) + 0.5);
        image(ind2rgb(heatImg, heatColormap));
        
        set(gca, 'XTick', heatXTickPos, ...
                 'XTickLabel', {sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(:, t).start]) / runs, mean([speed(:, t).start])), ...
                                sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(:, t).middle]) / runs, mean([speed(:, t).middle])), ...
                                sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(:, t).end]) / runs, mean([speed(:, t).end]))}, ...
                 'TickLength', [0 0], ...
                 'YTick', [], ...
                 'XAxisLocation', 'bottom');
        axis tight;
        set(gca, 'YDir', 'normal');
        my_xticklabels(gca, heatXTickPos, get(gca, 'XTickLabel'));
        set(gca, 'YDir', 'reverse');
        if t == 1
            set(get(gca, 'YLabel'), 'String', 'Tunnel Averages', ...
                                    'FontSize', get(get(gca, 'Title'), 'FontSize'));
        end
    end
    
    %% Draw the average heat map for all runs and all tunnels in the bottom-right plot.
    axesHandles(runs + 1, tunnels + 1) = axes('Units', 'normalized', ...
         'Position', [1.0 - axesWidth + axesPad, 0.0 + axesPad, axesWidth - axesPad * 2, axesHeight - axesPad * 2]);
    
    % Calculate the mean heatmap of all runs/tunnels then scale and saturate the result to fit the colormap.
    heatImg = flipud(squeeze(mean(mean(heatmap(:, :, :, :), 1), 2)));
    heatImg = uint16((heatImg - heatMin) / (heatMax - heatMin)  * length(heatColormap) + 0.5);
    image(ind2rgb(heatImg, heatColormap));

    set(gca, 'XTick', heatXTickPos, ...
             'XTickLabel', {sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(:, :).start]) / runs, mean([speed(:, :).start])), ...
                            sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(:, :).middle]) / runs, mean([speed(:, :).middle])), ...
                            sprintf('%.0f s\n%.1f mm/s', sum([timeSpent(:, :).end]) / runs, mean([speed(:, :).end]))}, ...
             'TickLength', [0 0], ...
             'YTick', []);
    axis tight;
    set(gca, 'YDir', 'normal');
    my_xticklabels(gca, heatXTickPos, get(gca, 'XTickLabel'));
    set(gca, 'YDir', 'reverse');
    
    %% Make the first axes the active one so that the color bar shows up there.
    axes(axesHandles(1, 1));
end
