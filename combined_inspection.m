function varargout = combined_inspection(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @combined_inspection_OpeningFcn, ...
                   'gui_OutputFcn',  @combined_inspection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function [handles] = initialize_data(hObject, eventdata, handles)
    
    saved_path = pwd;
    channel_no = pwd;
    disp(channel_no);
    handles.channel_no = strcat('channel',extractAfter(channel_no, string('channel')));
    handles.modification_tracker = [];
    handles.just_started = 1;
    
    handles.start_times = load('start_times.mat');
    handles.start_times = handles.start_times.start_indices;
    handles.hp_trace = NaN(2*handles.start_times{1,length(handles.start_times(1,:))},1);
    current = 1;
    
    for i = 1:length(handles.start_times(1,:))
        pieces = strsplit(saved_path, '/');
        pieces{1,length(pieces)-2} = handles.start_times{2,i};
        cd(strjoin(pieces, '/'));
        temp = rplhighpass('auto');
        handles.hp_trace(current:current-1+size(temp.data.analogData,1),1) = temp.data.analogData(:,1);
        current = current + size(temp.data.analogData,1);  
        
    end
    cd(saved_path);
    
    handles.hp_trace = handles.hp_trace(1:current-1,1);

    disp(length(handles.hp_trace));
    
%     load_highpass = rplhighpass('auto');
%     handles.hp_trace = load_highpass.data.analogData;
%     handles.time_trace = load_highpass.data.analogTime;
    
    channel_path = cd('./oSort/');
    folder_name = dir();
    cd(strcat(folder_name(length(folder_name)).name, '/sort/'));
    folder_name = dir();
    cd(folder_name(length(folder_name)).name);
    mat_name_d = dir('*_sorted_new.mat');
    mat_name = {mat_name_d.name};
    handles.spikes_data = load(mat_name{1});
    if size(handles.spikes_data.nrAssigned, 1) ~= 1
        handles.spikes_data.nrAssigned = flip(handles.spikes_data.nrAssigned);
    end
    cd(channel_path);

    disp(handles.spikes_data.nrAssigned);

    handles.cluster_names = handles.spikes_data.nrAssigned(:,1);
    disp('checkmark');

    handles.snap_nrAssigned = handles.spikes_data.nrAssigned;
    handles.noise_status = zeros(1,length(handles.spikes_data.nrAssigned(:,1)));
%     handles.noise_status = ones(1,length(handles.spikes_data.nrAssigned(:,1)));
    handles.noise_status2 = zeros(1,length(handles.spikes_data.nrAssigned(:,1)));    
    handles.noise_status3 = zeros(1,length(handles.spikes_data.nrAssigned(:,1))); 

    temp_arr = handles.spikes_data.nrAssigned(:,1);
    temp_arr = unique(temp_arr, 'stable');

    handles.distinct_plots = temp_arr;
    handles.distinct_plot_count = length(temp_arr);

    handles.autocorr_noise = load('./pngs/autocorr_noise.mat');
    handles.autocorr_noise = handles.autocorr_noise.tosave;
    handles.meandata = load('./pngs/meandata.mat');
    handles.meandata = handles.meandata.meandata;

    handles.counter_list = cell2mat(handles.meandata(2,:));
    handles.counter_list = flip(handles.counter_list);
    
    handles.sieved_means = cell2mat(handles.meandata(1,:)');
    handles.sieved_means = flip(handles.sieved_means);
    handles.sieved_means = handles.sieved_means';
%     size(handles.sieved_means)
    handles.min_min = handles.meandata{3,1};
    handles.max_max = handles.meandata{4,1};
    handles.contents_top = [0 0];
    
    handles.sieved_means_preserved = handles.sieved_means;
    handles.counter_list_preserved = handles.counter_list;
    
    % suggestions start
    
    for i = 1:length(handles.noise_status)
        average_data = handles.sieved_means(:,i);
        diff2 = abs(mean(average_data(1:5)) - mean(average_data(60:64)));
        if diff2 > 0.3*(max(average_data) - min(average_data))
            handles.noise_status(i) = 1;
        elseif diff2 > 3*handles.spikes_data.stdEstimateOrig
            handles.noise_status(i) = 1;
        end
    end

    disp('first round noise check:');
    disp(handles.noise_status);
    
    disp('autocorr round noise check:');
    handles.noise_status3 = flip(handles.autocorr_noise);
    disp(handles.noise_status3);
    
%     indices = handles.spikes_data.allSpikeInds;
% 
%     flagged = zeros(1, length(indices));
% 
%     for i = 1:length(indices)
%         central_peak = handles.hp_trace(indices(i));
%         left = max(1, indices(i) - 500);
%         right = min(length(handles.hp_trace), indices(i) + 500);
% 
%         window = handles.hp_trace(left:right);
%         if central_peak < 0
%             window = -window;
%         end
% 
%         thres = max(window);
% 
%         peaks = [];
%         warning('off','signal:findpeaks:largeMinPeakHeight');
%         while length(peaks) < 5
%             thres = thres * 0.95;
%             [peaks, pos] = findpeaks(window, 'MinPeakHeight', thres, 'MinPeakDistance', 20);
%         end
%         time_diff = diff(pos);
%         median_val = median(time_diff);
%         counter = 0;
%         for j = 1:length(time_diff)
%             if abs(median_val - time_diff(j)) < 10
%                 counter = counter + 1;
%             end
%         end
%         if counter/length(time_diff) > 0.75
%             for j = 1:length(indices)
%                 if abs(indices(j) - indices(i)) < 500
%                     flagged(j) = 1;
%                 end
%             end
%         end
%     end
%     sum(flagged)
%     length(flagged)
%     flagged_count = zeros(1,length(handles.counter_list));
%     for i = 1:length(flagged)
%         if flagged(i) == 1
%             for j = 1:length(handles.distinct_plots)
%                 if handles.distinct_plots(j) == handles.spikes_data.assignedNegative(i)
%                     flagged_count(j) = flagged_count(j) + 1;
%                 end
%             end
%         end
%     end
%     disp(flagged_count);
%     for i = 1:length(flagged_count)
%         ratio = flagged_count(i)/handles.spikes_data.nrAssigned(i,2);
%         disp(ratio);
%         if ratio >= 0.2
%             handles.noise_status2(i) = 1;
%         end
%     end
        
    disp('raw trace round noise check:');
    disp(handles.noise_status2);  
    handles.noise_status = handles.noise_status | handles.noise_status2 | handles.noise_status3;
    
    disp(handles.noise_status);
    
    % suggestions end
    
    
    plot(handles.bg1, 1, 1);
    cla(handles.bg1);
    set(handles.bg1,'visible','off');
    set(handles.bg1,'xtick',[]);
    plot(handles.bg2, 1, 1);
    cla(handles.bg2);
    set(handles.bg2,'visible','off');
    set(handles.bg2,'xtick',[]);    
    plot(handles.bg3, 1, 1);
    cla(handles.bg3);
    set(handles.bg3,'visible','off');
    set(handles.bg3,'xtick',[]);    
    
    handles.candidates = zeros(1,handles.distinct_plot_count);
    
    handles.png_stuff = cell(length(handles.counter_list), 3);

    for i = 1:length(handles.counter_list)
        file_name = char(strcat('./pngs/cluster',string(handles.spikes_data.nrAssigned(i,1)),'.png'));
        [handles.png_stuff{i,1}, handles.png_stuff{i,2}, handles.png_stuff{i,3}] = imread(file_name);
    end
   
    handles.total_pages = ceil(handles.distinct_plot_count/3);
    handles.current_page = 1;
    handles.current_display = 1;
    
%     handles.long1_plot = plot(handles.long1, handles.sieved_long(1,:,1), 'visible', 'off');
%     handles.long2_plot = plot(handles.long1, handles.sieved_long(2,:,1), 'visible', 'off');
    
    set(handles.long1,'visible','off');
    set(handles.long2,'visible','off');
    
    set(handles.long1_left, 'visible', 'off');
    set(handles.long1_right, 'visible', 'off');
    set(handles.long2_left, 'visible', 'off');
    set(handles.long2_right, 'visible', 'off');
    
    handles.long1_stat = 0;
    handles.long2_stat = 0;
    
    %%%%% initialize pca %%%%%
    
    handles.removing_noise_index = zeros(size(handles.spikes_data.assignedNegative));
    handles.noiseless_spikes = zeros(size(handles.spikes_data.newSpikesNegative));

    index = 1;

    for i = 1:length(handles.spikes_data.assignedNegative)
        if handles.spikes_data.assignedNegative(i) ~= 99999999
            handles.removing_noise_index(index) = handles.spikes_data.assignedNegative(i);
            handles.noiseless_spikes(index,:) = handles.spikes_data.newSpikesNegative(i,:);
            index = index + 1;
        end
    end

    handles.removing_noise_index = handles.removing_noise_index(1:index - 1);
    handles.noiseless_spikes = handles.noiseless_spikes(1:index - 1,:);

    disp('pca start');
    
    handles.peak_values = max(abs(handles.noiseless_spikes), [], 2);
    [~, handles.pca_score, ~] = pca(handles.noiseless_spikes);
    
%     [handles.L_R, handles.L, handles.IsolDist, handles.pca_score] = computeLratioALLclusters(handles.noiseless_spikes,handles.removing_noise_index,handles.spikes_data.nrAssigned(:,1));
%     disp(handles.L_R);
%     disp(handles.IsolDist);
    
    handles.pca_base_1 = scatter(handles.pca1, handles.pca_score(:,1), handles.pca_score(:,2));
    handles.pca_base_2 = scatter(handles.pca2, handles.pca_score(:,1), handles.pca_score(:,2));
    handles.pca_base_3 = scatter(handles.pca3, handles.pca_score(:,1), handles.pca_score(:,2));
    handles.pca_base_1.CData = [1 1 1];
    handles.pca_base_2.CData = [1 1 1];
    handles.pca_base_3.CData = [1 1 1];
    handles.pca_base_1.MarkerFaceAlpha = 0.2;
    handles.pca_base_2.MarkerFaceAlpha = 0.2;
    handles.pca_base_3.MarkerFaceAlpha = 0.2;
    handles.pca_overlay_data_1 = NaN(size(handles.pca_score,1),2);
    handles.pca_overlay_data_2 = NaN(size(handles.pca_score,1),2);
    handles.pca_overlay_data_3 = NaN(size(handles.pca_score,1),2);
    hold(handles.pca1, 'on');
        handles.pca_overlay_1 = scatter(handles.pca1, handles.pca_overlay_data_1(:,1), handles.pca_overlay_data_1(:,2), 'x', 'r');
    hold(handles.pca1, 'off');
    hold(handles.pca2, 'on');
        handles.pca_overlay_2 = scatter(handles.pca2, handles.pca_overlay_data_2(:,1), handles.pca_overlay_data_2(:,2), 'x', 'r');
    hold(handles.pca2, 'off');
    hold(handles.pca3, 'on');
        handles.pca_overlay_3 = scatter(handles.pca3, handles.pca_overlay_data_3(:,1), handles.pca_overlay_data_3(:,2), 'x', 'r');
    hold(handles.pca3, 'off');    
    
    handles.pca_base_1t = scatter(handles.pca1t, handles.pca_score(:,1), handles.peak_values);
    handles.pca_base_2t = scatter(handles.pca2t, handles.pca_score(:,1), handles.peak_values);
    handles.pca_base_3t = scatter(handles.pca3t, handles.pca_score(:,1), handles.peak_values);
    handles.pca_base_1t.CData = [1 1 1];
    handles.pca_base_2t.CData = [1 1 1];
    handles.pca_base_3t.CData = [1 1 1];
    handles.pca_base_1t.MarkerFaceAlpha = 0.2;
    handles.pca_base_2t.MarkerFaceAlpha = 0.2;
    handles.pca_base_3t.MarkerFaceAlpha = 0.2;
    handles.pca_overlay_data_1t = NaN(size(handles.pca_score,1),2);
    handles.pca_overlay_data_2t = NaN(size(handles.pca_score,1),2);
    handles.pca_overlay_data_3t = NaN(size(handles.pca_score,1),2);
    hold(handles.pca1t, 'on');
        handles.pca_overlay_1t = scatter(handles.pca1t, handles.pca_overlay_data_1t(:,1), handles.pca_overlay_data_1t(:,2), 'x', 'r');
    hold(handles.pca1t, 'off');
    hold(handles.pca2t, 'on');
        handles.pca_overlay_2t = scatter(handles.pca2t, handles.pca_overlay_data_2t(:,1), handles.pca_overlay_data_2t(:,2), 'x', 'r');
    hold(handles.pca2t, 'off');
    hold(handles.pca3t, 'on');
        handles.pca_overlay_3t = scatter(handles.pca3t, handles.pca_overlay_data_3t(:,1), handles.pca_overlay_data_3t(:,2), 'x', 'r');
    hold(handles.pca3t, 'off');     
    
    
    disp('pca end');

    handles.small_arr = [handles.small1, handles.small2, handles.small3, handles.small4, handles.small5, handles.small6, handles.small7, handles.small8, handles.small9, handles.small10, handles.small11, handles.small12, handles.small13, handles.small14, handles.small15];
    handles.current_highlight_details = 0;
    handles.color_purple = [0.611, 0.086, 0.584];
    handles.color_blue = [0.3 0.75 0.93];
    handles.color_beige = [0.901, 0.658, 0.733];
    handles.color_dark = [0.188, 0.152, 0.349];
    
    handles.data_bot = 0;
    
guidata(hObject, handles);



function [handles] = page_plots(hObject, eventdata, handles)

    
    for i = 1:length(handles.small_arr)
        cla(handles.small_arr(i));
        set(handles.small_arr(i), 'visible', 'off');
    end
    
    count = 1;
    
    start = 1 + 15*(handles.current_page - 1);

    for i = start:min(start+14, handles.distinct_plot_count)
        
        set(handles.small_arr(count), 'visible', 'on');
        
        plot(handles.small_arr(count), [1 64], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2,'Hittest','off');
        hold(handles.small_arr(count), 'on');
        plot(handles.small_arr(count), [1 64], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2,'Hittest','off');
        plot(handles.small_arr(count), 1:64, handles.sieved_means(:,i), 'r', 'LineWidth', 5,'Hittest','off');
        
        
        [children, indices, children_length] = find_children(handles.distinct_plots(i), hObject, eventdata, handles);

        for j = 1:length(indices)
            plot(handles.small_arr(count), 1:64, handles.sieved_means_preserved(:,indices(j)), 'LineWidth', 1,'Hittest','off');
        end
        
        title(handles.small_arr(count), strjoin(cellstr(num2str(children)), ' '));
        
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == handles.distinct_plots(i)
                if handles.noise_status(j) == 1
                    title(handles.small_arr(count), strcat(strjoin(cellstr(num2str(children)), ' '), ' (noise)'));
                end
            end
        end
        hold(handles.small_arr(count), 'off');

        
        ylim(handles.small_arr(count), [handles.min_min, handles.max_max]);
        xlim(handles.small_arr(count), [0 64]);
        set(handles.small_arr(count), 'ButtonDownFcn', {@small_click, count});
        
        count = count + 1;

    end
    
guidata(hObject, handles);


function [handles] = update_list(target_index, hObject, eventdata, handles)

    set(handles.list1, 'Value', 1);
    set(handles.list1,'String',[]);


        [children, ~, count] = find_children(handles.distinct_plots(target_index), hObject, eventdata, handles);

        if length(children) > 1
            set(handles.noise1, 'visible', 'off');
        else
            set(handles.noise1, 'visible', 'on');
            for i = 1:length(handles.snap_nrAssigned(:,1))
                if handles.snap_nrAssigned(i,1) == children(1)
                    set(handles.noise1,'Value',handles.noise_status(i));
                    break;
                end
            end
        end
        
        children = mat2cell(children, 1, count);
        children = [{'View All'},children];
        set(handles.list1,'String', children);


guidata(hObject, handles);


function [children, indices, length] = find_children(target, hObject, eventdata, handles)

    children = zeros(1, size(handles.spikes_data.nrAssigned,1));
    indices = zeros(1, size(handles.spikes_data.nrAssigned,1));
    index = 1;
    for i = 1:size(handles.spikes_data.nrAssigned,1)
        if handles.spikes_data.nrAssigned(i,1) == target
            indices(index) = i;
            children(index) = handles.snap_nrAssigned(i,1);
            index = index + 1;

        end
    end
    
    indices = indices(1:index - 1);
    children = children(1:index - 1);
    length = index - 1;

function small_click(hObject, eventdata, count)

    handles = guidata(hObject);
    index = count + 15*(handles.current_page - 1);
    
    figHandle = ancestor(hObject, 'figure');
    clickType = get(figHandle, 'SelectionType');
    
    if strcmp(clickType, 'alt')
        disp('right click');
        if handles.noise_status(index) == 0
            if handles.candidates(index) == 1
                handles.candidates(index) = 0;
            else
                handles.candidates(index) = 1;    
            end

            cla(handles.bg3);
            cla(handles.png6); 
            set(handles.pca_base_3, 'visible', 'off');
            set(handles.pca_overlay_3, 'visible', 'off');
            cla(handles.ac3, 'reset');
            cla(handles.sr3, 'reset');        

            if sum(handles.candidates) > 0
                handles = merge_plot1(hObject, eventdata, handles);
                handles = merge_plot2(hObject, eventdata, handles);
                handles = merge_plot34(hObject, eventdata, handles);
            end
        end
    else
        handles.current_highlight_details = index;
        handles = update_list(index, hObject, eventdata, handles);
    end
    
    handles = update_colors(hObject, eventdata, handles);
    
guidata(hObject, handles);
    
function [handles] = update_colors(hObject, eventdata, handles)

count = 1;

for i = 1 + 15*(handles.current_page - 1):min(15 + 15*(handles.current_page - 1), handles.distinct_plot_count)
    
    if handles.candidates(i) == 1
        if handles.current_highlight_details ~= i
            set(handles.small_arr(count), 'Color', handles.color_blue);
            
        else
            set(handles.small_arr(count), 'Color', handles.color_purple);
        end
    else
        if handles.current_highlight_details == i
            set(handles.small_arr(count), 'Color', handles.color_beige);
        else
            set(handles.small_arr(count), 'Color', [1 1 1]);
        end
    end
    count = count + 1;
end

guidata(hObject, handles);


function [handles] = merge_plot1(hObject, eventdata, handles)

        indices = zeros(1, length(handles.spikes_data.nrAssigned(:,1)));
        count = 1;
        for i = 1:handles.distinct_plot_count
            if handles.candidates(i) == 1
                for j = 1:length(handles.spikes_data.nrAssigned)
                    if handles.spikes_data.nrAssigned(j) == handles.distinct_plots(i)
                        indices(count) = j;
                        count = count + 1;
                    end
                end
            end
        end
        indices = indices(1:count-1);
        
            
        handles.box_m_handles = cell(1,length(indices));
        
        for j = 1:length(indices)
            if j == 1
                handles.box_m_handles{1,j} = image(handles.png6, handles.png_stuff{indices(j),1}, 'AlphaData', handles.png_stuff{indices(j),3});
                set(handles.box_m_handles{1,j}, 'AlphaData', handles.png_stuff{indices(j),3});
                hold(handles.png6, 'on');
            else
                handles.box_m_handles{1,j} = image(handles.png6, handles.png_stuff{indices(j),1}, 'AlphaData', handles.png_stuff{indices(j),3});
                set(handles.box_m_handles{1,j}, 'AlphaData', handles.png_stuff{indices(j),3});
            end
        end
           
        hold(handles.png6, 'off');

        set(handles.png6,'visible','off');
        set(handles.png6, 'XLimSpec', 'Tight');
        set(handles.png6, 'YLimSpec', 'Tight');

        count = 1;
        indices = zeros(1, handles.distinct_plot_count);
        for i = 1:handles.distinct_plot_count
            if handles.candidates(i) == 1
                indices(count) = i;
                count = count + 1;
            end
        end
        indices = indices(1:count-1);
        
        
        for j = 1:length(indices)
            if j == 1
                plot(handles.bg3, 1:64, handles.sieved_means(:,indices(j)), 'r', 'LineWidth', 2);
                hold(handles.bg3, 'on');
            else
                plot(handles.bg3, 1:64, handles.sieved_means(:,indices(j)), 'r', 'LineWidth', 2);
            end
        end
        plot(handles.bg3, [1 64], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.bg3, [1 64], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        hold(handles.bg3, 'off');
                
        ylim(handles.bg3, [handles.min_min, handles.max_max]);
        xlim(handles.bg3, [0 64]);
        set(handles.bg3, 'Color', 'None');

guidata(hObject, handles);

function [handles] = merge_plot2(hObject, eventdata, handles)

    
        contents = zeros(1, length(handles.spikes_data.nrAssigned(:,1)));
        count = 1;
        for i = 1:handles.distinct_plot_count
            if handles.candidates(i) == 1
                for j = 1:size(handles.spikes_data.nrAssigned,1)
                    if handles.spikes_data.nrAssigned(j,1) == handles.distinct_plots(i)
                        contents(count) = handles.snap_nrAssigned(j,1);
                        count = count + 1;
                    end
                end
            end
        end
        contents = contents(1:count-1);


    [out, outt, IsolDist, L_R] = pca_siever(handles, contents, handles.pca_score);

    set(handles.pca_base_3, 'visible', 'on');
    set(handles.pca_overlay_3, 'visible', 'on');
    handles.pca_base_3.CData = [0 1 0];

    handles.pca_overlay_3.XData = out(:,1);
    handles.pca_overlay_3.YData = out(:,2);
    
    set(handles.pca_base_3t, 'visible', 'on');
    set(handles.pca_overlay_3t, 'visible', 'on');
    handles.pca_base_3t.CData = [0 1 0];

    handles.pca_overlay_3t.XData = outt(:,1);
    handles.pca_overlay_3t.YData = outt(:,2);  
    
    title(handles.pca3t, strcat('IsolDist: ', num2str(IsolDist), ' L_R: ', num2str(L_R)));


guidata(hObject, handles);
    
function [handles] = merge_plot34(hObject, eventdata, handles)

        set(handles.ac3, 'visible', 'on');
        set(handles.sr3, 'visible', 'on');

        contents = zeros(1, length(handles.spikes_data.nrAssigned(:,1)));
        count = 1;
        for i = 1:handles.distinct_plot_count
            if handles.candidates(i) == 1
                for j = 1:size(handles.spikes_data.nrAssigned,1)
                    if handles.spikes_data.nrAssigned(j,1) == handles.distinct_plots(i)
                        contents(count) = handles.snap_nrAssigned(j,1);
                        count = count + 1;
                    end
                end
            end
        end
        contents = contents(1:count-1);    

                %%%% start autocorr/spikerates part %%%%

        
        curr_spikes = NaN(length(handles.spikes_data.assignedNegative), 64);
        curr_times = NaN(1,length(handles.spikes_data.assignedNegative));
        curr_times_spike_train = NaN(1,length(handles.spikes_data.assignedNegative));
        count = 1;
        for i = 1:length(handles.spikes_data.assignedNegative)
            for j = 1:length(contents)
                if handles.spikes_data.assignedNegative(i) == contents(j)
                    curr_spikes(count,:) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
                    curr_times(1,count) = handles.spikes_data.allSpikeInds(i);
                    curr_times_spike_train(1,count) = handles.spikes_data.newTimestampsNegative(i);
                    count = count + 1;
                    break;
                end
            end
        end
        curr_spikes = curr_spikes(1:count-1,:);
        curr_times = curr_times(1,1:count-1);
        curr_times_spike_train = curr_times_spike_train(1,1:count-1);
       
         
        spike_train = convertToSpiketrain(curr_times_spike_train, 1);
        
        [~,~,tvect,Cxx] = psautospk(spike_train, 1);
%         [tvect, Cxx, ~] = autocorr_2(spike_train);

        time_tracker = 1;
        counter = 0;
        amp_summer = 0;
        index = 1;
        total_length = length(handles.hp_trace);
        total_curr_spikes = length(curr_times);
        amp_data = zeros(floor(total_length/(30000*100)),1);
        sr_data = zeros(floor(total_length/(30000*100)),1);
        data_index = 1;

        while time_tracker < total_length
            if index > total_curr_spikes
                sr_data(data_index) = 0;
                amp_data(data_index) = 0;
                data_index = data_index + 1;
                time_tracker = time_tracker + 30000*100;
            elseif curr_times(index) < time_tracker + 30000*100
                amp_summer = amp_summer + curr_spikes(index,24);
                counter = counter + 1;
            else
                sr_data(data_index) = counter/100;
                if counter ~= 0
                    amp_data(data_index) = double(amp_summer)/counter;
                else
                    amp_data(data_index) = 0;
                end
                time_tracker = time_tracker + 30000*100;
                data_index = data_index + 1;
                amp_summer = 0;
                counter = 0;
            end 
            index = index + 1;
        end

        time_axis = linspace(0,total_length/30000,length(sr_data)); 
         
        axes(handles.sr3);
        cla reset;
        yyaxis left
        plot(time_axis, sr_data);
        ylim([min(sr_data)-0.005, max(sr_data)+0.005]);
        hold on
        yyaxis right
        plot(time_axis, amp_data, '--');
        ylim([min(amp_data)-5, max(amp_data)+5]);
        yyaxis left
        hold off
        xlim([0,(time_axis(end)+5)]);
        
        Cxx(1) = [];

        tvect = tvect(floor(length(tvect)/2):length(tvect));
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        tvect = tvect(2:100);
        Cxx = Cxx(2:100);
        
        plot(handles.ac3t, tvect(1:10), Cxx(1:10).');
        plot(handles.ac3, tvect, Cxx.');
        
        %%%%% end autocorr/spikerates part %%%%%

guidata(hObject, handles);

% --- Executes just before combined_inspection is made visible.
function combined_inspection_OpeningFcn(hObject, eventdata, handles, varargin)

    handles = initialize_data(hObject, eventdata, handles);
    disp('starting automerge');
    disp(handles.spikes_data.nrAssigned);
    handles = auto_merge(hObject, eventdata, handles);
    disp(handles.spikes_data.nrAssigned);
    handles = page_plots(hObject, eventdata, handles);


    % Choose default command line output for combined_inspection
    handles.output = hObject;

guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = combined_inspection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function next_Callback(hObject, eventdata, handles)

if handles.current_page < handles.total_pages
    handles.current_page = handles.current_page + 1;
    handles = page_plots(hObject, eventdata, handles);
    handles = update_colors(hObject, eventdata, handles);
    
end

guidata(hObject, handles);


function back_Callback(hObject, eventdata, handles)

if handles.current_page > 1
    handles.current_page = handles.current_page - 1;
    handles = page_plots(hObject, eventdata, handles);
    handles = update_colors(hObject, eventdata, handles);
    
end

guidata(hObject, handles);



function list1_Callback(hObject, eventdata, handles)

function list1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [handles] = details1_bot(hObject, eventdata, handles)

    
    contents = handles.contents_bot;
      
    index_for_png = zeros(1,length(contents));
    for i = 1:length(contents)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                index_for_png(i) = j;
            end
        end
    end
    
    axes_handles1 = [handles.png4, handles.png5];
    axes_handles2 = [handles.back4, handles.back5];
    
    pca_invis = [handles.pca_base_1, handles.pca_base_2];
    pca_overlay = [handles.pca_overlay_1, handles.pca_overlay_2];
    
    for j = 1:size(index_for_png, 2)
        if j == 1
            temp_handle = image(axes_handles1(2), handles.png_stuff{index_for_png(1,j),1}, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
            set(temp_handle, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
            hold(axes_handles1(2), 'on');
        else
            temp_handle = image(axes_handles1(2), handles.png_stuff{index_for_png(1,j),1}, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
            set(temp_handle, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
        end
    end
    
    hold(axes_handles1(2), 'off');
    
        set(axes_handles1(2),'visible','off');
        set(axes_handles1(2), 'XLimSpec', 'Tight');
        set(axes_handles1(2), 'YLimSpec', 'Tight');
    
        for j = 1:size(index_for_png, 2)
            if j == 1
                plot(axes_handles2(2), 1:64, handles.sieved_means_preserved(:,index_for_png(1,j)), 'r', 'LineWidth', 2);
                hold(axes_handles2(2), 'on');
            else
                plot(axes_handles2(2), 1:64, handles.sieved_means_preserved(:,index_for_png(1,j)), 'r', 'LineWidth', 2);
            end            
        end
        
        plot(axes_handles2(2), [1 64], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(axes_handles2(2), [1 64], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        hold(axes_handles2(2), 'off');
        
        ylim(axes_handles2(2), [handles.min_min, handles.max_max]);
        xlim(axes_handles2(2), [0 64]);
        set(axes_handles2(2), 'Color', 'None');   
        
            
        %%%%%%% pca part %%%%%%%%
        
        [out, outt, IsolDist, L_R] = pca_siever(handles, contents, handles.pca_score);
        
        title(handles.pca2t, strcat('IsolDist: ', num2str(IsolDist), ' L_R: ', num2str(L_R))); %strcat('IsolDist: ', num2str(IsolDist), ' L_R: ', num2str(L_R))
        pca_invis(2).CData = [0 1 0];
                
        pca_overlay(2).XData = out(:,1);
        pca_overlay(2).YData = out(:,2);
        
        handles.pca_base_2t.CData = [0 1 0];

        handles.pca_overlay_2t.XData = outt(:,1);
        handles.pca_overlay_2t.YData = outt(:,2);        
        
        %%%%% end pca part %%%%%%
        
        %%%% start autocorr/spikerates part %%%%
        
        axes_sr = [handles.sr1, handles.sr2];
        axes_ac = [handles.ac1, handles.ac2];
         
        
        curr_spikes = NaN(length(handles.spikes_data.assignedNegative), 64);
        curr_times = NaN(1,length(handles.spikes_data.assignedNegative));
        curr_times_spike_train = NaN(1,length(handles.spikes_data.assignedNegative));
        count = 1;
        for i = 1:length(handles.spikes_data.assignedNegative)
            for j = 1:length(contents)
                if handles.spikes_data.assignedNegative(i) == contents(j)
                    curr_spikes(count,:) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
                    curr_times(1,count) = handles.spikes_data.allSpikeInds(i);
                    curr_times_spike_train(1,count) = handles.spikes_data.newTimestampsNegative(i);
                    count = count + 1;
                    break;
                end
            end
        end
        curr_spikes = curr_spikes(1:count-1,:);
        curr_times = curr_times(1,1:count-1);
        curr_times_spike_train = curr_times_spike_train(1,1:count-1);
       
         
        spike_train = convertToSpiketrain(curr_times_spike_train, 1);
        
        [~,~,tvect,Cxx] = psautospk(spike_train, 1);
%         [tvect, Cxx, ~] = autocorr_2(spike_train);

        time_tracker = 1;
        counter = 0;
        amp_summer = 0;
        index = 1;
        total_length = length(handles.hp_trace);
        total_curr_spikes = length(curr_times);
        amp_data = zeros(floor(total_length/(30000*100)),1);
        sr_data = zeros(floor(total_length/(30000*100)),1);
        data_index = 1;

        while time_tracker < total_length
            if index > total_curr_spikes
                sr_data(data_index) = 0;
                amp_data(data_index) = 0;
                data_index = data_index + 1;
                time_tracker = time_tracker + 30000*100;
            elseif curr_times(index) < time_tracker + 30000*100
                amp_summer = amp_summer + curr_spikes(index,24);
                counter = counter + 1;
            else
                sr_data(data_index) = counter/100;
                if counter ~= 0
                    amp_data(data_index) = double(amp_summer)/counter;
                else
                    amp_data(data_index) = 0;
                end
                time_tracker = time_tracker + 30000*100;
                data_index = data_index + 1;
                amp_summer = 0;
                counter = 0;
            end 
            index = index + 1;
        end

        time_axis = linspace(0,total_length/1800000,length(sr_data)); 
         
        axes(axes_sr(2));
        cla reset;
        yyaxis left
        plot(time_axis, sr_data);
        ylim([min(sr_data)-0.005, max(sr_data)+0.005]);
        hold on
        yyaxis right
        plot(time_axis, amp_data, '--');
        ylim([min(amp_data)-5, max(amp_data)+5]);
        yyaxis left
        hold off
        xlim([0,(time_axis(end)+0.05)]);
        
        Cxx(1) = [];

        tvect = tvect(floor(length(tvect)/2):length(tvect));
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        tvect = tvect(2:100);
        Cxx = Cxx(2:100);
        
        plot(axes_ac(2), tvect, Cxx.');
        plot(handles.ac2t, tvect(1:10), Cxx(1:10).');
        
        %%%%% end autocorr/spikerates part %%%%%
           
    long_arr = [handles.long1, handles.long2];
    long_arr_l = [handles.long1_left, handles.long2_left];
    long_arr_r = [handles.long1_right, handles.long2_right];

    cla(long_arr(2));
    set(long_arr(2),'visible','off');

    set(long_arr_l(2), 'visible', 'off');
    set(long_arr_r(2), 'visible', 'off');

    handles.long2_stat = 0;

    
    if length(contents) == 1

            set(handles.uipanel2, 'Title', num2str(contents(1)));

    else

            set(handles.uipanel2, 'Title', strcat(num2str(contents(1)), ' - multiple'));

    end
    
    [ts, handles.updated_index_bot] = search_ts_forward(handles, 0, contents);
    handles.data_bot = ts;
    
guidata(hObject, handles);




function details1_Callback(hObject, eventdata, handles)

    
    contents = cellstr(get(handles.list1,'String'));
    contents = str2double(contents{get(handles.list1,'Value')});
    
    if isnan(contents) == 1
        contents = cellstr(get(handles.list1,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
    
    handles.contents_bot = handles.contents_top;
    handles.contents_top = contents;
    
         
    index_for_png = zeros(1,length(contents));
    for i = 1:length(contents)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                index_for_png(i) = j;
            end
        end
    end
    
    axes_handles1 = [handles.png4, handles.png5];
    axes_handles2 = [handles.back4, handles.back5];
    
    pca_invis = [handles.pca_base_1, handles.pca_base_2];
    pca_overlay = [handles.pca_overlay_1, handles.pca_overlay_2];
    
    for j = 1:size(index_for_png, 2)
        if j == 1
            temp_handle = image(axes_handles1(handles.current_display), handles.png_stuff{index_for_png(1,j),1}, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
            set(temp_handle, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
            hold(axes_handles1(handles.current_display), 'on');
        else
            temp_handle = image(axes_handles1(handles.current_display), handles.png_stuff{index_for_png(1,j),1}, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
            set(temp_handle, 'AlphaData', handles.png_stuff{index_for_png(1,j),3});
        end
    end
    
    hold(axes_handles1(handles.current_display), 'off');
    
        set(axes_handles1(handles.current_display),'visible','off');
        set(axes_handles1(handles.current_display), 'XLimSpec', 'Tight');
        set(axes_handles1(handles.current_display), 'YLimSpec', 'Tight');
    
        for j = 1:size(index_for_png, 2)
            if j == 1
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(:,index_for_png(1,j)), 'r', 'LineWidth', 2);
                hold(axes_handles2(handles.current_display), 'on');
            else
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(:,index_for_png(1,j)), 'r', 'LineWidth', 2);
            end            
        end
        
        plot(axes_handles2(handles.current_display), [1 64], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(axes_handles2(handles.current_display), [1 64], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        hold(axes_handles2(handles.current_display), 'off');
        
        ylim(axes_handles2(handles.current_display), [handles.min_min, handles.max_max]);
        xlim(axes_handles2(handles.current_display), [0 64]);
        set(axes_handles2(handles.current_display), 'Color', 'None');   
        
            
        %%%%%%% pca part %%%%%%%%
        
        [out, outt, IsolDist, L_R] = pca_siever(handles, contents, handles.pca_score);
        title(handles.pca1t, strcat('IsolDist: ', num2str(IsolDist), ' L_R: ', num2str(L_R)));
        
        pca_invis(handles.current_display).CData = [0 1 0];

        handles.pca_overlay_1.XData = out(:,1);
        handles.pca_overlay_1.YData = out(:,2);
        
        handles.pca_base_1t.CData = [0 1 0];
        
        handles.pca_overlay_1t.XData = outt(:,1);
        handles.pca_overlay_1t.YData = outt(:,2);
        
        %%%%% end pca part %%%%%%
        
        %%%% start autocorr/spikerates part %%%%
        
        axes_sr = [handles.sr1, handles.sr2];
        axes_ac = [handles.ac1, handles.ac2];
         
        
        curr_spikes = NaN(length(handles.spikes_data.assignedNegative), 64);
        curr_times = NaN(1,length(handles.spikes_data.assignedNegative));
        curr_times_spike_train = NaN(1,length(handles.spikes_data.assignedNegative));
        count = 1;
        for i = 1:length(handles.spikes_data.assignedNegative)
            for j = 1:length(contents)
                if handles.spikes_data.assignedNegative(i) == contents(j)
                    curr_spikes(count,:) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
                    curr_times(1,count) = handles.spikes_data.allSpikeInds(i);
                    curr_times_spike_train(1,count) = handles.spikes_data.newTimestampsNegative(i);
                    count = count + 1;
                    break;
                end
            end
        end
        curr_spikes = curr_spikes(1:count-1,:);
        curr_times = curr_times(1,1:count-1);
        curr_times_spike_train = curr_times_spike_train(1,1:count-1);
       
         
        spike_train = convertToSpiketrain(curr_times_spike_train, 1);
        
        [~,~,tvect,Cxx] = psautospk(spike_train, 1);
%         [tvect, Cxx, ~] = autocorr_2(spike_train); 

        time_tracker = 1;
        counter = 0;
        amp_summer = 0;
        index = 1;
        total_length = length(handles.hp_trace);
        total_curr_spikes = length(curr_times);
        amp_data = zeros(floor(total_length/(30000*100)),1);
        sr_data = zeros(floor(total_length/(30000*100)),1);
        data_index = 1;

        while time_tracker < total_length
            if index > total_curr_spikes
                sr_data(data_index) = 0;
                amp_data(data_index) = 0;
                data_index = data_index + 1;
                time_tracker = time_tracker + 30000*100;
            elseif curr_times(index) < time_tracker + 30000*100
                amp_summer = amp_summer + curr_spikes(index,24);
                counter = counter + 1;
            else
                sr_data(data_index) = counter/100;
                if counter ~= 0
                    amp_data(data_index) = double(amp_summer)/counter;
                else
                    amp_data(data_index) = 0;
                end
                time_tracker = time_tracker + 30000*100;
                data_index = data_index + 1;
                amp_summer = 0;
                counter = 0;
            end 
            index = index + 1;
        end

        time_axis = linspace(0,total_length/1800000,length(sr_data)); 
         
        axes(axes_sr(handles.current_display));
        cla reset;
        yyaxis left
        plot(time_axis, sr_data);
        ylim([min(sr_data)-0.005, max(sr_data)+0.005]);
        hold on
        yyaxis right
        plot(time_axis, amp_data, '--');
        ylim([min(amp_data)-5, max(amp_data)+5]);
        yyaxis left
        hold off
        xlim([0,(time_axis(end)+0.05)]);
        
        Cxx(1) = [];

        tvect = tvect(floor(length(tvect)/2):length(tvect));
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        tvect = tvect(2:100);
        Cxx = Cxx(2:100);
        
        plot(axes_ac(handles.current_display), tvect, Cxx.');
        plot(handles.ac1t, tvect(1:10), Cxx(1:10).');
        
        %%%%% end autocorr/spikerates part %%%%%
           
    long_arr = [handles.long1, handles.long2];
    long_arr_l = [handles.long1_left, handles.long2_left];
    long_arr_r = [handles.long1_right, handles.long2_right];

    cla(long_arr(handles.current_display));
    set(long_arr(handles.current_display),'visible','off');

    set(long_arr_l(handles.current_display), 'visible', 'off');
    set(long_arr_r(handles.current_display), 'visible', 'off');

    if handles.current_display == 1
        handles.long1_stat = 0;
    else
        handles.long2_stat = 0;
    end
    
    if length(contents) == 1
        if handles.current_display == 1
            set(handles.uipanel1, 'Title', num2str(contents(1)));
        else
            set(handles.uipanel2, 'Title', num2str(contents(1)));
        end
    else
        if handles.current_display == 1
            set(handles.uipanel1, 'Title', strcat(num2str(contents(1)), ' - multiple'));
        else
            set(handles.uipanel2, 'Title', strcat(num2str(contents(1)), ' - multiple'));
        end
    end
    
    if handles.just_started == 0
        handles.updated_index_bot = handles.updated_index_top;
    end
    
    [ts, handles.updated_index_top] = search_ts_forward(handles, 0, contents);
    
    if handles.just_started == 0
        handles.data_bot = handles.data_top;
    end
    
    handles.data_top = ts;
    
    if handles.just_started == 0
        handles = details1_bot(hObject, eventdata, handles);
    end
    
    handles.just_started = 0;
    
guidata(hObject, handles);




function [ts, updated_index] = search_ts_forward(handles, current_index, targets)

    start = 1;
    for i = current_index:length(handles.spikes_data.assignedNegative)
        if start ~= 1
            for j = 1:length(targets)
                if handles.spikes_data.assignedNegative(i) == targets(j)
                    updated_index = i;
                    ts = handles.hp_trace(handles.spikes_data.allSpikeInds(i)-499:handles.spikes_data.allSpikeInds(i)+500)';
                    return
                end
            end
        end
        start = 0;
    end

function [ts, updated_index] = search_ts_backward(handles, current_index, targets)

    start = 1;
    for i = current_index:-1:1
        if start ~= 1
            for j = 1:length(targets)
                if handles.spikes_data.assignedNegative(i) == targets(j)
                    updated_index = i;
                    ts = handles.hp_trace(handles.spikes_data.allSpikeInds(i)-499:handles.spikes_data.allSpikeInds(i)+500)';
                    return
                end
            end
        end
        start = 0;
    end        
      
    
    


function [sieved_pca, sieved_pcat, IsolDist, L_R] = pca_siever(handles, target_arr, full_pca)

sieved_pca = NaN(size(full_pca,1),2);    
sieved_pcat = NaN(size(full_pca,1),2);
pca3 = full_pca(:,1:3);
idx = NaN(1,size(full_pca,1));

counter = 1;
for i = 1:length(handles.removing_noise_index)
    for j = 1:length(target_arr)
        if handles.removing_noise_index(i) == target_arr(j)            
            sieved_pca(counter,:) = full_pca(i,1:2);
            sieved_pcat(counter,1) = full_pca(i,1);
            sieved_pcat(counter,2) = handles.peak_values(i);
            idx(counter) = i;
            counter = counter + 1;
            break;
        end
    end
end
idx = idx(1:counter-1);

IsolDist = IsolationDistance(pca3, idx);
[L,L_R]=L_Ratio(pca3,idx); %Mclust function 


function kick1_Callback(hObject, eventdata, handles)

    disp('kicking');
    contents = cellstr(get(handles.list1,'String'));
    contents = str2double(contents{get(handles.list1,'Value')});

    if isnan(contents) == 1
        contents = cellstr(get(handles.list1,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
    
    disp(contents);
    disp(handles.spikes_data.nrAssigned);
    
    found = 0;
    regroup = 0;
    
    for i = 1:length(contents)
        for j = 1:length(handles.spikes_data.nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.spikes_data.nrAssigned(j,1) = handles.snap_nrAssigned(j,1);
                found = 1;
            elseif found == 1
                if handles.spikes_data.nrAssigned(j,1) == contents(i)
                    if regroup == 0
                        handles.spikes_data.nrAssigned(j,1) = handles.snap_nrAssigned(j,1);
                        regroup = handles.snap_nrAssigned(j,1);
                    else
                        handles.spikes_data.nrAssigned(j,1) = regroup;
                    end
                end
            end
        end
    end
    
    disp(handles.spikes_data.nrAssigned);

    handles = refresher(hObject, eventdata, handles);

guidata(hObject, handles);


function [handles] = auto_merge(hObject, eventdata, handles)

    used = zeros(1,length(handles.snap_nrAssigned(:,1)));
    for i = 1:length(handles.snap_nrAssigned(:,1))
        handles.candidates = zeros(1,length(handles.snap_nrAssigned(:,1)));
        if used(i) == 1
            continue
        end
        if handles.noise_status(i) == 1
            continue
        end
        basis = NaN(64,length(used));
        basis(:,1) = handles.sieved_means_preserved(:,i);
        basis_counter = 1;
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if j == i || used(j) == 1
                continue;
            end
            if handles.noise_status(j) == 1
                continue;
            end
            pass = 1;
            for k = 1:basis_counter
                if mean(abs(handles.sieved_means_preserved(:,j)-basis(:,k))) > 3 || std(abs(handles.sieved_means_preserved(:,j)-basis(:,k))) > 20
                    pass = 0;
                end
                [val1a, index1a] = max(handles.sieved_means_preserved(:,j));
                [val2a, index2a] = max(basis(:,k));
                
                if abs(val1a - val2a) > 4 || abs(index1a - index2a) > 2
                    pass = 0;
                end

                [val1i, index1i] = min(handles.sieved_means_preserved(:,j));
                [val2i, index2i] = min(basis(:,k));
                
                if abs(val1i - val2i) > 4 || abs(index1i - index2i) > 1
                    pass = 0;
                end
                
                thres_dist = abs(index2a - index2i);
                if abs(thres_dist - abs(index1a - index1i)) > 2
                    pass = 0;
                end
                   
            end
            if pass == 1
                basis_counter = basis_counter + 1;
                basis(:,basis_counter) = handles.sieved_means_preserved(:,j);
                handles.candidates(j) = 1;
                handles.candidates(i) = 1;
            end
        end

        parent = 0;
        for k = 1:length(handles.candidates)
            if handles.candidates(k) == 1
                if parent == 0
                    parent = handles.spikes_data.nrAssigned(k,1); 
                end
                for j = 1:length(handles.spikes_data.nrAssigned(:,1))
                    if handles.spikes_data.nrAssigned(j,1) == handles.distinct_plots(k)
                        handles.spikes_data.nrAssigned(j,1) = parent;
                    end 
                end      
            end
        end
    end
    
handles = refresher(hObject, eventdata, handles);
guidata(hObject, handles);



% % --- Executes on button press in merge.
function [handles] = merge_Callback(hObject, eventdata, handles)

    disp('executing merger');

    parent = 0;
    for i = 1:length(handles.candidates)
        if handles.candidates(i) == 1
            if parent == 0
                parent = handles.distinct_plots(i);  
            end
            for j = 1:length(handles.spikes_data.nrAssigned(:,1))
                if handles.spikes_data.nrAssigned(j,1) == handles.distinct_plots(i)
                    handles.spikes_data.nrAssigned(j,1) = parent;
                end 
            end      
        end
    end

    handles = refresher(hObject, eventdata, handles);

guidata(hObject, handles);


function [handles] = refresher(hObject, eventdata, handles)

    disp('refreshing');

    set(handles.uipanel1, 'Title', 'Cluster Name');
    set(handles.uipanel2, 'Title', 'Cluster Name');
    
    temp_arr = handles.spikes_data.nrAssigned(:,1);
    temp_arr = unique(temp_arr, 'stable');

    handles.distinct_plots = temp_arr;
    handles.distinct_plot_count = length(temp_arr);
    
    new_counter = zeros(1,handles.distinct_plot_count);
    new_means = zeros(64,handles.distinct_plot_count);
    for i = 1:length(handles.distinct_plots)
        sum = zeros(64,1);
        weight = 0;
        for j = 1:length(handles.spikes_data.nrAssigned(:,1))
            if handles.spikes_data.nrAssigned(j,1) == handles.distinct_plots(i)
                sum = sum + handles.counter_list_preserved(j)*handles.sieved_means_preserved(:,j);
                weight = weight + handles.counter_list_preserved(j);
            end
        end
        
        new_means(:,i) = sum/weight;
        new_counter(i) = weight;
    end
    
    handles.sieved_means = new_means;
    handles.counter_list = new_counter;
    
    
    handles.candidates = zeros(1,handles.distinct_plot_count);
    handles.total_pages = ceil(handles.distinct_plot_count/3);
    handles.current_page = 1;
    
    set(handles.list1, 'Value', 1);
    set(handles.list1,'String',[]);
    
    handles = page_plots(hObject, eventdata, handles);
    
    handles.pca_base_1.CData = [1 1 1];
    handles.pca_base_2.CData = [1 1 1];
    handles.pca_base_3.CData = [1 1 1];

    handles.pca_overlay_1.XData = NaN(size(handles.pca_score,1),1);
    handles.pca_overlay_1.YData = NaN(size(handles.pca_score,1),1);    
    handles.pca_overlay_2.XData = NaN(size(handles.pca_score,1),1);
    handles.pca_overlay_2.YData = NaN(size(handles.pca_score,1),1);    
    handles.pca_overlay_3.XData = NaN(size(handles.pca_score,1),1);
    handles.pca_overlay_3.YData = NaN(size(handles.pca_score,1),1);
    
    right_plots = [handles.back4, handles.back5, handles.bg3, handles.png4, handles.png5, handles.png6, handles.ac1, handles.ac2, handles.ac3, handles.ac1t, handles.ac2t, handles.ac3t];
    for i = 1:length(right_plots)
        cla(right_plots(i));
    end
    cla(handles.sr1,'reset');
    cla(handles.sr2,'reset');
    cla(handles.sr3,'reset');
    
    handles.current_highlight_details = 0;

guidata(hObject, handles);
    


function [f,Pxxn,tvect,Cxx] = psautospk(spk,tstep,nfft,window,noverlap,dflag)


    if ( (nargin ~= 6) & (nargin ~= 2) )
      disp(' ');
      disp('usage1: psautospk(spk,tstep) ');
      disp(' ');
      disp('usage2: psautospk(spk,tstep,nfft,window,noverlap,dflag) ');
      disp(' ');
      disp('       for more information type "help psautospk" in the main');
      disp('       matlab window');
      disp(' ');
      return;
    end;

    if ( nargin == 2 )
      nfft = 2048;
      window = bartlett(nfft);
      noverlap = 1024;
      dflag = 'none';
    end;

    %computes the sampling frequency in Hz
    tstep_s = tstep*1e-3;  %converts to sec
    Fs = 1/tstep_s; %in Hz

    %computes and subtracts the mean firing rate
    spk = spk(:); %convertes to column vector if necessary
    spk = spk*Fs; %converts to units of spikes/sec
    l_spk = length(spk);
    s_spk = sum(spk);
    m_spk = s_spk/l_spk;
    spk = spk - m_spk;

    % [Pxx_,f_] = psd(spk,nfft,Fs,window,noverlap,dflag);
    % %converts to units of (spk/Hz)^2
    % Pxxn_ = Pxx_ * tstep_s;

    % JD/aug16 : update to use pwelch instead of deprecated psd
    [Pxxn,f] = pwelch(spk,window,noverlap,nfft,Fs);
    %converts to units of (spk/Hz)^2
    Pxx = Pxxn / tstep_s;

    %prepares the data to compute the autocorrelation
    Pxxx = zeros(nfft,1);
    Pxxx(1:nfft/2+1,1) = Pxx(1:nfft/2+1,1);
    for k = 2:nfft/2
      Pxxx(nfft+2-k,1) = Pxx(k,1);
    end;

    %computes the autocorrelation function
    Cxxx = fft(Pxxx,nfft);
    %normalizes to get the usual definition of autocorrelation
    Cxxx = Cxxx/nfft;

    tvect = -(nfft/2)*tstep:tstep:(nfft/2)*tstep;
    Cxx = zeros(nfft+1,1);
    for k = 1:nfft/2
      Cxx(k,1) = real(Cxxx(nfft/2 + k,1));
    end;
    Cxx(nfft/2+1,1) = real(Cxxx(1,1));
    for k = nfft/2+2:nfft+2
      Cxx(k,1) = real(Cxxx(k-nfft/2,1));
    end;

function n = convertToSpiketrain(timestamps, binsize)
    if nargin<2
        binsize=1;
    end

    spiketrain=(timestamps/1000);  %now in ms
    spiketrain=spiketrain-spiketrain(1); %offset gone
    roundedSpiketrain = round(spiketrain);
    if binsize==1
        n=zeros(1,roundedSpiketrain(end));
        n( roundedSpiketrain(find(roundedSpiketrain>0)) )=1; 
    else
       n = histc( roundedSpiketrain, [0:binsize:roundedSpiketrain(end)] );
    end


% --- Executes on button press in export1.
function export1_Callback(hObject, eventdata, handles)
    
    new = handles.spikes_data.nrAssigned;
    orig = handles.snap_nrAssigned;
    noise_labels = handles.noise_status;
    
    disp(new);
    disp(noise_labels);
%     save('changes.mat', 'new', 'orig', 'noise_labels');

    assignedNegative_final = handles.spikes_data.assignedNegative;
    for i = 1:length(handles.spikes_data.assignedNegative)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == assignedNegative_final(i)
                assignedNegative_final(i) = handles.spikes_data.nrAssigned(j,1);
                break;
            end
        end
    end
    
    grouping_arr = NaN(length(handles.start_times(1,:)),length(handles.spikes_data.allSpikeInds));
    spikes_arr = NaN(length(handles.start_times(1,:)),length(handles.spikes_data.allSpikeInds),64);
    temp_arr = NaN(length(handles.start_times(1,:)),length(handles.spikes_data.allSpikeInds));
    counters = ones(1,length(handles.start_times(1,:)));
    for i = 1:length(handles.spikes_data.allSpikeInds)
        for j = 1:length(handles.start_times(1,:))
            if handles.spikes_data.allSpikeInds(i) > handles.start_times{1,j}
                target = j;
            end
        end
        
        noise = 0;
        cluster_name = assignedNegative_final(i);
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == cluster_name
                if handles.noise_status(j) == 1
                    noise = 1;
                end
            end
        end
        
        if noise == 0
            temp_arr(target,counters(target)) = handles.spikes_data.allSpikeInds(i) - handles.start_times{1,target};
            grouping_arr(target,counters(target)) = assignedNegative_final(i);
            spikes_arr(target,counters(target),:) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
            counters(target) = counters(target) + 1;
        end
    end
    
    curr_dir = pwd;
    
    for i = 1:length(handles.start_times(1,:))
        to_create = unique(grouping_arr(i,1:counters(i)-1), 'stable');
        for j = 1:length(to_create)
            if to_create(j) == 99999999
                to_create(j) = [];
                break;
            end
        end
        
        curr_dir_split = strsplit(curr_dir, '/');
        curr_dir_split{1,length(curr_dir_split)-2} = handles.start_times{2,i};

        cd(strjoin(curr_dir_split, '/'));
        temp2 = pwd;
        
        files = dir();
        dirFlags1 = [files.isdir];
        folders = files(dirFlags1);
        
        for m = 1:length(folders)
            disp('comparing');
            disp(folders(m).name);
            disp(class(folders(m).name));

            if strncmpi(folders(m).name,'old_',4) == 1
                folders(m).name
                disp('found');
                rmdir(folders(m).name, 's');
            end
        end

        files = dir();
        dirFlags1 = [files.isdir];
        folders = files(dirFlags1);        
        
        for m = 1:length(folders)
            disp('comparing');
            disp(folders(m).name);
            disp(class(folders(m).name));

            if strncmpi(folders(m).name,'cell',4) == 1
                disp('found');
                movefile(folders(m).name,strcat('old_',folders(m).name));
            end
        end
       
        
        
        for j = 1:length(to_create)
            
            if j < 10
                string_temp = strcat('cell0',num2str(j));
            else
                string_temp = strcat('cell',num2str(j));
            end               
            
            disp(string_temp); %num2str(to_create(j))
            mkdir(string_temp);
            cd(string_temp);
            
            timestamps = zeros(1,counters(i)-1);
            waves = NaN(counters(i)-1,64);

            counter2 = 1;
            for k = 1:length(grouping_arr(i,1:counters(i)-1))
                if grouping_arr(i,k) == to_create(j)
                    timestamps(1,counter2) = double(temp_arr(i,k)/0.03);
                    waves(counter2,:) = spikes_arr(i,k,:);
                    counter2 = counter2 + 1;
                end
            end
            timestamps = timestamps(1,1:counter2-1);
            
            disp(length(timestamps));
            strain.timestamps = timestamps/1000;
            strain.spikeForm = nanmean(waves, 1);
            save('spiketrain.mat', '-struct', 'strain');            
            writetable(array2table(timestamps), 'spiketrain.csv');
            cd(temp2);
        end
    end
    
    cd(curr_dir);

guidata(hObject, handles);


% --- Executes on button press in toggle2.
function toggle2_Callback(hObject, eventdata, handles)

    if handles.long2_stat == 0
        handles.long2_stat = 1;
        set(handles.long2,'visible','on');
        set(handles.long2_left,'visible','on');
        set(handles.long2_right,'visible','on');
        
        handles.long2_plot = plot(handles.long2, handles.data_bot, 'visible', 'on');
        hold(handles.long2, 'on');
        text(handles.long2, 5,min(handles.data_bot) + 0.05*(max(handles.data_bot) - min(handles.data_bot)),strcat('Minute: ', num2str(handles.spikes_data.allSpikeInds(handles.updated_index_bot)/1800000)),'FontSize',20);
        
        for j = 1:length(handles.spikes_data.allSpikeInds)
            if handles.spikes_data.allSpikeInds(j) > handles.spikes_data.allSpikeInds(handles.updated_index_bot) - 500
                if handles.spikes_data.allSpikeInds(j) < handles.spikes_data.allSpikeInds(handles.updated_index_bot) + 501
                    
                    x_pos = handles.spikes_data.allSpikeInds(j) - handles.spikes_data.allSpikeInds(handles.updated_index_bot) + 500;
                    plot(handles.long2, [x_pos x_pos], [max(handles.data_bot), min(handles.data_bot)], '-.');
                    text(handles.long2, x_pos,min(handles.data_bot) + 0.25*(max(handles.data_bot) - min(handles.data_bot)),num2str(handles.spikes_data.assignedNegative(j)),'FontSize',15);
                else
                    break;
                end
            end
        end        

        plot(handles.long2, [0 1000], [0 0]);
        plot(handles.long2, [0 1000], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long2, [0 1000], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        hold(handles.long2, 'off');
    else
        handles.long2_stat = 0;
        set(handles.long2,'visible','off');
        cla(handles.long2);
        set(handles.long2_left,'visible','off');
        set(handles.long2_right,'visible','off');
    end

guidata(hObject, handles);



% --- Executes on button press in toggle1.
function toggle1_Callback(hObject, eventdata, handles)

    if handles.long1_stat == 0
        handles.long1_stat = 1;
        set(handles.long1,'visible','on');
        set(handles.long1_left,'visible','on');
        set(handles.long1_right,'visible','on');
        
        handles.long1_plot = plot(handles.long1, handles.data_top, 'visible', 'on');
        hold(handles.long1, 'on');
        text(handles.long1, 5,min(handles.data_top) + 0.05*(max(handles.data_top) - min(handles.data_top)),strcat('Minute: ', num2str(handles.spikes_data.allSpikeInds(handles.updated_index_top)/1800000)),'FontSize',20);
        
        for j = 1:length(handles.spikes_data.allSpikeInds)
            if handles.spikes_data.allSpikeInds(j) > handles.spikes_data.allSpikeInds(handles.updated_index_top) - 500
                if handles.spikes_data.allSpikeInds(j) < handles.spikes_data.allSpikeInds(handles.updated_index_top) + 501
                    
                    x_pos = handles.spikes_data.allSpikeInds(j) - handles.spikes_data.allSpikeInds(handles.updated_index_top) + 500;
                    plot(handles.long1, [x_pos x_pos], [max(handles.data_top), min(handles.data_top)], '-.');
                    text(handles.long1, x_pos,min(handles.data_top) + 0.25*(max(handles.data_top) - min(handles.data_top)),num2str(handles.spikes_data.assignedNegative(j)),'FontSize',15);
                else
                    break;
                end
            end
        end
        
        plot(handles.long1, [0 1000], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long1, [0 1000], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long1, [0 1000], [0 0]);
        hold(handles.long1, 'off');
    else
        handles.long1_stat = 0;
        set(handles.long1,'visible','off');
        cla(handles.long1);
        set(handles.long1_left,'visible','off');
        set(handles.long1_right,'visible','off');
    end

guidata(hObject, handles);


% --- Executes on button press in long2_left.
function long2_left_Callback(hObject, eventdata, handles)

    [ts, handles.updated_index_bot] = search_ts_backward(handles, handles.updated_index_bot, handles.contents_bot);
       
    handles.data_bot = ts;
        handles.long2_plot = plot(handles.long2, handles.data_bot, 'visible', 'on');
        hold(handles.long2, 'on');
        text(handles.long2, 5,min(handles.data_bot) + 0.05*(max(handles.data_bot) - min(handles.data_bot)),strcat('Minute: ', num2str(handles.spikes_data.allSpikeInds(handles.updated_index_bot)/1800000)),'FontSize',20);
        
        for j = 1:length(handles.spikes_data.allSpikeInds)
            if handles.spikes_data.allSpikeInds(j) > handles.spikes_data.allSpikeInds(handles.updated_index_bot) - 500
                if handles.spikes_data.allSpikeInds(j) < handles.spikes_data.allSpikeInds(handles.updated_index_bot) + 501
                    x_pos = handles.spikes_data.allSpikeInds(j) - handles.spikes_data.allSpikeInds(handles.updated_index_bot) + 500;
                    plot(handles.long2, [x_pos x_pos], [max(handles.data_bot), min(handles.data_bot)], '-.');
                    text(handles.long2, x_pos, min(handles.data_bot) + 0.1*(max(handles.data_bot) - min(handles.data_bot)),num2str(handles.spikes_data.assignedNegative(j)),'FontSize',15);
                else
                    break;
                end
            end
        end

        plot(handles.long2, [0 1000], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long2, [0 1000], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long2, [0 1000], [0 0]);
        hold(handles.long2, 'off');    
    

guidata(hObject, handles);

% --- Executes on button press in long2_right.
function long2_right_Callback(hObject, eventdata, handles)

    [ts, handles.updated_index_bot] = search_ts_forward(handles, handles.updated_index_bot, handles.contents_bot);
       
    handles.data_bot = ts;
        handles.long2_plot = plot(handles.long2, handles.data_bot, 'visible', 'on');
        hold(handles.long2, 'on');
        text(handles.long2, 5,min(handles.data_bot) + 0.05*(max(handles.data_bot) - min(handles.data_bot)),strcat('Minute: ', num2str(handles.spikes_data.allSpikeInds(handles.updated_index_bot)/1800000)),'FontSize',20);
        
        for j = 1:length(handles.spikes_data.allSpikeInds)
            if handles.spikes_data.allSpikeInds(j) > handles.spikes_data.allSpikeInds(handles.updated_index_bot) - 500
                if handles.spikes_data.allSpikeInds(j) < handles.spikes_data.allSpikeInds(handles.updated_index_bot) + 501
                    x_pos = handles.spikes_data.allSpikeInds(j) - handles.spikes_data.allSpikeInds(handles.updated_index_bot) + 500;
                    plot(handles.long2, [x_pos x_pos], [max(handles.data_bot), min(handles.data_bot)], '-.');
                    text(handles.long2, x_pos, min(handles.data_bot) + 0.1*(max(handles.data_bot) - min(handles.data_bot)),num2str(handles.spikes_data.assignedNegative(j)),'FontSize',15);
                else
                    break;
                end
            end
        end
        
        plot(handles.long2, [0 1000], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long2, [0 1000], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long2, [0 1000], [0 0]);
        hold(handles.long2, 'off');
    
guidata(hObject, handles);

% --- Executes on button press in long1_left.
function long1_left_Callback(hObject, eventdata, handles)
    
    
    [ts, handles.updated_index_top] = search_ts_backward(handles, handles.updated_index_top, handles.contents_top);
       
    handles.data_top = ts;
        handles.long1_plot = plot(handles.long1, handles.data_top, 'visible', 'on');
        hold(handles.long1, 'on');
        text(handles.long1, 5,min(handles.data_top) + 0.05*(max(handles.data_top) - min(handles.data_top)),strcat('Minute: ', num2str(handles.spikes_data.allSpikeInds(handles.updated_index_top)/1800000)),'FontSize',20);
        
        for j = 1:length(handles.spikes_data.allSpikeInds)
            if handles.spikes_data.allSpikeInds(j) > handles.spikes_data.allSpikeInds(handles.updated_index_top) - 500
                if handles.spikes_data.allSpikeInds(j) < handles.spikes_data.allSpikeInds(handles.updated_index_top) + 501
                    x_pos = handles.spikes_data.allSpikeInds(j) - handles.spikes_data.allSpikeInds(handles.updated_index_top) + 500;
                    plot(handles.long1, [x_pos x_pos], [max(handles.data_top), min(handles.data_top)], '-.');
                    text(handles.long1, x_pos, min(handles.data_top) + 0.1*(max(handles.data_top) - min(handles.data_top)),num2str(handles.spikes_data.assignedNegative(j)),'FontSize',15);
                else
                    break;
                end
            end
        end

        plot(handles.long1, [0 1000], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long1, [0 1000], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long1, [0 1000], [0 0]);
        hold(handles.long1, 'off');    
    

guidata(hObject, handles);


% --- Executes on button press in long1_right.
function long1_right_Callback(hObject, eventdata, handles)

    [ts, handles.updated_index_top] = search_ts_forward(handles, handles.updated_index_top, handles.contents_top);
       
    handles.data_top = ts;
        handles.long1_plot = plot(handles.long1, handles.data_top, 'visible', 'on');
        hold(handles.long1, 'on');
        text(handles.long1, 5,min(handles.data_top) + 0.05*(max(handles.data_top) - min(handles.data_top)),strcat('Minute: ', num2str(handles.spikes_data.allSpikeInds(handles.updated_index_top)/1800000)),'FontSize',20);
        
        for j = 1:length(handles.spikes_data.allSpikeInds)
            if handles.spikes_data.allSpikeInds(j) > handles.spikes_data.allSpikeInds(handles.updated_index_top) - 500
                if handles.spikes_data.allSpikeInds(j) < handles.spikes_data.allSpikeInds(handles.updated_index_top) + 501
                    x_pos = handles.spikes_data.allSpikeInds(j) - handles.spikes_data.allSpikeInds(handles.updated_index_top) + 500;
                    plot(handles.long1, [x_pos x_pos], [max(handles.data_top), min(handles.data_top)], '-.');
                    text(handles.long1, x_pos, min(handles.data_top) + 0.1*(max(handles.data_top) - min(handles.data_top)),num2str(handles.spikes_data.assignedNegative(j)),'FontSize',15);
                else
                    break;
                end
            end
        end

        plot(handles.long1, [0 1000], [5*handles.spikes_data.stdEstimateOrig 5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);
        plot(handles.long1, [0 1000], [-5*handles.spikes_data.stdEstimateOrig -5*handles.spikes_data.stdEstimateOrig], '--', 'LineWidth', 2);        
        plot(handles.long1, [0 1000], [0 0]);
        hold(handles.long1, 'off');

guidata(hObject, handles);


function noise1_Callback(hObject, eventdata, handles)

    contents = cellstr(get(handles.list1,'String'));
    contents = str2double(contents);
    contents = contents(2:length(contents));
    
    for i = 1:length(contents)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.noise_status(j) = get(handles.noise1,'Value');
%                 if handles.noise_status(j) == 1
%                     if handles.candidates(j) == 1
%                         handles.candidates(j) = 0;
%                         cla(handles.bg3);
%                         cla(handles.png6); 
%                         set(handles.pca_base_3, 'visible', 'off');
%                         set(handles.pca_overlay_3, 'visible', 'off');
%                         cla(handles.ac3, 'reset');
%                         cla(handles.sr3, 'reset');        
% 
%                         if sum(handles.candidates) > 0
%                             handles = merge_plot1(hObject, eventdata, handles);
%                             handles = merge_plot2(hObject, eventdata, handles);
%                             handles = merge_plot34(hObject, eventdata, handles);
%                         end
%                     end
%                 end
            end
        end
    end
    
    disp('noise: ');
    disp(handles.noise_status);
    
    handles = page_plots(hObject, eventdata, handles);
    handles = update_colors(hObject, eventdata, handles);
        
guidata(hObject, handles);











