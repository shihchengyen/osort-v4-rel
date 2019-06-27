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
    disp('asdfawfasfd');
    disp(length(handles.hp_trace));
    
%     load_highpass = rplhighpass('auto');
%     handles.hp_trace = load_highpass.data.analogData;
%     handles.time_trace = load_highpass.data.analogTime;
    
    channel_path = cd('./oSort/');
    folder_name = dir();
    cd(strcat(folder_name(length(folder_name)).name, '/sort/'));
    folder_name = dir();
    cd(folder_name(length(folder_name)).name);
    handles.spikes_data = load('P1_sorted_new.mat');
    cd(channel_path);

    disp(handles.spikes_data.nrAssigned);

    handles.cluster_names = handles.spikes_data.nrAssigned(:,1);
    handles.sieved = NaN(length(handles.spikes_data.assignedNegative), 64, length(handles.spikes_data.nrAssigned));
    handles.sieved_long = NaN(length(handles.spikes_data.assignedNegative), 1000, length(handles.spikes_data.nrAssigned));
    handles.counter_list = ones(1,length(handles.spikes_data.nrAssigned));
    handles.snap_nrAssigned = handles.spikes_data.nrAssigned;
    handles.noise_status = zeros(1,length(handles.spikes_data.nrAssigned));
    
    temp_arr = handles.spikes_data.nrAssigned(:,1);
    temp_arr = unique(temp_arr, 'stable');

    handles.distinct_plots = temp_arr;
    handles.distinct_plot_count = length(temp_arr);

    handles.max_max = 0;
    handles.min_min = 0;

    for i = 1:length(handles.spikes_data.assignedNegative)
        for j = 1:length(handles.spikes_data.nrAssigned)
            if handles.spikes_data.assignedNegative(i) == handles.snap_nrAssigned(j,1)
                for k = 1:length(handles.cluster_names)
                    if handles.cluster_names(k) == handles.spikes_data.nrAssigned(j,1)
                        handles.sieved(handles.counter_list(k),:,k) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
                        handles.sieved_long(handles.counter_list(k),:,k) = handles.hp_trace(handles.spikes_data.allSpikeInds(i)-499:handles.spikes_data.allSpikeInds(i)+500);

                        if max(handles.sieved(handles.counter_list(k),:,k)) > handles.max_max
                            handles.max_max = max(handles.sieved(handles.counter_list(k),:,k));
                        end
                        if min(handles.sieved(handles.counter_list(k),:,k)) < handles.min_min
                            handles.min_min = min(handles.sieved(handles.counter_list(k),:,k));
                        end

                        handles.counter_list(k) = handles.counter_list(k) + 1;
                        break;
                    end
                end
            end
        end
    end
    for i = 1:length(handles.counter_list)
        handles.counter_list(i) = handles.counter_list(i) - 1;
    end

    handles.sieved_means = nanmean(handles.sieved, 1);
    handles.sieved_means_preserved = handles.sieved_means;
    handles.sieved_long_preserved = handles.sieved_long;
    handles.counter_list_preserved = handles.counter_list;
    
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

    [~, handles.pca_score, ~] = pca(handles.noiseless_spikes);
    
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
    

guidata(hObject, handles);


function [handles] = page_plots(hObject, eventdata, handles)

    back_arr = [handles.back1, handles.back2, handles.back3];
    png_arr = [handles.png1, handles.png2, handles.png3];
    
    for i = 1:3
        cla(back_arr(i));
        cla(png_arr(i));
    end
    
    count = 1;
    handles.box_handles = cell(1,3,length(handles.spikes_data.nrAssigned(:,1)));

    start = 1 + 3*(handles.current_page - 1);

    for i = start:min(start+2, handles.distinct_plot_count)
        
        [children, indices, children_length] = find_children(handles.distinct_plots(i), hObject, eventdata, handles);
        
        
        for j = 1:length(indices)
            if j == 1
                handles.box_handles{1,count,j} = image(png_arr(count), handles.png_stuff{indices(j),1}, 'AlphaData', handles.png_stuff{indices(j),3});
                set(handles.box_handles{1,count,j}, 'AlphaData', handles.png_stuff{indices(j),3});
                hold(png_arr(count), 'on');
            else
                handles.box_handles{1,count,j} = image(png_arr(count), handles.png_stuff{indices(j),1}, 'AlphaData', handles.png_stuff{indices(j),3});
                set(handles.box_handles{1,count,j}, 'AlphaData', handles.png_stuff{indices(j),3});
            end
        end
           
        hold(png_arr(count), 'off');

        set(png_arr(count),'visible','off');
        set(png_arr(count), 'XLimSpec', 'Tight');
        set(png_arr(count), 'YLimSpec', 'Tight');

        plot(back_arr(count), 1:64, handles.sieved_means(1,:,i), 'r', 'LineWidth', 2);
        ylim(back_arr(count), [handles.min_min, handles.max_max]);
        xlim(back_arr(count), [0 64]);
        set(back_arr(count), 'Color', 'None');
     
        
        count = count + 1;

    end
    
    handles = page_lists(hObject, eventdata, handles);
    handles = page_noise(hObject, eventdata, handles);

guidata(hObject, handles);

function [handles] = page_noise(hObject, eventdata, handles)

    set(handles.noise1,'visible','off');
    set(handles.noise2,'visible','off');
    set(handles.noise3,'visible','off');

    noise_list = [handles.noise1, handles.noise2, handles.noise3];
    lists_list = [handles.list1, handles.list2, handles.list3];
    
    for i = 1:length(lists_list)

        if i + 3*(handles.current_page - 1) > handles.distinct_plot_count
            break;
        end
        
        contents = cellstr(get(lists_list(i),'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
        
        if length(contents) == 1
            set(noise_list(i),'visible','on');
            for j = 1:length(handles.noise_status)
                if handles.snap_nrAssigned(j,1) == contents(1)
                    set(noise_list(i),'value',handles.noise_status(j));
                end
            end
        end
        
    end

guidata(hObject, handles);


function [handles] = page_lists(hObject, eventdata, handles)

list_of_lists = [handles.list1, handles.list2, handles.list3];
for counter = 1:length(list_of_lists)
    set(list_of_lists(counter), 'Value', 1);
    set(list_of_lists(counter),'String',[]);
end

index_to_display = 1 + 3*(handles.current_page - 1);

for counter = 1:length(list_of_lists)
    if index_to_display > handles.distinct_plot_count
        break;
    end
    [children, ~, count] = find_children(handles.distinct_plots(index_to_display), hObject, eventdata, handles);
    children = mat2cell(children, 1, count);
    children = [{'View All'},children];
    set(list_of_lists(counter),'String', children);
    index_to_display = index_to_display + 1;
end


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

function [handles] = merge_plot1(hObject, eventdata, handles)

        cla(handles.back6);
        cla(handles.png6);
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
                plot(handles.back6, 1:64, handles.sieved_means(1,:,indices(j)), 'r', 'LineWidth', 2);
                hold(handles.back6, 'on');
            else
                plot(handles.back6, 1:64, handles.sieved_means(1,:,indices(j)), 'r', 'LineWidth', 2);
            end
        end
        hold(handles.back6, 'off');
                
        ylim(handles.back6, [handles.min_min, handles.max_max]);
        xlim(handles.back6, [0 64]);
        set(handles.back6, 'Color', 'None');

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


    out = pca_siever(handles, contents, handles.pca_score);

    handles.pca_base_3.CData = [0 1 0];

    handles.pca_overlay_3.XData = out(:,1);
    handles.pca_overlay_3.YData = out(:,2);


guidata(hObject, handles);
    
function [handles] = merge_plot34(hObject, eventdata, handles)

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
        
        plot(handles.ac3, tvect, Cxx.');
        
        %%%%% end autocorr/spikerates part %%%%%

guidata(hObject, handles);

% --- Executes just before combined_inspection is made visible.
function combined_inspection_OpeningFcn(hObject, eventdata, handles, varargin)

    handles = initialize_data(hObject, eventdata, handles);
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
    
    merge_list = [handles.merge1, handles.merge2, handles.merge3];
    
    for i = 1:3
        if i + 3*(handles.current_page - 1) > handles.distinct_plot_count
            break;
        end
        if handles.noise_status(i + 3*(handles.current_page - 1)) == 1
            set(merge_list(i),'visible','off');
            set(merge_list(i),'value',0);
        else
            set(merge_list(i),'visible','on');
            set(merge_list(i),'value',handles.candidates(i + 3*(handles.current_page - 1)));
        end
    end
end

guidata(hObject, handles);


function back_Callback(hObject, eventdata, handles)

if handles.current_page > 1
    handles.current_page = handles.current_page - 1;
    handles = page_plots(hObject, eventdata, handles);
    
    merge_list = [handles.merge1, handles.merge2, handles.merge3];
    
    for i = 1:3
        if i + 3*(handles.current_page - 1) > handles.distinct_plot_count
            break;
        end
        if handles.noise_status(i + 3*(handles.current_page - 1)) == 1
            set(merge_list(i),'visible','off');
            set(merge_list(i),'value',0);
        else
            set(merge_list(i),'visible','on');
            set(merge_list(i),'value',handles.candidates(i + 3*(handles.current_page - 1)));
        end
    end
end

guidata(hObject, handles);



function list1_Callback(hObject, eventdata, handles)

function list1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function list2_Callback(hObject, eventdata, handles)

function list2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function list3_Callback(hObject, eventdata, handles)

function list3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function details1_Callback(hObject, eventdata, handles)

    
    contents = cellstr(get(handles.list1,'String'));
    contents = str2double(contents{get(handles.list1,'Value')});
    
    singular = 1;
    if isnan(contents) == 1
        singular = 0;
        contents = cellstr(get(handles.list1,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
    
      
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
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(1,:,index_for_png(1,j)), 'r', 'LineWidth', 2);
                hold(axes_handles2(handles.current_display), 'on');
            else
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(1,:,index_for_png(1,j)), 'r', 'LineWidth', 2);
            end            
        end
        
        hold(axes_handles2(handles.current_display), 'off');
        
        ylim(axes_handles2(handles.current_display), [handles.min_min, handles.max_max]);
        xlim(axes_handles2(handles.current_display), [0 64]);
        set(axes_handles2(handles.current_display), 'Color', 'None');   
        
            
        %%%%%%% pca part %%%%%%%%
        
        out = pca_siever(handles, contents, handles.pca_score);
        
        pca_invis(handles.current_display).CData = [0 1 0];
        
        pca_overlay(handles.current_display).XData = out(:,1);
        pca_overlay(handles.current_display).YData = out(:,2);
        
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
        xlim([0,(time_axis(end)+5)]);
        
        Cxx(1) = [];

        tvect = tvect(floor(length(tvect)/2):length(tvect));
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        tvect = tvect(2:100);
        Cxx = Cxx(2:100);
        
        plot(axes_ac(handles.current_display), tvect, Cxx.');
        
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
    
    if handles.current_display == 1
        if singular == 1
            
          
            handles.data_top = handles.sieved_long_preserved(1:handles.counter_list_preserved(index_for_png(1)),:,index_for_png(1));
        else

            
            handles.data_top = handles.sieved_long(1:handles.counter_list(3*(handles.current_page-1)+1),:,3*(handles.current_page-1)+1); % dummy
        end
        handles.current_display = 2;
    else
        if singular == 1
            
           
            handles.data_btm = handles.sieved_long_preserved(1:handles.counter_list_preserved(index_for_png(1)),:,index_for_png(1));
        else
            
          
            handles.data_btm = handles.sieved_long(1:handles.counter_list(3*(handles.current_page-1)+1),:,3*(handles.current_page-1)+1); % dummy
        end        
        handles.current_display = 1;
    end
        

guidata(hObject, handles);


function [sieved_pca] = pca_siever(handles, target_arr, full_pca)

sieved_pca = NaN(size(full_pca,1),2);    
    
counter = 1;
for i = 1:length(handles.removing_noise_index)
    for j = 1:length(target_arr)
        if handles.removing_noise_index(i) == target_arr(j)            
            sieved_pca(counter,:) = full_pca(i,1:2);
            counter = counter + 1;
            break;
        end
    end
end


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
    
    for i = 1:length(contents)
        for j = 1:length(handles.spikes_data.nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.spikes_data.nrAssigned(j,1) = handles.snap_nrAssigned(j,1);
                break;
            end
        end
    end
    
    disp(handles.spikes_data.nrAssigned);

    handles = refresher(hObject, eventdata, handles);

guidata(hObject, handles);

function kick2_Callback(hObject, eventdata, handles)

    disp('kicking');
    contents = cellstr(get(handles.list2,'String'));
    contents = str2double(contents{get(handles.list2,'Value')});

    if isnan(contents) == 1
        contents = cellstr(get(handles.list2,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
    
    disp(contents);
    disp(handles.spikes_data.nrAssigned);
    
    for i = 1:length(contents)
        for j = 1:length(handles.spikes_data.nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.spikes_data.nrAssigned(j,1) = handles.snap_nrAssigned(j,1);
                break;
            end
        end
    end
    
    disp(handles.spikes_data.nrAssigned);

    handles = refresher(hObject, eventdata, handles);

guidata(hObject, handles);

function kick3_Callback(hObject, eventdata, handles)

    disp('kicking');
    contents = cellstr(get(handles.list3,'String'));
    contents = str2double(contents{get(handles.list3,'Value')});

    if isnan(contents) == 1
        contents = cellstr(get(handles.list3,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
    
    disp(contents);
    disp(handles.spikes_data.nrAssigned);
    
    for i = 1:length(contents)
        for j = 1:length(handles.spikes_data.nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.spikes_data.nrAssigned(j,1) = handles.snap_nrAssigned(j,1);
                break;
            end
        end
    end
    
    disp(handles.spikes_data.nrAssigned);

    handles = refresher(hObject, eventdata, handles);

guidata(hObject, handles);

function merge1_Callback(hObject, eventdata, handles)

    handles.candidates(1 + 3*(handles.current_page - 1)) = get(handles.merge1,'Value');
    handles = merge_plot1(hObject, eventdata, handles);
    handles = merge_plot2(hObject, eventdata, handles);
    handles = merge_plot34(hObject, eventdata, handles);
    
guidata(hObject, handles);

function noise1_Callback(hObject, eventdata, handles)

    contents = cellstr(get(handles.list1,'String'));
    contents = str2double(contents);
    contents = contents(2:length(contents));
    
    for i = 1:length(contents)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.noise_status(j) = get(handles.noise1,'Value');
                if handles.noise_status(j) == 1
                    set(handles.merge1,'visible','off');
                else
                    set(handles.merge1,'visible','on');
                end               
            end
        end
    end
    
    refresher(hObject, eventdata, handles);

guidata(hObject, handles);


function details2_Callback(hObject, eventdata, handles)


    contents = cellstr(get(handles.list2,'String'));
    contents = str2double(contents{get(handles.list2,'Value')});

    singular = 1;
    if isnan(contents) == 1
        singular = 0;
        contents = cellstr(get(handles.list2,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
  
    disp(handles.spikes_data.nrAssigned);
    disp(contents);
    
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
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(1,:,index_for_png(1,j)), 'r', 'LineWidth', 2);
                hold(axes_handles2(handles.current_display), 'on');
            else
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(1,:,index_for_png(1,j)), 'r', 'LineWidth', 2);
            end            
        end
        
        hold(axes_handles2(handles.current_display), 'off');
        
        ylim(axes_handles2(handles.current_display), [handles.min_min, handles.max_max]);
        xlim(axes_handles2(handles.current_display), [0 64]);
        set(axes_handles2(handles.current_display), 'Color', 'None');    
        
        %%%%%%% pca part %%%%%%%%
        
        out = pca_siever(handles, contents, handles.pca_score);
        
        pca_invis(handles.current_display).CData = [0 1 0];
        
        pca_overlay(handles.current_display).XData = out(:,1);
        pca_overlay(handles.current_display).YData = out(:,2);
        
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
        xlim([0,(time_axis(end)+5)]);
        
        Cxx(1) = [];

        tvect = tvect(floor(length(tvect)/2):length(tvect));
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        tvect = tvect(2:100);
        Cxx = Cxx(2:100);
        
        plot(axes_ac(handles.current_display), tvect, Cxx.');
        
        %%%%% end autocorr/spikerates part %%%%%        

    long_arr = [handles.long1, handles.long2];
    long_arr_l = [handles.long1_left, handles.long2_left];
    long_arr_r = [handles.long1_right, handles.long2_right];
    long_arr_stat = [handles.long1_stat, handles.long2_stat];

    cla(long_arr(handles.current_display));
    set(long_arr(handles.current_display),'visible','off');

    set(long_arr_l(handles.current_display), 'visible', 'off');
    set(long_arr_r(handles.current_display), 'visible', 'off');

    if handles.current_display == 1
        handles.long1_stat = 0;
    else
        handles.long2_stat = 0;
    end
    
    if handles.current_display == 1
        if singular == 1
            handles.data_top = handles.sieved_long_preserved(1:handles.counter_list_preserved(index_for_png(1)),:,index_for_png(1));
        else
            handles.data_top = handles.sieved_long(1:handles.counter_list(3*(handles.current_page-1)+2),:,3*(handles.current_page-1)+2); % dummy
        end
        handles.current_display = 2;
    else
        if singular == 1
            handles.data_btm = handles.sieved_long_preserved(1:handles.counter_list_preserved(index_for_png(1)),:,index_for_png(1));
        else
            handles.data_btm = handles.sieved_long(1:handles.counter_list(3*(handles.current_page-1)+2),:,3*(handles.current_page-1)+2); % dummy
        end        
        handles.current_display = 1;
    end
        
    
guidata(hObject, handles);



function details3_Callback(hObject, eventdata, handles)


    contents = cellstr(get(handles.list3,'String'));
    contents = str2double(contents{get(handles.list3,'Value')});

    singular = 1;
    if isnan(contents) == 1
        singular = 0;
        contents = cellstr(get(handles.list3,'String'));
        contents = str2double(contents);
        contents = contents(2:length(contents));
    end
    
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
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(1,:,index_for_png(1,j)), 'r', 'LineWidth', 2);
                hold(axes_handles2(handles.current_display), 'on');
            else
                plot(axes_handles2(handles.current_display), 1:64, handles.sieved_means_preserved(1,:,index_for_png(1,j)), 'r', 'LineWidth', 2);
            end            
        end
        
        hold(axes_handles2(handles.current_display), 'off');
        
        ylim(axes_handles2(handles.current_display), [handles.min_min, handles.max_max]);
        xlim(axes_handles2(handles.current_display), [0 64]);
        set(axes_handles2(handles.current_display), 'Color', 'None');    
        
        
        %%%%%%% pca part %%%%%%%%
        
        out = pca_siever(handles, contents, handles.pca_score);
        
        pca_invis(handles.current_display).CData = [0 1 0];
        
        pca_overlay(handles.current_display).XData = out(:,1);
        pca_overlay(handles.current_display).YData = out(:,2);
        
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
        xlim([0,(time_axis(end)+5)]);
        
        Cxx(1) = [];

        tvect = tvect(floor(length(tvect)/2):length(tvect));
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        tvect = tvect(2:100);
        Cxx = Cxx(2:100);
        
        plot(axes_ac(handles.current_display), tvect, Cxx.');
        
        %%%%% end autocorr/spikerates part %%%%%        

    long_arr = [handles.long1, handles.long2];
    long_arr_l = [handles.long1_left, handles.long2_left];
    long_arr_r = [handles.long1_right, handles.long2_right];
    long_arr_stat = [handles.long1_stat, handles.long2_stat];
           
    cla(long_arr(handles.current_display));
    set(long_arr(handles.current_display),'visible','off');

    set(long_arr_l(handles.current_display), 'visible', 'off');
    set(long_arr_r(handles.current_display), 'visible', 'off');
    
    if handles.current_display == 1
        handles.long1_stat = 0;
    else
        handles.long2_stat = 0;
    end     
    
    if handles.current_display == 1
        if singular == 1
            handles.data_top = handles.sieved_long_preserved(1:handles.counter_list_preserved(index_for_png(1)),:,index_for_png(1));
        else
            handles.data_top = handles.sieved_long(1:handles.counter_list(3*(handles.current_page-1)+3),:,3*(handles.current_page-1)+3); % dummy
        end
        handles.current_display = 2;
    else
        if singular == 1
            handles.data_btm = handles.sieved_long_preserved(1:handles.counter_list_preserved(index_for_png(1)),:,index_for_png(1));
        else
            handles.data_btm = handles.sieved_long(1:handles.counter_list(3*(handles.current_page-1)+3),:,3*(handles.current_page-1)+3); % dummy
        end        
        handles.current_display = 1;
    end        
        
    
guidata(hObject, handles);



% --- Executes on button press in merge2.
function merge2_Callback(hObject, eventdata, handles)

    handles.candidates(2 + 3*(handles.current_page - 1)) = get(handles.merge2,'Value');
    handles = merge_plot1(hObject, eventdata, handles);
    handles = merge_plot2(hObject, eventdata, handles);
    handles = merge_plot34(hObject, eventdata, handles);
    
guidata(hObject, handles);


function noise2_Callback(hObject, eventdata, handles)

    contents = cellstr(get(handles.list2,'String'));
    contents = str2double(contents);
    contents = contents(2:length(contents));
    
    for i = 1:length(contents)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.noise_status(j) = get(handles.noise2,'Value');
                if handles.noise_status(j) == 1
                    set(handles.merge2,'visible','off');
                else
                    set(handles.merge2,'visible','on');
                end
            end
        end
    end

%     refresher(hObject, eventdata, handles);
    
guidata(hObject, handles);

function noise3_Callback(hObject, eventdata, handles)

    contents = cellstr(get(handles.list3,'String'));
    contents = str2double(contents);
    contents = contents(2:length(contents));
    
    for i = 1:length(contents)
        for j = 1:length(handles.snap_nrAssigned(:,1))
            if handles.snap_nrAssigned(j,1) == contents(i)
                handles.noise_status(j) = get(handles.noise3,'Value');
                if handles.noise_status(j) == 1
                    set(handles.merge3,'visible','off');
                else
                    set(handles.merge2,'visible','on');
                end
            end
        end
    end
    
    
%     refresher(hObject, eventdata, handles);

guidata(hObject, handles);


% --- Executes on button press in merge3.
function merge3_Callback(hObject, eventdata, handles)

    handles.candidates(3 + 3*(handles.current_page - 1)) = get(handles.merge3,'Value');
    
    handles = merge_plot1(hObject, eventdata, handles);
    handles = merge_plot2(hObject, eventdata, handles);
    handles = merge_plot34(hObject, eventdata, handles);
    
guidata(hObject, handles);


% --- Executes on button press in merge.
function merge_Callback(hObject, eventdata, handles)

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
    
    sieved = NaN(length(handles.spikes_data.assignedNegative), 64, handles.distinct_plot_count);
    counter_list = ones(1,length(handles.distinct_plots));
    
    for i = 1:length(handles.spikes_data.assignedNegative)
        for j = 1:length(handles.snap_nrAssigned)
            if handles.spikes_data.assignedNegative(i) == handles.snap_nrAssigned(j,1)
                for k = 1:length(handles.distinct_plots)
                    if handles.distinct_plots(k) == handles.spikes_data.nrAssigned(j,1)
                        sieved(counter_list(k),:,k) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
                        handles.sieved_long(counter_list(k),:,k) = handles.hp_trace(handles.spikes_data.allSpikeInds(i)-499:handles.spikes_data.allSpikeInds(i)+500);
                        if max(sieved(counter_list(k),:,k)) > handles.max_max
                            handles.max_max = max(sieved(counter_list(k),:,k));
                        end
                        if min(sieved(counter_list(k),:,k)) < handles.min_min
                            handles.min_min = min(sieved(counter_list(k),:,k));
                        end

                        counter_list(k) = counter_list(k) + 1;
                        break;
                    end
                end
            end
        end
    end
    for i = 1:length(counter_list)
        counter_list(i) = counter_list(i) - 1;
    end
    handles.counter_list = counter_list;
    
    handles.sieved = sieved;
    handles.sieved_means = nanmean(handles.sieved, 1);
    
    handles.candidates = zeros(1,handles.distinct_plot_count);
    handles.total_pages = ceil(handles.distinct_plot_count/3);
    handles.current_page = 1;
    
    handles = page_plots(hObject, eventdata, handles);
    set(handles.merge1, 'value', 0);
    set(handles.merge2, 'value', 0);
    set(handles.merge3, 'value', 0);
    
    handles.pca_base_1.CData = [1 1 1];
    handles.pca_base_2.CData = [1 1 1];
    handles.pca_base_3.CData = [1 1 1];

    handles.pca_overlay_1.XData = NaN(size(handles.pca_score,1),1);
    handles.pca_overlay_1.YData = NaN(size(handles.pca_score,1),1);    
    handles.pca_overlay_2.XData = NaN(size(handles.pca_score,1),1);
    handles.pca_overlay_2.YData = NaN(size(handles.pca_score,1),1);    
    handles.pca_overlay_3.XData = NaN(size(handles.pca_score,1),1);
    handles.pca_overlay_3.YData = NaN(size(handles.pca_score,1),1);
    
    right_plots = [handles.back4, handles.back5, handles.back6, handles.png4, handles.png5, handles.png6, handles.ac1, handles.ac2, handles.ac3];
    for i = 1:length(right_plots)
        cla(right_plots(i));
    end
    cla(handles.sr1,'reset');
    cla(handles.sr2,'reset');
    cla(handles.sr3,'reset');

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
            
            disp(strcat('cell',num2str(to_create(j))));
            mkdir(strcat('cell',num2str(to_create(j))));
            cd(strcat('cell',num2str(to_create(j))));
            
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
            strain.timestamps = timestamps;
            strain.spikeForm = nanmean(waves, 1);
            save('spiketrain.mat', '-struct', 'strain');            
            
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
        handles.number_btm = 1;
        handles.long2_plot = plot(handles.long2, handles.data_btm(handles.number_btm,:), 'visible', 'on');
        hold(handles.long2, 'on');
        plot(handles.long2, [350 350], [max(handles.data_btm(handles.number_btm,:)), min(handles.data_btm(handles.number_btm,:))], '--');
        plot(handles.long2, [650 650], [max(handles.data_btm(handles.number_btm,:)), min(handles.data_btm(handles.number_btm,:))], '--');
        plot(handles.long2, [0 1000], [0 0]);
        hold(handles.long2, 'off');
    else
        handles.long2_stat = 0;
        cla(handles.long2);
        set(handles.long2,'visible','off');
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
        handles.number_top = 1;
        handles.long1_plot = plot(handles.long1, handles.data_top(handles.number_top,:), 'visible', 'on');
        hold(handles.long1, 'on');
        plot(handles.long1, [350 350], [max(handles.data_top(handles.number_top,:)), min(handles.data_top(handles.number_top,:))], '--');
        plot(handles.long1, [650 650], [max(handles.data_top(handles.number_top,:)), min(handles.data_top(handles.number_top,:))], '--');
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

    if handles.number_btm > 1
        handles.number_btm = handles.number_btm - 1;
    end

    set(handles.long2_plot, 'YData', handles.data_btm(handles.number_btm,:));

guidata(hObject, handles);

% --- Executes on button press in long2_right.
function long2_right_Callback(hObject, eventdata, handles)

    if handles.number_btm < 100
        handles.number_btm = handles.number_btm + 1;
    end

    set(handles.long2_plot, 'YData', handles.data_btm(handles.number_btm,:));
    
guidata(hObject, handles);

% --- Executes on button press in long1_left.
function long1_left_Callback(hObject, eventdata, handles)
    
    if handles.number_top > 1
        handles.number_top = handles.number_top - 1;
    end

    set(handles.long1_plot, 'YData', handles.data_top(handles.number_top,:));

guidata(hObject, handles);


% --- Executes on button press in long1_right.
function long1_right_Callback(hObject, eventdata, handles)

    if handles.number_top < 100
        handles.number_top = handles.number_top + 1;
    end

    set(handles.long1_plot, 'YData', handles.data_top(handles.number_top,:));

guidata(hObject, handles);





