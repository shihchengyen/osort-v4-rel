function create_pngs(thres_string)

%     channel_no = current_path;
%     disp(channel_no);
%     handles.channel_no = strcat('channel',extractAfter(channel_no, string('channel')));
    handles.modification_tracker = [];

    channel_path = cd('./oSort/');
    folder_name = dir();
    % saved_identifier = folder_name(length(folder_name)).name;
    cd(strcat(thres_string, '/sort/')); % cd(strcat('detect1Thresh6kern18', '/sort/'));
    folder_name = dir();
    cd(folder_name(length(folder_name)).name);
    mat_name_d = dir('*_sorted_new.mat');
    mat_name = {mat_name_d.name};
    handles.spikes_data = load(mat_name{1}); % 'P1_sorted_new.mat'
    cd(channel_path);
    
    title = ['pngs_', thres_string];
    mkdir(title);
    cd(title);

    disp(handles.spikes_data.nrAssigned);

    sieved = NaN(length(handles.spikes_data.assignedNegative), 64, length(handles.spikes_data.nrAssigned(:,1)));
    counter_list = ones(1,length(handles.spikes_data.nrAssigned(:,1)));
    
    max_arr = NaN(1, length(handles.spikes_data.assignedNegative));
    min_arr = NaN(1, length(handles.spikes_data.assignedNegative));
    new_counter = 1;
    
    for i = 1:length(handles.spikes_data.assignedNegative)
        for j = 1:length(handles.spikes_data.nrAssigned(:,1))
            if handles.spikes_data.assignedNegative(i) == handles.spikes_data.nrAssigned(j,1)
                for k = 1:length(handles.spikes_data.nrAssigned(:,1))
                    if handles.spikes_data.nrAssigned(k,1) == handles.spikes_data.nrAssigned(j,1)
                        sieved(counter_list(k),:,k) = downsample(handles.spikes_data.newSpikesNegative(i,:),4,2);
                        
                        max_arr(new_counter) = max(sieved(counter_list(k),:,k));
                        min_arr(new_counter) = min(sieved(counter_list(k),:,k));
                        
%                         if max(sieved(counter_list(k),:,k)) > max_max
%                             max_max = max(sieved(counter_list(k),:,k));
%                         end
%                         if min(sieved(counter_list(k),:,k)) < min_min
%                             min_min = min(sieved(counter_list(k),:,k));
%                         end
                        new_counter = new_counter + 1;
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
    
    max_arr = max_arr(1:new_counter-1);
    min_arr = min_arr(1:new_counter-1);
    
    max_max = prctile(max_arr, 97.5) + 50;
    min_min = prctile(min_arr, 2.5) - 50;
    
    sieved_means = nanmean(sieved, 1);
    meandata = cell(4, length(handles.spikes_data.nrAssigned(:,1)));
    for i = 1:length(handles.spikes_data.nrAssigned(:,1))
        meandata{1,i} = sieved_means(1,:,i);
        meandata{2,i} = counter_list(i);
    end
    meandata{3,1} = min_min;
    meandata{4,1} = max_max;
    
    save('meandata.mat', 'meandata');

    fig1 = figure();
    ax1 = axes('Parent', fig1);

    for i = 1:length(counter_list)
        
        disp(size(counter_list));
        disp(size(sieved));
        disp(i);
        cluster1 = plot(ax1, 1:64, sieved(1:counter_list(i),:,i), '-', 'Color', [0 0 0 0.05]);
        hold(ax1,'on');
        cluster2 = plot(ax1, [0 64], [0 0], '-.', 'Color', [0 0 0], 'LineWidth', 3);
        hold(ax1,'off');
        xlim([0 64]);
        ylim([min_min max_max]);
        grid on
        set(ax1,'XTickLabel',[], 'YTickLabel', []);
        set(ax1, 'Color', 'None');
        ax = gca;

        if i == 1
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end

        ax.LineWidth = 3;
        ax.YTick = round(min_min/50)*50:50:max_max;
        ax.XTick = [10 20 30 40 50 60];
        set(ax, 'Color', 'None');

        file_name = char(strcat('cluster',string(handles.spikes_data.nrAssigned(i,1)),'.png'));
        file_name2 = char(strcat('cluster',string(handles.spikes_data.nrAssigned(i,1)),'.fig'));
        file_name2

        %saveas(gcf,file_name);
        export_fig(file_name, '-transparent', '-png');
        %savefig(gcf, file_name2);


    end
    
    % stitched on autocorr-noise processing here %
    handles.distinct_plots = handles.spikes_data.nrAssigned(:,1);
    handles.noise_status3 = zeros(1, length(handles.spikes_data.nrAssigned(:,1)));
    
    for k = 1:length(handles.noise_status3)
        
        curr_times_spike_train = NaN(1,length(handles.spikes_data.assignedNegative));
        count = 1;
        for i = 1:length(handles.spikes_data.assignedNegative)
            
                if handles.spikes_data.assignedNegative(i) == handles.distinct_plots(k)
                    curr_times_spike_train(1,count) = handles.spikes_data.newTimestampsNegative(i);
                    count = count + 1;
                end
            
        end

        curr_times_spike_train = curr_times_spike_train(1,1:count-1);
        spike_train = convertToSpiketrain(curr_times_spike_train, 1);
        [~,~,~,Cxx] = psautospk(spike_train, 1);
    
        Cxx(1) = [];
        Cxx = Cxx(floor(length(Cxx)/2):length(Cxx));
        Cxx = Cxx(15:100);
        
        [~, ~, ~, p] = findpeaks(Cxx);
        thres = median(p) + 3*std(p);
        [~, identified] = findpeaks(Cxx, 'MinPeakProminence', thres);
        
        diff_c = diff(identified);
        median_c = median(diff_c);
        counter_c = 0;
        for i = 1:length(diff_c)
            if abs(diff_c(i) - median_c) < 2
                counter_c = counter_c + 1;
            end
        end
        
        if counter_c/length(diff_c) >= 0.7 && length(diff_c) > 1
            handles.noise_status3(k) = 1;
        end
    end
    
    tosave = handles.noise_status3;
    save('autocorr_noise.mat', 'tosave');
    
    
    cd(channel_path);
end



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
    end
end

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
end

function [imageData, alpha] = export_fig(varargin) %#ok<*STRCL1>
%EXPORT_FIG  Exports figures in a publication-quality format
%
% Examples:
%   imageData = export_fig
%   [imageData, alpha] = export_fig
%   export_fig filename
%   export_fig filename -format1 -format2
%   export_fig ... -nocrop
%   export_fig ... -c[<val>,<val>,<val>,<val>]
%   export_fig ... -transparent
%   export_fig ... -native
%   export_fig ... -m<val>
%   export_fig ... -r<val>
%   export_fig ... -a<val>
%   export_fig ... -q<val>
%   export_fig ... -p<val>
%   export_fig ... -d<gs_option>
%   export_fig ... -depsc
%   export_fig ... -<renderer>
%   export_fig ... -<colorspace>
%   export_fig ... -append
%   export_fig ... -bookmark
%   export_fig ... -clipboard
%   export_fig ... -update
%   export_fig ... -nofontswap
%   export_fig ... -font_space <char>
%   export_fig ... -linecaps
%   export_fig ... -noinvert
%   export_fig ... -preserve_size
%   export_fig ... -options <optionsStruct>
%   export_fig(..., handle)
%
% This function saves a figure or single axes to one or more vector and/or
% bitmap file formats, and/or outputs a rasterized version to the workspace,
% with the following properties:
%   - Figure/axes reproduced as it appears on screen
%   - Cropped borders (optional)
%   - Embedded fonts (vector formats)
%   - Improved line and grid line styles
%   - Anti-aliased graphics (bitmap formats)
%   - Render images at native resolution (optional for bitmap formats)
%   - Transparent background supported (pdf, eps, png, tif)
%   - Semi-transparent patch objects supported (png, tif)
%   - RGB, CMYK or grayscale output (CMYK only with pdf, eps, tif)
%   - Variable image compression, including lossless (pdf, eps, jpg)
%   - Optional rounded line-caps (pdf, eps)
%   - Optionally append to file (pdf, tif)
%   - Vector formats: pdf, eps, svg
%   - Bitmap formats: png, tif, jpg, bmp, export to workspace
%
% This function is especially suited to exporting figures for use in
% publications and presentations, because of the high quality and
% portability of media produced.
%
% Note that the background color and figure dimensions are reproduced
% (the latter approximately, and ignoring cropping & magnification) in the
% output file. For transparent background (and semi-transparent patch
% objects), use the -transparent option or set the figure 'Color' property
% to 'none'. To make axes transparent set the axes 'Color' property to
% 'none'. PDF, EPS, TIF & PNG are the only formats that support a transparent
% background; only TIF & PNG formats support transparency of patch objects.
%
% The choice of renderer (opengl/zbuffer/painters) has a large impact on the
% output quality. The default value (opengl for bitmaps, painters for vector
% formats) generally gives good results, but if you aren't satisfied
% then try another renderer.  Notes:
%   1) For vector formats (EPS,PDF), only painters generates vector graphics
%   2) For bitmap formats, only opengl correctly renders transparent patches
%   3) For bitmap formats, only painters correctly scales line dash and dot
%      lengths when magnifying or anti-aliasing
%   4) Fonts may be substitued with Courier when using painters
%
% When exporting to vector format (PDF & EPS) and bitmap format using the
% painters renderer, this function requires that ghostscript is installed
% on your system. You can download this from:
%   http://www.ghostscript.com
% When exporting to EPS it additionally requires pdftops, from the Xpdf
% suite of functions. You can download this from: http://xpdfreader.com
%
% SVG output uses the fig2svg (https://github.com/kupiqu/fig2svg) or plot2svg
% (https://github.com/jschwizer99/plot2svg) utilities, or Matlab's built-in
% SVG export if neither of these utilities are available on Matlab's path.
% Note: cropping/padding are not supported in export_fig's SVG output.
%
% Inputs:
%   filename - string containing the name (optionally including full or
%             relative path) of the file the figure is to be saved as. If
%             a path is not specified, the figure is saved in the current
%             directory. If no name and no output arguments are specified,
%             the default name, 'export_fig_out', is used. If neither a
%             file extension nor a format are specified, a ".png" is added
%             and the figure saved in that format.
%   -<format> - string(s) containing the output file extension(s). Options:
%             '-pdf', '-eps', '-svg', '-png', '-tif', '-jpg' and '-bmp'.
%             Multiple formats can be specified, without restriction.
%             For example: export_fig('-jpg', '-pdf', '-png', ...)
%             Either '-tif','-tiff' can be specified, and either '-jpg','-jpeg'.
%   -nocrop - option indicating that empty margins should not be cropped.
%   -c[<val>,<val>,<val>,<val>] - option indicating crop amounts. Must be
%             a 4-element vector of numeric values: [top,right,bottom,left]
%             where NaN/Inf indicate auto-cropping, 0 means no cropping,
%             and any other value mean cropping in pixel amounts.
%   -transparent - option indicating that the figure background is to be made
%             transparent (PNG,PDF,TIF,EPS formats only). Implies -noinvert.
%   -m<val> - option where val indicates the factor to magnify the
%             on-screen figure pixel dimensions by when generating bitmap
%             outputs (does not affect vector formats). Default: '-m1'.
%   -r<val> - option val indicates the resolution (in pixels per inch) to
%             export bitmap and vector outputs at, keeping the dimensions
%             of the on-screen figure. Default: '-r864' (for vector output
%             only). Note that the -m option overides the -r option for
%             bitmap outputs only.
%   -native - option indicating that the output resolution (when outputting
%             a bitmap format) should be such that the vertical resolution
%             of the first suitable image found in the figure is at the
%             native resolution of that image. To specify a particular
%             image to use, give it the tag 'export_fig_native'. Notes:
%             This overrides any value set with the -m and -r options. It
%             also assumes that the image is displayed front-to-parallel
%             with the screen. The output resolution is approximate and
%             should not be relied upon. Anti-aliasing can have adverse
%             effects on image quality (disable with the -a1 option).
%   -a1, -a2, -a3, -a4 - option indicating the amount of anti-aliasing to
%             use for bitmap outputs. '-a1' means no anti-aliasing;
%             '-a4' is the maximum amount (default).
%   -<renderer> - option to force a particular renderer (painters, opengl or
%             zbuffer). Default value: opengl for bitmap formats or
%             figures with patches and/or transparent annotations;
%             painters for vector formats without patches/transparencies.
%   -<colorspace> - option indicating which colorspace color figures should
%             be saved in: RGB (default), CMYK or gray. Usage example: '-gray'.
%             Note: CMYK is only supported in PDF, EPS and TIF formats.
%   -q<val> - option to vary bitmap image quality (PDF, EPS, JPG formats only).
%             A larger val, in the range 0-100, produces higher quality and
%             lower compression. val > 100 results in lossless compression.
%             Default: '-q95' for JPG, ghostscript prepress default for PDF,EPS.
%             Note: lossless compression can sometimes give a smaller file size
%             than the default lossy compression, depending on the image type.
%   -p<val> - option to pad a border of width val to exported files, where
%             val is either a relative size with respect to cropped image
%             size (i.e. p=0.01 adds a 1% border). For EPS & PDF formats,
%             val can also be integer in units of 1/72" points (abs(val)>1).
%             val can be positive (padding) or negative (extra cropping).
%             If used, the -nocrop flag will be ignored, i.e. the image will
%             always be cropped and then padded. Default: 0 (i.e. no padding).
%   -append - option indicating that if the file already exists the figure is to
%             be appended as a new page, instead of being overwritten (default).
%             PDF & TIF output formats only.
%   -bookmark - option to indicate that a bookmark with the name of the
%             figure is to be created in the output file (PDF format only).
%   -clipboard - option to save output as an image on the system clipboard.
%             Note: background transparency is not preserved in clipboard
%   -d<gs_option> - option to indicate a ghostscript setting. For example,
%             -dMaxBitmap=0 or -dNoOutputFonts (Ghostscript 9.15+).
%   -depsc -  option to use EPS level-3 rather than the default level-2 print
%             device. This solves some bugs with Matlab's default -depsc2 device
%             such as discolored subplot lines on images (vector formats only).
%   -update - option to download and install the latest version of export_fig
%   -nofontswap - option to avoid font swapping. Font swapping is automatically
%             done in vector formats (only): 11 standard Matlab fonts are
%             replaced by the original figure fonts. This option prevents this.
%   -font_space <char> - option to set a spacer character for font-names that
%             contain spaces, used by EPS/PDF. Default: ''
%   -linecaps - option to create rounded line-caps (vector formats only).
%   -noinvert - option to avoid setting figure's InvertHardcopy property to
%             'off' during output (this solves some problems of empty outputs).
%   -preserve_size - option to preserve the figure's PaperSize property in output
%             file (PDF/EPS formats only; default is to not preserve it).
%   -options <optionsStruct> - format-specific parameters as defined in Matlab's
%             documentation of the imwrite function, contained in a struct under
%             the format name. For example to specify the JPG Comment parameter,
%             pass a struct such as this: options.JPG.Comment='abc'. Similarly,
%             options.PNG.BitDepth=4. Valid only for PNG,TIF,JPG output formats.
%   handle -  The handle of the figure, axes or uipanels (can be an array of
%             handles, but the objects must be in the same figure) which is
%             to be saved. Default: gcf (handle of current figure).
%
% Outputs:
%   imageData - MxNxC uint8 image array of the exported image.
%   alpha     - MxN single array of alphamatte values in the range [0,1],
%               for the case when the background is transparent.
%
%   Some helpful examples and tips can be found at:
%      https://github.com/altmany/export_fig
%
%   See also PRINT, SAVEAS, ScreenCapture (on the Matlab File Exchange)

%{
% Copyright (C) Oliver Woodford 2008-2014, Yair Altman 2015-

% The idea of using ghostscript is inspired by Peder Axensten's SAVEFIG
% (fex id: 10889) which is itself inspired by EPS2PDF (fex id: 5782).
% The idea for using pdftops came from the MATLAB newsgroup (id: 168171).
% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928).
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)
% The idea of anti-aliasing bitmaps came from Anders Brun's MYAA (fex id:
% 20979).
% The idea of appending figures in pdfs came from Matt C in comments on the
% FEX (id: 23629)

% Thanks to Roland Martin for pointing out the colour MATLAB
% bug/feature with colorbar axes and transparent backgrounds.
% Thanks also to Andrew Matthews for describing a bug to do with the figure
% size changing in -nodisplay mode. I couldn't reproduce it, but included a
% fix anyway.
% Thanks to Tammy Threadgill for reporting a bug where an axes is not
% isolated from gui objects.
%}
%{
% 23/02/12: Ensure that axes limits don't change during printing
% 14/03/12: Fix bug in fixing the axes limits (thanks to Tobias Lamour for reporting it).
% 02/05/12: Incorporate patch of Petr Nechaev (many thanks), enabling bookmarking of figures in pdf files.
% 09/05/12: Incorporate patch of Arcelia Arrieta (many thanks), to keep tick marks fixed.
% 12/12/12: Add support for isolating uipanels. Thanks to michael for suggesting it.
% 25/09/13: Add support for changing resolution in vector formats. Thanks to Jan Jaap Meijer for suggesting it.
% 07/05/14: Add support for '~' at start of path. Thanks to Sally Warner for suggesting it.
% 24/02/15: Fix Matlab R2014b bug (issue #34): plot markers are not displayed when ZLimMode='manual'
% 25/02/15: Fix issue #4 (using HG2 on R2014a and earlier)
% 25/02/15: Fix issue #21 (bold TeX axes labels/titles in R2014b)
% 26/02/15: If temp dir is not writable, use the user-specified folder for temporary EPS/PDF files (Javier Paredes)
% 27/02/15: Modified repository URL from github.com/ojwoodford to /altmany
%           Indented main function
%           Added top-level try-catch block to display useful workarounds
% 28/02/15: Enable users to specify optional ghostscript options (issue #36)
% 06/03/15: Improved image padding & cropping thanks to Oscar Hartogensis
% 26/03/15: Fixed issue #49 (bug with transparent grayscale images); fixed out-of-memory issue
% 26/03/15: Fixed issue #42: non-normalized annotations on HG1
% 26/03/15: Fixed issue #46: Ghostscript crash if figure units <> pixels
% 27/03/15: Fixed issue #39: bad export of transparent annotations/patches
% 28/03/15: Fixed issue #50: error on some Matlab versions with the fix for issue #42
% 29/03/15: Fixed issue #33: bugs in Matlab's print() function with -cmyk
% 29/03/15: Improved processing of input args (accept space between param name & value, related to issue #51)
% 30/03/15: When exporting *.fig files, then saveas *.fig if figure is open, otherwise export the specified fig file
% 30/03/15: Fixed edge case bug introduced yesterday (commit #ae1755bd2e11dc4e99b95a7681f6e211b3fa9358)
% 09/04/15: Consolidated header comment sections; initialize output vars only if requested (nargout>0)
% 14/04/15: Workaround for issue #45: lines in image subplots are exported in invalid color
% 15/04/15: Fixed edge-case in parsing input parameters; fixed help section to show the -depsc option (issue #45)
% 21/04/15: Bug fix: Ghostscript croaks on % chars in output PDF file (reported by Sven on FEX page, 15-Jul-2014)
% 22/04/15: Bug fix: Pdftops croaks on relative paths (reported by Tintin Milou on FEX page, 19-Jan-2015)
% 04/05/15: Merged fix #63 (Kevin Mattheus Moerman): prevent tick-label changes during export
% 07/05/15: Partial fix for issue #65: PDF export used painters rather than opengl renderer (thanks Nguyenr)
% 08/05/15: Fixed issue #65: bad PDF append since commit #e9f3cdf 21/04/15 (thanks Robert Nguyen)
% 12/05/15: Fixed issue #67: exponent labels cropped in export, since fix #63 (04/05/15)
% 28/05/15: Fixed issue #69: set non-bold label font only if the string contains symbols (\beta etc.), followup to issue #21
% 29/05/15: Added informative error message in case user requested SVG output (issue #72)
% 09/06/15: Fixed issue #58: -transparent removed anti-aliasing when exporting to PNG
% 19/06/15: Added -update option to download and install the latest version of export_fig
% 07/07/15: Added -nofontswap option to avoid font-swapping in EPS/PDF
% 16/07/15: Fixed problem with anti-aliasing on old Matlab releases
% 11/09/15: Fixed issue #103: magnification must never become negative; also fixed reported error msg in parsing input params
% 26/09/15: Alert if trying to export transparent patches/areas to non-PNG outputs (issue #108)
% 04/10/15: Do not suggest workarounds for certain errors that have already been handled previously
% 01/11/15: Fixed issue #112: use same renderer in print2eps as export_fig (thanks to Jesús Pestana Puerta)
% 10/11/15: Custom GS installation webpage for MacOS. Thanks to Andy Hueni via FEX
% 19/11/15: Fixed clipboard export in R2015b (thanks to Dan K via FEX)
% 21/02/16: Added -c option for indicating specific crop amounts (idea by Cedric Noordam on FEX)
% 08/05/16: Added message about possible error reason when groot.Units~=pixels (issue #149)
% 17/05/16: Fixed case of image YData containing more than 2 elements (issue #151)
% 08/08/16: Enabled exporting transparency to TIF, in addition to PNG/PDF (issue #168)
% 11/12/16: Added alert in case of error creating output PDF/EPS file (issue #179)
% 13/12/16: Minor fix to the commit for issue #179 from 2 days ago
% 22/03/17: Fixed issue #187: only set manual ticks when no exponent is present
% 09/04/17: Added -linecaps option (idea by Baron Finer, issue #192)
% 15/09/17: Fixed issue #205: incorrect tick-labels when Ticks number don't match the TickLabels number
% 15/09/17: Fixed issue #210: initialize alpha map to ones instead of zeros when -transparent is not used
% 18/09/17: Added -font_space option to replace font-name spaces in EPS/PDF (workaround for issue #194)
% 18/09/17: Added -noinvert option to solve some export problems with some graphic cards (workaround for issue #197)
% 08/11/17: Fixed issue #220: axes exponent is removed in HG1 when TickMode is 'manual' (internal Matlab bug)
% 08/11/17: Fixed issue #221: alert if the requested folder does not exist
% 19/11/17: Workaround for issue #207: alert when trying to use transparent bgcolor with -opengl
% 29/11/17: Workaround for issue #206: warn if exporting PDF/EPS for a figure that contains an image
% 11/12/17: Fixed issue #230: use OpenGL renderer when exported image contains transparency (also see issue #206)
% 30/01/18: Updated SVG message to point to https://github.com/kupiqu/plot2svg and display user-selected filename if available
% 27/02/18: Fixed issue #236: axes exponent cropped from output if on right-hand axes
% 29/05/18: Fixed issue #245: process "string" inputs just like 'char' inputs
% 13/08/18: Fixed issue #249: correct black axes color to off-black to avoid extra cropping with -transparent
% 27/08/18: Added a possible file-open reason in EPS/PDF write-error message (suggested by "craq" on FEX page)
% 22/09/18: Xpdf website changed to xpdfreader.com
% 23/09/18: Fixed issue #243: only set non-bold font (workaround for issue #69) in R2015b or earlier; warn if changing font
% 23/09/18: Workaround for issue #241: don't use -r864 in EPS/PDF outputs when -native is requested (solves black lines problem)
% 18/11/18: Issue #261: Added informative alert when trying to export a uifigure (which is not currently supported)
% 13/12/18: Issue #261: Fixed last commit for cases of specifying axes/panel handle as input, rather than a figure handle
% 13/01/19: Issue #72: Added basic SVG output support
% 04/02/19: Workaround for issues #207 and #267: -transparent implies -noinvert
% 08/03/19: Issue #269: Added ability to specify format-specific options for PNG,TIF,JPG outputs; fixed help section
% 21/03/19: Fixed the workaround for issues #207 and #267 from 4/2/19 (-transparent now does *NOT* imply -noinvert; -transparent output should now be ok in all formats)
% 12/06/19: Issue #277: Enabled preservation of figure's PaperSize in output PDF/EPS file
%}

    if nargout
        [imageData, alpha] = deal([]);
    end
    hadError = false;
    displaySuggestedWorkarounds = true;

    % Ensure the figure is rendered correctly _now_ so that properties like axes limits are up-to-date
    drawnow;
    pause(0.05);  % this solves timing issues with Java Swing's EDT (http://undocumentedmatlab.com/blog/solving-a-matlab-hang-problem)

    % Parse the input arguments
    fig = get(0, 'CurrentFigure');
    [fig, options] = parse_args(nargout, fig, varargin{:});

    % Ensure that we have a figure handle
    if isequal(fig,-1)
        return;  % silent bail-out
    elseif isempty(fig)
        error('No figure found');
    else
        oldWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        try jf = get(handle(ancestor(fig,'figure')),'JavaFrame'); catch, jf=1; end
        warning(oldWarn);
        if isempty(jf)
            error('Figures created using the uifigure command or App Designer are not supported by export_fig. See <a href="https://github.com/altmany/export_fig/issues/261">issue #261</a> for details.');
        end
    end

    % Isolate the subplot, if it is one
    cls = all(ismember(get(fig, 'Type'), {'axes', 'uipanel'}));
    if cls
        % Given handles of one or more axes, so isolate them from the rest
        fig = isolate_axes(fig);
    else
        % Check we have a figure
        if ~isequal(get(fig, 'Type'), 'figure')
            error('Handle must be that of a figure, axes or uipanel');
        end
        % Get the old InvertHardcopy mode
        old_mode = get(fig, 'InvertHardcopy');
    end

    % Hack the font units where necessary (due to a font rendering bug in print?).
    % This may not work perfectly in all cases.
    % Also it can change the figure layout if reverted, so use a copy.
    magnify = options.magnify * options.aa_factor;
    if isbitmap(options) && magnify ~= 1
        fontu = findall(fig, 'FontUnits', 'normalized');
        if ~isempty(fontu)
            % Some normalized font units found
            if ~cls
                fig = copyfig(fig);
                set(fig, 'Visible', 'off');
                fontu = findall(fig, 'FontUnits', 'normalized');
                cls = true;
            end
            set(fontu, 'FontUnits', 'points');
        end
    end

    try
        % MATLAB "feature": axes limits and tick marks can change when printing
        Hlims = findall(fig, 'Type', 'axes');
        if ~cls
            % Record the old axes limit and tick modes
            Xlims = make_cell(get(Hlims, 'XLimMode'));
            Ylims = make_cell(get(Hlims, 'YLimMode'));
            Zlims = make_cell(get(Hlims, 'ZLimMode'));
            Xtick = make_cell(get(Hlims, 'XTickMode'));
            Ytick = make_cell(get(Hlims, 'YTickMode'));
            Ztick = make_cell(get(Hlims, 'ZTickMode'));
            Xlabel = make_cell(get(Hlims, 'XTickLabelMode')); 
            Ylabel = make_cell(get(Hlims, 'YTickLabelMode')); 
            Zlabel = make_cell(get(Hlims, 'ZTickLabelMode')); 
        end

        % Set all axes limit and tick modes to manual, so the limits and ticks can't change
        % Fix Matlab R2014b bug (issue #34): plot markers are not displayed when ZLimMode='manual'
        set(Hlims, 'XLimMode', 'manual', 'YLimMode', 'manual');
        set_tick_mode(Hlims, 'X');
        set_tick_mode(Hlims, 'Y');
        if ~using_hg2(fig)
            set(Hlims,'ZLimMode', 'manual');
            set_tick_mode(Hlims, 'Z');
        end
    catch
        % ignore - fix issue #4 (using HG2 on R2014a and earlier)
    end

    % Fix issue #21 (bold TeX axes labels/titles in R2014b when exporting to EPS/PDF)
    try
        if using_hg2(fig) && isvector(options)
            % Set the FontWeight of axes labels/titles to 'normal'
            % Fix issue #69: set non-bold font only if the string contains symbols (\beta etc.)
            % Issue #243: only set non-bold font (workaround for issue #69) in R2015b or earlier
            try isPreR2016a = verLessThan('matlab','8.7'); catch, isPreR2016a = true; end
            if isPreR2016a
                texLabels = findall(fig, 'type','text', 'FontWeight','bold');
                symbolIdx = ~cellfun('isempty',strfind({texLabels.String},'\'));
                if ~isempty(symbolIdx)
                    set(texLabels(symbolIdx), 'FontWeight','normal');
                    warning('export_fig:BoldTexLabels', 'Bold labels with Tex symbols converted into non-bold in export_fig (fix for issue #69)');
                end
            end
        end
    catch
        % ignore
    end

    % Fix issue #42: non-normalized annotations on HG1 (internal Matlab bug)
    annotationHandles = [];
    try
        if ~using_hg2(fig)
            annotationHandles = findall(fig,'Type','hggroup','-and','-property','Units','-and','-not','Units','norm');
            try  % suggested by Jesús Pestana Puerta (jespestana) 30/9/2015
                originalUnits = get(annotationHandles,'Units');
                set(annotationHandles,'Units','norm');
            catch
            end
        end
    catch
        % should never happen, but ignore in any case - issue #50
    end

    % Fix issue #46: Ghostscript crash if figure units <> pixels
    oldFigUnits = get(fig,'Units');
    set(fig,'Units','pixels');

    % Set to print exactly what is there
    if options.invert_hardcopy
        try set(fig, 'InvertHardcopy', 'off'); catch, end  % fail silently in uifigures
    end

    % Set the renderer
    switch options.renderer
        case 1
            renderer = '-opengl';
        case 2
            renderer = '-zbuffer';
        case 3
            renderer = '-painters';
        otherwise
            renderer = '-opengl'; % Default for bitmaps
    end

    hImages = findall(fig,'type','image');

    % Handle transparent patches
    hasTransparency = ~isempty(findall(fig,'-property','FaceAlpha','-and','-not','FaceAlpha',1));
    hasPatches      = ~isempty(findall(fig,'type','patch'));
    if hasTransparency
        % Alert if trying to export transparent patches/areas to non-supported outputs (issue #108)
        % http://www.mathworks.com/matlabcentral/answers/265265-can-export_fig-or-else-draw-vector-graphics-with-transparent-surfaces
        % TODO - use transparency when exporting to PDF by not passing via print2eps
        msg = 'export_fig currently supports transparent patches/areas only in PNG output. ';
        if options.pdf
            warning('export_fig:transparency', '%s\nTo export transparent patches/areas to PDF, use the print command:\n print(gcf, ''-dpdf'', ''%s.pdf'');', msg, options.name);
        elseif ~options.png && ~options.tif  % issue #168
            warning('export_fig:transparency', '%s\nTo export the transparency correctly, try using the ScreenCapture utility on the Matlab File Exchange: http://bit.ly/1QFrBip', msg);
        end
    elseif ~isempty(hImages)
        % Fix for issue #230: use OpenGL renderer when exported image contains transparency
        for idx = 1 : numel(hImages)
            cdata = get(hImages(idx),'CData');
            if any(isnan(cdata(:)))
                hasTransparency = true;
                break
            end
        end
    end

    try
        % Do the bitmap formats first
        if isbitmap(options)
            if abs(options.bb_padding) > 1
                displaySuggestedWorkarounds = false;
                error('For bitmap output (png,jpg,tif,bmp) the padding value (-p) must be between -1<p<1')
            end
            % Get the background colour
            if options.transparent && (options.png || options.alpha)
                % Get out an alpha channel
                % MATLAB "feature": black colorbar axes can change to white and vice versa!
                hCB = findall(fig, 'Type','axes', 'Tag','Colorbar');
                if isempty(hCB)
                    yCol = [];
                    xCol = [];
                else
                    yCol = get(hCB, 'YColor');
                    xCol = get(hCB, 'XColor');
                    if iscell(yCol)
                        yCol = cell2mat(yCol);
                        xCol = cell2mat(xCol);
                    end
                    yCol = sum(yCol, 2);
                    xCol = sum(xCol, 2);
                end
                % MATLAB "feature": apparently figure size can change when changing
                % colour in -nodisplay mode
                pos = get(fig, 'Position');
                % Set the background colour to black, and set size in case it was
                % changed internally
                tcol = get(fig, 'Color');
                set(fig, 'Color', 'k', 'Position', pos);
                % Correct the colorbar axes colours
                set(hCB(yCol==0), 'YColor', [0 0 0]);
                set(hCB(xCol==0), 'XColor', [0 0 0]);
                % Correct black axes color to off-black (issue #249)
                hAxes = findall(fig, 'Type','axes');
                hXs = fixBlackAxle(hAxes, 'XColor');
                hYs = fixBlackAxle(hAxes, 'YColor');
                hZs = fixBlackAxle(hAxes, 'ZColor');

                % The following code might cause out-of-memory errors
                try
                    % Print large version to array
                    B = print2array(fig, magnify, renderer);
                    % Downscale the image
                    B = downsize(single(B), options.aa_factor);
                catch
                    % This is more conservative in memory, but kills transparency (issue #58)
                    B = single(print2array(fig, magnify/options.aa_factor, renderer));
                end

                % Set background to white (and set size)
                set(fig, 'Color', 'w', 'Position', pos);
                % Correct the colorbar axes colours
                set(hCB(yCol==3), 'YColor', [1 1 1]);
                set(hCB(xCol==3), 'XColor', [1 1 1]);
                % Revert the black axes colors
                set(hXs, 'XColor', [0,0,0]);
                set(hYs, 'YColor', [0,0,0]);
                set(hZs, 'ZColor', [0,0,0]);

                % The following code might cause out-of-memory errors
                try
                    % Print large version to array
                    A = print2array(fig, magnify, renderer);
                    % Downscale the image
                    A = downsize(single(A), options.aa_factor);
                catch
                    % This is more conservative in memory, but kills transparency (issue #58)
                    A = single(print2array(fig, magnify/options.aa_factor, renderer));
                end

                % Set the background colour (and size) back to normal
                set(fig, 'Color', tcol, 'Position', pos);
                % Compute the alpha map
                alpha = round(sum(B - A, 3)) / (255 * 3) + 1;
                A = alpha;
                A(A==0) = 1;
                A = B ./ A(:,:,[1 1 1]);
                clear B
                % Convert to greyscale
                if options.colourspace == 2
                    A = rgb2grey(A);
                end
                A = uint8(A);
                % Crop the background
                if options.crop
                    %[alpha, v] = crop_borders(alpha, 0, 1, options.crop_amounts);
                    %A = A(v(1):v(2),v(3):v(4),:);
                    [alpha, vA, vB] = crop_borders(alpha, 0, options.bb_padding, options.crop_amounts);
                    if ~any(isnan(vB)) % positive padding
                        B = repmat(uint8(zeros(1,1,size(A,3))),size(alpha));
                        B(vB(1):vB(2), vB(3):vB(4), :) = A(vA(1):vA(2), vA(3):vA(4), :); % ADDED BY OH
                        A = B;
                    else  % negative padding
                        A = A(vA(1):vA(2), vA(3):vA(4), :);
                    end
                end
                if options.png
                    % Compute the resolution
                    res = options.magnify * get(0, 'ScreenPixelsPerInch') / 25.4e-3;
                    % Save the png
                    [format_options, bitDepth] = getFormatOptions(options, 'png');  %Issue #269
                    if ~isempty(bitDepth) && bitDepth < 16 && size(A,3) == 3
                        % BitDepth specification requires using a color-map
                        [A, map] = rgb2ind(A, 256);
                        imwrite(A, map, [options.name '.png'], 'Alpha',double(alpha), 'ResolutionUnit','meter', 'XResolution',res, 'YResolution',res, format_options{:});
                    else
                        imwrite(A, [options.name '.png'], 'Alpha',double(alpha), 'ResolutionUnit','meter', 'XResolution',res, 'YResolution',res, format_options{:});
                    end
                    % Clear the png bit
                    options.png = false;
                end
                % Return only one channel for greyscale
                if isbitmap(options)
                    A = check_greyscale(A);
                end
                if options.alpha
                    % Store the image
                    imageData = A;
                    % Clear the alpha bit
                    options.alpha = false;
                end
                % Get the non-alpha image
                if isbitmap(options)
                    alph = alpha(:,:,ones(1, size(A, 3)));
                    A = uint8(single(A) .* alph + 255 * (1 - alph));
                    clear alph
                end
                if options.im
                    % Store the new image
                    imageData = A;
                end
            else
                % Print large version to array
                if options.transparent
                    % MATLAB "feature": apparently figure size can change when changing
                    % colour in -nodisplay mode
                    pos = get(fig, 'Position');
                    tcol = get(fig, 'Color');
                    set(fig, 'Color', 'w', 'Position', pos);
                    A = print2array(fig, magnify, renderer);
                    set(fig, 'Color', tcol, 'Position', pos);
                    tcol = 255;
                else
                    [A, tcol] = print2array(fig, magnify, renderer);
                end
                % Crop the background
                if options.crop
                    A = crop_borders(A, tcol, options.bb_padding, options.crop_amounts);
                end
                % Downscale the image
                A = downsize(A, options.aa_factor);
                if options.colourspace == 2
                    % Convert to greyscale
                    A = rgb2grey(A);
                else
                    % Return only one channel for greyscale
                    A = check_greyscale(A);
                end
                % Outputs
                if options.im
                    imageData = A;
                end
                if options.alpha
                    imageData = A;
                    alpha = ones(size(A, 1), size(A, 2), 'single');
                end
            end
            % Save the images
            if options.png
                res = options.magnify * get(0, 'ScreenPixelsPerInch') / 25.4e-3;
                [format_options, bitDepth] = getFormatOptions(options, 'png');  %Issue #269
                if ~isempty(bitDepth) && bitDepth < 16 && size(A,3) == 3
                    % BitDepth specification requires using a color-map
                    [A, map] = rgb2ind(A, 256);
                    imwrite(A, map, [options.name '.png'], 'ResolutionUnit','meter', 'XResolution',res, 'YResolution',res, format_options{:});
                else
                    imwrite(A, [options.name '.png'], 'ResolutionUnit','meter', 'XResolution',res, 'YResolution',res, format_options{:});
                end
            end
            if options.bmp
                imwrite(A, [options.name '.bmp']);
            end
            % Save jpeg with given quality
            if options.jpg
                quality = options.quality;
                if isempty(quality)
                    quality = 95;
                end
                format_options = getFormatOptions(options, 'jpg');  %Issue #269
                if quality > 100
                    imwrite(A, [options.name '.jpg'], 'Mode','lossless', format_options{:});
                else
                    imwrite(A, [options.name '.jpg'], 'Quality',quality, format_options{:});
                end
            end
            % Save tif images in cmyk if wanted (and possible)
            if options.tif
                if options.colourspace == 1 && size(A, 3) == 3
                    A = double(255 - A);
                    K = min(A, [], 3);
                    K_ = 255 ./ max(255 - K, 1);
                    C = (A(:,:,1) - K) .* K_;
                    M = (A(:,:,2) - K) .* K_;
                    Y = (A(:,:,3) - K) .* K_;
                    A = uint8(cat(3, C, M, Y, K));
                    clear C M Y K K_
                end
                append_mode = {'overwrite', 'append'};
                format_options = getFormatOptions(options, 'tif');  %Issue #269
                imwrite(A, [options.name '.tif'], 'Resolution',options.magnify*get(0,'ScreenPixelsPerInch'), 'WriteMode',append_mode{options.append+1}, format_options{:});
            end
        end

        % Now do the vector formats
        if isvector(options)
            % Set the default renderer to painters
            if ~options.renderer
                if hasTransparency || hasPatches
                    % This is *MUCH* slower, but more accurate for patches and transparent annotations (issue #39)
                    renderer = '-opengl';
                else
                    renderer = '-painters';
                end
            end
            options.rendererStr = renderer;  % fix for issue #112
            % Generate some filenames
            tmp_nam = [tempname '.eps'];
            try
                % Ensure that the temp dir is writable (Javier Paredes 30/1/15)
                fid = fopen(tmp_nam,'w');
                fwrite(fid,1);
                fclose(fid);
                delete(tmp_nam);
                isTempDirOk = true;
            catch
                % Temp dir is not writable, so use the user-specified folder
                [dummy,fname,fext] = fileparts(tmp_nam); %#ok<ASGLU>
                fpath = fileparts(options.name);
                tmp_nam = fullfile(fpath,[fname fext]);
                isTempDirOk = false;
            end
            if isTempDirOk
                pdf_nam_tmp = [tempname '.pdf'];
            else
                pdf_nam_tmp = fullfile(fpath,[fname '.pdf']);
            end
            if options.pdf
                pdf_nam = [options.name '.pdf'];
                try copyfile(pdf_nam, pdf_nam_tmp, 'f'); catch, end  % fix for issue #65
            else
                pdf_nam = pdf_nam_tmp;
            end
            % Generate the options for print
            printArgs = {renderer};
            if ~isempty(options.resolution)  % issue #241
                printArgs{end+1} = sprintf('-r%d', options.resolution);
            end
            if options.colourspace == 1  % CMYK
                % Issue #33: due to internal bugs in Matlab's print() function, we can't use its -cmyk option
                %printArgs{end+1} = '-cmyk';
            end
            if ~options.crop
                % Issue #56: due to internal bugs in Matlab's print() function, we can't use its internal cropping mechanism,
                % therefore we always use '-loose' (in print2eps.m) and do our own cropping (in crop_borders)
                %printArgs{end+1} = '-loose';
            end
            if any(strcmpi(varargin,'-depsc'))
                % Issue #45: lines in image subplots are exported in invalid color.
                % The workaround is to use the -depsc parameter instead of the default -depsc2
                printArgs{end+1} = '-depsc';
            end
            try
                % Remove background if requested (issue #207)
                originalBgColor = get(fig, 'Color');
                [hXs, hYs, hZs] = deal([]);
                if options.transparent %&& ~isequal(get(fig, 'Color'), 'none')
                    if options.renderer == 1  % OpenGL
                        warning('export_fig:openglTransparentBG', '-opengl sometimes fails to produce transparent backgrounds; in such a case, try to use -painters instead');
                    end

                    % Fix for issue #207, #267 (corrected)
                    set(fig,'Color','none');

                    % Correct black axes color to off-black (issue #249)
                    hAxes = findall(fig, 'Type','axes');
                    hXs = fixBlackAxle(hAxes, 'XColor');
                    hYs = fixBlackAxle(hAxes, 'YColor');
                    hZs = fixBlackAxle(hAxes, 'ZColor');
                end
                % Generate an eps
                print2eps(tmp_nam, fig, options, printArgs{:});
                % {
                % Remove the background, if desired
                if options.transparent %&& ~isequal(get(fig, 'Color'), 'none')
                    eps_remove_background(tmp_nam, 1 + using_hg2(fig));

                    % Revert the black axes colors
                    set(hXs, 'XColor', [0,0,0]);
                    set(hYs, 'YColor', [0,0,0]);
                    set(hZs, 'ZColor', [0,0,0]);
                end
                %}
                % Restore the figure's previous background color (if modified)
                try set(fig,'Color',originalBgColor); drawnow; catch, end
                % Fix colorspace to CMYK, if requested (workaround for issue #33)
                if options.colourspace == 1  % CMYK
                    % Issue #33: due to internal bugs in Matlab's print() function, we can't use its -cmyk option
                    change_rgb_to_cmyk(tmp_nam);
                end
                % Add a bookmark to the PDF if desired
                if options.bookmark
                    fig_nam = get(fig, 'Name');
                    if isempty(fig_nam)
                        warning('export_fig:EmptyBookmark', 'Bookmark requested for figure with no name. Bookmark will be empty.');
                    end
                    add_bookmark(tmp_nam, fig_nam);
                end
                % Generate a pdf
                eps2pdf(tmp_nam, pdf_nam_tmp, 1, options.append, options.colourspace==2, options.quality, options.gs_options);
                % Ghostscript croaks on % chars in the output PDF file, so use tempname and then rename the file
                try
                    % Rename the file (except if it is already the same)
                    % Abbie K's comment on the commit for issue #179 (#commitcomment-20173476)
                    if ~isequal(pdf_nam_tmp, pdf_nam)
                        movefile(pdf_nam_tmp, pdf_nam, 'f');
                    end
                catch
                    % Alert in case of error creating output PDF/EPS file (issue #179)
                    if exist(pdf_nam_tmp, 'file')
                        errMsg = ['Could not create ' pdf_nam ' - perhaps the folder does not exist, or you do not have write permissions, or the file is open in another application'];
                        error(errMsg);
                    else
                        error('Could not generate the intermediary EPS file.');
                    end
                end
            catch ex
                % Restore the figure's previous background color (in case it was not already restored)
                try set(fig,'Color',originalBgColor); drawnow; catch, end
                % Delete the eps
                delete(tmp_nam);
                % Rethrow the EPS/PDF-generation error
                rethrow(ex);
            end
            % Delete the eps
            delete(tmp_nam);
            if options.eps || options.linecaps
                try
                    % Generate an eps from the pdf
                    % since pdftops can't handle relative paths (e.g., '..\'), use a temp file
                    eps_nam_tmp = strrep(pdf_nam_tmp,'.pdf','.eps');
                    pdf2eps(pdf_nam, eps_nam_tmp);

                    % Issue #192: enable rounded line-caps
                    if options.linecaps
                        fstrm = read_write_entire_textfile(eps_nam_tmp);
                        fstrm = regexprep(fstrm, '[02] J', '1 J');
                        read_write_entire_textfile(eps_nam_tmp, fstrm);
                        if options.pdf
                            eps2pdf(eps_nam_tmp, pdf_nam, 1, options.append, options.colourspace==2, options.quality, options.gs_options);
                        end
                    end

                    if options.eps
                        movefile(eps_nam_tmp, [options.name '.eps'], 'f');
                    else  % if options.pdf
                        try delete(eps_nam_tmp); catch, end
                    end
                catch ex
                    if ~options.pdf
                        % Delete the pdf
                        delete(pdf_nam);
                    end
                    try delete(eps_nam_tmp); catch, end
                    rethrow(ex);
                end
                if ~options.pdf
                    % Delete the pdf
                    delete(pdf_nam);
                end
            end
            % Issue #206: warn if the figure contains an image
            if ~isempty(hImages) && strcmpi(renderer,'-opengl')  % see addendum to issue #206
                warnMsg = ['exporting images to PDF/EPS may result in blurry images on some viewers. ' ...
                           'If so, try to change viewer, or increase the image''s CData resolution, or use -opengl renderer, or export via the print function. ' ...
                           'See <a href="matlab:web(''https://github.com/altmany/export_fig/issues/206'',''-browser'');">issue #206</a> for details.'];
                warning('export_fig:pdf_eps:blurry_image', warnMsg);
            end
        end

        % SVG format
        if options.svg
            oldUnits = get(fig,'Units');
            filename = [options.name '.svg'];
            % Adapted from Dan Joshea's https://github.com/djoshea/matlab-save-figure :
            try %if verLessThan('matlab', '8.4')
                % Try using the fig2svg/plot2svg utilities
                try
                    fig2svg(filename, fig);  %https://github.com/kupiqu/fig2svg
                catch
                    plot2svg(filename, fig); %https://github.com/jschwizer99/plot2svg
                    warning('export_fig:SVG:plot2svg', 'export_fig used the plot2svg utility for SVG output. Better results may be gotten via the fig2svg utility (https://github.com/kupiqu/fig2svg).');
                end
            catch %else  % (neither fig2svg nor plot2svg are available)
                % Try Matlab's built-in svg engine (from Batik Graphics2D for java)
                try
                    set(fig,'Units','pixels');   % All data in the svg-file is saved in pixels
                    printArgs = {renderer};
                    if ~isempty(options.resolution)
                        printArgs{end+1} = sprintf('-r%d', options.resolution);
                    end
                    print(fig, '-dsvg', printArgs{:}, filename);
                    warning('export_fig:SVG:print', 'export_fig used Matlab''s built-in SVG output engine. Better results may be gotten via the fig2svg utility (https://github.com/kupiqu/fig2svg).');
                catch err  % built-in print() failed - maybe an old Matlab release (no -dsvg)
                    set(fig,'Units',oldUnits);
                    filename = strrep(filename,'export_fig_out','filename');
                    msg = ['SVG output is not supported for your figure: ' err.message '\n' ...
                        'Try one of the following alternatives:\n' ...
                        '  1. saveas(gcf,''' filename ''')\n' ...
                        '  2. fig2svg utility: https://github.com/kupiqu/fig2svg\n' ...  % Note: replaced defunct https://github.com/jschwizer99/plot2svg with up-to-date fork on https://github.com/kupiqu/fig2svg
                        '  3. export_fig to EPS/PDF, then convert to SVG using non-Matlab tools\n'];
                    error(sprintf(msg)); %#ok<SPERR>
                end
            end
            % SVG output was successful if we reached this point
            % Restore original figure units
            set(fig,'Units',oldUnits);
            % Add warning about unsupported export_fig options with SVG output
            if any(~isnan(options.crop_amounts)) || any(options.bb_padding)
                warning('export_fig:SVG:options', 'export_fig''s SVG output does not [currently] support cropping/padding.');
            end
        end

        % Revert the figure or close it (if requested)
        if cls || options.closeFig
            % Close the created figure
            close(fig);
        else
            % Reset the hardcopy mode
            try set(fig, 'InvertHardcopy', old_mode); catch, end  % fail silently in uifigures
            % Reset the axes limit and tick modes
            for a = 1:numel(Hlims)
                try
                    set(Hlims(a), 'XLimMode', Xlims{a}, 'YLimMode', Ylims{a}, 'ZLimMode', Zlims{a},... 
                                  'XTickMode', Xtick{a}, 'YTickMode', Ytick{a}, 'ZTickMode', Ztick{a},...
                                  'XTickLabelMode', Xlabel{a}, 'YTickLabelMode', Ylabel{a}, 'ZTickLabelMode', Zlabel{a}); 
                catch
                    % ignore - fix issue #4 (using HG2 on R2014a and earlier)
                end
            end
            % Revert the tex-labels font weights
            try set(texLabels, 'FontWeight','bold'); catch, end
            % Revert annotation units
            for handleIdx = 1 : numel(annotationHandles)
                try
                    oldUnits = originalUnits{handleIdx};
                catch
                    oldUnits = originalUnits;
                end
                try set(annotationHandles(handleIdx),'Units',oldUnits); catch, end
            end
            % Revert figure units
            set(fig,'Units',oldFigUnits);
        end

        % Output to clipboard (if requested)
        if options.clipboard
            % Delete the output file if unchanged from the default name ('export_fig_out.png')
            if strcmpi(options.name,'export_fig_out')
                try
                    fileInfo = dir('export_fig_out.png');
                    if ~isempty(fileInfo)
                        timediff = now - fileInfo.datenum;
                        ONE_SEC = 1/24/60/60;
                        if timediff < ONE_SEC
                            delete('export_fig_out.png');
                        end
                    end
                catch
                    % never mind...
                end
            end

            % Save the image in the system clipboard
            % credit: Jiro Doke's IMCLIPBOARD: http://www.mathworks.com/matlabcentral/fileexchange/28708-imclipboard
            try
                error(javachk('awt', 'export_fig -clipboard output'));
            catch
                warning('export_fig:clipboardJava', 'export_fig -clipboard output failed: requires Java to work');
                return;
            end
            try
                % Import necessary Java classes
                import java.awt.Toolkit
                import java.awt.image.BufferedImage
                import java.awt.datatransfer.DataFlavor

                % Get System Clipboard object (java.awt.Toolkit)
                cb = Toolkit.getDefaultToolkit.getSystemClipboard();

                % Add java class (ImageSelection) to the path
                if ~exist('ImageSelection', 'class')
                    javaaddpath(fileparts(which(mfilename)), '-end');
                end

                % Get image size
                ht = size(imageData, 1);
                wd = size(imageData, 2);

                % Convert to Blue-Green-Red format
                try
                    imageData2 = imageData(:, :, [3 2 1]);
                catch
                    % Probably gray-scaled image (2D, without the 3rd [RGB] dimension)
                    imageData2 = imageData(:, :, [1 1 1]);
                end

                % Convert to 3xWxH format
                imageData2 = permute(imageData2, [3, 2, 1]);

                % Append Alpha data (unused - transparency is not supported in clipboard copy)
                alphaData2 = uint8(permute(255*alpha,[3,2,1])); %=255*ones(1,wd,ht,'uint8')
                imageData2 = cat(1, imageData2, alphaData2);

                % Create image buffer
                imBuffer = BufferedImage(wd, ht, BufferedImage.TYPE_INT_RGB);
                imBuffer.setRGB(0, 0, wd, ht, typecast(imageData2(:), 'int32'), 0, wd);

                % Create ImageSelection object from the image buffer
                imSelection = ImageSelection(imBuffer);

                % Set clipboard content to the image
                cb.setContents(imSelection, []);
            catch
                warning('export_fig:clipboardFailed', 'export_fig -clipboard output failed: %s', lasterr); %#ok<LERR>
            end
        end

        % Don't output the data to console unless requested
        if ~nargout
            clear imageData alpha
        end
    catch err
        % Display possible workarounds before the error message
        if displaySuggestedWorkarounds && ~strcmpi(err.message,'export_fig error')
            if ~hadError,  fprintf(2, 'export_fig error. ');  end
            fprintf(2, 'Please ensure:\n');
            fprintf(2, '  that you are using the <a href="https://github.com/altmany/export_fig/archive/master.zip">latest version</a> of export_fig\n');
            if isvector(options)
                if ismac
                    fprintf(2, '  and that you have <a href="http://pages.uoregon.edu/koch">Ghostscript</a> installed\n');
                else
                    fprintf(2, '  and that you have <a href="http://www.ghostscript.com">Ghostscript</a> installed\n');
                end
            end
            try
                if options.eps
                    fprintf(2, '  and that you have <a href="http://xpdfreader.com/download.html">pdftops</a> installed\n');
                end
            catch
                % ignore - probably an error in parse_args
            end
            fprintf(2, '  and that you do not have <a href="matlab:which export_fig -all">multiple versions</a> of export_fig installed by mistake\n');
            fprintf(2, '  and that you did not made a mistake in the <a href="matlab:help export_fig">expected input arguments</a>\n');
            try
                % Alert per issue #149
                if ~strncmpi(get(0,'Units'),'pixel',5)
                    fprintf(2, '  or try to set groot''s Units property back to its default value of ''pixels'' (<a href="matlab:web(''https://github.com/altmany/export_fig/issues/149'',''-browser'');">details</a>)\n');
                end
            catch
                % ignore - maybe an old MAtlab release
            end
            fprintf(2, '\nIf the problem persists, then please <a href="https://github.com/altmany/export_fig/issues">report a new issue</a>.\n\n');
        end
        rethrow(err)
    end
end

function options = default_options()
    % Default options used by export_fig
    options = struct(...
        'name',            'export_fig_out', ...
        'crop',            true, ...
        'crop_amounts',    nan(1,4), ...  % auto-crop all 4 image sides
        'transparent',     false, ...
        'renderer',        0, ...         % 0: default, 1: OpenGL, 2: ZBuffer, 3: Painters
        'pdf',             false, ...
        'eps',             false, ...
        'svg',             false, ...
        'png',             false, ...
        'tif',             false, ...
        'jpg',             false, ...
        'bmp',             false, ...
        'clipboard',       false, ...
        'colourspace',     0, ...         % 0: RGB/gray, 1: CMYK, 2: gray
        'append',          false, ...
        'im',              false, ...
        'alpha',           false, ...
        'aa_factor',       0, ...
        'bb_padding',      0, ...
        'magnify',         [], ...
        'resolution',      [], ...
        'bookmark',        false, ...
        'closeFig',        false, ...
        'quality',         [], ...
        'update',          false, ...
        'fontswap',        true, ...
        'font_space',      '', ...
        'linecaps',        false, ...
        'invert_hardcopy', true, ...
        'format_options',  struct, ...
        'preserve_size',   false, ...
        'gs_options',      {{}});
end

function [fig, options] = parse_args(nout, fig, varargin)
    % Parse the input arguments

    % Convert strings => chars
    varargin = cellfun(@str2char,varargin,'un',false);

    % Set the defaults
    native = false; % Set resolution to native of an image
    options = default_options();
    options.im =    (nout == 1);  % user requested imageData output
    options.alpha = (nout == 2);  % user requested alpha output

    % Go through the other arguments
    skipNext = false;
    for a = 1:nargin-2
        if skipNext
            skipNext = false;
            continue;
        end
        if all(ishandle(varargin{a}))
            fig = varargin{a};
        elseif ischar(varargin{a}) && ~isempty(varargin{a})
            if varargin{a}(1) == '-'
                switch lower(varargin{a}(2:end))
                    case 'nocrop'
                        options.crop = false;
                        options.crop_amounts = [0,0,0,0];
                    case {'trans', 'transparent'}
                        options.transparent = true;
                    case 'opengl'
                        options.renderer = 1;
                    case 'zbuffer'
                        options.renderer = 2;
                    case 'painters'
                        options.renderer = 3;
                    case 'pdf'
                        options.pdf = true;
                    case 'eps'
                        options.eps = true;
                    case 'svg'
                        options.svg = true;
                    case 'png'
                        options.png = true;
                    case {'tif', 'tiff'}
                        options.tif = true;
                    case {'jpg', 'jpeg'}
                        options.jpg = true;
                    case 'bmp'
                        options.bmp = true;
                    case 'rgb'
                        options.colourspace = 0;
                    case 'cmyk'
                        options.colourspace = 1;
                    case {'gray', 'grey'}
                        options.colourspace = 2;
                    case {'a1', 'a2', 'a3', 'a4'}
                        options.aa_factor = str2double(varargin{a}(3));
                    case 'append'
                        options.append = true;
                    case 'bookmark'
                        options.bookmark = true;
                    case 'native'
                        native = true;
                    case 'clipboard'
                        options.clipboard = true;
                        options.im = true;
                        options.alpha = true;
                    case 'update'
                        % Download the latest version of export_fig into the export_fig folder
                        try
                            zipFileName = 'https://github.com/altmany/export_fig/archive/master.zip';
                            folderName = fileparts(which(mfilename('fullpath')));
                            targetFileName = fullfile(folderName, datestr(now,'yyyy-mm-dd.zip'));
                            urlwrite(zipFileName,targetFileName); %#ok<URLWR>
                        catch
                            error('Could not download %s into %s\n',zipFileName,targetFileName);
                        end

                        % Unzip the downloaded zip file in the export_fig folder
                        try
                            unzip(targetFileName,folderName);
                        catch
                            error('Could not unzip %s\n',targetFileName);
                        end
                    case 'nofontswap'
                        options.fontswap = false;
                    case 'font_space'
                        options.font_space = varargin{a+1};
                        skipNext = true;
                    case 'linecaps'
                        options.linecaps = true;
                    case 'noinvert'
                        options.invert_hardcopy = false;
                    case 'preserve_size'
                        options.preserve_size = true;
                    case 'options'
                        % Issue #269: format-specific options
                        inputOptions = varargin{a+1};
                        %options.format_options  = inputOptions;
                        if isempty(inputOptions), continue, end
                        formats = fieldnames(inputOptions(1));
                        for idx = 1 : numel(formats)
                            optionsStruct = inputOptions.(formats{idx});
                            %optionsCells = [fieldnames(optionsStruct) struct2cell(optionsStruct)]';
                            formatName = regexprep(lower(formats{idx}),{'tiff','jpeg'},{'tif','jpg'});
                            options.format_options.(formatName) = optionsStruct; %=optionsCells(:)';
                        end
                        skipNext = true;
                    otherwise
                        try
                            wasError = false;
                            if strcmpi(varargin{a}(1:2),'-d')
                                varargin{a}(2) = 'd';  % ensure lowercase 'd'
                                options.gs_options{end+1} = varargin{a};
                            elseif strcmpi(varargin{a}(1:2),'-c')
                                if numel(varargin{a})==2
                                    skipNext = true;
                                    vals = str2num(varargin{a+1}); %#ok<ST2NM>
                                else
                                    vals = str2num(varargin{a}(3:end)); %#ok<ST2NM>
                                end
                                if numel(vals)~=4
                                    wasError = true;
                                    error('option -c cannot be parsed: must be a 4-element numeric vector');
                                end
                                options.crop_amounts = vals;
                                options.crop = true;
                            else  % scalar parameter value
                                val = str2double(regexp(varargin{a}, '(?<=-(m|M|r|R|q|Q|p|P))-?\d*.?\d+', 'match'));
                                if isempty(val) || isnan(val)
                                    % Issue #51: improved processing of input args (accept space between param name & value)
                                    val = str2double(varargin{a+1});
                                    if isscalar(val) && ~isnan(val)
                                        skipNext = true;
                                    end
                                end
                                if ~isscalar(val) || isnan(val)
                                    wasError = true;
                                    error('option %s is not recognised or cannot be parsed', varargin{a});
                                end
                                switch lower(varargin{a}(2))
                                    case 'm'
                                        % Magnification may never be negative
                                        if val <= 0
                                            wasError = true;
                                            error('Bad magnification value: %g (must be positive)', val);
                                        end
                                        options.magnify = val;
                                    case 'r'
                                        options.resolution = val;
                                    case 'q'
                                        options.quality = max(val, 0);
                                    case 'p'
                                        options.bb_padding = val;
                                end
                            end
                        catch err
                            % We might have reached here by raising an intentional error
                            if wasError  % intentional raise
                                rethrow(err)
                            else  % unintentional
                                error(['Unrecognized export_fig input option: ''' varargin{a} '''']);
                            end
                        end
                end
            else
                [p, options.name, ext] = fileparts(varargin{a});
                if ~isempty(p)
                    % Issue #221: alert if the requested folder does not exist
                    if ~exist(p,'dir'),  error(['Folder ' p ' does not exist!']);  end
                    options.name = [p filesep options.name];
                end
                switch lower(ext)
                    case {'.tif', '.tiff'}
                        options.tif = true;
                    case {'.jpg', '.jpeg'}
                        options.jpg = true;
                    case '.png'
                        options.png = true;
                    case '.bmp'
                        options.bmp = true;
                    case '.eps'
                        options.eps = true;
                    case '.pdf'
                        options.pdf = true;
                    case '.fig'
                        % If no open figure, then load the specified .fig file and continue
                        if isempty(fig)
                            fig = openfig(varargin{a},'invisible');
                            varargin{a} = fig;
                            options.closeFig = true;
                        else
                            % save the current figure as the specified .fig file and exit
                            saveas(fig(1),varargin{a});
                            fig = -1;
                            return
                        end
                    case '.svg'
                        options.svg = true;
                    otherwise
                        options.name = varargin{a};
                end
            end
        end
    end

    % Quick bail-out if no figure found
    if isempty(fig),  return;  end

    % Do border padding with repsect to a cropped image
    if options.bb_padding
        options.crop = true;
    end

    % Set default anti-aliasing now we know the renderer
    if options.aa_factor == 0
        try isAA = strcmp(get(ancestor(fig, 'figure'), 'GraphicsSmoothing'), 'on'); catch, isAA = false; end
        options.aa_factor = 1 + 2 * (~(using_hg2(fig) && isAA) | (options.renderer == 3));
    end

    % Convert user dir '~' to full path
    if numel(options.name) > 2 && options.name(1) == '~' && (options.name(2) == '/' || options.name(2) == '\')
        options.name = fullfile(char(java.lang.System.getProperty('user.home')), options.name(2:end));
    end

    % Compute the magnification and resolution
    if isempty(options.magnify)
        if isempty(options.resolution)
            options.magnify = 1;
            options.resolution = 864;
        else
            options.magnify = options.resolution ./ get(0, 'ScreenPixelsPerInch');
        end
    elseif isempty(options.resolution)
        options.resolution = 864;
    end

    % Set the default format
    if ~isvector(options) && ~isbitmap(options)
        options.png = true;
    end

    % Check whether transparent background is wanted (old way)
    if isequal(get(ancestor(fig(1), 'figure'), 'Color'), 'none')
        options.transparent = true;
    end

    % If requested, set the resolution to the native vertical resolution of the
    % first suitable image found
    if native
        if isbitmap(options)
            % Find a suitable image
            list = findall(fig, 'Type','image', 'Tag','export_fig_native');
            if isempty(list)
                list = findall(fig, 'Type','image', 'Visible','on');
            end
            for hIm = list(:)'
                % Check height is >= 2
                height = size(get(hIm, 'CData'), 1);
                if height < 2
                    continue
                end
                % Account for the image filling only part of the axes, or vice versa
                yl = get(hIm, 'YData');
                if isscalar(yl)
                    yl = [yl(1)-0.5 yl(1)+height+0.5];
                else
                    yl = [min(yl), max(yl)];  % fix issue #151 (case of yl containing more than 2 elements)
                    if ~diff(yl)
                        continue
                    end
                    yl = yl + [-0.5 0.5] * (diff(yl) / (height - 1));
                end
                hAx = get(hIm, 'Parent');
                yl2 = get(hAx, 'YLim');
                % Find the pixel height of the axes
                oldUnits = get(hAx, 'Units');
                set(hAx, 'Units', 'pixels');
                pos = get(hAx, 'Position');
                set(hAx, 'Units', oldUnits);
                if ~pos(4)
                    continue
                end
                % Found a suitable image
                % Account for stretch-to-fill being disabled
                pbar = get(hAx, 'PlotBoxAspectRatio');
                pos = min(pos(4), pbar(2)*pos(3)/pbar(1));
                % Set the magnification to give native resolution
                options.magnify = abs((height * diff(yl2)) / (pos * diff(yl)));  % magnification must never be negative: issue #103
                break
            end
        elseif options.resolution == 864  % don't use -r864 in vector mode if user asked for -native
            options.resolution = []; % issue #241 (internal Matlab bug produces black lines with -r864)
        end
    end
end

% Convert a possible string => char (issue #245)
function value = str2char(value)
    if isa(value,'string')
        value = char(value);
    end
end

function A = downsize(A, factor)
    % Downsample an image
    if factor == 1
        % Nothing to do
        return
    end
    try
        % Faster, but requires image processing toolbox
        A = imresize(A, 1/factor, 'bilinear');
    catch
        % No image processing toolbox - resize manually
        % Lowpass filter - use Gaussian as is separable, so faster
        % Compute the 1d Gaussian filter
        filt = (-factor-1:factor+1) / (factor * 0.6);
        filt = exp(-filt .* filt);
        % Normalize the filter
        filt = single(filt / sum(filt));
        % Filter the image
        padding = floor(numel(filt) / 2);
        for a = 1:size(A, 3)
            A(:,:,a) = conv2(filt, filt', single(A([ones(1, padding) 1:end repmat(end, 1, padding)],[ones(1, padding) 1:end repmat(end, 1, padding)],a)), 'valid');
        end
        % Subsample
        A = A(1+floor(mod(end-1, factor)/2):factor:end,1+floor(mod(end-1, factor)/2):factor:end,:);
    end
end

function A = rgb2grey(A)
    A = cast(reshape(reshape(single(A), [], 3) * single([0.299; 0.587; 0.114]), size(A, 1), size(A, 2)), class(A)); % #ok<ZEROLIKE>
end

function A = check_greyscale(A)
    % Check if the image is greyscale
    if size(A, 3) == 3 && ...
            all(reshape(A(:,:,1) == A(:,:,2), [], 1)) && ...
            all(reshape(A(:,:,2) == A(:,:,3), [], 1))
        A = A(:,:,1); % Save only one channel for 8-bit output
    end
end

function eps_remove_background(fname, count)
    % Remove the background of an eps file
    % Open the file
    fh = fopen(fname, 'r+');
    if fh == -1
        error('Not able to open file %s.', fname);
    end
    % Read the file line by line
    while count
        % Get the next line
        l = fgets(fh);
        if isequal(l, -1)
            break; % Quit, no rectangle found
        end
        % Check if the line contains the background rectangle
        if isequal(regexp(l, ' *0 +0 +\d+ +\d+ +r[fe] *[\n\r]+', 'start'), 1)
            % Set the line to whitespace and quit
            l(1:regexp(l, '[\n\r]', 'start', 'once')-1) = ' ';
            fseek(fh, -numel(l), 0);
            fprintf(fh, l);
            % Reduce the count
            count = count - 1;
        end
    end
    % Close the file
    fclose(fh);
end

function b = isvector(options)
    b = options.pdf || options.eps;
end

function b = isbitmap(options)
    b = options.png || options.tif || options.jpg || options.bmp || options.im || options.alpha;
end

% Helper function
function A = make_cell(A)
    if ~iscell(A)
        A = {A};
    end
end

function add_bookmark(fname, bookmark_text)
    % Adds a bookmark to the temporary EPS file after %%EndPageSetup
    % Read in the file
    fh = fopen(fname, 'r');
    if fh == -1
        error('File %s not found.', fname);
    end
    try
        fstrm = fread(fh, '*char')';
    catch ex
        fclose(fh);
        rethrow(ex);
    end
    fclose(fh);

    % Include standard pdfmark prolog to maximize compatibility
    fstrm = strrep(fstrm, '%%BeginProlog', sprintf('%%%%BeginProlog\n/pdfmark where {pop} {userdict /pdfmark /cleartomark load put} ifelse'));
    % Add page bookmark
    fstrm = strrep(fstrm, '%%EndPageSetup', sprintf('%%%%EndPageSetup\n[ /Title (%s) /OUT pdfmark',bookmark_text));

    % Write out the updated file
    fh = fopen(fname, 'w');
    if fh == -1
        error('Unable to open %s for writing.', fname);
    end
    try
        fwrite(fh, fstrm, 'char*1');
    catch ex
        fclose(fh);
        rethrow(ex);
    end
    fclose(fh);
end

function set_tick_mode(Hlims, ax)
    % Set the tick mode of linear axes to manual
    % Leave log axes alone as these are tricky
    M = get(Hlims, [ax 'Scale']);
    if ~iscell(M)
        M = {M};
    end
    %idx = cellfun(@(c) strcmp(c, 'linear'), M);
    idx = find(strcmp(M,'linear'));
    %set(Hlims(idx), [ax 'TickMode'], 'manual');  % issue #187
    %set(Hlims(idx), [ax 'TickLabelMode'], 'manual');  % this hides exponent label in HG2!
    for idx2 = 1 : numel(idx)
        try
            % Fix for issue #187 - only set manual ticks when no exponent is present
            hAxes = Hlims(idx(idx2));
            props = {[ax 'TickMode'],'manual', [ax 'TickLabelMode'],'manual'};
            tickVals = get(hAxes,[ax 'Tick']);
            tickStrs = get(hAxes,[ax 'TickLabel']);
            try % Fix issue #236
                exponents = [hAxes.([ax 'Axis']).SecondaryLabel];
            catch
                exponents = [hAxes.([ax 'Ruler']).SecondaryLabel];
            end
            if isempty([exponents.String])
                % Fix for issue #205 - only set manual ticks when the Ticks number match the TickLabels number
                if numel(tickVals) == numel(tickStrs)
                    set(hAxes, props{:});  % no exponent and matching ticks, so update both ticks and tick labels to manual
                end
            end
        catch  % probably HG1
            % Fix for issue #220 - exponent is removed in HG1 when TickMode is 'manual' (internal Matlab bug)
            if isequal(tickVals, str2num(tickStrs)') %#ok<ST2NM>
                set(hAxes, props{:});  % revert back to old behavior
            end
        end
    end
end

function change_rgb_to_cmyk(fname)  % convert RGB => CMYK within an EPS file
    % Do post-processing on the eps file
    try
        % Read the EPS file into memory
        fstrm = read_write_entire_textfile(fname);

        % Replace all gray-scale colors
        fstrm = regexprep(fstrm, '\n([\d.]+) +GC\n', '\n0 0 0 ${num2str(1-str2num($1))} CC\n');
        
        % Replace all RGB colors
        fstrm = regexprep(fstrm, '\n[0.]+ +[0.]+ +[0.]+ +RC\n', '\n0 0 0 1 CC\n');  % pure black
        fstrm = regexprep(fstrm, '\n([\d.]+) +([\d.]+) +([\d.]+) +RC\n', '\n${sprintf(''%.4g '',[1-[str2num($1),str2num($2),str2num($3)]/max([str2num($1),str2num($2),str2num($3)]),1-max([str2num($1),str2num($2),str2num($3)])])} CC\n');

        % Overwrite the file with the modified contents
        read_write_entire_textfile(fname, fstrm);
    catch
        % never mind - leave as is...
    end
end

function hBlackAxles = fixBlackAxle(hAxes, axleName)
    hBlackAxles = [];
    for idx = 1 : numel(hAxes)
        ax = hAxes(idx);
        axleColor = get(ax, axleName);
        if isequal(axleColor,[0,0,0]) || isequal(axleColor,'k')
            hBlackAxles(end+1) = ax; %#ok<AGROW>
        end
    end
    set(hBlackAxles, axleName, [0,0,0.01]);  % off-black
end

% Issue #269: format-specific options
function [optionsCells, bitDepth] = getFormatOptions(options, formatName)
    bitDepth = [];
    try
        optionsStruct = options.format_options.(lower(formatName));
    catch
        % User did not specify any extra parameters for this format
        optionsCells = {};
        return
    end
    optionNames = fieldnames(optionsStruct);
    optionVals  = struct2cell(optionsStruct);
    optionsCells = [optionNames, optionVals]';
    if nargout < 2, return, end  % bail out if BitDepth is not required
    try
        idx = find(strcmpi(optionNames,'BitDepth'), 1, 'last');
        if ~isempty(idx)
            bitDepth = optionVals{idx};
        end
    catch
        % never mind - ignore
    end
end

function fh = copyfig(fh)
%COPYFIG Create a copy of a figure, without changing the figure
%
% Examples:
%   fh_new = copyfig(fh_old)
%
% This function will create a copy of a figure, but not change the figure,
% as copyobj sometimes does, e.g. by changing legends.
%
% IN:
%    fh_old - The handle of the figure to be copied. Default: gcf.
%
% OUT:
%    fh_new - The handle of the created figure.

% Copyright (C) Oliver Woodford 2012, Yair Altman 2015

% 26/02/15: If temp dir is not writable, use the dest folder for temp
%           destination files (Javier Paredes)
% 15/04/15: Suppress warnings during copyobj (Dun Kirk comment on FEX page 2013-10-02)
% 09/09/18: Fix issue #252: Workaround for cases where copyobj() fails for any reason

    % Set the default
    if nargin == 0
        fh = gcf;
    end
    % Is there a legend?
    useCopyobj = isempty(findall(fh, 'Type', 'axes', 'Tag', 'legend'));
    if useCopyobj
        % Safe to copy using copyobj
        oldWarn = warning('off'); %Suppress warnings during copyobj (Dun Kirk comment on FEX page 2013-10-02)
        try
            fh = copyobj(fh, 0);
        catch
            % Fix issue #252: Workaround for cases where copyobj() fails for any reason
            useCopyobj = false;  % if copyobj() croaks, use file save/load below
        end
        warning(oldWarn);
    end
    if ~useCopyobj
        % copyobj will change the figure, so save and then load it instead
        tmp_nam = [tempname '.fig'];
        try
            % Ensure that the temp dir is writable (Javier Paredes 26/2/15)
            fid = fopen(tmp_nam,'w');
            fwrite(fid,1);
            fclose(fid);
            delete(tmp_nam);  % cleanup
        catch
            % Temp dir is not writable, so use the current folder
            [dummy,fname,fext] = fileparts(tmp_nam); %#ok<ASGLU>
            fpath = pwd;
            tmp_nam = fullfile(fpath,[fname fext]);
        end
        hgsave(fh, tmp_nam);
        fh = hgload(tmp_nam);
        delete(tmp_nam);
    end
end

function [A, vA, vB, bb_rel] = crop_borders(A, bcol, padding, crop_amounts)
%CROP_BORDERS Crop the borders of an image or stack of images
%
%   [B, vA, vB, bb_rel] = crop_borders(A, bcol, [padding])
%
%IN:
%   A - HxWxCxN stack of images.
%   bcol - Cx1 background colour vector.
%   padding - scalar indicating how much padding to have in relation to
%             the cropped-image-size (0<=padding<=1). Default: 0
%   crop_amounts - 4-element vector of crop amounts: [top,right,bottom,left]
%             where NaN/Inf indicate auto-cropping, 0 means no cropping,
%             and any other value mean cropping in pixel amounts.
%
%OUT:
%   B - JxKxCxN cropped stack of images.
%   vA     - coordinates in A that contain the cropped image
%   vB     - coordinates in B where the cropped version of A is placed
%   bb_rel - relative bounding box (used for eps-cropping)

%{
% 06/03/15: Improved image cropping thanks to Oscar Hartogensis
% 08/06/15: Fixed issue #76: case of transparent figure bgcolor
% 21/02/16: Enabled specifying non-automated crop amounts
% 04/04/16: Fix per Luiz Carvalho for old Matlab releases
% 23/10/16: Fixed issue #175: there used to be a 1px minimal padding in case of crop, now removed
%}

    if nargin < 3
        padding = 0;
    end
    if nargin < 4
        crop_amounts = nan(1,4);  % =auto-cropping
    end
    crop_amounts(end+1:4) = NaN;  % fill missing values with NaN

    [h, w, c, n] = size(A);
    if isempty(bcol)  % case of transparent bgcolor
        bcol = A(ceil(end/2),1,:,1);
    end
    if isscalar(bcol)
        bcol = bcol(ones(c, 1));
    end

    % Crop margin from left
    if ~isfinite(crop_amounts(4))
        bail = false;
        for l = 1:w
            for a = 1:c
                if ~all(col(A(:,l,a,:)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        l = 1 + abs(crop_amounts(4));
    end

    % Crop margin from right
    if ~isfinite(crop_amounts(2))
        bcol = A(ceil(end/2),w,:,1);
        bail = false;
        for r = w:-1:l
            for a = 1:c
                if ~all(col(A(:,r,a,:)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        r = w - abs(crop_amounts(2));
    end

    % Crop margin from top
    if ~isfinite(crop_amounts(1))
        bcol = A(1,ceil(end/2),:,1);
        bail = false;
        for t = 1:h
            for a = 1:c
                if ~all(col(A(t,:,a,:)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        t = 1 + abs(crop_amounts(1));
    end

    % Crop margin from bottom
    bcol = A(h,ceil(end/2),:,1);
    if ~isfinite(crop_amounts(3))
        bail = false;
        for b = h:-1:t
            for a = 1:c
                if ~all(col(A(b,:,a,:)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        b = h - abs(crop_amounts(3));
    end

    if padding == 0  % no padding
        % Issue #175: there used to be a 1px minimal padding in case of crop, now removed
        %{
        if ~isequal([t b l r], [1 h 1 w]) % Check if we're actually croppping
            padding = 1; % Leave one boundary pixel to avoid bleeding on resize
            bcol(:) = nan;  % make the 1px padding transparent
        end
        %}
    elseif abs(padding) < 1  % pad value is a relative fraction of image size
        padding = sign(padding)*round(mean([b-t r-l])*abs(padding)); % ADJUST PADDING
    else  % pad value is in units of 1/72" points
        padding = round(padding);  % fix cases of non-integer pad value
    end

    if padding > 0  % extra padding
        % Create an empty image, containing the background color, that has the
        % cropped image size plus the padded border
        B = repmat(bcol,[(b-t)+1+padding*2,(r-l)+1+padding*2,1,n]);  % Fix per Luiz Carvalho
        % vA - coordinates in A that contain the cropped image
        vA = [t b l r];
        % vB - coordinates in B where the cropped version of A will be placed
        vB = [padding+1, (b-t)+1+padding, padding+1, (r-l)+1+padding];
        % Place the original image in the empty image
        B(vB(1):vB(2), vB(3):vB(4), :, :) = A(vA(1):vA(2), vA(3):vA(4), :, :);
        A = B;
    else  % extra cropping
        vA = [t-padding b+padding l-padding r+padding];
        A = A(vA(1):vA(2), vA(3):vA(4), :, :);
        vB = [NaN NaN NaN NaN];
    end

    % For EPS cropping, determine the relative BoundingBox - bb_rel
    bb_rel = [l-1 h-b-1 r+1 h-t+1]./[w h w h];
end

function A = col(A)
    A = A(:);
end

%FIX_LINES  Improves the line style of eps files generated by print
%
% Examples:
%   fix_lines fname
%   fix_lines fname fname2
%   fstrm_out = fixlines(fstrm_in)
%
% This function improves the style of lines in eps files generated by
% MATLAB's print function, making them more similar to those seen on
% screen. Grid lines are also changed from a dashed style to a dotted
% style, for greater differentiation from dashed lines.
% 
% The function also places embedded fonts after the postscript header, in
% versions of MATLAB which place the fonts first (R2006b and earlier), in
% order to allow programs such as Ghostscript to find the bounding box
% information.
%
%IN:
%   fname - Name or path of source eps file.
%   fname2 - Name or path of destination eps file. Default: same as fname.
%   fstrm_in - File contents of a MATLAB-generated eps file.
%
%OUT:
%   fstrm_out - Contents of the eps file with line styles fixed.

% Copyright: (C) Oliver Woodford, 2008-2014

% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928)
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)

% Thank you to Sylvain Favrot for bringing the embedded font/bounding box
% interaction in older versions of MATLAB to my attention.
% Thank you to D Ko for bringing an error with eps files with tiff previews
% to my attention.
% Thank you to Laurence K for suggesting the check to see if the file was
% opened.

% 01/03/15: Issue #20: warn users if using this function in HG2 (R2014b+)
% 27/03/15: Fixed out of memory issue with enormous EPS files (generated by print() with OpenGL renderer), related to issue #39

function fstrm = fix_lines(fstrm, fname2)

% Issue #20: warn users if using this function in HG2 (R2014b+)
if using_hg2
    warning('export_fig:hg2','The fix_lines function should not be used in this Matlab version.');
end
    
if nargout == 0 || nargin > 1
    if nargin < 2
        % Overwrite the input file
        fname2 = fstrm;
    end
    % Read in the file
    fstrm = read_write_entire_textfile(fstrm);
end

% Move any embedded fonts after the postscript header
if strcmp(fstrm(1:15), '%!PS-AdobeFont-')
    % Find the start and end of the header
    ind = regexp(fstrm, '[\n\r]%!PS-Adobe-');
    [ind2, ind2] = regexp(fstrm, '[\n\r]%%EndComments[\n\r]+');
    % Put the header first
    if ~isempty(ind) && ~isempty(ind2) && ind(1) < ind2(1)
        fstrm = fstrm([ind(1)+1:ind2(1) 1:ind(1) ind2(1)+1:end]);
    end
end

% Make sure all line width commands come before the line style definitions,
% so that dash lengths can be based on the correct widths
% Find all line style sections
ind = [regexp(fstrm, '[\n\r]SO[\n\r]'),... % This needs to be here even though it doesn't have dots/dashes!
       regexp(fstrm, '[\n\r]DO[\n\r]'),...
       regexp(fstrm, '[\n\r]DA[\n\r]'),...
       regexp(fstrm, '[\n\r]DD[\n\r]')];
ind = sort(ind);
% Find line width commands
[ind2, ind3] = regexp(fstrm, '[\n\r]\d* w[\n\r]');
% Go through each line style section and swap with any line width commands
% near by
b = 1;
m = numel(ind);
n = numel(ind2);
for a = 1:m
    % Go forwards width commands until we pass the current line style
    while b <= n && ind2(b) < ind(a)
        b = b + 1;
    end
    if b > n
        % No more width commands
        break;
    end
    % Check we haven't gone past another line style (including SO!)
    if a < m && ind2(b) > ind(a+1)
        continue;
    end
    % Are the commands close enough to be confident we can swap them?
    if (ind2(b) - ind(a)) > 8
        continue;
    end
    % Move the line style command below the line width command
    fstrm(ind(a)+1:ind3(b)) = [fstrm(ind(a)+4:ind3(b)) fstrm(ind(a)+1:ind(a)+3)];
    b = b + 1;
end

% Find any grid line definitions and change to GR format
% Find the DO sections again as they may have moved
ind = int32(regexp(fstrm, '[\n\r]DO[\n\r]'));
if ~isempty(ind)
    % Find all occurrences of what are believed to be axes and grid lines
    ind2 = int32(regexp(fstrm, '[\n\r] *\d* *\d* *mt *\d* *\d* *L[\n\r]'));
    if ~isempty(ind2)
        % Now see which DO sections come just before axes and grid lines
        ind2 = repmat(ind2', [1 numel(ind)]) - repmat(ind, [numel(ind2) 1]);
        ind2 = any(ind2 > 0 & ind2 < 12); % 12 chars seems about right
        ind = ind(ind2);
        % Change any regions we believe to be grid lines to GR
        fstrm(ind+1) = 'G';
        fstrm(ind+2) = 'R';
    end
end

% Define the new styles, including the new GR format
% Dot and dash lengths have two parts: a constant amount plus a line width
% variable amount. The constant amount comes after dpi2point, and the
% variable amount comes after currentlinewidth. If you want to change
% dot/dash lengths for a one particular line style only, edit the numbers
% in the /DO (dotted lines), /DA (dashed lines), /DD (dot dash lines) and
% /GR (grid lines) lines for the style you want to change.
new_style = {'/dom { dpi2point 1 currentlinewidth 0.08 mul add mul mul } bdef',... % Dot length macro based on line width
             '/dam { dpi2point 2 currentlinewidth 0.04 mul add mul mul } bdef',... % Dash length macro based on line width
             '/SO { [] 0 setdash 0 setlinecap } bdef',... % Solid lines
             '/DO { [1 dom 1.2 dom] 0 setdash 0 setlinecap } bdef',... % Dotted lines
             '/DA { [4 dam 1.5 dam] 0 setdash 0 setlinecap } bdef',... % Dashed lines
             '/DD { [1 dom 1.2 dom 4 dam 1.2 dom] 0 setdash 0 setlinecap } bdef',... % Dot dash lines
             '/GR { [0 dpi2point mul 4 dpi2point mul] 0 setdash 1 setlinecap } bdef'}; % Grid lines - dot spacing remains constant

% Construct the output
% This is the original (memory-intensive) code:
%first_sec = strfind(fstrm, '% line types:'); % Isolate line style definition section
%[second_sec, remaining] = strtok(fstrm(first_sec+1:end), '/');
%[remaining, remaining] = strtok(remaining, '%');
%fstrm = [fstrm(1:first_sec) second_sec sprintf('%s\r', new_style{:}) remaining];
fstrm = regexprep(fstrm,'(% line types:.+?)/.+?%',['$1',sprintf('%s\r',new_style{:}),'%']);

% Write the output file
if nargout == 0 || nargin > 1
    read_write_entire_textfile(fname2, fstrm);
end
end


function varargout = ghostscript(cmd)
%GHOSTSCRIPT  Calls a local GhostScript executable with the input command
%
% Example:
%   [status result] = ghostscript(cmd)
%
% Attempts to locate a ghostscript executable, finally asking the user to
% specify the directory ghostcript was installed into. The resulting path
% is stored for future reference.
% 
% Once found, the executable is called with the input command string.
%
% This function requires that you have Ghostscript installed on your
% system. You can download this from: http://www.ghostscript.com
%
% IN:
%   cmd - Command string to be passed into ghostscript.
%
% OUT:
%   status - 0 iff command ran without problem.
%   result - Output from ghostscript.

% Copyright: Oliver Woodford, 2009-2015, Yair Altman 2015-
%{
% Thanks to Jonas Dorn for the fix for the title of the uigetdir window on Mac OS.
% Thanks to Nathan Childress for the fix to default location on 64-bit Windows systems.
% 27/04/11 - Find 64-bit Ghostscript on Windows. Thanks to Paul Durack and
%            Shaun Kline for pointing out the issue
% 04/05/11 - Thanks to David Chorlian for pointing out an alternative
%            location for gs on linux.
% 12/12/12 - Add extra executable name on Windows. Thanks to Ratish
%            Punnoose for highlighting the issue.
% 28/06/13 - Fix error using GS 9.07 in Linux. Many thanks to Jannick
%            Steinbring for proposing the fix.
% 24/10/13 - Fix error using GS 9.07 in Linux. Many thanks to Johannes 
%            for the fix.
% 23/01/14 - Add full path to ghostscript.txt in warning. Thanks to Koen
%            Vermeer for raising the issue.
% 27/02/15 - If Ghostscript croaks, display suggested workarounds
% 30/03/15 - Improved performance by caching status of GS path check, if ok
% 14/05/15 - Clarified warning message in case GS path could not be saved
% 29/05/15 - Avoid cryptic error in case the ghostscipt path cannot be saved (issue #74)
% 10/11/15 - Custom GS installation webpage for MacOS. Thanks to Andy Hueni via FEX
%}

    try
        % Call ghostscript
        [varargout{1:nargout}] = system([gs_command(gs_path()) cmd]);
    catch err
        % Display possible workarounds for Ghostscript croaks
        url1 = 'https://github.com/altmany/export_fig/issues/12#issuecomment-61467998';  % issue #12
        url2 = 'https://github.com/altmany/export_fig/issues/20#issuecomment-63826270';  % issue #20
        hg2_str = ''; if using_hg2, hg2_str = ' or Matlab R2014a'; end
        fprintf(2, 'Ghostscript error. Rolling back to GS 9.10%s may possibly solve this:\n * <a href="%s">%s</a> ',hg2_str,url1,url1);
        if using_hg2
            fprintf(2, '(GS 9.10)\n * <a href="%s">%s</a> (R2014a)',url2,url2);
        end
        fprintf('\n\n');
        if ismac || isunix
            url3 = 'https://github.com/altmany/export_fig/issues/27';  % issue #27
            fprintf(2, 'Alternatively, this may possibly be due to a font path issue:\n * <a href="%s">%s</a>\n\n',url3,url3);
            % issue #20
            fpath = which(mfilename);
            if isempty(fpath), fpath = [mfilename('fullpath') '.m']; end
            fprintf(2, 'Alternatively, if you are using csh, modify shell_cmd from "export..." to "setenv ..."\nat the bottom of <a href="matlab:opentoline(''%s'',174)">%s</a>\n\n',fpath,fpath);
        end
        rethrow(err);
    end
end

function path_ = gs_path
    % Return a valid path
    % Start with the currently set path
    path_ = user_string('ghostscript');
    % Check the path works
    if check_gs_path(path_)
        return
    end
    % Check whether the binary is on the path
    if ispc
        bin = {'gswin32c.exe', 'gswin64c.exe', 'gs'};
    else
        bin = {'gs'};
    end
    for a = 1:numel(bin)
        path_ = bin{a};
        if check_store_gs_path(path_)
            return
        end
    end
    % Search the obvious places
    if ispc
        default_location = 'C:\Program Files\gs\';
        dir_list = dir(default_location);
        if isempty(dir_list)
            default_location = 'C:\Program Files (x86)\gs\'; % Possible location on 64-bit systems
            dir_list = dir(default_location);
        end
        executable = {'\bin\gswin32c.exe', '\bin\gswin64c.exe'};
        ver_num = 0;
        % If there are multiple versions, use the newest
        for a = 1:numel(dir_list)
            ver_num2 = sscanf(dir_list(a).name, 'gs%g');
            if ~isempty(ver_num2) && ver_num2 > ver_num
                for b = 1:numel(executable)
                    path2 = [default_location dir_list(a).name executable{b}];
                    if exist(path2, 'file') == 2
                        path_ = path2;
                        ver_num = ver_num2;
                    end
                end
            end
        end
        if check_store_gs_path(path_)
            return
        end
    else
        executable = {'/usr/bin/gs', '/usr/local/bin/gs'};
        for a = 1:numel(executable)
            path_ = executable{a};
            if check_store_gs_path(path_)
                return
            end
        end
    end
    % Ask the user to enter the path
    while true
        if strncmp(computer, 'MAC', 3) % Is a Mac
            % Give separate warning as the uigetdir dialogue box doesn't have a
            % title
            uiwait(warndlg('Ghostscript not found. Please locate the program.'))
        end
        base = uigetdir('/', 'Ghostcript not found. Please locate the program.');
        if isequal(base, 0)
            % User hit cancel or closed window
            break;
        end
        base = [base filesep]; %#ok<AGROW>
        bin_dir = {'', ['bin' filesep], ['lib' filesep]};
        for a = 1:numel(bin_dir)
            for b = 1:numel(bin)
                path_ = [base bin_dir{a} bin{b}];
                if exist(path_, 'file') == 2
                    if check_store_gs_path(path_)
                        return
                    end
                end
            end
        end
    end
    if ismac
        error('Ghostscript not found. Have you installed it (http://pages.uoregon.edu/koch)?');
    else
        error('Ghostscript not found. Have you installed it from www.ghostscript.com?');
    end
end

function good = check_store_gs_path(path_)
    % Check the path is valid
    good = check_gs_path(path_);
    if ~good
        return
    end
    % Update the current default path to the path found
    if ~user_string('ghostscript', path_)
        filename = fullfile(fileparts(which('user_string.m')), '.ignore', 'ghostscript.txt');
        warning('Path to ghostscript installation could not be saved in %s (perhaps a permissions issue). You can manually create this file and set its contents to %s, to improve performance in future invocations (this warning is safe to ignore).', filename, path_);
        return
    end
end

function good = check_gs_path(path_)
    persistent isOk
    if isempty(path_)
        isOk = false;
    elseif ~isequal(isOk,true)
        % Check whether the path is valid
        [status, message] = system([gs_command(path_) '-h']); %#ok<ASGLU>
        isOk = status == 0;
    end
    good = isOk;
end

function cmd = gs_command(path_)
    % Initialize any required system calls before calling ghostscript
    % TODO: in Unix/Mac, find a way to determine whether to use "export" (bash) or "setenv" (csh/tcsh)
    shell_cmd = '';
    if isunix
        shell_cmd = 'export LD_LIBRARY_PATH=""; '; % Avoids an error on Linux with GS 9.07
    end
    if ismac
        shell_cmd = 'export DYLD_LIBRARY_PATH=""; ';  % Avoids an error on Mac with GS 9.07
    end
    % Construct the command string
    cmd = sprintf('%s"%s" ', shell_cmd, path_);
end


function fh = isolate_axes(ah, vis)
%ISOLATE_AXES Isolate the specified axes in a figure on their own
%
% Examples:
%   fh = isolate_axes(ah)
%   fh = isolate_axes(ah, vis)
%
% This function will create a new figure containing the axes/uipanels
% specified, and also their associated legends and colorbars. The objects
% specified must all be in the same figure, but they will generally only be
% a subset of the objects in the figure.
%
% IN:
%    ah - An array of axes and uipanel handles, which must come from the
%         same figure.
%    vis - A boolean indicating whether the new figure should be visible.
%          Default: false.
%
% OUT:
%    fh - The handle of the created figure.

% Copyright (C) Oliver Woodford 2011-2013

% Thank you to Rosella Blatt for reporting a bug to do with axes in GUIs
% 16/03/12: Moved copyfig to its own function. Thanks to Bob Fratantonio
%           for pointing out that the function is also used in export_fig.m
% 12/12/12: Add support for isolating uipanels. Thanks to michael for suggesting it
% 08/10/13: Bug fix to allchildren suggested by Will Grant (many thanks!)
% 05/12/13: Bug fix to axes having different units. Thanks to Remington Reid for reporting
% 21/04/15: Bug fix for exporting uipanels with legend/colorbar on HG1 (reported by Alvaro
%           on FEX page as a comment on 24-Apr-2014); standardized indentation & help section
% 22/04/15: Bug fix: legends and colorbars were not exported when exporting axes handle in HG2

    % Make sure we have an array of handles
    if ~all(ishandle(ah))
        error('ah must be an array of handles');
    end
    % Check that the handles are all for axes or uipanels, and are all in the same figure
    fh = ancestor(ah(1), 'figure');
    nAx = numel(ah);
    for a = 1:nAx
        if ~ismember(get(ah(a), 'Type'), {'axes', 'uipanel'})
            error('All handles must be axes or uipanel handles.');
        end
        if ~isequal(ancestor(ah(a), 'figure'), fh)
            error('Axes must all come from the same figure.');
        end
    end
    % Tag the objects so we can find them in the copy
    old_tag = get(ah, 'Tag');
    if nAx == 1
        old_tag = {old_tag};
    end
    set(ah, 'Tag', 'ObjectToCopy');
    % Create a new figure exactly the same as the old one
    fh = copyfig(fh); %copyobj(fh, 0);
    if nargin < 2 || ~vis
        set(fh, 'Visible', 'off');
    end
    % Reset the object tags
    for a = 1:nAx
        set(ah(a), 'Tag', old_tag{a});
    end
    % Find the objects to save
    ah = findall(fh, 'Tag', 'ObjectToCopy');
    if numel(ah) ~= nAx
        close(fh);
        error('Incorrect number of objects found.');
    end
    % Set the axes tags to what they should be
    for a = 1:nAx
        set(ah(a), 'Tag', old_tag{a});
    end
    % Keep any legends and colorbars which overlap the subplots
    % Note: in HG1 these are axes objects; in HG2 they are separate objects, therefore we
    %       don't test for the type, only the tag (hopefully nobody but Matlab uses them!)
    lh = findall(fh, 'Tag', 'legend', '-or', 'Tag', 'Colorbar');
    nLeg = numel(lh);
    if nLeg > 0
        set([ah(:); lh(:)], 'Units', 'normalized');
        try
            ax_pos = get(ah, 'OuterPosition'); % axes and figures have the OuterPosition property
        catch
            ax_pos = get(ah, 'Position'); % uipanels only have Position, not OuterPosition
        end
        if nAx > 1
            ax_pos = cell2mat(ax_pos(:));
        end
        ax_pos(:,3:4) = ax_pos(:,3:4) + ax_pos(:,1:2);
        try
            leg_pos = get(lh, 'OuterPosition');
        catch
            leg_pos = get(lh, 'Position');  % No OuterPosition in HG2, only in HG1
        end
        if nLeg > 1;
            leg_pos = cell2mat(leg_pos);
        end
        leg_pos(:,3:4) = leg_pos(:,3:4) + leg_pos(:,1:2);
        ax_pos = shiftdim(ax_pos, -1);
        % Overlap test
        M = bsxfun(@lt, leg_pos(:,1), ax_pos(:,:,3)) & ...
            bsxfun(@lt, leg_pos(:,2), ax_pos(:,:,4)) & ...
            bsxfun(@gt, leg_pos(:,3), ax_pos(:,:,1)) & ...
            bsxfun(@gt, leg_pos(:,4), ax_pos(:,:,2));
        ah = [ah; lh(any(M, 2))];
    end
    % Get all the objects in the figure
    axs = findall(fh);
    % Delete everything except for the input objects and associated items
    delete(axs(~ismember(axs, [ah; allchildren(ah); allancestors(ah)])));
end

function ah = allchildren(ah)
    ah = findall(ah);
    if iscell(ah)
        ah = cell2mat(ah);
    end
    ah = ah(:);
end

function ph = allancestors(ah)
    ph = [];
    for a = 1:numel(ah)
        h = get(ah(a), 'parent');
        while h ~= 0
            ph = [ph; h];
            h = get(h, 'parent');
        end
    end
end


function [A, bcol] = print2array(fig, res, renderer, gs_options)
%PRINT2ARRAY  Exports a figure to an image array
%
% Examples:
%   A = print2array
%   A = print2array(figure_handle)
%   A = print2array(figure_handle, resolution)
%   A = print2array(figure_handle, resolution, renderer)
%   A = print2array(figure_handle, resolution, renderer, gs_options)
%   [A bcol] = print2array(...)
%
% This function outputs a bitmap image of the given figure, at the desired
% resolution.
%
% If renderer is '-painters' then ghostcript needs to be installed. This
% can be downloaded from: http://www.ghostscript.com
%
% IN:
%   figure_handle - The handle of the figure to be exported. Default: gcf.
%   resolution - Resolution of the output, as a factor of screen
%                resolution. Default: 1.
%   renderer - string containing the renderer paramater to be passed to
%              print. Default: '-opengl'.
%   gs_options - optional ghostscript options (e.g.: '-dNoOutputFonts'). If
%                multiple options are needed, enclose in call array: {'-a','-b'}
%
% OUT:
%   A - MxNx3 uint8 image of the figure.
%   bcol - 1x3 uint8 vector of the background color

% Copyright (C) Oliver Woodford 2008-2014, Yair Altman 2015-
%{
% 05/09/11: Set EraseModes to normal when using opengl or zbuffer
%           renderers. Thanks to Pawel Kocieniewski for reporting the issue.
% 21/09/11: Bug fix: unit8 -> uint8! Thanks to Tobias Lamour for reporting it.
% 14/11/11: Bug fix: stop using hardcopy(), as it interfered with figure size
%           and erasemode settings. Makes it a bit slower, but more reliable.
%           Thanks to Phil Trinh and Meelis Lootus for reporting the issues.
% 09/12/11: Pass font path to ghostscript.
% 27/01/12: Bug fix affecting painters rendering tall figures. Thanks to
%           Ken Campbell for reporting it.
% 03/04/12: Bug fix to median input. Thanks to Andy Matthews for reporting it.
% 26/10/12: Set PaperOrientation to portrait. Thanks to Michael Watts for
%           reporting the issue.
% 26/02/15: If temp dir is not writable, use the current folder for temp
%           EPS/TIF files (Javier Paredes)
% 27/02/15: Display suggested workarounds to internal print() error (issue #16)
% 28/02/15: Enable users to specify optional ghostscript options (issue #36)
% 10/03/15: Fixed minor warning reported by Paul Soderlind; fixed code indentation
% 28/05/15: Fixed issue #69: patches with LineWidth==0.75 appear wide (internal bug in Matlab's print() func)
% 07/07/15: Fixed issue #83: use numeric handles in HG1
% 11/12/16: Fixed cropping issue reported by Harry D.
% 29/09/18: Fixed issue #254: error in print2array>read_tif_img
%}

    % Generate default input arguments, if needed
    if nargin < 2
        res = 1;
        if nargin < 1
            fig = gcf;
        end
    end
    % Warn if output is large
    old_mode = get(fig, 'Units');
    set(fig, 'Units', 'pixels');
    px = get(fig, 'Position');
    set(fig, 'Units', old_mode);
    npx = prod(px(3:4)*res)/1e6;
    if npx > 30
        % 30M pixels or larger!
        warning('MATLAB:LargeImage', 'print2array generating a %.1fM pixel image. This could be slow and might also cause memory problems.', npx);
    end
    % Retrieve the background colour
    bcol = get(fig, 'Color');
    % Set the resolution parameter
    res_str = ['-r' num2str(ceil(get(0, 'ScreenPixelsPerInch')*res))];
    % Generate temporary file name
    tmp_nam = [tempname '.tif'];
    try
        % Ensure that the temp dir is writable (Javier Paredes 26/2/15)
        fid = fopen(tmp_nam,'w');
        fwrite(fid,1);
        fclose(fid);
        delete(tmp_nam);  % cleanup
        isTempDirOk = true;
    catch
        % Temp dir is not writable, so use the current folder
        [dummy,fname,fext] = fileparts(tmp_nam); %#ok<ASGLU>
        fpath = pwd;
        tmp_nam = fullfile(fpath,[fname fext]);
        isTempDirOk = false;
    end
    % Enable users to specify optional ghostscript options (issue #36)
    if nargin > 3 && ~isempty(gs_options)
        if iscell(gs_options)
            gs_options = sprintf(' %s',gs_options{:});
        elseif ~ischar(gs_options)
            error('gs_options input argument must be a string or cell-array of strings');
        else
            gs_options = [' ' gs_options];
        end
    else
        gs_options = '';
    end
    if nargin > 2 && strcmp(renderer, '-painters')
        % First try to print directly to tif file
        try
            % Print the file into a temporary TIF file and read it into array A
            [A, err, ex] = read_tif_img(fig, res_str, renderer, tmp_nam);
            if err, rethrow(ex); end
        catch  % error - try to print to EPS and then using Ghostscript to TIF
            % Print to eps file
            if isTempDirOk
                tmp_eps = [tempname '.eps'];
            else
                tmp_eps = fullfile(fpath,[fname '.eps']);
            end
            print2eps(tmp_eps, fig, 0, renderer, '-loose');
            try
                % Initialize the command to export to tiff using ghostscript
                cmd_str = ['-dEPSCrop -q -dNOPAUSE -dBATCH ' res_str ' -sDEVICE=tiff24nc'];
                % Set the font path
                fp = font_path();
                if ~isempty(fp)
                    cmd_str = [cmd_str ' -sFONTPATH="' fp '"'];
                end
                % Add the filenames
                cmd_str = [cmd_str ' -sOutputFile="' tmp_nam '" "' tmp_eps '"' gs_options];
                % Execute the ghostscript command
                ghostscript(cmd_str);
            catch me
                % Delete the intermediate file
                delete(tmp_eps);
                rethrow(me);
            end
            % Delete the intermediate file
            delete(tmp_eps);
            % Read in the generated bitmap
            A = imread(tmp_nam);
            % Delete the temporary bitmap file
            delete(tmp_nam);
        end
        % Set border pixels to the correct colour
        if isequal(bcol, 'none')
            bcol = [];
        elseif isequal(bcol, [1 1 1])
            bcol = uint8([255 255 255]);
        else
            for l = 1:size(A, 2)
                if ~all(reshape(A(:,l,:) == 255, [], 1))
                    break;
                end
            end
            for r = size(A, 2):-1:l
                if ~all(reshape(A(:,r,:) == 255, [], 1))
                    break;
                end
            end
            for t = 1:size(A, 1)
                if ~all(reshape(A(t,:,:) == 255, [], 1))
                    break;
                end
            end
            for b = size(A, 1):-1:t
                if ~all(reshape(A(b,:,:) == 255, [], 1))
                    break;
                end
            end
            bcol = uint8(median(single([reshape(A(:,[l r],:), [], size(A, 3)); reshape(A([t b],:,:), [], size(A, 3))]), 1));
            for c = 1:size(A, 3)
                A(:,[1:l-1, r+1:end],c) = bcol(c);
                A([1:t-1, b+1:end],:,c) = bcol(c);
            end
        end
    else
        if nargin < 3
            renderer = '-opengl';
        end
        % Print the file into a temporary TIF file and read it into array A
        [A, err, ex] = read_tif_img(fig, res_str, renderer, tmp_nam);
        % Throw any error that occurred
        if err
            % Display suggested workarounds to internal print() error (issue #16)
            fprintf(2, 'An error occured with Matlab''s builtin print function.\nTry setting the figure Renderer to ''painters'' or use opengl(''software'').\n\n');
            rethrow(ex);
        end
        % Set the background color
        if isequal(bcol, 'none')
            bcol = [];
        else
            bcol = bcol * 255;
            if isequal(bcol, round(bcol))
                bcol = uint8(bcol);
            else
                bcol = squeeze(A(1,1,:));
            end
        end
    end
    % Check the output size is correct
    if isequal(res, round(res))
        px = round([px([4 3])*res 3]);  % round() to avoid an indexing warning below
        if ~isequal(size(A), px)
            % Correct the output size
            A = A(1:min(end,px(1)),1:min(end,px(2)),:);
        end
    end
end

% Function to create a TIF image of the figure and read it into an array
function [A, err, ex] = read_tif_img(fig, res_str, renderer, tmp_nam)
    A =  [];  % fix for issue #254
    err = false;
    ex = [];
    % Temporarily set the paper size
    old_pos_mode    = get(fig, 'PaperPositionMode');
    old_orientation = get(fig, 'PaperOrientation');
    set(fig, 'PaperPositionMode','auto', 'PaperOrientation','portrait');
    try
        % Workaround for issue #69: patches with LineWidth==0.75 appear wide (internal bug in Matlab's print() function)
        fp = [];  % in case we get an error below
        fp = findall(fig, 'Type','patch', 'LineWidth',0.75);
        set(fp, 'LineWidth',0.5);
        % Fix issue #83: use numeric handles in HG1
        if ~using_hg2(fig),  fig = double(fig);  end
        % Print to tiff file
        print(fig, renderer, res_str, '-dtiff', tmp_nam);
        % Read in the printed file
        A = imread(tmp_nam);
        % Delete the temporary file
        delete(tmp_nam);
    catch ex
        err = true;
    end
    set(fp, 'LineWidth',0.75);  % restore original figure appearance
    % Reset the paper size
    set(fig, 'PaperPositionMode',old_pos_mode, 'PaperOrientation',old_orientation);
end

% Function to return (and create, where necessary) the font path
function fp = font_path()
    fp = user_string('gs_font_path');
    if ~isempty(fp)
        return
    end
    % Create the path
    % Start with the default path
    fp = getenv('GS_FONTPATH');
    % Add on the typical directories for a given OS
    if ispc
        if ~isempty(fp)
            fp = [fp ';'];
        end
        fp = [fp getenv('WINDIR') filesep 'Fonts'];
    else
        if ~isempty(fp)
            fp = [fp ':'];
        end
        fp = [fp '/usr/share/fonts:/usr/local/share/fonts:/usr/share/fonts/X11:/usr/local/share/fonts/X11:/usr/share/fonts/truetype:/usr/local/share/fonts/truetype'];
    end
    user_string('gs_font_path', fp);
end


function print2eps(name, fig, export_options, varargin)
%PRINT2EPS  Prints figures to eps with improved line styles
%
% Examples:
%   print2eps filename
%   print2eps(filename, fig_handle)
%   print2eps(filename, fig_handle, export_options)
%   print2eps(filename, fig_handle, export_options, print_options)
%
% This function saves a figure as an eps file, with two improvements over
% MATLAB's print command. First, it improves the line style, making dashed
% lines more like those on screen and giving grid lines a dotted line style.
% Secondly, it substitutes original font names back into the eps file,
% where these have been changed by MATLAB, for up to 11 different fonts.
%
% Inputs:
%   filename - string containing the name (optionally including full or
%              relative path) of the file the figure is to be saved as. A
%              ".eps" extension is added if not there already. If a path is
%              not specified, the figure is saved in the current directory.
%   fig_handle - The handle of the figure to be saved. Default: gcf().
%   export_options - array or struct of optional scalar values:
%       bb_padding - Scalar value of amount of padding to add to border around
%                    the cropped image, in points (if >1) or percent (if <1).
%                    Can be negative as well as positive; Default: 0
%       crop       - Cropping flag. Deafult: 0
%       fontswap   - Whether to swap non-default fonts in figure. Default: true
%       preserve_size - Whether to preserve the figure's PaperSize. Default: false
%       font_space - Character used to separate font-name terms in the EPS output
%                    e.g. "Courier New" => "Courier-New". Default: ''
%                    (available only via the struct alternative)
%       renderer   - Renderer used to generate bounding-box. Default: 'opengl'
%                    (available only via the struct alternative)
%       crop_amounts - 4-element vector of crop amounts: [top,right,bottom,left]
%                    (available only via the struct alternative)
%   print_options - Additional parameter strings to be passed to the print command

%{
% Copyright (C) Oliver Woodford 2008-2014, Yair Altman 2015-

% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928)
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)
%}
%{
% 14/11/11: Fix a MATLAB bug rendering black or white text incorrectly.
%           Thanks to Mathieu Morlighem for reporting the issue and
%           obtaining a fix from TMW.
% 08/12/11: Added ability to correct fonts. Several people have requested
%           this at one time or another, and also pointed me to printeps
%           (fex id: 7501), so thank you to them. My implementation (which
%           was not inspired by printeps - I'd already had the idea for my
%           approach) goes slightly further in that it allows multiple
%           fonts to be swapped.
% 14/12/11: Fix bug affecting font names containing spaces. Thanks to David
%           Szwer for reporting the issue.
% 25/01/12: Add a font not to be swapped. Thanks to Anna Rafferty and Adam
%           Jackson for reporting the issue. Also fix a bug whereby using a
%           font alias can lead to another font being swapped in.
% 10/04/12: Make the font swapping case insensitive.
% 26/10/12: Set PaperOrientation to portrait. Thanks to Michael Watts for
%           reporting the issue.
% 26/10/12: Fix issue to do with swapping fonts changing other fonts and
%           sizes we don't want, due to listeners. Thanks to Malcolm Hudson
%           for reporting the issue.
% 22/03/13: Extend font swapping to axes labels. Thanks to Rasmus Ischebeck
%           for reporting the issue.
% 23/07/13: Bug fix to font swapping. Thanks to George for reporting the
%           issue.
% 13/08/13: Fix MATLAB feature of not exporting white lines correctly.
%           Thanks to Sebastian Hesslinger for reporting it.
% 24/02/15: Fix for Matlab R2014b bug (issue #31): LineWidths<0.75 are not
%           set in the EPS (default line width is used)
% 25/02/15: Fixed issue #32: BoundingBox problem caused uncropped EPS/PDF files
% 05/03/15: Fixed issue #43: Inability to perform EPS file post-processing
% 06/03/15: Improved image padding & cropping thanks to Oscar Hartogensis
% 21/03/15: Fixed edge-case of missing handles having a 'FontName' property
% 26/03/15: Attempt to fix issue #45: white lines in subplots do not print correctly
% 27/03/15: Attempt to fix issue #44: white artifact lines appearing in patch exports
% 30/03/15: Fixed issue #52: improved performance on HG2 (R2014b+)
% 09/04/15: Comment blocks consolidation and minor code cleanup (no real code change)
% 12/04/15: Fixed issue #56: bad cropping
% 14/04/15: Workaround for issue #45: lines in image subplots are exported in invalid color
% 07/07/15: Added option to avoid font-swapping in EPS/PDF
% 07/07/15: Fixed issue #83: use numeric handles in HG1
% 22/07/15: Fixed issue #91 (thanks to Carlos Moffat)
% 28/09/15: Fixed issue #108 (thanks to JacobD10)
% 01/11/15: Fixed issue #112: optional renderer for bounding-box computation (thanks to Jes�s Pestana Puerta)
% 21/02/16: Enabled specifying non-automated crop amounts
% 22/02/16: Better support + backward compatibility for transparency (issue #108)
% 10/06/16: Fixed issue #159: text handles get cleared by Matlab in the print() command
% 12/06/16: Improved the fix for issue #159 (in the previous commit)
% 12/06/16: Fixed issue #158: transparent patch color in PDF/EPS
% 18/09/17: Fixed issue #194: incorrect fonts in EPS/PDF output
% 18/09/17: Fixed issue #195: relaxed too-tight cropping in EPS/PDF
% 14/11/17: Workaround for issue #211: dashed/dotted lines in 3D axes appear solid
% 15/11/17: Updated issue #211: only set SortMethod='ChildOrder' in HG2, and when it looks the same onscreen; support multiple figure axes
% 18/11/17: Fixed issue #225: transparent/translucent dashed/dotted lines appear solid in EPS/PDF
% 24/03/18: Fixed issue #239: black title meshes with temporary black background figure bgcolor, causing bad cropping
% 21/03/19: Improvement for issue #258: missing fonts in output EPS/PDF (still *NOT* fully solved)
% 21/03/19: Fixed issues #166,#251: Arial font is no longer replaced with Helvetica but rather treated as a non-standard user font
% 14/05/19: Made Helvetica the top default font-swap, replacing Courier
% 12/06/19: Issue #277: Enabled preservation of figure's PaperSize in output PDF/EPS file
%}

    options = {'-loose'};
    if nargin > 3
        options = [options varargin];
    elseif nargin < 3
        export_options = 0;
        if nargin < 2
            fig = gcf();
        end
    end

    % Retrieve padding, crop & font-swap values
    crop_amounts = nan(1,4);  % auto-crop all 4 sides by default
    if isstruct(export_options)
        try preserve_size = export_options.preserve_size; catch, preserve_size = false; end
        try fontswap      = export_options.fontswap;      catch, fontswap = true;       end
        try font_space    = export_options.font_space;    catch, font_space = '';       end
        font_space(2:end) = '';
        try bb_crop       = export_options.crop;          catch, bb_crop = 0;           end
        try crop_amounts  = export_options.crop_amounts;  catch,                        end
        try bb_padding    = export_options.bb_padding;    catch, bb_padding = 0;        end
        try renderer      = export_options.rendererStr;   catch, renderer = 'opengl';   end  % fix for issue #110
        if renderer(1)~='-',  renderer = ['-' renderer];  end
    else
        if numel(export_options) > 3  % preserve_size
            preserve_size = export_options(4);
        else
            preserve_size = false;
        end
        if numel(export_options) > 2  % font-swapping
            fontswap = export_options(3);
        else
            fontswap = true;
        end
        if numel(export_options) > 1  % cropping
            bb_crop = export_options(2);
        else
            bb_crop = 0;  % scalar value, so use default bb_crop value of 0
        end
        if numel(export_options) > 0  % padding
            bb_padding = export_options(1);
        else
            bb_padding = 0;
        end
        renderer = '-opengl';
        font_space = '';
    end

    % Construct the filename
    if numel(name) < 5 || ~strcmpi(name(end-3:end), '.eps')
        name = [name '.eps']; % Add the missing extension
    end

    % Set paper size
    old_pos_mode    = get(fig, 'PaperPositionMode');
    old_orientation = get(fig, 'PaperOrientation');
    old_paper_units = get(fig, 'PaperUnits');
    set(fig, 'PaperPositionMode','auto', 'PaperOrientation','portrait', 'PaperUnits','points');

    % Find all the used fonts in the figure
    font_handles = findall(fig, '-property', 'FontName');
    fonts = get(font_handles, 'FontName');
    if isempty(fonts)
        fonts = {};
    elseif ~iscell(fonts)
        fonts = {fonts};
    end

    % Map supported font aliases onto the correct name
    fontsl = lower(fonts);
    for a = 1:numel(fonts)
        f = fontsl{a};
        f(f==' ') = [];
        switch f
            case {'times', 'timesnewroman', 'times-roman'}
                fontsl{a} = 'times';
            %case {'arial', 'helvetica'}  % issues #166, #251
            %    fontsl{a} = 'helvetica';
            case {'newcenturyschoolbook', 'newcenturyschlbk'}
                fontsl{a} = 'newcenturyschlbk';
            otherwise
        end
    end
    fontslu = unique(fontsl);

    % Determine the font swap table
    if fontswap
        % Issue #258: Rearrange standard fonts list based on decending "problematicness"
        % The issue is still *NOT* fully solved because I cannot figure out how to force
        % the EPS postscript engine to look for the user's font on disk
        % Also see: https://stat.ethz.ch/pipermail/r-help/2005-January/064374.html
        matlab_fonts = {'Helvetica', 'Times', 'Courier', 'Symbol', 'ZapfDingbats', ...
                        'Palatino', 'Bookman', 'ZapfChancery', 'AvantGarde', ...
                        'NewCenturySchlbk', 'Helvetica-Narrow'};
        matlab_fontsl = lower(matlab_fonts);
        require_swap = find(~ismember(fontslu, matlab_fontsl));
        unused_fonts = find(~ismember(matlab_fontsl, fontslu));
        font_swap = cell(3, min(numel(require_swap), numel(unused_fonts)));
        fonts_new = fonts;
        for a = 1:size(font_swap, 2)
            font_swap{1,a} = find(strcmp(fontslu{require_swap(a)}, fontsl));
            font_swap{2,a} = matlab_fonts{unused_fonts(a)};
            font_swap{3,a} = fonts{font_swap{1,a}(1)};
            fonts_new(font_swap{1,a}) = font_swap(2,a);
        end
    else
        font_swap = [];
    end

    % Swap the fonts
    if ~isempty(font_swap)
        fonts_size = get(font_handles, 'FontSize');
        if iscell(fonts_size)
            fonts_size = cell2mat(fonts_size);
        end
        M = false(size(font_handles));

        % Loop because some changes may not stick first time, due to listeners
        c = 0;
        update = zeros(1000, 1);
        for b = 1:10 % Limit number of loops to avoid infinite loop case
            for a = 1:numel(M)
                M(a) = ~isequal(get(font_handles(a), 'FontName'), fonts_new{a}) || ~isequal(get(font_handles(a), 'FontSize'), fonts_size(a));
                if M(a)
                    set(font_handles(a), 'FontName', fonts_new{a}, 'FontSize', fonts_size(a));
                    c = c + 1;
                    update(c) = a;
                end
            end
            if ~any(M)
                break;
            end
        end

        % Compute the order to revert fonts later, without the need of a loop
        [update, M] = unique(update(1:c));
        [dummy, M] = sort(M); %#ok<ASGLU>
        update = reshape(update(M), 1, []);
    end

    % MATLAB bug fix - black and white text can come out inverted sometimes
    % Find the white and black text
    black_text_handles = findall(fig, 'Type', 'text', 'Color', [0 0 0]);
    white_text_handles = findall(fig, 'Type', 'text', 'Color', [1 1 1]);
    % Set the font colors slightly off their correct values
    set(black_text_handles, 'Color', [0 0 0] + eps);
    set(white_text_handles, 'Color', [1 1 1] - eps);

    % MATLAB bug fix - white lines can come out funny sometimes
    % Find the white lines
    white_line_handles = findall(fig, 'Type', 'line', 'Color', [1 1 1]);
    % Set the line color slightly off white
    set(white_line_handles, 'Color', [1 1 1] - 0.00001);

    % MATLAB bug fix (issue #211): dashed/dotted lines in 3D axes appear solid
    % Note: this "may limit other functionality in plotting such as hidden line/surface removal"
    % reference: Technical Support Case #02838114, https://mail.google.com/mail/u/0/#inbox/15fb7659f70e7bd8
    hAxes = findall(fig, 'Type', 'axes');
    if using_hg2 && ~isempty(hAxes)  % issue #211 presumably happens only in HG2, not HG1
        try
            % If there are any axes using SortMethod~='ChildOrder'
            oldSortMethods = get(hAxes,{'SortMethod'});  % use {'SortMethod'} to ensure we get a cell array, even for single axes
            if any(~strcmpi('ChildOrder',oldSortMethods))  % i.e., any oldSortMethods=='depth'
                % Check if the axes look visually different onscreen when SortMethod='ChildOrder'
                imgBefore = print2array(fig);
                set(hAxes,'SortMethod','ChildOrder');
                imgAfter  = print2array(fig);
                if isequal(imgBefore, imgAfter)
                    % They look the same, so use SortMethod='ChildOrder' when generating the EPS
                else
                    % They look different, so revert SortMethod and issue a warning message
                    warning('YMA:export_fig:issue211', ...
                            ['You seem to be using axes that have overlapping/hidden graphic elements. ' 10 ...
                             'Setting axes.SortMethod=''ChildOrder'' may solve potential problems in EPS/PDF export. ' 10 ...
                             'Additional info: https://github.com/altmany/export_fig/issues/211'])
                    set(hAxes,{'SortMethod'},oldSortMethods);
                end
            end
        catch err
            % ignore
            a=err;  %#ok<NASGU> % debug breakpoint
        end
    end

    % Workaround for issue #45: lines in image subplots are exported in invalid color
    % In this case the -depsc driver solves the problem, but then all the other workarounds
    % below (for all the other issues) will fail, so it's better to let the user decide by
    % just issuing a warning and accepting the '-depsc' input parameter
    epsLevel2 = ~any(strcmpi(options,'-depsc'));
    if epsLevel2
        % Use -depsc2 (EPS color level-2) if -depsc (EPS color level-3) was not specifically requested
        options{end+1} = '-depsc2';
        % Issue a warning if multiple images & lines were found in the figure, and HG1 with painters renderer is used
        isPainters = any(strcmpi(options,'-painters'));
        if isPainters && ~using_hg2 && numel(findall(fig,'Type','image'))>1 && ~isempty(findall(fig,'Type','line'))
            warning('YMA:export_fig:issue45', ...
                    ['Multiple images & lines detected. In such cases, the lines might \n' ...
                     'appear with an invalid color due to an internal MATLAB bug (fixed in R2014b). \n' ...
                     'Possible workaround: add a ''-depsc'' or ''-opengl'' parameter to the export_fig command.']);
        end
    end

    % Fix issue #83: use numeric handles in HG1
    if ~using_hg2(fig),  fig = double(fig);  end

    % Workaround for when transparency is lost through conversion fig>EPS>PDF (issue #108)
    % Replace transparent patch RGB values with an ID value (rare chance that ID color is being used already)
    if using_hg2
        origAlphaColors = eps_maintainAlpha(fig);
    end

    % Print to eps file
    print(fig, options{:}, name);

    % Restore the original axes SortMethods (if updated)
    try set(hAxes,{'SortMethod'},oldSortMethods); catch, end

    % Do post-processing on the eps file
    try
        % Read the EPS file into memory
        fstrm = read_write_entire_textfile(name);
    catch
        fstrm = '';
    end

    % Restore colors for transparent patches/lines and apply the
    % setopacityalpha setting in the EPS file (issue #108)
    if using_hg2
        [~,fstrm,foundFlags] = eps_maintainAlpha(fig, fstrm, origAlphaColors);

        % If some of the transparencies were not found in the EPS file, then rerun the
        % export with only the found transparencies modified (backward compatibility)
        if ~isempty(fstrm) && ~all(foundFlags)
            foundIdx = find(foundFlags);
            for objIdx = 1 : sum(foundFlags)
                colorsIdx = foundIdx(objIdx);
                colorsData = origAlphaColors{colorsIdx};
                hObj     = colorsData{1};
                propName = colorsData{2};
                newColor = colorsData{4};
                hObj.(propName).ColorData = newColor;
            end
            delete(name);
            print(fig, options{:}, name);
            fstrm = read_write_entire_textfile(name);
            [~,fstrm] = eps_maintainAlpha(fig, fstrm, origAlphaColors(foundFlags));
        end
    end

    % Bail out if EPS post-processing is not possible
    if isempty(fstrm)
        warning('Loading EPS file failed, so unable to perform post-processing. This is usually because the figure contains a large number of patch objects. Consider exporting to a bitmap format in this case.');
        return
    end

    % Fix for Matlab R2014b bug (issue #31): LineWidths<0.75 are not set in the EPS (default line width is used)
    try
        if using_hg2(fig)
            % Convert miter joins to line joins
            %fstrm = regexprep(fstrm, '\n10.0 ML\n', '\n1 LJ\n');
            % This is faster (the original regexprep could take many seconds when the axes contains many lines):
            fstrm = strrep(fstrm, sprintf('\n10.0 ML\n'), sprintf('\n1 LJ\n'));

            % In HG2, grid lines and axes Ruler Axles have a default LineWidth of 0.5 => replace en-bulk (assume that 1.0 LineWidth = 1.333 LW)
            %   hAxes=gca; hAxes.YGridHandle.LineWidth, hAxes.YRuler.Axle.LineWidth
            %fstrm = regexprep(fstrm, '(GC\n2 setlinecap\n1 LJ)\nN', '$1\n0.667 LW\nN');
            % This is faster:
            fstrm = strrep(fstrm, sprintf('GC\n2 setlinecap\n1 LJ\nN'), sprintf('GC\n2 setlinecap\n1 LJ\n0.667 LW\nN'));

            % This is more accurate but *MUCH* slower (issue #52)
            %{
            % Modify all thin lines in the figure to have 10x LineWidths
            hLines = findall(fig,'Type','line');
            hThinLines = [];
            for lineIdx = 1 : numel(hLines)
                thisLine = hLines(lineIdx);
                if thisLine.LineWidth < 0.75 && strcmpi(thisLine.Visible,'on')
                    hThinLines(end+1) = thisLine; %#ok<AGROW>
                    thisLine.LineWidth = thisLine.LineWidth * 10;
                end
            end

            % If any thin lines were found
            if ~isempty(hThinLines)
                % Prepare an EPS with large-enough line widths
                print(fig, options{:}, name);
                % Restore the original LineWidths in the figure
                for lineIdx = 1 : numel(hThinLines)
                    thisLine = handle(hThinLines(lineIdx));
                    thisLine.LineWidth = thisLine.LineWidth / 10;
                end

                % Compare the original and the new EPS files and correct the original stream's LineWidths
                fstrm_new = read_write_entire_textfile(name);
                idx = 500;  % skip heading with its possibly-different timestamp
                markerStr = sprintf('10.0 ML\nN');
                markerLen = length(markerStr);
                while ~isempty(idx) && idx < length(fstrm)
                    lastIdx = min(length(fstrm), length(fstrm_new));
                    delta = fstrm(idx+1:lastIdx) - fstrm_new(idx+1:lastIdx);
                    idx = idx + find(delta,1);
                    if ~isempty(idx) && ...
                            isequal(fstrm(idx-markerLen+1:idx), markerStr) && ...
                            ~isempty(regexp(fstrm_new(idx-markerLen+1:idx+12),'10.0 ML\n[\d\.]+ LW\nN')) %#ok<RGXP1>
                        value = str2double(regexprep(fstrm_new(idx:idx+12),' .*',''));
                        if isnan(value), break; end  % something's wrong... - bail out
                        newStr = sprintf('%0.3f LW\n',value/10);
                        fstrm = [fstrm(1:idx-1) newStr fstrm(idx:end)];
                        idx = idx + 12;
                    else
                        break;
                    end
                end
            end
            %}

            % This is much faster although less accurate: fix all non-gray lines to have a LineWidth of 0.75 (=1 LW)
            % Note: This will give incorrect LineWidth of 075 for lines having LineWidth<0.75, as well as for non-gray grid-lines (if present)
            %       However, in practice these edge-cases are very rare indeed, and the difference in LineWidth should not be noticeable
            %fstrm = regexprep(fstrm, '([CR]C\n2 setlinecap\n1 LJ)\nN', '$1\n1 LW\nN');
            % This is faster (the original regexprep could take many seconds when the axes contains many lines):
            fstrm = strrep(fstrm, sprintf('\n2 setlinecap\n1 LJ\nN'), sprintf('\n2 setlinecap\n1 LJ\n1 LW\nN'));
        end
    catch err
        fprintf(2, 'Error fixing LineWidths in EPS file: %s\n at %s:%d\n', err.message, err.stack(1).file, err.stack(1).line);
    end

    % Reset the font and line colors
    try
        set(black_text_handles, 'Color', [0 0 0]);
        set(white_text_handles, 'Color', [1 1 1]);
    catch
        % Fix issue #159: redo findall() '*text_handles'
        black_text_handles = findall(fig, 'Type', 'text', 'Color', [0 0 0]+eps);
        white_text_handles = findall(fig, 'Type', 'text', 'Color', [1 1 1]-eps);
        set(black_text_handles, 'Color', [0 0 0]);
        set(white_text_handles, 'Color', [1 1 1]);
    end
    set(white_line_handles, 'Color', [1 1 1]);

    % Preserve the figure's PaperSize in the output file, if requested (issue #277)
    if preserve_size
        paper_size = get(fig, 'PaperSize');  % in [points]
        fstrm = sprintf('<< /PageSize [%d %d] >> setpagedevice\n%s', paper_size, fstrm);
    end

    % Reset paper size
    set(fig, 'PaperPositionMode',old_pos_mode, 'PaperOrientation',old_orientation, 'PaperUnits',old_paper_units);

    % Reset the font names in the figure
    if ~isempty(font_swap)
        for a = update
            set(font_handles(a), 'FontName', fonts{a}, 'FontSize', fonts_size(a));
        end
    end

    % Replace the font names
    if ~isempty(font_swap)
        for a = 1:size(font_swap, 2)
            fontName = font_swap{3,a};
            %fontName = fontName(~isspace(font_swap{3,a}));
            if length(fontName) > 29
                warning('YMA:export_fig:font_name','Font name ''%s'' is longer than 29 characters. This might cause problems in some EPS/PDF readers. Consider using a different font.',fontName);
            end
            if isempty(font_space)
                fontName(fontName==' ') = '';
            else
                fontName(fontName==' ') = char(font_space);
            end

            % Replace all instances of the standard Matlab fonts with the original user's font names
            %fstrm = regexprep(fstrm, [font_swap{1,a} '-?[a-zA-Z]*\>'], fontName);
            %fstrm = regexprep(fstrm, [font_swap{2,a} '([ \n])'], [fontName '$1']);
            %fstrm = regexprep(fstrm, font_swap{2,a}, fontName);  % also replace -Bold, -Italic, -BoldItalic

            % Times-Roman's Bold/Italic fontnames don't include '-Roman'
            fstrm = regexprep(fstrm, [font_swap{2,a} '(\-Roman)?'], fontName);
        end
    end

    % Move the bounding box to the top of the file (HG2 only), or fix the line styles (HG1 only)
    if using_hg2(fig)
        % Move the bounding box to the top of the file (HG2 only)
        [s, e] = regexp(fstrm, '%%BoundingBox: [^%]*%%');
        if numel(s) == 2
            fstrm = fstrm([1:s(1)-1 s(2):e(2)-2 e(1)-1:s(2)-1 e(2)-1:end]);
        end
    else
        % Fix the line styles (HG1 only)
        fstrm = fix_lines(fstrm);
    end

    % Apply the bounding box padding & cropping, replacing Matlab's print()'s bounding box
    if bb_crop
        % Calculate a new bounding box based on a bitmap print using crop_border.m
        % 1. Determine the Matlab BoundingBox and PageBoundingBox
        [s,e] = regexp(fstrm, '%%BoundingBox: [^%]*%%'); % location BB in eps file
        if numel(s)==2, s=s(2); e=e(2); end
        aa = fstrm(s+15:e-3); % dimensions bb - STEP1
        bb_matlab = cell2mat(textscan(aa,'%f32%f32%f32%f32'));  % dimensions bb - STEP2

        [s,e] = regexp(fstrm, '%%PageBoundingBox: [^%]*%%'); % location bb in eps file
        if numel(s)==2, s=s(2); e=e(2); end
        aa = fstrm(s+19:e-3); % dimensions bb - STEP1
        pagebb_matlab = cell2mat(textscan(aa,'%f32%f32%f32%f32'));  % dimensions bb - STEP2

        % 1b. Fix issue #239: black title meshes with temporary black background figure bgcolor, causing bad cropping
        hTitles = [];
        if isequal(get(fig,'Color'),'none')
            hAxes = findall(fig,'type','axes');
            for idx = 1 : numel(hAxes)
                hAx = hAxes(idx);
                try
                    hTitle = hAx.Title;
                    oldColor = hTitle.Color;
                    if all(oldColor < 5*eps) || (ischar(oldColor) && lower(oldColor(1))=='k')
                        hTitles(end+1) = hTitle; %#ok<AGROW>
                        hTitle.Color = [0,0,.01];
                    end
                catch
                end
            end
        end

        % 2. Create a bitmap image and use crop_borders to create the relative
        %    bb with respect to the PageBoundingBox
        [A, bcol] = print2array(fig, 1, renderer);
        [aa, aa, aa, bb_rel] = crop_borders(A, bcol, bb_padding, crop_amounts); %#ok<ASGLU>

        try set(hTitles,'Color','k'); catch, end

        % 3. Calculate the new Bounding Box
        pagew = pagebb_matlab(3)-pagebb_matlab(1);
        pageh = pagebb_matlab(4)-pagebb_matlab(2);
        %bb_new = [pagebb_matlab(1)+pagew*bb_rel(1) pagebb_matlab(2)+pageh*bb_rel(2) ...
        %          pagebb_matlab(1)+pagew*bb_rel(3) pagebb_matlab(2)+pageh*bb_rel(4)];
        bb_new = pagebb_matlab([1,2,1,2]) + [pagew,pageh,pagew,pageh].*bb_rel;  % clearer
        bb_offset = (bb_new-bb_matlab) + [-2,-2,2,2];  % 2px margin so that cropping is not TOO tight (issue #195)

        % Apply the bounding box padding
        if bb_padding
            if abs(bb_padding)<1
                bb_padding = round((mean([bb_new(3)-bb_new(1) bb_new(4)-bb_new(2)])*bb_padding)/0.5)*0.5; % ADJUST BB_PADDING
            end
            add_padding = @(n1, n2, n3, n4) sprintf(' %.0f', str2double({n1, n2, n3, n4}) + bb_offset + bb_padding*[-1,-1,1,1]); %#ok<NASGU>
        else
            add_padding = @(n1, n2, n3, n4) sprintf(' %.0f', str2double({n1, n2, n3, n4}) + bb_offset); %#ok<NASGU> % fix small but noticeable bounding box shift
        end
        fstrm = regexprep(fstrm, '%%BoundingBox:[ ]+([-]?\d+)[ ]+([-]?\d+)[ ]+([-]?\d+)[ ]+([-]?\d+)', '%%BoundingBox:${add_padding($1, $2, $3, $4)}');
    end

    % Fix issue #44: white artifact lines appearing in patch exports
    % Note: the problem is due to the fact that Matlab's print() function exports patches
    %       as a combination of filled triangles, and a white line appears where the triangles touch
    % In the workaround below, we will modify such dual-triangles into a filled rectangle.
    % We are careful to only modify regexps that exactly match specific patterns - it's better to not
    % correct some white-line artifacts than to change the geometry of a patch, or to corrupt the EPS.
    %   e.g.: '0 -450 937 0 0 450 3 MP PP 937 0 0 -450 0 450 3 MP PP' => '0 -450 937 0 0 450 0 0 4 MP'
    fstrm = regexprep(fstrm, '\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\2 \1 \3 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n');
    fstrm = regexprep(fstrm, '\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\2 \3 \1 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n');
    fstrm = regexprep(fstrm, '\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\3 \1 \2 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n');
    fstrm = regexprep(fstrm, '\n([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) ([-\d.]+ [-\d.]+) 3 MP\nPP\n\3 \2 \1 3 MP\nPP\n','\n$1 $2 $3 0 0 4 MP\nPP\n');

    % Write out the fixed eps file
    read_write_entire_textfile(name, fstrm);
end

function [StoredColors, fstrm, foundFlags] = eps_maintainAlpha(fig, fstrm, StoredColors)
    if nargin == 1  % in: convert transparency in Matlab figure into unique RGB colors
        hObjs = findall(fig); %findobj(fig,'Type','Area');
        StoredColors = {};
        propNames = {'Face','Edge'};
        for objIdx = 1:length(hObjs)
            hObj = hObjs(objIdx);
            for propIdx = 1 : numel(propNames)
                try
                    propName = propNames{propIdx};
                    if strcmp(hObj.(propName).ColorType, 'truecoloralpha')
                        nColors = length(StoredColors);
                        oldColor = hObj.(propName).ColorData;
                        newColor = uint8([101; 102+floor(nColors/255); mod(nColors,255); 255]);
                        StoredColors{end+1} = {hObj, propName, oldColor, newColor}; %#ok<AGROW>
                        hObj.(propName).ColorData = newColor;
                    end
                catch
                    % Never mind - ignore (either doesn't have the property or cannot change it)
                end
            end
        end
    else  % restore transparency in Matlab figure by converting back from the unique RGBs
        %Find the transparent patches
        wasError = false;
        nColors = length(StoredColors);
        foundFlags = false(1,nColors);
        for objIdx = 1 : nColors
            colorsData = StoredColors{objIdx};
            hObj      = colorsData{1};
            propName  = colorsData{2};
            origColor = colorsData{3};
            newColor  = colorsData{4};
            try
                %Restore the EPS files patch color
                colorID   = num2str(round(double(newColor(1:3)') /255,3),'%.3g %.3g %.3g'); %ID for searching
                origRGB   = num2str(round(double(origColor(1:3)')/255,3),'%.3g %.3g %.3g'); %Replace with original color
                origAlpha = num2str(round(double(origColor(end)) /255,3),'%.3g'); %Convert alpha value for EPS

                %Find and replace the RGBA values within the EPS text fstrm
                if strcmpi(propName,'Face')
                    oldStr = sprintf(['\n' colorID ' RC\n']);  % ...N\n (removed to fix issue #225)
                    newStr = sprintf(['\n' origRGB ' RC\n' origAlpha ' .setopacityalpha true\n']);  % ...N\n
                else  %'Edge'
                    oldStr = sprintf(['\n' colorID ' RC\n']);  % ...1 LJ\n (removed to fix issue #225)
                    newStr = sprintf(['\n' origRGB ' RC\n' origAlpha ' .setopacityalpha true\n']);
                end
                foundFlags(objIdx) = ~isempty(strfind(fstrm, oldStr)); %#ok<STREMP>
                fstrm = strrep(fstrm, oldStr, newStr);

                %Restore the figure object's original color
                hObj.(propName).ColorData = origColor;
            catch err
                % something is wrong - cannot restore transparent color...
                if ~wasError
                    fprintf(2, 'Error maintaining transparency in EPS file: %s\n at %s:%d\n', err.message, err.stack(1).file, err.stack(1).line);
                    wasError = true;
                end
            end
        end
    end
end



%READ_WRITE_ENTIRE_TEXTFILE Read or write a whole text file to/from memory
%
% Read or write an entire text file to/from memory, without leaving the
% file open if an error occurs.
%
% Reading:
%   fstrm = read_write_entire_textfile(fname)
% Writing:
%   read_write_entire_textfile(fname, fstrm)
%
%IN:
%   fname - Pathname of text file to be read in.
%   fstrm - String to be written to the file, including carriage returns.
%
%OUT:
%   fstrm - String read from the file. If an fstrm input is given the
%           output is the same as that input. 

function fstrm = read_write_entire_textfile(fname, fstrm)
modes = {'rt', 'wt'};
writing = nargin > 1;
fh = fopen(fname, modes{1+writing});
if fh == -1
    error('Unable to open file %s.', fname);
end
try
    if writing
        fwrite(fh, fstrm, 'char*1');
    else
        fstrm = fread(fh, '*char')';
    end
catch ex
    fclose(fh);
    rethrow(ex);
end
fclose(fh);
end


function string = user_string(string_name, string)
%USER_STRING  Get/set a user specific string
%
% Examples:
%   string  = user_string(string_name)
%   isSaved = user_string(string_name, new_string)
%
% Function to get and set a string in a system or user specific file. This
% enables, for example, system specific paths to binaries to be saved.
%
% The specified string will be saved in a file named <string_name>.txt,
% either in a subfolder named .ignore under this file's folder, or in the
% user's prefdir folder (in case this file's folder is non-writable).
%
% IN:
%   string_name - String containing the name of the string required, which
%                 sets the filename storing the string: <string_name>.txt
%   new_string  - The new string to be saved in the <string_name>.txt file
%
% OUT:
%   string  - The currently saved string. Default: ''
%   isSaved - Boolean indicating whether the save was succesful

% Copyright (C) Oliver Woodford 2011-2014, Yair Altman 2015-

% This method of saving paths avoids changing .m files which might be in a
% version control system. Instead it saves the user dependent paths in
% separate files with a .txt extension, which need not be checked in to
% the version control system. Thank you to Jonas Dorn for suggesting this
% approach.

% 10/01/2013 - Access files in text, not binary mode, as latter can cause
%              errors. Thanks to Christian for pointing this out.
% 29/05/2015 - Save file in prefdir if current folder is non-writable (issue #74)
% 09/01/2018 - Fix issue #232: if the string looks like a file/folder path, ensure it actually exists

    if ~ischar(string_name)
        error('string_name must be a string.');
    end
    % Create the full filename
    fname = [string_name '.txt'];
    dname = fullfile(fileparts(mfilename('fullpath')), '.ignore');
    file_name = fullfile(dname, fname);
    if nargin > 1
        % Set string
        if ~ischar(string)
            error('new_string must be a string.');
        end
        % Make sure the save directory exists
        %dname = fileparts(file_name);
        if ~exist(dname, 'dir')
            % Create the directory
            try
                if ~mkdir(dname)
                    string = false;
                    return
                end
            catch
                string = false;
                return
            end
            % Make it hidden
            try
                fileattrib(dname, '+h');
            catch
            end
        end
        % Write the file
        fid = fopen(file_name, 'wt');
        if fid == -1
            % file cannot be created/updated - use prefdir if file does not already exist
            % (if file exists but is simply not writable, don't create a duplicate in prefdir)
            if ~exist(file_name,'file')
                file_name = fullfile(prefdir, fname);
                fid = fopen(file_name, 'wt');
            end
            if fid == -1
                string = false;
                return;
            end
        end
        try
            fprintf(fid, '%s', string);
        catch
            fclose(fid);
            string = false;
            return
        end
        fclose(fid);
        string = true;
    else
        % Get string
        fid = fopen(file_name, 'rt');
        if fid == -1
            % file cannot be read, try to read the file in prefdir
            file_name = fullfile(prefdir, fname);
            fid = fopen(file_name, 'rt');
            if fid == -1
                string = '';
                return
            end
        end
        string = fgetl(fid);
        fclose(fid);

        % Fix issue #232: if the string looks like a file/folder path, ensure it actually exists
        if ~isempty(string) && any(string=='\' | string=='/') && ~exist(string) %#ok<EXIST>
            string = '';
        end
    end
end


%USING_HG2 Determine if the HG2 graphics engine is used
%
%   tf = using_hg2(fig)
%
%IN:
%   fig - handle to the figure in question.
%
%OUT:
%   tf - boolean indicating whether the HG2 graphics engine is being used
%        (true) or not (false).

% 19/06/2015 - Suppress warning in R2015b; cache result for improved performance
% 06/06/2016 - Fixed issue #156 (bad return value in R2016b)

function tf = using_hg2(fig)
    persistent tf_cached
    if isempty(tf_cached)
        try
            if nargin < 1,  fig = figure('visible','off');  end
            oldWarn = warning('off','MATLAB:graphicsversion:GraphicsVersionRemoval');
            try
                % This generates a [supressed] warning in R2015b:
                tf = ~graphicsversion(fig, 'handlegraphics');
            catch
                tf = ~verLessThan('matlab','8.4');  % =R2014b
            end
            warning(oldWarn);
        catch
            tf = false;
        end
        if nargin < 1,  delete(fig);  end
        tf_cached = tf;
    else
        tf = tf_cached;
    end
end
