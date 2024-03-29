% % AD9361 Filter Wizard
%
% Name-Value Pair Arguments
%
% Specify optional comma-separated pairs of Name,Value arguments. Name is
% the argument name and Value is the corresponding value. Name must
% appear inside single quotes (' '). You can specify several name and
% value pair arguments in any order as Name1,Value1,...,NameN,ValueN.
% Example: 'remote','192.168.0.1','MarkerFaceColor','red'
%
%   Name       Value
%   'remote'          'hide' | 'IP_number' | 'IP_number:port'
%   'PathConfig'      'rx' | 'tx' | 'either'
%   'ApplyString'     'Save String'
%   'ApplyCallback'   'Callback for Apply button'
%   'HelpCallback'    'Callback for Help button'
%   'DefaultRxVals'   Structure of default values (in Hz)
%   'DefaultTxVals'   Structure of default values (in Hz)
%
%%
% Copyright 2014(c) Analog Devices, Inc.
%
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without modification,
%  are permitted provided that the following conditions are met:
%      - Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      - Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the
%        distribution.
%      - Neither the name of Analog Devices, Inc. nor the names of its
%        contributors may be used to endorse or promote products derived
%        from this software without specific prior written permission.
%      - The use of this software may or may not infringe the patent rights
%        of one or more patent holders.  This license does not release you
%        from the requirement that you obtain separate licenses from these
%        patent holders to use this software.
%      - Use of the software either in source or binary form or filter designs
%        resulting from the use of this software, must be connected to, run
%        on or loaded to an Analog Devices Inc. component.
%
%  THIS SOFTWARE IS PROVIDED BY ANALOG DEVICES "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
%  INCLUDING, BUT NOT LIMITED TO, NON-INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A
%  PARTICULAR PURPOSE ARE DISCLAIMED.
%
%  IN NO EVENT SHALL ANALOG DEVICES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, INTELLECTUAL PROPERTY
%  RIGHTS, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
%  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
%  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function varargout = AD9361_Filter_Wizard(varargin)
% AD9361_FILTER_WIZARD MATLAB code for AD9361_Filter_Wizard.fig
%      AD9361_FILTER_WIZARD, by itself, creates a new AD9361_FILTER_WIZARD or raises the existing
%      singleton*.
%
%      H = AD9361_FILTER_WIZARD returns the handle to a new AD9361_FILTER_WIZARD or the handle to
%      the existing singleton*.
%
%      AD9361_FILTER_WIZARD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AD9361_FILTER_WIZARD.M with the given input arguments.
%
%      AD9361_FILTER_WIZARD('Property','Value',...) creates a new AD9361_FILTER_WIZARD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AD9361_Filter_Wizard_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AD9361_Filter_Wizard_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AD9361_Filter_Wizard

% Last Modified by GUIDE v2.5 19-Jan-2017 14:03:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AD9361_Filter_Wizard_OpeningFcn, ...
    'gui_OutputFcn',  @AD9361_Filter_Wizard_OutputFcn, ...
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

% --- Executes just before AD9361_Filter_Wizard is made visible.
function AD9361_Filter_Wizard_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AD9361_Filter_Wizard (see VARARGIN)

% Choose default command line output for AD9361_Filter_Wizard
handles.output = hObject;

% check MATLAB version compatibility
if ~check_compat()
    warndlg('AD9361 Filter Wizard requires MATLAB 2015b or later, older versions are unsupported.')
end

% import various value bounds
bounds9361 = rate_bounds;
handles.bounds = bounds9361;

new = 0;
handles.freq_units = 3;

handles.applycallback = {};

% inputs need to be name/value _pairs_
if rem(length(varargin),2)
    error('myApp:argChk', 'Wrong number of input arguments')
end

for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'remote')
        % 'remote'      'ip_number' | 'none'
        if (strcmpi(varargin{i + 1}, 'hide')) || isempty(varargin{i + 1})
            set(handles.target_remote, 'Visible', 'off');
        else
            set(handles.IP_num, 'String', varargin{i+ 1});
        end
    elseif strcmpi(varargin{i}, 'PathConfig')
        % 'PathConfig'     'rx' | 'tx' | 'either'
        if (strcmpi(varargin{i + 1}, 'rx'))
            set(handles.filter_type, 'Value', 1);
            set(handles.filter_type, 'Enable', 'off');
        elseif (strcmpi(varargin{i + 1}, 'tx'))
            set(handles.filter_type, 'Value', 2);
            set(handles.filter_type, 'Enable', 'off');
        elseif (strcmpi(varargin{i + 1}, 'either'))
            set(handles.filter_type, 'Value', 1);
            set(handles.filter_type, 'Enable', 'on');
        else
            error('Unknown value to "PathConfig"');
        end
        %filter_type_Callback(handles.filter_type, eventdata, handles);
    elseif strcmpi(varargin{i}, 'ApplyString')
        % 'ApplyString'  'Save String'
        set(handles.save2workspace, 'String', varargin{i + 1});
    elseif strcmpi(varargin{i}, 'DefaultRxVals')
        % 'DefaultRxVals'     'structure (in Hz)'
        input = varargin{i + 1};
        input.RxTx = 'Rx';
        handles.rx = cook_input(input);
        new = 1;
    elseif strcmpi(varargin{i}, 'DefaultTxVals')
        % 'DefaultTxVals'     'structure (in Hz)'
        input = varargin{i + 1};
        input.RxTx = 'Tx';
        handles.tx = cook_input(input);
        new = 1;
    elseif strcmpi(varargin{i}, 'ApplyCallback')
        handles.applycallback = str2func(varargin{i + 1});
    elseif strcmpi(varargin{i}, 'HelpCallback')
        handles.helpcallback = str2func(varargin{i + 1});
    elseif strcmpi(varargin{i}, 'CallbackObj')
        handles.callbackObj = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'SetDevice')
        % Set default device
        switch lower(varargin{i + 1})
            case 'ad9361'
                handles.device.Children(15).Value = 1;                        
            case 'ad9364'
                handles.device.Children(15).Value = 2;        
            case 'pluto'
                handles.device.Children(15).Value = 3; 
            otherwise
                error('Unknown default device');
        end
    else
        error('Unknown input to function');
    end
end

handles = autoselect_rates(handles);

if isfield(handles, 'tx') || isfield(handles, 'rx')
    set(handles.store_filter, 'Visible', 'off');
    guidata(hObject, handles);
    data2gui(hObject, handles);
else
    handles.rx = 0;
    handles.tx = 0;
end

handles.Original_Size = get(handles.AD9361_Filter_app, 'Position');

% set the window name
window_name = sprintf('AD9361 Filter Wizard %s', get_version);
set(handles.AD9361_Filter_app, 'Name', window_name);

guidata(hObject, handles);

axes(handles.ADI_logo);
pict = imread('Analog_Devices_Logo.png');
image(pict);
axis image;
box off;
set(handles.ADI_logo, 'XTickLabel', []);
set(handles.ADI_logo, 'YTickLabel', []);
set(handles.ADI_logo, 'XTick', []);
set(handles.ADI_logo, 'YTick', []);
set(handles.ADI_logo, 'Box', 'off');
set(handles.ADI_logo, 'HandleVisibility', 'off');

% Set the defaults
set(handles.Use_FIR, 'Value', 1);
set(handles.use_FPGAfilter, 'Value', 0);
hide_advanced(handles);

% restore previously used IP address
[pathstr, ~, ~] = fileparts(mfilename('fullpath'));
cached_ip_file = fullfile(pathstr, '.previous_ip_addr');
if exist(cached_ip_file, 'file')
    fd = fopen(cached_ip_file, 'rt');
    set(handles.IP_num, 'String', fgets(fd));
    fclose(fd);
end

axes(handles.magnitude_plot);

handles.libiio_ctrl_dev = {};

for i = 1:4
    handles.arrows{i} = annotation('arrow');
    set(handles.arrows{i}, 'Visible', 'off');
end

% Update handles structure
guidata(hObject, handles);

if ~new
    reset_input(hObject, handles);
    load_settings(hObject, handles);
else
    reset_input(hObject, handles);
    dirty(hObject, handles);
    handles.active_plot = 0;
    display_default_image(hObject);
end
handles = guidata(hObject);

handles.freq_units = get(handles.Freq_units, 'Value');
handles.active_plot = 0;

% set(zoom(gca),'ActionPostCallback',@(x,y) zoom_axis(gca));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AD9361_Filter_Wizard wait for user response (see UIRESUME)
uiwait(handles.AD9361_Filter_app);

% --- Outputs from this function are returned to the command line.
function varargout = AD9361_Filter_Wizard_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

channels = {'rx' 'tx'};
numarg = 1;

% Return empty structs if the GUI was closed before a filter was designed.
[varargout{1:nargout}] = deal(struct);

for i = 1:length(channels)
    output = getfield(handles, char(channels(i)));
    if isfield(output, 'firtaps')
        varargout{numarg} = clean_filter_fields(output);
        numarg = numarg + 1;
    end
end

% The main figure can be deleted now
delete(handles.AD9361_Filter_app);


% --- Executes on selection change in Freq_units.
function Freq_units_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
units = get(hObject, 'Value');

if (handles.freq_units ~= units)
    fstop = value2Hz(handles, handles.freq_units, str2double(get(handles.Fstop, 'String')));
    fpass = value2Hz(handles, handles.freq_units, str2double(get(handles.Fpass, 'String')));
    fcutoff = value2Hz(handles, handles.freq_units, str2double(get(handles.Fcutoff, 'String')));
    data_rate = value2Hz(handles, handles.freq_units, str2double(get(handles.data_clk, 'String')));
    rf_bandwidth = value2Hz(handles, handles.freq_units, str2double(get(handles.RFbw, 'String')));
    fir_rate = value2Hz(handles, handles.freq_units, str2double(get(handles.FIR_rate, 'String')));
    hb1_rate = value2Hz(handles, handles.freq_units, str2double(get(handles.HB1_rate, 'String')));
    hb2_rate = value2Hz(handles, handles.freq_units, str2double(get(handles.HB2_rate, 'String')));
    hb3_rate = value2Hz(handles, handles.freq_units, str2double(get(handles.HB3_rate, 'String')));
    fpga_rate = value2Hz(handles, handles.freq_units, str2double(get(handles.FPGA_rate, 'String')));
    
    handles.freq_units = units;
    set(handles.Fstop, 'String', num2str(Hz2value(handles, handles.freq_units, fstop)));
    set(handles.Fpass, 'String', num2str(Hz2value(handles, handles.freq_units, fpass)));
    set(handles.Fcutoff, 'String', num2str(Hz2value(handles, handles.freq_units, fcutoff)));
    set(handles.data_clk, 'String', num2str(Hz2value(handles, handles.freq_units, data_rate)));
    set(handles.RFbw, 'String', num2str(Hz2value(handles, handles.freq_units, rf_bandwidth)));
    set(handles.FPGA_rate, 'String', num2str(Hz2value(handles, handles.freq_units, fpga_rate)));
    set(handles.FIR_rate, 'String', num2str(Hz2value(handles, handles.freq_units, fir_rate)));
    set(handles.HB1_rate, 'String', num2str(Hz2value(handles, handles.freq_units, hb1_rate)));
    set(handles.HB2_rate, 'String', num2str(Hz2value(handles, handles.freq_units, hb2_rate)));
    set(handles.HB3_rate, 'String', num2str(Hz2value(handles, handles.freq_units, hb3_rate)));
    set(handles.HBs_rate, 'String', get(handles.HB3_rate, 'String'));
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function Freq_units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dirty(hObject, handles)
handles.active_plot = 0;
plot_buttons_off(handles);
Filters = get(handles.Saved_Filters, 'String');
filter_selection = 0;
for i = 1:length(Filters)
    if strcmp(Filters{i}, 'Custom')
        filter_selection = i;
        break;
    end
end

if filter_selection == 0
    filter_selection = length(Filters) + 1;
    Filters{filter_selection} = 'Custom';
    set(handles.Saved_Filters, 'String', Filters);
end
set(handles.Saved_Filters, 'Value', filter_selection);
set(handles.store_filter, 'Visible', 'on');

function Fpass_Callback(hObject, eventdata, handles)
% hObject    handle to Fpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

Fpass = value2Hz(handles, handles.freq_units, str2double(get(hObject,'String')));
if get(handles.filter_type, 'Value') == 1
    handles.rx.Fpass = Fpass;
else
    handles.tx.Fpass = Fpass;
end

data2gui(hObject, handles);
handles = guidata(hObject);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Fpass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Fstop_Callback(hObject, eventdata, handles)
% hObject    handle to Fstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

Fstop = value2Hz(handles, handles.freq_units, str2double(get(hObject,'String')));
if get(handles.filter_type, 'Value') == 1
    handles.rx.Fstop = Fstop;
else
    handles.tx.Fstop = Fstop;
end

set(handles.Fcutoff, 'String', '0');

% reset caldiv so the advanced mode setting is respected
caldiv = get_caldiv(handles);
if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

% recalculate fcutoff for advanced mode
if get(handles.Advanced_options, 'Value')
    set_fcutoff(handles, caldiv);
end

data2gui(hObject, handles);
handles = guidata(hObject);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Fstop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Pll_rate_Callback(hObject, eventdata, handles)
% hObject    handle to Pll_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.active_plot ~= 0)
    handles.active_plot = 0;
end
plot_buttons_off(handles);

% --- Executes during object creation, after setting all properties.
function Pll_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pll_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interpolate = converter_interp(handles)
interpolate = cellstr(get(handles.converter2PLL, 'String'));
if isempty(interpolate{1})
    error('For some reason converter2PLL is empty. This will cause me to crash');
    dbstack('-completenames');
end
interpolate = char(interpolate(get(handles.converter2PLL, 'Value')));
interpolate = str2double(interpolate(1:2));

function handles = autoselect_rates(handles)
% sanity check the PLL rate and DAC divider values and alter them if necessary
if isfield(handles, 'tx') && isfield(handles, 'rx')
    if (handles.rx.PLL_rate ~= handles.tx.PLL_rate)
        % If dec3 is used by Rx or Tx, both must use it in order for the
        % PLL rates to match.
        if handles.tx.HB3 == 3 || handles.rx.HB3 == 3
            handles.rx.HB3 = 3;
            handles.tx.HB3 = 3;
        end
        
        ADC_rate = handles.rx.Rdata * handles.rx.FIR * ...
            handles.rx.HB1 * handles.rx.HB2 * handles.rx.HB3;
        DAC_rate = handles.tx.Rdata * handles.tx.FIR * ...
            handles.tx.HB1 * handles.tx.HB2 * handles.tx.HB3;
        DAC_div = ADC_rate / DAC_rate;
        if (handles.tx.DAC_div ~= DAC_div)
            if (DAC_div == 1 || DAC_div == 2)
                handles.tx.DAC_div = DAC_div;
                handles.tx.PLL_mult = handles.rx.PLL_mult;
                filter_type = get(handles.filter_type, 'Value');
                set(handles.filter_type, 'Value', 0);
                handles.tx.caldiv = default_caldiv(handles);
                set(handles.filter_type, 'Value', filter_type);
            end
        end
        
        handles.rx.PLL_mult = fastest_FIR([64 32 16 8 4 2 1], handles.bounds.MAX_BBPLL_FREQ, handles.bounds.MIN_BBPLL_FREQ, ...
            handles.rx.Rdata * handles.rx.FIR * handles.rx.HB1 * handles.rx.HB2 * handles.rx.HB3 * handles.rx.DAC_div);
        handles.tx.PLL_mult = handles.rx.PLL_mult;
        
        if handles.rx.PLL_mult > 64
            X = ['Date rate = ', num2str(tohwTx.TXSAMP), ' Hz. Tx BBPLL is too high for Rx to match.'];
            disp(X);
        end
        
        handles.rx.PLL_rate = handles.rx.Rdata * handles.rx.FIR * handles.rx.HB1 * ...
            handles.rx.HB2 * handles.rx.HB3 * handles.rx.PLL_mult;
    else
        ADC_rate = handles.rx.Rdata * handles.rx.FIR * ...
            handles.rx.HB1 * handles.rx.HB2 * handles.rx.HB3;
        DAC_rate = handles.tx.Rdata * handles.tx.FIR * ...
            handles.tx.HB1 * handles.tx.HB2 * handles.tx.HB3;
        DAC_div = ADC_rate / DAC_rate;
        if (handles.tx.DAC_div ~= DAC_div)
            if (DAC_div == 1 || DAC_div == 2)
                handles.tx.DAC_div = DAC_div;
                handles.tx.PLL_mult = handles.rx.PLL_mult;
                filter_type = get(handles.filter_type, 'Value');
                set(handles.filter_type, 'Value', 0);
                handles.tx.caldiv = default_caldiv(handles);
                set(handles.filter_type, 'Value', filter_type);
            end
        end
    end
end

function sel = get_current_rxtx(handles)
if (get(handles.filter_type, 'Value') == 1)
    % receive
    sel = handles.rx;
else
    % transmit
    sel = handles.tx;
end

function data_clk_Callback(hObject, eventdata, handles)
% hObject    handle to data_clk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

data_rate = value2Hz(handles, handles.freq_units, str2double(get(hObject,'String')));

if get(handles.which_device,'Value')==3 && data_rate<=520.83e3
    set(handles.use_FPGAfilter, 'Value', 1);
    set(handles.FPGA_label, 'Visible', 'on');
    set(handles.use_FPGAfilter, 'Visible', 'on');
    set(handles.FPGA_rate, 'Visible', 'on');
    freq_off(handles);
    data_rate = data_rate*8;
else
    set(handles.use_FPGAfilter, 'Value', 0);
    set(handles.FPGA_label, 'Visible', 'off');
    set(handles.use_FPGAfilter, 'Visible', 'off');
    set(handles.FPGA_rate, 'Visible', 'off');
    freq_on(handles);
end

input = {};
input.Rdata = data_rate;

input.RxTx = 'Rx';
handles.rx = cook_input(input);
input.RxTx = 'Tx';
handles.tx = cook_input(input);
handles = autoselect_rates(handles);

% reset fcutoff to the default value for advanced mode
if get(handles.Advanced_options, 'Value')
    set(handles.Fcutoff, 'String', '0');
    caldiv = get_caldiv(handles);
    set_fcutoff(handles, caldiv);
end

data2gui(hObject, handles);
handles = guidata(hObject);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function data_clk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_clk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Astop_Callback(hObject, eventdata, handles)
% hObject    handle to Astop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dirty(hObject, handles);
handles = guidata(hObject);

Astop = str2double(get(hObject, 'String'));
if get(handles.filter_type, 'Value') == 1
    handles.rx.Astop = Astop;
else
    handles.tx.Astop = Astop;
end

if get(handles.FIR_Astop, 'Value') >= str2double(get(hObject,'String'))
    set(handles.FIR_Astop, 'Value', str2double(get(hObject,'String')));
end

data2gui(hObject, handles);
handles = guidata(hObject);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Astop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Astop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Apass_Callback(hObject, eventdata, handles)
% hObject    handle to Apass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if str2double(get(hObject,'String')) == 0
%     set(hObject,'String', '0.00001');
% end

dirty(hObject, handles);
handles = guidata(hObject);

Apass = str2double(get(hObject, 'String'));
if get(handles.filter_type, 'Value') == 1
    handles.rx.Apass = Apass;
else
    handles.tx.Apass = Apass;
end

data2gui(hObject, handles);
handles = guidata(hObject);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Apass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Apass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in filter_type.
function filter_type_Callback(hObject, eventdata, handles)
% hObject    handle to filter_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data2gui(hObject, handles);


% --- Executes during object creation, after setting all properties.
function filter_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Determine frequencies at different points in the filter path.
function [rx_FIR_rate, rx_HB1_rate, rx_HB2_rate, rx_HB3_rate, rx_RFbw_hw, ...
    tx_FIR_rate, tx_HB1_rate, tx_HB2_rate, tx_HB3_rate, tx_RFbw_hw] = get_path_rates(handles)
rx_FIR_rate = handles.rx.Rdata * handles.rx.FIR;
rx_HB1_rate = rx_FIR_rate * handles.rx.HB1;
rx_HB2_rate = rx_HB1_rate * handles.rx.HB2;
rx_HB3_rate = rx_HB2_rate * handles.rx.HB3;
tx_FIR_rate = handles.tx.Rdata * handles.tx.FIR;
tx_HB1_rate = tx_FIR_rate * handles.tx.HB1;
tx_HB2_rate = tx_HB1_rate * handles.tx.HB2;
tx_HB3_rate = tx_HB2_rate * handles.tx.HB3;

filter_type = get(handles.filter_type, 'Value');
set(handles.filter_type, 'Value', 1);
sel = get_current_rxtx(handles);
rx_RFbw_hw = get_rfbw_hw(handles, sel.caldiv);
set(handles.filter_type, 'Value', 2);
sel = get_current_rxtx(handles);
tx_RFbw_hw = get_rfbw_hw(handles, sel.caldiv);
set(handles.filter_type, 'Value', filter_type);

% --- Executes on button press in save2noos.
function save2noos_Callback(hObject, eventdata, handles)
% hObject    handle to save2noos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,path] = uiputfile('*.c', 'Save coefficients as');
if filename == 0
    return;
else
    newpath = strcat(path, filename);
end

fid = fopen(newpath,'w');

PLL_rate = str2double(get(handles.Pll_rate, 'String')) * 1e6;
[rx_FIR_rate, rx_HB1_rate, rx_HB2_rate, rx_HB3_rate, rx_RFbw_hw, ...
    tx_FIR_rate, tx_HB1_rate, tx_HB2_rate, tx_HB3_rate, tx_RFbw_hw] = get_path_rates(handles);

fprintf(fid, '// Generated with AD9361 Filter Design Wizard %s\n', get_version);
fprintf(fid, '// MATLAB %s, %s\n', version(), datestr(now()));
fprintf(fid, '// Inputs:\n');

fprintf(fid, '// Data Sample Frequency = %.0f Hz\n', handles.rx.Rdata);
if get(handles.phase_eq, 'Value')
    fprintf(fid, '// RX Phase equalization = %f ns\n', handles.rx.phEQ);
    fprintf(fid, '// TX Phase equalization = %f ns\n', handles.tx.phEQ);
end

% Rx
fprintf(fid, '\nAD9361_RXFIRConfig rx_fir_config = {\n');
fprintf(fid, '\t3, // rx\n');
fprintf(fid, '\t%d, // rx_gain\n', handles.rx.gain);
fprintf(fid, '\t%d, // rx_dec\n', handles.rx.FIR);
coefficients = sprintf('%.0f,', flip(rot90(handles.rfirtaps)));
coefficients = coefficients(1:end-1); % strip final comma
fprintf(fid, '\t{%s}, // rx_coef[128]\n', coefficients);
fprintf(fid, '\t%d, // rx_coef_size\n', handles.nfirtaps);
fprintf(fid, '\t{%.0f,%.0f,%.0f,%.0f,%.0f,%.0f}, // rx_path_clks[6]\n', ...
    PLL_rate, rx_HB3_rate, rx_HB2_rate, rx_HB1_rate, rx_FIR_rate, handles.rx.Rdata);
fprintf(fid, '\t%.0f // rx_bandwidth\n', rx_RFbw_hw);
fprintf(fid, '};\n\n');

% Tx
fprintf(fid, 'AD9361_TXFIRConfig tx_fir_config = {\n');
fprintf(fid, '\t3, // tx\n');
fprintf(fid, '\t%d, // tx_gain\n', handles.tx.gain);
fprintf(fid, '\t%d, // tx_int\n', handles.tx.FIR);
coefficients = sprintf('%.0f,', flip(rot90(handles.tfirtaps)));
coefficients = coefficients(1:end-1); % strip final comma
fprintf(fid, '\t{%s}, // tx_coef[128]\n', coefficients);
fprintf(fid, '\t%d, // tx_coef_size\n', handles.nfirtaps);
fprintf(fid, '\t{%.0f,%.0f,%.0f,%.0f,%.0f,%.0f}, // tx_path_clks[6]\n', ...
    PLL_rate, tx_HB3_rate, tx_HB2_rate, tx_HB1_rate, tx_FIR_rate, handles.tx.Rdata);
fprintf(fid, '\t%.0f // tx_bandwidth\n', tx_RFbw_hw);
fprintf(fid, '};\n');

fclose(fid);


% --- Executes on button press in save2coefficients.
function save2coefficients_Callback(hObject, eventdata, handles)
% hObject    handle to save2coefficients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,path] = uiputfile('*.ftr', 'Save coefficients as');
if filename == 0
    return;
else
    newpath = strcat(path,filename);
end

fid = fopen(newpath,'w');

fprintf(fid, '# Generated with AD9361 Filter Design Wizard %s\r\n', get_version);
fprintf(fid, '# MATLAB %s, %s\r\n', version(), datestr(now()));
fprintf(fid, '# Inputs:\r\n');

%FIXME
%converter_rate = get_converter_clk(handles);
%converter_rate = get_ADC_clk(handles);

PLL_rate = str2double(get(handles.Pll_rate, 'String')) * 1e6;
[rx_FIR_rate, rx_HB1_rate, rx_HB2_rate, rx_HB3_rate, rx_RFbw_hw, ...
    tx_FIR_rate, tx_HB1_rate, tx_HB2_rate, tx_HB3_rate, tx_RFbw_hw] = get_path_rates(handles);

%fprintf(fid, '# PLL CLK Frequency = %f Hz\r\n', pll_rate);
%fprintf(fid, '# Converter Sample Frequency = %f Hz\r\n', converter_rate);
fprintf(fid, '# Data Sample Frequency = %.0f Hz\r\n', handles.rx.Rdata);
if get(handles.phase_eq, 'Value')
    fprintf(fid, '# RX Phase equalization = %f ns\r\n', handles.rx.phEQ);
    fprintf(fid, '# TX Phase equalization = %f ns\r\n', handles.tx.phEQ);
end
fprintf(fid, 'TX 3 GAIN %d INT %d\r\n', handles.tx.gain, handles.tx.FIR);
fprintf(fid, 'RX 3 GAIN %d DEC %d\r\n', handles.rx.gain, handles.rx.FIR);
fprintf(fid, 'RTX %.0f %.0f %.0f %.0f %.0f %.0f\r\n', PLL_rate, tx_HB3_rate, tx_HB2_rate, tx_HB1_rate, tx_FIR_rate, handles.tx.Rdata);
fprintf(fid, 'RRX %.0f %.0f %.0f %.0f %.0f %.0f\r\n', PLL_rate, rx_HB3_rate, rx_HB2_rate, rx_HB1_rate, rx_FIR_rate, handles.rx.Rdata);
fprintf(fid, 'BWTX %.0f\r\n', tx_RFbw_hw);
fprintf(fid, 'BWRX %.0f\r\n', rx_RFbw_hw);

% concat and transform Rx and Tx coefficient matrices for output
coefficients = flip(rot90(vertcat(handles.tfirtaps, handles.rfirtaps)));

% output all non-zero coefficients since they're padded to 128 with zeros
for i = 1:handles.nfirtaps
    fprintf(fid, '%d,%d\r\n', coefficients(i,:));
end

fclose(fid);


% --- Executes on button press in save2target.
function save2target_Callback(hObject, eventdata, handles)
% hObject    handle to save2target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PLL_rate = str2double(get(handles.Pll_rate, 'String')) * 1e6;
[rx_FIR_rate, rx_HB1_rate, rx_HB2_rate, rx_HB3_rate, rx_RFbw_hw, ...
    tx_FIR_rate, tx_HB1_rate, tx_HB2_rate, tx_HB3_rate, tx_RFbw_hw] = get_path_rates(handles);

fir_filter_str = sprintf('TX 3 GAIN %d INT %d', handles.tx.gain, handles.tx.FIR);
fir_filter_str = strcat(fir_filter_str, sprintf('\nRX 3 GAIN %d DEC %d', handles.rx.gain, handles.rx.FIR));
fir_filter_str = strcat(fir_filter_str, sprintf('\nRTX %.0f %.0f %.0f %.0f %.0f %.0f', PLL_rate, tx_HB3_rate, tx_HB2_rate, tx_HB1_rate, tx_FIR_rate, handles.tx.Rdata));
fir_filter_str = strcat(fir_filter_str, sprintf('\nRRX %.0f %.0f %.0f %.0f %.0f %.0f', PLL_rate, rx_HB3_rate, rx_HB2_rate, rx_HB1_rate, rx_FIR_rate, handles.rx.Rdata));
fir_filter_str = strcat(fir_filter_str, sprintf('\nBWTX %.0f', tx_RFbw_hw));
fir_filter_str = strcat(fir_filter_str, sprintf('\nBWRX %.0f', rx_RFbw_hw));

% concat and transform Rx and Tx coefficient matrices for outputting
coefficients = flip(rot90(vertcat(handles.tfirtaps, handles.rfirtaps)));

% output all non-zero coefficients since they're padded to 128 with zeros
for i = 1:handles.nfirtaps
    fir_filter_str = strcat(fir_filter_str, sprintf('\n%d,%d', coefficients(i,:)));
end

% write FIR filter to target
ret = writeAttributeString(handles.libiio_ctrl_dev, 'filter_fir_config', fir_filter_str);
if (ret < 0)
    msgbox('Could not write FIR filter to target!', 'Error', 'error');
    return;
end

% write Rx/Tx data rates to target
ret = writeAttributeString(handles.libiio_ctrl_dev, 'in_voltage_sampling_frequency', num2str(handles.rx.Rdata));
if (ret < 0)
    msgbox('Could not write Rx data rate to target!', 'Error', 'error');
    return;
end
ret = writeAttributeString(handles.libiio_ctrl_dev, 'out_voltage_sampling_frequency', num2str(handles.tx.Rdata));
if (ret < 0)
    msgbox('Could not write Tx data rate to target!', 'Error', 'error');
    return;
end

% explicitly write Rx/Tx RF bandwidth to target
ret = writeAttributeString(handles.libiio_ctrl_dev, 'in_voltage_rf_bandwidth', num2str(rx_RFbw_hw));
if (ret < 0)
    msgbox('Could not write Rx RF bandwidth to target!', 'Error', 'error');
    return;
end
ret = writeAttributeString(handles.libiio_ctrl_dev, 'out_voltage_rf_bandwidth', num2str(tx_RFbw_hw));
if (ret < 0)
    msgbox('Could not write Tx RF bandwidth to target!', 'Error', 'error');
    return;
end

% enable both Rx/Tx FIR filters on the target
ret = writeAttributeString(handles.libiio_ctrl_dev, 'in_out_voltage_filter_fir_en', '1');
if (ret < 0)
    msgbox('Could not enable Rx/Tx FIR filters on target!', 'Error', 'error');
    return;
end

msgbox('FIR coefficients successfully written and enabled on target.', 'Filter enabled on target');


function IP_num_Callback(hObject, eventdata, handles)
% hObject    handle to IP_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.connect2target, 'Enable', 'on');
set(handles.connect2target, 'String', 'Connect to Target');


% --- Executes during object creation, after setting all properties.
function IP_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IP_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function FIR_Astop_Callback(hObject, eventdata, handles)
% hObject    handle to FIR_Astop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.active_plot = 0;
plot_buttons_off(handles);
data2gui(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FIR_Astop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FIR_Astop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function [data_clk, bbpll, converter_rate] = get_target_path_rates(libiio, path)
if strcmp(path, 'rx')
    path = 'rx_path_rates';
    sampling_freq = 'in_voltage_sampling_frequency';
    scan_str = 'BBPLL:%d ADC:%d';
else
    path = 'tx_path_rates';
    sampling_freq = 'out_voltage_sampling_frequency';
    scan_str = 'BBPLL:%d DAC:%d';
end

% Read the data clock
[ret, data_clk] = readAttributeDouble(libiio, sampling_freq);
if (ret < 0)
    msgbox('Could not read clocks!', 'Error', 'error');
    return;
end

% Read clocks
[ret, rbuf] = readAttributeString(libiio, path);
if (ret < 0)
    msgbox('Could not read clocks!', 'Error', 'error');
    return;
end

clocks = num2cell(sscanf(rbuf, scan_str));
[bbpll, converter_rate] = clocks{:};

% --- Executes on button press in target_get_clock.
function target_get_clock_Callback(hObject, eventdata, handles)
% hObject    handle to target_get_clock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~ isempty(handles.libiio_ctrl_dev)
    % Read rx/tx clocks
    if (get(handles.filter_type, 'Value') == 1)
        [data_clk, bbpll, converter_rate] = get_target_path_rates(handles.libiio_ctrl_dev, 'rx');
    else
        [data_clk, bbpll, converter_rate] = get_target_path_rates(handles.libiio_ctrl_dev, 'tx');
    end
    
    div = num2str(converter_rate / data_clk);
    decimate = cellstr(get(handles.HB1, 'String'))';
    idx = find(strncmp(decimate, div, length(div)) == 1);
    if (~isempty(idx))
        set(handles.HB1, 'Value', idx(1));
    end
    
    % Set the BPLL div
    opts = get(handles.converter2PLL, 'String');
    for i = 1:length(opts)
        j = char(opts(i));
        j = str2double(j(1:2));
        if j == bbpll / converter_rate
            set(handles.converter2PLL, 'Value', i);
            break;
        end
    end
    
    % Update the data clock
    put_data_clk(handles, data_clk);
    data_clk_Callback(handles.data_clk, eventdata, handles);
end


% --- Executes on button press in FVTool_deeper.
function FVTool_deeper_Callback(hObject, eventdata, handles)
% hObject    handle to FVTool_deeper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel = get_current_rxtx(handles);
converter_rate = sel.Rdata * sel.FIR * sel.HB1 * sel.HB2 * sel.HB3;

N = 500;
Fs = converter_rate; % sampling frequency
F = linspace(0,converter_rate/2,2048);

if (get(handles.filter_type, 'Value') == 1)
    A = sinc(F/Fs).^3;
else
    A = sinc(F/Fs);
end

d = fdesign.arbmag('N,F,A',N,F,A,Fs);
Hcon = design(d,'SystemObject',false);

apass = str2double(get(handles.Apass, 'String'));
astop = str2double(get(handles.Astop, 'String'));

str = sprintf('%s Filter\nFpass = %g MHz; Fstop = %g MHz\nApass = %g dB; Astop = %g dB', sel.RxTx, sel.Fpass/1e6, sel.Fstop/1e6, apass, astop);

if get(handles.use_FPGAfilter, 'Value')== 1
    hfvt1 = fvtool(Hcon,handles.analogfilter,handles.Hmiddle,handles.grpdelaycal,handles.plutofilter,...
        'FrequencyRange','Specify freq. vector', ...
        'FrequencyVector',linspace(0,converter_rate/2,2048),'Fs',...
        converter_rate, ...
        'ShowReference','off','Color','White');
    set(hfvt1, 'Color', [1 1 1]);
    set(hfvt1.CurrentAxes, 'YLim', [-100 20]);
    legend(hfvt1, 'Converter','Analog','Analog + Half Band','Analog + HB + FIR','Analog + HB + FIR + Pluto');
else
    hfvt1 = fvtool(Hcon,handles.analogfilter,handles.Hmiddle,handles.grpdelaycal,...
        'FrequencyRange','Specify freq. vector', ...
        'FrequencyVector',linspace(0,converter_rate/2,2048),'Fs',...
        converter_rate, ...
        'ShowReference','off','Color','White');
    set(hfvt1, 'Color', [1 1 1]);
    set(hfvt1.CurrentAxes, 'YLim', [-100 20]);
    legend(hfvt1, 'Converter','Analog','Analog + Half Band','Analog + HB + FIR');
    
end
text(1, 10,...
    str,...
    'BackgroundColor','white',...
    'EdgeColor','red');

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over FIR_Astop.
function FIR_Astop_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to FIR_Astop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function zoom_axis(ax)
A = get(ax, 'XTick') / 1e6;
arrayfun(@num2str, A, 'unif', 0);
set(ax, 'XTickLabel', A);

function plot_buttons_off(handles)
set(handles.FVTool_deeper, 'Visible', 'off');
set(handles.FVTool_datarate, 'Visible', 'off');
set(handles.save2workspace, 'Visible', 'on');
set(handles.save2workspace, 'Enable', 'off')
set(handles.save2coefficients, 'Visible', 'on');
set(handles.save2coefficients, 'Enable', 'off');
set(handles.save2noos, 'Visible', 'on');
set(handles.save2noos, 'Enable', 'off');
set(handles.save2target, 'Visible', 'on');
set(handles.save2HDL, 'Visible', 'off');

set(handles.results_Apass, 'Visible', 'off');
set(handles.results_Astop, 'Visible', 'off');
set(handles.results_maxinput, 'Visible', 'off');
set(handles.results_maxinputdB, 'Visible', 'off');
set(handles.results_taps, 'Visible', 'off');
set(handles.results_group_delay, 'Visible', 'off');

function freq_off(handles)
set(handles.Fpass, 'Enable', 'off');
set(handles.Fstop, 'Enable', 'off');
set(handles.Fcenter, 'Enable', 'off');
% set(handles.Fcutoff, 'Enable', 'off');
% set(handles.RFbw, 'Enable', 'off');


function freq_on(handles)
set(handles.Fpass, 'Enable', 'on');
set(handles.Fstop, 'Enable', 'on');
set(handles.Fcenter, 'Enable', 'on');
% set(handles.Fcutoff, 'Enable', 'on');
% set(handles.RFbw, 'Enable', 'on');


function create_filter(hObject, handles)
handles = guidata(hObject);
v = version('-release');
v = str2double(v(1:4));
if (v < 2014)
    choice = questdlg('Sorry. The AD9361/AD9364 Filter Design Wizard requires at least the R2014 version of MATLAB. You do not seem to have it installed.', ...
        'Error Message', ...
        'More Information','OK','OK');
    switch choice
        case 'More Information'
            web('http://www.mathworks.com/products/matlab/');
        case 'OK'
    end
    return
end

error_msg = '';
if ~ license('test','signal_blocks') || ~ license('checkout','signal_blocks')
    error_msg = 'You are missing the license for it.';
elseif isempty(ver('dsp'))
    error_msg = 'You do not seem to have it installed.';
end

if error_msg
    choice = questdlg(['Sorry. The AD9361/AD9364 Filter Design Wizard requires the DSP System Toolbox. You do not seem to have it installed. ' error_msg], ...
        'Error Message', ...
        'More Information','OK','OK');
    switch choice
        case 'More Information'
            web('http://www.mathworks.com/products/dsp-system/');
        case 'OK'
    end
    return
end

set(handles.design_filter, 'Enable', 'off');

sel = get_current_rxtx(handles);
converter_rate = sel.Rdata * sel.FIR * sel.HB1 * sel.HB2 * sel.HB3;
PLL_rate = str2double(get(handles.Pll_rate, 'String')) * 1e6;
FIR_rate = sel.Rdata * sel.FIR;
HB1_rate = FIR_rate * sel.HB1;
HB2_rate = HB1_rate * sel.HB2;
HB3_rate = HB2_rate * sel.HB3;

% determine the RF bandwidth from the current caldiv
RFbw = get_rfbw(handles, sel.caldiv);

% filter design input structure
filter_input = sel;
filter_input.converter_rate = converter_rate;
filter_input.PLL_rate = PLL_rate;
filter_input.dBstop_FIR = sel.FIRdBmin;
filter_input.wnom = value2Hz(handles, handles.freq_units, str2double(get(handles.Fcutoff, 'String')));
filter_input.int_FIR = get(handles.Use_FIR, 'Value');
filter_input.FPGAfilter = get(handles.use_FPGAfilter, 'Value');
filter_input.RFbw = RFbw;

plot_buttons_off(handles);

% make sure things are sane before drawing
if sel.Fpass >= sel.Fstop || sel.Fpass <= 0 || sel.Fstop <= 0
    display_default_image(hObject);
    plot_buttons_off(handles);
    handles.active_plot = 0;
    set(handles.design_filter, 'Enable', 'on');
    guidata(hObject, handles);
    return;
end

oldpointer = get(gcf, 'pointer');
set(gcf,'Pointer','watch');
drawnow;

if (filter_input.phEQ == 0)
    filter_input.phEQ = minimize_group_delay(handles, filter_input);
end

filter_result = internal_design_filter(filter_input);

handles.filter = filter_result.filter;
handles.analogfilter = filter_result.Hanalog;
handles.grpdelayvar = filter_result.grpdelayvar;
handles.grpdelaycal = clone(filter_result.filter);
handles.plotpluto = clone(filter_result.filter);

if get(handles.filter_type, 'Value') == 1  % Rx
    addStage(handles.grpdelaycal,filter_result.Hd2,1);
    addStage(handles.grpdelaycal,filter_result.Hd1,2);
    handles.plutofilter = clone(handles.grpdelaycal);
    ast = 80;
    tw = 0.08;
    f = fdesign.decimator(8, 'Nyquist', 8, 'TW,Ast', tw, ast);
    hf = design(f,'SystemObject',true);
    hf.FullPrecisionOverride = false;
    hf.OutputDataType='Custom';
    hf.CustomOutputDataType=numerictype([],16,15);
    hf.CoefficientsDataType='Custom';
    hf.CustomCoefficientsDataType=numerictype([],16,15);
    hf.ProductDataType='Custom';
    hf.CustomProductDataType=numerictype([],16,15);
    hf.AccumulatorDataType = 'Custom';
    hf.CustomAccumulatorDataType=numerictype([],16,15);
    
    addStage(handles.plutofilter, hf);
    addStage(handles.plotpluto, hf);
    handles.rfirtaps = int16(filter_result.firtaps);
    
    % values used for saving to a filter file or pushing to the target directly
    handles.rx = filter_result;
    
else  % Tx
    addStage(handles.grpdelaycal,filter_result.Hd1);
    addStage(handles.grpdelaycal,filter_result.Hd2);
    handles.plutofilter = clone(handles.grpdelaycal);
    ast = 80;
    tw = 0.08;
    f = fdesign.interpolator(8, 'Nyquist', 8, 'TW,Ast', tw, ast);
    hf = design(f,'kaiserwin','SystemObject',true);
    hf.Numerator = hf.Numerator./8;
    hf.FullPrecisionOverride = false;
    hf.OutputDataType='Custom';
    hf.CustomOutputDataType=numerictype([],16,15);
    hf.CoefficientsDataType='Custom';
    hf.CustomCoefficientsDataType=numerictype([],16,15);
    hf.ProductDataType='Custom';
    hf.CustomProductDataType=numerictype([],16,15);
    hf.AccumulatorDataType = 'Custom';
    hf.CustomAccumulatorDataType=numerictype([],16,15);
    
    addStage(handles.plutofilter, hf, 1);
    addStage(handles.plotpluto, hf, 1);
    handles.tfirtaps = int16(filter_result.firtaps);
    
    % values used for saving to a filter file or pushing to the target directly
    handles.tx = filter_result;
end

set(gcf, 'Pointer', oldpointer);

if get(handles.phase_eq, 'Value')
    set(handles.target_delay, 'String', num2str(filter_input.phEQ, 8));
end

handles.nfirtaps = filter_result.nfirtaps;
handles.Hmiddle = filter_result.Hmiddle;

set(handles.FVTool_deeper, 'Visible', 'on');
set(handles.FVTool_datarate, 'Visible', 'on');

set(handles.save2workspace, 'Enable', 'on')
if isfield(handles, 'rfirtaps') && isfield(handles, 'tfirtaps')
    set(handles.save2coefficients, 'Enable', 'on');
    set(handles.save2noos, 'Enable', 'on');
    if ~ isempty(handles.libiio_ctrl_dev)
        set(handles.save2target, 'Enable', 'on');
    end
end

units = cellstr(get(handles.Freq_units, 'String'));
units = char(units(get(handles.Freq_units, 'Value')));
if get(handles.use_FPGAfilter, 'Value')== 0
    set(handles.FVTool_datarate, 'String', sprintf('FVTool to %g %s', str2double(get(handles.data_clk, 'String'))/2, units));
else
    set(handles.FVTool_datarate, 'String', sprintf('FVTool to %g %s', str2double(get(handles.data_clk, 'String'))*8/2, units));
end

if ~ get(handles.Use_FIR, 'Value')
    set(handles.save2HDL, 'Visible', 'on');
end

set(handles.results_Apass, 'Visible', 'on');
set(handles.results_Astop, 'Visible', 'on');
set(handles.results_maxinput, 'Visible', 'on');
set(handles.results_maxinputdB, 'Visible', 'on');
set(handles.results_taps, 'Visible', 'on');
set(handles.results_group_delay, 'Visible', 'on');

set(handles.results_taps, 'String', [num2str(handles.nfirtaps) ' ']);
set(handles.RFbw, 'String', num2str(Hz2value(handles, handles.freq_units, RFbw)));

G = 8192;
axes(handles.magnitude_plot);
cla(handles.magnitude_plot);

if get(handles.filter_type, 'Value') == 1
    channel = 'Rx';
else
    channel = 'Tx';
end

if get(handles.use_FPGAfilter, 'Value')== 1
    handles.active_plot = plot(handles.magnitude_plot, linspace(0,sel.Rdata/2,G),mag2db(...
        abs(analogresp(channel,linspace(0,sel.Rdata/2,G),converter_rate,filter_result.b1,filter_result.a1,filter_result.b2,filter_result.a2).*freqz(...
        handles.plotpluto,linspace(0,sel.Rdata/2,G),converter_rate))));
    xlim([0 sel.Rdata/8]);
else
    handles.active_plot = plot(handles.magnitude_plot, linspace(0,sel.Rdata/2,G),mag2db(...
        abs(analogresp(channel,linspace(0,sel.Rdata/2,G),converter_rate,filter_result.b1,filter_result.a1,filter_result.b2,filter_result.a2).*freqz(...
        filter_result.filter,linspace(0,sel.Rdata/2,G),converter_rate))));
    xlim([0 sel.Rdata/2]);
end

ylim([-100 10]);
zoom_axis(gca);
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');

% plot the mask that we are interested in
if get(handles.use_FPGAfilter, 'Value')== 1
    sel.Fpass = sel.Rdata/2*0.085;
    sel.Fstop = sel.Rdata/2*0.165;
end
line([sel.Fpass sel.Fpass], [-(sel.Apass/2) -100], 'Color', 'Red');
line([0 sel.Fpass], [-(sel.Apass/2) -(sel.Apass/2)], 'Color', 'Red');
line([0 sel.Fstop], [sel.Apass/2 sel.Apass/2], 'Color', 'Red');
line([sel.Fstop sel.Fstop], [sel.Apass/2 -sel.Astop], 'Color', 'Red');
line([sel.Fstop sel.Rdata], [-sel.Astop -sel.Astop], 'Color', 'Red');

% add the quantitative values about actual passband, stopband, and group delay
if get(handles.use_FPGAfilter, 'Value')== 1
    if get(handles.filter_type, 'Value') == 1
        rg_stop = abs(freqz(handles.plotpluto,linspace(sel.Fstop,sel.Rdata/8,G),converter_rate).*analogresp('Rx',linspace(sel.Fstop,sel.Rdata/8,G),converter_rate,filter_result.b1,filter_result.a1,filter_result.b2,filter_result.a2));
    else
        rg_stop = abs(freqz(handles.plotpluto,linspace(sel.Fstop,sel.Rdata/8,G),converter_rate).*analogresp('Tx',linspace(sel.Fstop,sel.Rdata/8,G),converter_rate,filter_result.b1,filter_result.a1,filter_result.b2,filter_result.a2));
    end
    Astop_actual = -mag2db(max(rg_stop));
    set(handles.results_Astop, 'String', [num2str(Astop_actual) ' dB ']);
else
    set(handles.results_Astop, 'String', [num2str(filter_result.Astop_actual) ' dB ']);
end
set(handles.results_Apass, 'String', [num2str(filter_result.Apass_actual) ' dB ']);
set(handles.results_maxinput, 'String', [num2str(1.0 / max(abs(freqz(filter_result.Hmd)))) ' FS ']);
set(handles.results_maxinputdB, 'String', [num2str(20 * log10(1.0 / max(abs(freqz(filter_result.Hmd))))) ' dB ']);
set(handles.results_group_delay, 'String', [num2str(filter_result.grpdelayvar * 1e9, 3) ' ns ']);

if ~isempty(ver('fixedpoint')) && license('test','fixed_point_toolbox') && license('checkout','fixed_point_toolbox')
    set(handles.results_fixed, 'String', 'Fixed Point');
else
    set(handles.results_fixed, 'String', 'Floating point approx');
end

% check if the filter is compatible with available clock rates
clk_divs = [1 2 3 4 6 8 12 16 24 32 48];
clk_rates = converter_rate ./ clk_divs;
if ~any(4 * sel.Rdata == clk_rates)
    warndlg('This filter is not compatible with 2Rx2Tx mode.', 'Filter compatibility');
end
if ~any(2 * sel.Rdata == clk_rates)
    warndlg('This filter is not compatible with 1Rx1Tx mode.', 'Filter compatibility');
end

set(handles.design_filter, 'Visible', 'on');
guidata(hObject, handles);

function phEQ = minimize_group_delay(handles, filter_input)
set(handles.results_group_delay, 'Visible', 'on');
filter_input.phEQ = 0;

set(handles.target_delay, 'String', num2str(filter_input.phEQ, 4));
drawnow;

filter_result = internal_design_filter(filter_input);
nom = filter_result.delay * 1e9;

results = [nom (filter_result.grpdelayvar * 1e9);];
nom = round(nom);

box on;
axis on;
axes(handles.magnitude_plot);
cla(handles.magnitude_plot);
plot(results(:,1), results(:,2) ,'r.');
xlabel('Group Delay target (ns)');
ylabel('Group Delay variance (ns)');
drawnow;

span = 20;
initial_step = 1; % 1/(2^n)

for j = initial_step:10
    i = -(span * initial_step);
    i_start = i;
    i_end = (span * initial_step);
    i_min = (nom + i/(2^j));
    i_max = nom + (i_end/(2^j));
    while (i < i_end)
        filter_input.phEQ = nom + i/(2^j);
        i = i + 1;
        idx = find(abs(results( :, 1 ) - filter_input.phEQ) < 1/(2^(1+j)));
        if isempty(idx)
            set(handles.target_delay, 'String', num2str(filter_input.phEQ, 8));
            drawnow;
            
            filter_result = internal_design_filter(filter_input);
            results = [results ; (filter_result.delay * 1e9) (filter_result.grpdelayvar * 1e9);];
            results = sortrows(results, 1);
            
            [minval, minidx] = min(results(:,2));
            str = sprintf('%1.2f@%3.2f', minval, results(minidx));
            set(handles.results_group_delay, 'String', str);
            if (results(minidx) + span/(2^j) > nom + i_end/(2^j))
                i_end = ceil(abs(results(minidx) + span/(2^j) - nom) * (2^j));
                i_max = max(i_max, nom + (i_end/(2^j)));
            end
            
            plot(results(:,1), results(:,2) ,'r.');
            xlim([i_min i_max]);
            xlabel('Group Delay target (ns)');
            ylabel('Group Delay variance (ns)');
            drawnow;
        end
        if i == i_end
            [minval, minidx] = min(results(:,2));
            if minidx == 1
                i = i_start - (span * initial_step);
                i_end = i_start;
                i_start = i;
                i_min = (nom + i/(2^j));
            end
        end
    end
    
    [minval, minidx] = min(results(:,2));
    nom = round(results(minidx) * (2^j)) / (2^j);
    
    if (j > initial_step )
        sub = results(minidx - (span * initial_step):minidx + i_end,:);
        [sub_minval, sub_minidx] = min(sub(:,2));
        [sub_maxval, sub_maxidx] = max(sub(:,2));
        x = std(sub(find(sub(:,2) <= (sub_maxval + sub_minval)/2),2));
        y = std(sub(find(sub(:,2) >= (sub_maxval + sub_minval)/2),2));
        if (x < 0.05) && (y < 0.05)
            break;
        end
    end
end

% Best result found
phEQ = results(minidx);


function load_settings(hObject, handles)

filename = 'ad9361_settings.mat';

if ~ exist(filename, 'file')
    errordlg('I can not find the required files, must be some sort of installation error', ...
        'File Error');
    set(handles.Saved_Filters, 'Visible', 'off');
    return;
end

options = load(filename);
Tx = fieldnames(options.ad9361_settings.tx);
Rx = fieldnames(options.ad9361_settings.rx);

Tx_numRows = size(Tx, 1);
Rx_numRows = size(Rx, 1);

needle = cellstr(get(handles.Saved_Filters, 'String'));
needle = char(needle(get(handles.Saved_Filters, 'Value')));
needle = strsplit(needle, '(');
needle = deblank(needle{1});

for matchRx = 1:Rx_numRows;
    match = Rx{matchRx, end};
    if find(strncmp(needle, match, length(needle)))
        break
    end
end

for matchTx = 1:Tx_numRows;
    match = Tx{matchTx, end};
    if find(strncmp(needle, match, length(needle)))
        break
    end
end

% remove generated filter taps if they exist
if isfield(handles, 'rfirtaps')
    handles = rmfield(handles, 'rfirtaps');
end
if isfield(handles, 'tfirtaps')
    handles = rmfield(handles, 'tfirtaps');
end

handles.rx = cook_input(getfield(options.ad9361_settings.rx, Rx{matchRx}));
handles.tx = cook_input(getfield(options.ad9361_settings.tx, Tx{matchTx}));

set(handles.store_filter, 'Visible', 'off');
guidata(hObject, handles);
data2gui(hObject, handles);

function data2gui(hObject, handles)

OK = 1;

rates_uipanel = findall(gcf, 'type', 'uipanel', 'Tag', 'rates_uipanel');
if get(handles.filter_type, 'Value') == 1
    % Receive
    sel = handles.rx;
    max_HB = handles.bounds.MAX_RX;
    set(rates_uipanel, 'Tag', 'rates_uipanel', 'Title', 'AD936x Decimation Rates');
else
    % Transmit
    sel = handles.tx;
    max_HB = handles.bounds.MAX_TX;
    set(rates_uipanel, 'Tag', 'rates_uipanel', 'Title', 'AD936x Interpolation Rates');
end

% set things from the file.
if sel.phEQ == -1
    set(handles.phase_eq, 'Value', 0);
else
    show_advanced(handles);
    set(handles.phase_eq, 'Value', 1);
    set(handles.target_delay, 'Value', sel.phEQ);
end

if sel.caldiv && sel.caldiv ~= default_caldiv(handles)
    show_advanced(handles);
    set_fcutoff(handles, sel.caldiv);
end

if get(handles.use_FPGAfilter, 'Value')==1
    sel.Rdata = sel.Rdata/8;
    FPGA_rate = sel.Rdata*8;
else
    FPGA_rate = sel.Rdata;
end
FIR_rate = FPGA_rate * sel.FIR;
HB1_rate = FIR_rate * sel.HB1;
HB2_rate = HB1_rate * sel.HB2;
HB3_rate = HB2_rate * sel.HB3;

set(handles.FPGA_rate, 'String', num2str(Hz2value(handles, handles.freq_units, FPGA_rate)));
set(handles.Fpass, 'String', num2str(Hz2value(handles, handles.freq_units, sel.Fpass)));
set(handles.Fstop, 'String', num2str(Hz2value(handles, handles.freq_units, sel.Fstop)));
set(handles.RFbw, 'String', num2str(Hz2value(handles, handles.freq_units, get_rfbw(handles, sel.caldiv))));

set(handles.Fcenter, 'String', num2str(Hz2value(handles, handles.freq_units, sel.Fcenter)));

set(handles.Apass, 'String', num2str(sel.Apass));
set(handles.Astop, 'String', num2str(sel.Astop));

set(handles.data_clk, 'String', num2str(Hz2value(handles, handles.freq_units, sel.Rdata)));

if handles.rx.Fstop <= handles.rx.Rdata / 2
    set(handles.Fstop, 'ForegroundColor', [0 0 0]);
else
    set(handles.Fstop, 'ForegroundColor', [1 0 0]);
    if OK
        warn = 'Data rate must be at least 2 times Fstop';
    end
    OK = 0;
end

if handles.rx.Rdata == handles.tx.Rdata
    set(handles.data_clk, 'ForegroundColor', [0 0 0]);
else
    set(handles.data_clk, 'ForegroundColor', [1 0 0]);
    if OK
        warn = 'Rx and Tx data rates need to be the same';
    end
    OK = 0;
end

opts = get(handles.FIR, 'String');
for i = 1:length(opts)
    j = char(opts(i));
    j = str2double(j(1:2));
    if j == sel.FIR
        set(handles.FIR, 'Value', i);
        break;
    end
end
set(handles.FIR_rate, 'String', num2str(Hz2value(handles, handles.freq_units, FIR_rate)));
if FIR_rate > handles.bounds.MAX_FIR
    set(handles.FIR_rate, 'ForegroundColor', [1 0 0]);
    if OK
        warn = 'FIR rate too high';
    end
    OK = 0;
else
    set(handles.FIR_rate, 'ForegroundColor', [0 0 0]);
end

opts = get(handles.HB1, 'String');
for i = 1:length(opts)
    j = char(opts(i));
    j = str2double(j(1:2));
    if j == sel.HB1
        set(handles.HB1, 'Value', i);
        break;
    end
end

set(handles.HB1_rate, 'String', num2str(Hz2value(handles, handles.freq_units, HB1_rate)));
if HB1_rate > max_HB.HB1
    set(handles.HB1_rate, 'ForegroundColor', [1 0 0]);
    if OK
        warn = 'HB1 rate too high';
    end
    OK = 0;
else
    set(handles.HB1_rate, 'ForegroundColor', [0 0 0]);
end

opts = get(handles.HB2, 'String');
for i = 1:length(opts)
    j = char(opts(i));
    j = str2double(j(1:2));
    if j == sel.HB2
        set(handles.HB2, 'Value', i);
        break;
    end
end
set(handles.HB2_rate, 'String', num2str(Hz2value(handles, handles.freq_units, HB2_rate)));
if HB2_rate > max_HB.HB2
    set(handles.HB2_rate, 'ForegroundColor', [1 0 0]);
    if OK
        warn = 'HB2 rate too high';
    end
    OK = 0;
else
    set(handles.HB2_rate, 'ForegroundColor', [0 0 0]);
end

opts = get(handles.HB3, 'String');
for i = 1:length(opts)
    j = char(opts(i));
    j = str2double(j(1:2));
    if j == sel.HB3
        set(handles.HB3, 'Value', i);
        break;
    end
end
set(handles.HB3_rate, 'String', num2str(Hz2value(handles, handles.freq_units, HB3_rate)));
if HB3_rate > max_HB.HB3
    set(handles.HB3_rate, 'ForegroundColor', [1 0 0]);
    if OK
        warn = 'HB3 rate too high';
    end
    OK = 0;
else
    set(handles.HB3_rate, 'ForegroundColor', [0 0 0]);
end

HBs = sel.HB1 * sel.HB2 * sel.HB3;
set(handles.HBs_rate, 'String', get(handles.HB3_rate, 'String'));
opts = get(handles.HBs, 'String');
for i = 1:length(opts)
    j = char(opts(i));
    j = str2double(j(1:2));
    if j == HBs
        set(handles.HBs, 'Value', i);
        break;
    end
end

% show correct PLL div option to match settings
opts = get(handles.converter2PLL, 'String');
for i = 1:length(opts)
    j = char(opts(i));
    j = str2double(j(1:2));
    if j == sel.PLL_mult
        set(handles.converter2PLL, 'Value', i);
        break;
    end
end

% PLL Settings
set(handles.DAC_by2, 'Value', handles.tx.DAC_div);
ADC_clk = handles.rx.Rdata * handles.rx.FIR * handles.rx.HB1 * handles.rx.HB2 * handles.rx.HB3;
set(handles.ADC_clk, 'String', num2str(ADC_clk / 1e6));
DAC_clk = handles.tx.Rdata * handles.tx.FIR * handles.tx.HB1 * handles.tx.HB2 * handles.tx.HB3;
set(handles.DAC_clk, 'String', num2str(DAC_clk / 1e6));

ADC_clk = round(ADC_clk);
DAC_clk = round(DAC_clk);

% Make sure ADC clocks are within bounds
if (ADC_clk <= handles.bounds.MAX_ADC_CLK) && (ADC_clk >= handles.bounds.MIN_ADC_CLK)
    if OK
        set(handles.ADC_clk, 'ForegroundColor', [0 0 0]);
    end
else
    set(handles.ADC_clk, 'ForegroundColor', [1 0 0]);
    if OK
        if (ADC_clk > handles.bounds.MAX_ADC_CLK)
            max_adc = num2str(Hz2value(handles, 3, handles.bounds.MAX_ADC_CLK));
            warn = sprintf('ADC Clock above maximum (%s %s)', max_adc, 'MHz');
        else
            min_adc = num2str(Hz2value(handles, 3, handles.bounds.MIN_ADC_CLK));
            warn = sprintf('ADC Clock below minimum (%s %s)', min_adc, 'MHz');
        end
        OK = 0;
    end
end

% Make sure DAC clocks are within bounds
if (DAC_clk <= handles.bounds.MAX_DAC_CLK) && (DAC_clk >= handles.bounds.MIN_DAC_CLK)
    if OK
        set(handles.DAC_clk, 'ForegroundColor', [0 0 0]);
    end
else
    set(handles.DAC_clk, 'ForegroundColor', [1 0 0]);
    if OK
        if (DAC_clk > handles.bounds.MAX_DAC_CLK)
            max_dac = num2str(Hz2value(handles, 3, handles.bounds.MAX_DAC_CLK));
            warn = sprintf('DAC Clock above maximum (%s %s)', max_dac, 'MHz');
        else
            min_dac = num2str(Hz2value(handles, 3, handles.bounds.MIN_DAC_CLK));
            warn = sprintf('DAC Clock below minimum (%s %s)', min_dac, 'MHz');
        end
        OK = 0;
    end
end

% Make sure Rx and Tx PLL rates are equal
if (handles.rx.FIR * handles.rx.HB1 * handles.rx.HB2 * handles.rx.HB3) == ...
        (handles.tx.FIR * handles.tx.HB1 * handles.tx.HB2 * handles.tx.HB3 * handles.tx.DAC_div)
    if OK
        set(handles.ADC_clk, 'ForegroundColor', [0 0 0]);
        set(handles.DAC_clk, 'ForegroundColor', [0 0 0]);
    end
else
    set(handles.ADC_clk, 'ForegroundColor', [1 0 0]);
    set(handles.DAC_clk, 'ForegroundColor', [1 0 0]);
    set(handles.Pll_rate, 'String', '?');
    if OK
        warn = '(DAC * multipler) and ADC rates do not match';
    end
    OK = 0;
end

% Make sure the PLL rate is within the allowed bounds
pll = get_pll_rate(handles);
set(handles.Pll_rate, 'String', num2str(pll / 1e6));

if (pll <= handles.bounds.MAX_BBPLL_FREQ) && (pll >= handles.bounds.MIN_BBPLL_FREQ)
    set(handles.Pll_rate, 'ForegroundColor', [0 0 0]);
else
    set(handles.Pll_rate, 'ForegroundColor', [1 0 0]);
    if OK
        if (pll > handles.bounds.MAX_BBPLL_FREQ)
            max_bbpll = num2str(Hz2value(handles, 3, handles.bounds.MAX_BBPLL_FREQ));
            warn = sprintf('PLL rate above maximum (%s %s)', max_bbpll, 'MHz');
        else
            min_bbpll = num2str(Hz2value(handles, 3, handles.bounds.MIN_BBPLL_FREQ));
            warn = sprintf('PLL rate below minimum (%s %s)', min_bbpll, 'MHz');
        end
    end
    OK = 0;
end

% don't have data - so don't display the FVTool button
set(handles.FVTool_deeper, 'Visible', 'off');
set(handles.FVTool_datarate, 'Visible', 'off');

set(handles.save2HDL, 'Visible', 'off');

set(handles.save2workspace, 'Enable', 'off');
set(handles.save2coefficients, 'Enable', 'off');
set(handles.save2noos, 'Enable', 'off');
set(handles.save2target, 'Enable', 'off');

% disable connect to target button if libiio bindings are missing
[pathstr, ~, ~] = fileparts(mfilename('fullpath'));
if ~exist(fullfile(pathstr, 'libiio', 'libiio_if.m'), 'file')
    set(handles.connect2target, 'Enable', 'off');
end

%set(handles.target_get_clock, 'Visible', 'off');

if OK
    set(handles.design_filter, 'Enable', 'on');
    set(handles.design_warning, 'String', '');
    set(handles.design_warning, 'Visible', 'off');
else
    set(handles.design_filter, 'Enable', 'off');
    set(handles.design_warning, 'String', warn);
    set(handles.design_warning, 'Visible', 'on');
end

guidata(hObject, handles);

function reset_input(hObject, handles)
handles.active_plot = 0;
display_default_image(hObject);
handles = guidata(hObject);
filename = 'ad9361_settings.mat';

if ~ exist(filename, 'file')
    errordlg('I can not find the required files, must be some sort of installation error', ...
        'File Error');
    set(handles.Saved_Filters, 'Visible', 'off');
    return;
end

options = load(filename);
Tx = fieldnames(options.ad9361_settings.tx);
Rx = fieldnames(options.ad9361_settings.rx);

Tx_numRows = size(Tx, 1);
Rx_numRows = size(Rx, 1);
choices = cell(1, Tx_numRows);

for rowTx = 1:Tx_numRows;
    needle = Tx{rowTx, end};
    if find(strncmp(needle, Rx, length(needle)))
        choices{rowTx} = strcat(needle, ' (Rx & Tx)');
    else
        choices{rowTx} = strcat(needle, ' (Tx only)');
    end
end
for rowRx = 1:Rx_numRows;
    needle = Rx{rowRx, end};
    if ~ find(strncmp(needle, choices, length(needle)))
        rowTx =+ 1;
        choices{rowTx} = strcat(needle, ' (Rx only)');
    end
end

choices = sort_nat(choices);
set(handles.Saved_Filters, 'String', choices);
set(handles.Saved_Filters, 'Value', 1);

guidata(hObject, handles);

%                               Rx           Tx
% FIR = [1 2 4];                0.23|61.44        61.44
% HB1 = [1 2];                  0.93|122.88       122.88
% HB2 = [1 2];                  1.86|245.76       160
% HB3 = [1 2 3];                3.72|320          320
% converter                    11.2 |640          5.58|320
% DAC_by2 [1 2]                                   11.2|640
% PLL = [1 2 4 8 16 32 64];     715|1430          715|1430


% --------------------------------------------------------------------
function new_tooltip_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to new_tooltip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

reset_input(hObject, handles);
load_settings(hObject, handles);
handles = guidata(hObject);
handles.freq_units = get(handles.Freq_units, 'Value');
handles.active_plot = 0;

guidata(hObject, handles);


function display_default_image(hObject)
handles = guidata(gcf);
set(handles.FVTool_deeper, 'Visible', 'off');
axes(handles.magnitude_plot);
handles.active_plot = 0;

cla(handles.magnitude_plot);
set(handles.magnitude_plot, 'Tag', 'magnitude_plot');

max_y = 20;
max_x = 200;
ripple = 10;
Fpass = 90;
Fcenter = 100;
Fstop = 110;
label_colour = 'Red';

box on;
axis on;

set(gca,'XTickLabel',{});
xlabel('');
xlim([0 max_x]);
text(-10, 0, '0dB');

ylabel('Mag (dB)');
set(gca,'YTickLabel',{});
ylim([-100 max_y]);

switch get(get(handles.Response_Type, 'SelectedObject'), 'String')
    case 'Lowpass'
        % Low part of the low pass
        line([0 Fpass], [-ripple -ripple], 'Color', 'Black');
        line([0 Fpass], [ripple ripple], 'Color', 'Black');
        line([Fpass Fpass], [max_y ripple], 'Color', 'Black');
        line([Fpass Fpass], [-ripple -100], 'Color', 'Black');
        line([0 Fstop], [0 0], 'Color', label_colour, 'LineStyle', ':');
        line([Fpass Fstop+30], [-ripple -ripple], 'Color', label_colour, 'LineStyle', ':');
        line([Fpass Fstop+30], [ripple ripple], 'Color', label_colour, 'LineStyle', ':');
        
        handles.arrows{1} = annotation('arrow', 'Y',[max_y ripple], 'X',[130 130]);
        set(handles.arrows{1}, 'Color', label_colour);
        handles.arrows{2} = annotation('arrow', 'Y',[-max_y -ripple], 'X',[130 130]);
        set(handles.arrows{2}, 'Color', label_colour);
        text(Fstop + 12, 0, 'A_{pass}', 'BackgroundColor','white', 'EdgeColor','white');
        
        % Stop band
        line([Fstop max_x-10], [-80 -80], 'Color', 'Black');
        line([max_x-10 max_x-10], [-80+ripple -80], 'Color', 'Black');
        line([Fstop Fstop], [-80+ripple -80], 'Color', 'Black');
        line([Fstop Fstop], [-80+ripple -80], 'Color', 'Black');
        line([Fstop Fstop], [-80 -100], 'Color', label_colour, 'LineStyle', ':');
        line([max_x-10 max_x-10], [-80 -100], 'Color', label_colour, 'LineStyle', ':');
        
        line([150 170], [0 0], 'Color', label_colour, 'LineStyle', ':');
        text(0, -108, '0');
        text(Fpass - 5, -108, 'F_{pass}');
        text(Fstop - 5, -108, 'F_{stop}');
        text(max_x - 15, -108, 'Fs_{/2}');
        
        % A(stop) label and arrows
        hTest = text(150, -40, 'A_{stop}');
        textExt = get(hTest,'Extent');
        w = textExt(1) + textExt(3)/2;
        handles.arrows{3} = annotation('arrow', 'Y',[-35 0], 'X',[w w]);
        set(handles.arrows{3}, 'Color', label_colour);
        handles.arrows{4} = annotation('arrow', 'Y',[-45 -80], 'X',[w w]);
        set(handles.arrows{4}, 'Color', label_colour);
        
        % reparent arrows within the filter plot so they're displayed properly
        plot = findall(gcf, 'type', 'axes', 'Tag', 'magnitude_plot');
        for i = 1:4
            set(handles.arrows{i}, 'Parent', plot);
        end
        
    case 'Root Raised Cosine'
        % Pass band
        line([0 Fpass], [-ripple -ripple], 'Color', 'Black');
        line([0 Fpass], [ripple ripple], 'Color', 'Black');
        line([Fpass Fpass], [max_y ripple], 'Color', 'Black');
        line([Fpass Fpass], [-ripple -100], 'Color', 'Black');
        
        % Stop band
        line([Fstop max_x-10], [-80 -80], 'Color', 'Black');
        line([max_x-10 max_x-10], [-80+ripple -80], 'Color', 'Black');
        line([Fstop Fstop], [-80+ripple -80], 'Color', 'Black');
        line([Fstop Fstop], [-80+ripple -80], 'Color', 'Black');
        line([Fstop Fstop], [-80 -100], 'Color', label_colour, 'LineStyle', ':');
        line([max_x-10 max_x-10], [-80 -100], 'Color', label_colour, 'LineStyle', ':');
        
        text(0, -108, '0');
        text(Fpass - 5, -108, 'F_{pass}');
        text(Fstop - 5, -108, 'F_{stop}');
        text(max_x - 15, -108, 'Fs_{/2}');
        
    case 'Bandpass'
        Fpass = 20;
        Fcenter = 80;
        Fstop = 30;
        % pass part of the bandpass
        line([Fcenter-Fpass Fcenter+Fpass], [-ripple -ripple], 'Color', 'Black');
        line([Fcenter-Fstop Fcenter+Fstop], [ripple ripple], 'Color', 'Black');
        line([Fcenter-Fstop Fcenter-Fstop], [max_y ripple], 'Color', 'Black');
        line([Fcenter+Fstop Fcenter+Fstop], [max_y ripple], 'Color', 'Black');
        line([Fcenter-Fpass Fcenter-Fpass], [-ripple -80], 'Color', 'Black');
        line([Fcenter+Fpass Fcenter+Fpass], [-ripple -80], 'Color', 'Black');
        
        line([0 Fcenter+Fstop], [0 0], 'Color', label_colour, 'LineStyle', ':');
        line([Fcenter+Fstop Fcenter+Fstop+30], [-ripple -ripple], 'Color', label_colour, 'LineStyle', ':');
        line([Fcenter+Fstop Fcenter+Fstop+30], [ripple ripple], 'Color', label_colour, 'LineStyle', ':');
        
        [x1, y1] = xy2norm(130, ripple, handles);
        [x2, y2] = xy2norm(130, max_y, handles);
        handles.arrows{1} = annotation('arrow', 'Y',[y2 y1], 'X',[x1 x2]);
        set(handles.arrows{1}, 'Color', label_colour);
        [x1, y1] = xy2norm(130, -ripple, handles);
        [x2, y2] = xy2norm(130, -max_y, handles);
        handles.arrows{2} = annotation('arrow', 'Y',[y2 y1], 'X',[x1 x2]);
        set(handles.arrows{2}, 'Color', label_colour);
        text(Fcenter + Fstop + 12, 0, 'A_{pass}', 'BackgroundColor','white', 'EdgeColor','white');
        
        % Stop band
        line([Fcenter+Fstop max_x-10], [-80 -80], 'Color', 'Black');
        line([0 Fcenter-Fstop], [-80 -80], 'Color', 'Black');
        
        line([max_x-10 max_x-10], [-80+ripple -80], 'Color', 'Black');
        line([Fcenter+Fstop Fcenter+Fstop], [-80+ripple -80], 'Color', 'Black');
        line([Fcenter-Fstop Fcenter-Fstop], [-80+ripple -80], 'Color', 'Black');
        
        line([Fcenter-Fstop Fcenter-Fstop], [ripple -100], 'Color', label_colour, 'LineStyle', ':');
        line([Fcenter+Fstop Fcenter+Fstop], [ripple -100], 'Color', label_colour, 'LineStyle', ':');
        line([max_x-10 max_x-10], [-80 -100], 'Color', label_colour, 'LineStyle', ':');
        line([Fcenter Fcenter], [0 -100], 'Color', label_colour, 'LineStyle', ':');
        %        line([150 170], [0 0], 'Color', label_colour, 'LineStyle', ':');
        
        % Labels "0" "Fcenter"
        text(0, -108, '0');
        text(Fcenter - 5, -108, 'F_{center}');
        line([Fcenter Fcenter+Fpass], [-40 -40], 'Color', label_colour, 'LineStyle', ':');
        text(Fcenter + Fstop +5, -40, 'F_{pass}');
        line([Fcenter Fcenter+Fstop], [-60 -60], 'Color', label_colour, 'LineStyle', ':');
        text(Fcenter + Fstop +5, -60, 'F_{stop}');
        text(max_x - 15, -108, 'Fs_{/2}');
        
        % A(stop) label and arrows
        hTest = text(Fcenter/4, -40, 'A_{stop}');
        textExt = get(hTest,'Extent');
        w = textExt(1) + textExt(3)/2;
        [x1, y1] = xy2norm(w, 0, handles);
        [x2, y2] = xy2norm(w, -35, handles);
        handles.arrows{3} = annotation('arrow', 'Y',[y2 y1], 'X',[x1 x2]);
        set(handles.arrows{3}, 'Color', label_colour);
        [x1, y1] = xy2norm(w, -80, handles);
        [x2, y2] = xy2norm(w, -45, handles);
        handles.arrows{4} = annotation('arrow', 'Y',[y2 y1], 'X',[x1 x2]);
        set(handles.arrows{4}, 'Color', label_colour);
        
    case 'Equalize'
        line([0 max_x-10], [0 0], 'Color', 'Black');
        line([max_x-10 max_x-10], [0 ripple], 'Color', 'Black');
        text(0, -108, '0');
        text(max_x - 15, -108, 'Fs_{/2}');
end

guidata(hObject, handles);


function [x1, y1] = xy2norm(x, y, handles)
x_pos = 1;
y_pos = 2;
width = 3;
height = 4;

y_limits = get(handles.magnitude_plot, 'YLim');
x_limits = get(handles.magnitude_plot, 'XLim');
% Position = [left bottom width height]
axesoffsets0 = get(handles.magnitude_plot, 'Position');
axesoffsets1 = get(handles.filter_specs, 'Position');

Figure_Size = get(gcf, 'Position');
y1 = (axesoffsets1(y_pos) + axesoffsets0(y_pos)) / Figure_Size(height);
y2 = axesoffsets0(height) / Figure_Size(height);
y3 = abs((y - y_limits(1)) / abs(y_limits(2) - y_limits(1)));

y1 = y1 + y2 * y3;

x1 = (axesoffsets1(x_pos) + axesoffsets0(x_pos)) / Figure_Size(width) + ...
    axesoffsets0(width)/Figure_Size(width) * ...
    abs((x - x_limits(1)) / abs(x_limits(2) - x_limits(1)));


% --------------------------------------------------------------------
function save_filter_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,path] = uiputfile('*.txt', 'Save filter setup as');
if filename == 0
    return;
else
    newpath = strcat(path,filename);
end

fp = fopen(newpath, 'wt');
fprintf(fp, 'Filter = %d\n', get(handles.filter_type, 'Value'));
fprintf(fp, 'Phase Equalization = %d\n', get(handles.phase_eq, 'Value'));
fprintf(fp, 'Use AD936x FIR = %d\n', get(handles.Use_FIR, 'Value'));

fstop = str2double(get(handles.Fstop, 'String'));
fpass = str2double(get(handles.Fpass, 'String'));
if (handles.freq_units == 2)
    fstop = fstop * 1e6;
    fpass = fpass * 1e6;
end
fprintf(fp, 'Fpass = %d\n', fpass);
fprintf(fp, 'Fstop = %d\n', fstop);
fprintf(fp, 'Apass = %d\n', str2double(get(handles.Apass, 'String')));
fprintf(fp, 'Astop = %d\n', str2double(get(handles.Astop, 'String')));
fprintf(fp, 'Param = %f\n', get(handles.FIR_Astop, 'Value'));

pll_rate = str2double(get(handles.Pll_rate, 'String'));
converter_rate = str2double(get(handles.ADC_clk, 'String'));
data_rate = str2double(get(handles.data_clk, 'String'));
if (handles.freq_units == 2)
    pll_rate = pll_rate * 1e6;
    converter_rate = converter_rate * 1e6;
    data_rate = data_rate * 1e6;
end

fprintf(fp, 'PLL rate  = %d\n', pll_rate);
fprintf(fp, 'Converter = %d\n', converter_rate);
fprintf(fp, 'Data rate = %d\n', data_rate);

fclose(fp);


% --------------------------------------------------------------------
function open_filter_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = uigetfile('*.txt', 'Open Filter as');
if filename == 0
    return;
end
%t = fgets(fp);

reset_input(hObject, handles);
load_settings(hObject, handles);

fp = fopen(filename, 'rt');

set(handles.filter_type, 'Value', fscanf(fp, 'Filter = %d\n'));
set(handles.phase_eq, 'Value', fscanf(fp, 'Phase Equalization = %d\n'));
set(handles.Use_FIR, 'Value', fscanf(fp, 'Use AD936x FIR = %d\n'));

handles.freq_units = 3;
set(handles.Freq_units, 'Value', handles.freq_units);
set(handles.Fpass, 'String', num2str((fscanf(fp, 'Fpass = %d\n') / 1e6)));
set(handles.Fstop, 'String', num2str((fscanf(fp, 'Fstop = %d\n') / 1e6)));
set(handles.Apass, 'String', num2str((fscanf(fp, 'Apass = %e\n'))));
set(handles.Astop, 'String', num2str((fscanf(fp, 'stop = %e\n'))));
set(handles.FIR_Astop, 'Value', fscanf(fp, 'Param = %f\n'));
set(handles.Pll_rate, 'String', num2str((fscanf(fp, 'LL rate = %d\n') / 1e6)));
% FIXME
set(handles.ADC_clk, 'String', num2str((fscanf(fp, 'Converter = %d\n') / 1e6)));
set(handles.data_clk, 'String', num2str((fscanf(fp, 'Data rate = %d\n') / 1e6)));

fclose(fp);
guidata(hObject, handles);


% --- Executes when AD9361_Filter_app is resized.
function AD9361_Filter_app_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to AD9361_Filter_app (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return;

Figure_Size = get(hObject, 'Position');

x_pos = 1;
y_pos = 2;
width = 3;
height = 4;

if ~exist('handles', 'var')
    return;
end

if ~isfield(handles, 'Original_Size')
    return;
end

if (Figure_Size(width) < handles.Original_Size(width)) || (Figure_Size(height) < handles.Original_Size(height))
    set(hObject, 'Position', handles.Original_Size);
    Figure_Size = handles.Original_Size;
end

pos_filter = get(handles.filter_specs,'Position');
pos_mag = get(handles.magnitude_plot,'Position');
pos_cont = get(handles.controls,'Position');
pos_current = get(handles.current_filter,'Position');
pos_logo = get(handles.ADI_logo, 'Position');
pos_help = get(handles.help_button, 'Position');
pos_deeper = get(handles.FVTool_deeper, 'Position');

% 10 from the top
pos_filter(height) = Figure_Size(height) - pos_cont(height) - 4.5;
pos_filter(width) = Figure_Size(width) - pos_current(width) - 5;
set(handles.filter_specs, 'Position', pos_filter);

% 10 on each side
pos_mag(width) = pos_filter(width) - 11;
pos_mag(height) = Figure_Size(height) - pos_cont(height) - 9;
set(handles.magnitude_plot, 'Position', pos_mag);

pos_cont(x_pos) = (Figure_Size(width) - pos_cont(width))/2;
set(handles.controls,'Position', pos_cont);

pos_logo(y_pos) = Figure_Size(height) - pos_logo(height) - .5 ;
set(handles.ADI_logo, 'Position', pos_logo);

pos_help(x_pos) = Figure_Size(width) - pos_help(width) - 2;
pos_help(y_pos) = Figure_Size(height) - pos_help(height) - 1 ;
set(handles.help_button, 'Position', pos_help);

pos_deeper(x_pos) = Figure_Size(width) - pos_deeper(width) - 10;
pos_deeper(y_pos) = pos_filter(x_pos) - 6;
set(handles.FVTool_deeper, 'Position', pos_deeper);

guidata(hObject, handles);
movegui(hObject, 'onscreen');


% --- Executes on button press in phase_eq.
function phase_eq_Callback(hObject, eventdata, handles)
% hObject    handle to phase_eq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    phEQ = 0;
    set(handles.target_delay_label, 'Visible', 'on');
    set(handles.target_delay, 'Visible', 'on');
    set(handles.target_delay, 'String', '0');
    set(handles.target_delay_units, 'Visible', 'on');
else
    phEQ = -1;
    set(handles.target_delay_label, 'Visible', 'off');
    set(handles.target_delay, 'Visible', 'off');
    set(handles.target_delay_units, 'Visible', 'off');
    set(handles.target_delay, 'String', '-1');
end

if (get(handles.filter_type, 'Value') == 1)
    handles.rx.phEQ = phEQ;
else
    handles.tx.phEQ = phEQ;
end

if (handles.active_plot ~= 0)
    handles.active_plot = 0;
end
plot_buttons_off(handles);

dirty(hObject, handles);
set(handles.design_filter, 'Enable', 'on');
guidata(hObject, handles);


% --- Executes on button press in Use_FIR.
function Use_FIR_Callback(hObject, eventdata, handles)
% hObject    handle to Use_FIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.active_plot ~= 0)
    handles.active_plot = 0;
end
plot_buttons_off(handles);

dirty(hObject, handles);
set(handles.design_filter, 'Enable', 'on');


% --- Executes on button press in save2workspace.
function save2workspace_Callback(hObject, eventdata, handles)
% hObject    handle to save2workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.save2workspace, 'Enable', 'off');
drawnow;

if get(handles.filter_type, 'Value') == 1
    filter = handles.rx;
else
    filter = handles.tx;
end

if (~isempty(handles.applycallback))
    handles.applycallback(handles.callbackObj, filter);
else
    if get(handles.filter_type, 'Value') == 1
        assignin('base', 'AD9361_Rx_Filter_object', filter);
    else
        assignin('base', 'AD9361_Tx_Filter_object', filter);
    end
end


% --- Executes on selection change in HB1.
function HB1_Callback(hObject, eventdata, handles)
% hObject    handle to HB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dirty(hObject, handles);
handles = guidata(hObject);

HB = cellstr(get(hObject,'String'));
HB = HB{get(hObject,'Value')};
HB = str2double(HB(1:2));

if get(handles.filter_type, 'Value') == 1
    handles.rx.HB1 = HB;
    handles.rx.HBs = HB * handles.rx.HB2 * handles.rx.HB3;
else
    handles.tx.HB1 = HB;
    handles.tx.HBs = HB * handles.tx.HB2 * handles.tx.HB3;
end

% Update handles structure
data2gui(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function HB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in HBs.
function HBs_Callback(hObject, eventdata, handles)
% hObject    handle to HB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dirty(hObject, handles);
handles = guidata(hObject);

HBs = cellstr(get(hObject,'String'));
HBs = HBs{get(hObject,'Value')};
HBs = str2double(HBs(1:2));
% handles.HBs = HBs;


HBs = factor(HBs);
for i = 1:3;
    if i > length(HBs);
        HBs(i) = 1;
    end
end

% The matrix is sorted since '3' isn't a valid value for HB1 or HB2 so if it
% exists in the set of factors it should be placed last in order to be assigned
% to HB3. Note that this should be changed to check for valid values per field
% if the options for each field changes a lot.
HBs = num2cell(sort(HBs));

[HB1, HB2, HB3] = HBs{:};
if (get(handles.filter_type, 'Value') == 1)
    handles.rx.HB1 = HB1;
    handles.rx.HB2 = HB2;
    handles.rx.HB3 = HB3;
else
    handles.tx.HB1 = HB1;
    handles.tx.HB2 = HB2;
    handles.tx.HB3 = HB3;
end

% reset caldiv so the advanced mode setting is respected
caldiv = get_caldiv(handles);
if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

% Update handles structure
data2gui(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function HBs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HBs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in converter2PLL.
function converter2PLL_Callback(hObject, eventdata, handles)
% hObject    handle to converter2PLL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

HB = cellstr(get(hObject,'String'));
HB = HB{get(hObject,'Value')};
HB = str2double(HB(1:2));

handles.tx.PLL_mult = HB;
handles.rx.PLL_mult = HB;

% reset caldiv so the advanced mode setting is respected
caldiv = get_caldiv(handles);
if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

data2gui(hObject, handles);
handles = guidata(hObject);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function converter2PLL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to converter2PLL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_button.
function help_button_Callback(hObject, eventdata, handles)
if ~isfield(handles, 'helpcallback')
    web('http://wiki.analog.com/resources/eval/user-guides/ad-fmcomms2-ebz/software/filters');
else
    handles.helpcallback(handles.callbackObj);
end


% --- Executes on button press in save2HDL.
function save2HDL_Callback(hObject, eventdata, handles)
% hObject    handle to save2HDL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

error_msg = '';
if ~ license('test','fixed_point_toolbox') || ~ license('checkout','fixed_point_toolbox')
    error_msg = 'You are missing the license for it.';
elseif isempty(ver('fixedpoint'))
    error_msg = 'You do not seem to have it installed.';
end

if error_msg
    choice = questdlg(['Sorry. The AD9361/AD9364 Filter Design Wizard requires the Fixed-Point Designer. You do not seem to have it installed. ' error_msg], ...
        'Error Message', ...
        'More Information','OK','OK');
    switch choice
        case 'More Information'
            web('http://www.mathworks.com/products/fixed-point-designer/');
        case 'OK'
    end
    return
    
end

sel = get_current_rxtx(handles);
Hd = sel.Hmd;
fdhdltool(Hd,numerictype(1,16,15));

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over HB1.
function HB1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to HB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ret = value2Hz(handles, pulldown, value)
if isnan(value)
    value = 0;
end

switch pulldown
    case 1
        ret = value;
    case 2
        ret = value * 1e3;
    case 3
        ret = value * 1e6;
end

function ret = Hz2value(handles, pulldown, value)
switch pulldown
    case 1
        ret = value;
    case 2
        ret = value / 1e3;
    case 3
        ret = value / 1e6;
end

function put_data_clk(handles, data_clk)
data_clk = Hz2value(handles, handles.freq_units, data_clk);
set(handles.data_clk, 'String', num2str(data_clk, 10));

function caldiv = default_caldiv(handles)
if (get(handles.filter_type, 'Value') == 1)
    wnom = 1.4 * handles.rx.Fstop;  % Rx
else
    wnom = 1.6 * handles.tx.Fstop;  % Tx
end

pll = get_pll_rate(handles);
div = ceil((pll/wnom)*(log(2)/(2*pi)));
caldiv = min(max(div,1),511);

function caldiv = get_caldiv(handles)
if (get(handles.filter_type, 'Value') == 1)
    wnom = 1.4 * handles.rx.Fstop;  % Rx
else
    wnom = 1.6 * handles.tx.Fstop;  % Tx
end

Fcutoff = str2double(get(handles.Fcutoff, 'String'));

if Fcutoff
    wnom = value2Hz(handles, handles.freq_units, Fcutoff);
end

pll = get_pll_rate(handles);
div = ceil((pll/wnom)*(log(2)/(2*pi)));
caldiv = min(max(div,1),511);

function set_fcutoff(handles, caldiv)
wc = (get_pll_rate(handles) / caldiv)*(log(2)/(2*pi));
set(handles.Fcutoff, 'String', num2str(Hz2value(handles, handles.freq_units, wc)));

function pll = get_pll_rate(handles)
if (get(handles.filter_type, 'Value') == 1)
    % Rx
    pll = handles.rx.Rdata * handles.rx.FIR * handles.rx.HB1 * ...
        handles.rx.HB2 * handles.rx.HB3 * handles.rx.PLL_mult;
else
    % Tx
    pll = handles.tx.Rdata * handles.tx.FIR * handles.tx.HB1 * ...
        handles.tx.HB2 * handles.tx.HB3 * handles.tx.DAC_div * ...
        handles.tx.PLL_mult;
end

% calculate a channel's complex bandwidth that matches 32 bit integer precision
% on the driver
function rfbw_hw = get_rfbw_hw(handles, caldiv)
sel = get_current_rxtx(handles);
pll_rate = get_pll_rate(handles);
[rfbw_hw, caldiv] = calculate_rfbw(pll_rate, caldiv, sel.RxTx, true);

% calculate a channel's full precision complex bandwidth
function rfbw = get_rfbw(handles, caldiv)
sel = get_current_rxtx(handles);
pll_rate = get_pll_rate(handles);
[rfbw, caldiv] = calculate_rfbw(pll_rate, caldiv, sel.RxTx, false);

if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

function Fcutoff_Callback(hObject, eventdata, handles)
% hObject    handle to Fcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

caldiv = get_caldiv(handles);
if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

set_fcutoff(handles, caldiv);

data2gui(hObject, handles);
handles = guidata(hObject);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Fcutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function target_delay_Callback(hObject, eventdata, handles)
% hObject    handle to target_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(handles.filter_type, 'Value') == 1)
    handles.rx.phEQ = str2double(get(hObject, 'String'));
else
    handles.tx.phEQ = str2double(get(hObject, 'String'));
end
data2gui(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function target_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FVTool_datarate.
function FVTool_datarate_Callback(hObject, eventdata, handles)
% hObject    handle to FVTool_datarate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel = get_current_rxtx(handles);
converter_rate = sel.Rdata * sel.FIR * sel.HB1 * sel.HB2 * sel.HB3;

str = sprintf('%s Filter\nFpass = %g MHz; Fstop = %g MHz\nApass = %g dB; Astop = %g dB', sel.RxTx, sel.Fpass/1e6, sel.Fstop/1e6, sel.Apass, sel.Astop);

if get(handles.use_FPGAfilter, 'Value')== 1
    hfvt3 = fvtool(handles.analogfilter,handles.Hmiddle,handles.grpdelaycal,handles.plutofilter,...
        'FrequencyRange','Specify freq. vector', ...
        'FrequencyVector',linspace(0,sel.Rdata/2,2048),'Fs',...
        converter_rate, ...
        'ShowReference','off','Color','White');
    set(hfvt3, 'Color', [1 1 1]);
    set(hfvt3.CurrentAxes, 'YLim', [-100 20]);
    legend(hfvt3, 'Analog','Analog + Half Band','Analog + HB + FIR','Analog + HB + FIR + Pluto');
else
    hfvt3 = fvtool(handles.analogfilter,handles.Hmiddle,handles.grpdelaycal,...
        'FrequencyRange','Specify freq. vector', ...
        'FrequencyVector',linspace(0,sel.Rdata/2,2048),'Fs',...
        converter_rate, ...
        'ShowReference','off','Color','White');
    set(hfvt3, 'Color', [1 1 1]);
    set(hfvt3.CurrentAxes, 'YLim', [-100 20]);
    legend(hfvt3, 'Analog','Analog + Half Band','Analog + HB + FIR');
end
text(0.5, 10,...
    str,...
    'BackgroundColor','white',...
    'EdgeColor','red');

hfvt4 = fvtool(...
    sel.Hmd,...
    'FrequencyRange','Specify freq. vector', ...
    'FrequencyVector',linspace(0,sel.Rdata/2,2048),'Fs',...
    sel.Rdata*sel.FIR, ...
    'ShowReference','off','Color','White');
set(hfvt4.CurrentAxes, 'YLim', [-100 20]);
legend(hfvt4, 'FIR Filter');

% add the quantitative values about FIR magnitude
[h,~] = freqz(sel.Hmd,1024);
maxmag = max(20*log10(abs(h)));

[gd,~] = grpdelay(handles.grpdelaycal,2048);
I = round(sel.Fpass/(converter_rate/2)*2048);
gd2 = gd(1:I).*(1/converter_rate);
str2 = sprintf('Delay Variance = %g ns', handles.grpdelayvar*1e9);

hfvt0 = fvtool(handles.grpdelaycal,...
    'FrequencyRange','Specify freq. vector', ...
    'FrequencyVector',linspace(0,sel.Fpass,2048),...
    'Fs',converter_rate,'Analysis','grpdelay');
hfvt0.GroupDelayUnits = 'Time';
text(0.1,(mean(gd2))*1e6,...
    str2,...
    'BackgroundColor','white',...
    'EdgeColor','red');

function show_advanced(handles)
set(handles.phase_eq, 'Visible', 'on');
if get(handles.phase_eq, 'Value')
    set(handles.target_delay_label, 'Visible', 'on');
    set(handles.target_delay, 'Visible', 'on');
    set(handles.target_delay_units, 'Visible', 'on');
end
set(handles.Advanced_options, 'Value', 1);
set(handles.Fcutoff_label, 'Visible', 'on');
set(handles.RFbw_label, 'Visible', 'on');
set(handles.Fcutoff, 'Visible', 'on');
set(handles.RFbw, 'Visible', 'on');
set(handles.Use_FIR, 'Visible', 'on');
set(handles.FIR_Astop, 'Visible', 'on');
set(handles.FIR_Astop_label, 'Visible', 'on');

set(handles.HBs_label, 'Visible', 'off');
set(handles.HBs, 'Visible', 'off');
set(handles.HBs_rate, 'Visible', 'off');

set(handles.HB1_label, 'Visible', 'on');
set(handles.HB1, 'Visible', 'on');
set(handles.HB1_rate, 'Visible', 'on');
set(handles.HB2_label, 'Visible', 'on');
set(handles.HB2, 'Visible', 'on');
set(handles.HB2_rate, 'Visible', 'on');
set(handles.HB3_label, 'Visible', 'on');
set(handles.HB3, 'Visible', 'on');
set(handles.HB3_rate, 'Visible', 'on');


function hide_advanced(handles)
set(handles.Advanced_options, 'Value', 0);
set(handles.phase_eq, 'Value', 0);
set(handles.Use_FIR, 'Value', 1);
set(handles.target_delay, 'String', '0');
set(handles.Fcutoff, 'String', '0');

set(handles.phase_eq, 'Visible', 'off');
set(handles.target_delay_label, 'Visible', 'off');
set(handles.target_delay, 'Visible', 'off');
set(handles.target_delay_units, 'Visible', 'off');
set(handles.Fcutoff_label, 'Visible', 'off');
set(handles.RFbw_label, 'Visible', 'off');
set(handles.Fcutoff, 'Visible', 'off');
set(handles.RFbw, 'Visible', 'off');
set(handles.Use_FIR, 'Visible', 'off');
set(handles.FIR_Astop, 'Visible', 'off');
set(handles.FIR_Astop_label, 'Visible', 'off');

set(handles.HBs_label, 'Visible', 'on');
set(handles.HBs, 'Visible', 'on');
set(handles.HBs_rate, 'Visible', 'on');

set(handles.HB1_label, 'Visible', 'off');
set(handles.HB1, 'Visible', 'off');
set(handles.HB1_rate, 'Visible', 'off');
set(handles.HB2_label, 'Visible', 'off');
set(handles.HB2, 'Visible', 'off');
set(handles.HB2_rate, 'Visible', 'off');
set(handles.HB3_label, 'Visible', 'off');
set(handles.HB3, 'Visible', 'off');
set(handles.HB3_rate, 'Visible', 'off');


% --- Executes on button press in Advanced_options.
function Advanced_options_Callback(hObject, eventdata, handles)
% hObject    handle to Advanced_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if get(hObject,'Value')
    show_advanced(handles);
    
    % reset fcutoff to the default value when advanced is re-enabled
    set(handles.Fcutoff, 'String', '0');
    caldiv = get_caldiv(handles);
    set_fcutoff(handles, caldiv);
else
    hide_advanced(handles);
end

% reset caldiv for Rx/Tx to defaults, forces advanced state to be respected for
% both channels when data2gui is run
if isstruct(get_current_rxtx(handles))
    filter_type = get(handles.filter_type, 'Value');
    
    set(handles.filter_type, 'Value', 1);
    caldiv = default_caldiv(handles);
    handles.rx.caldiv = caldiv;
    set(handles.filter_type, 'Value', 2);
    caldiv = default_caldiv(handles);
    handles.tx.caldiv = caldiv;
    
    set(handles.filter_type, 'Value', filter_type);
end

% remove generated filter taps if they exist, forces both Rx/Tx filters to be
% redesigned when advanced is toggled
if isfield(handles, 'rfirtaps')
    handles = rmfield(handles, 'rfirtaps');
end
if isfield(handles, 'tfirtaps')
    handles = rmfield(handles, 'tfirtaps');
end

dirty(hObject, handles);
guidata(hObject, handles);
data2gui(hObject, handles);

% --- Executes on button press in connect2target.
function connect2target_Callback(hObject, eventdata, handles)
% hObject    handle to connect2target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ip_address = get(handles.IP_num, 'String');
set(handles.connect2target, 'Enable', 'off');
set(handles.connect2target, 'String', 'Connecting to Target');
drawnow;

% If the libiio is already initialized delete the libiio_if object
if (~isempty(handles.libiio_ctrl_dev))
    delete(handles.libiio_ctrl_dev);
end

% Add libiio sys object library to search path
% (assumes we're running in the full repo checkout)
[pathstr, ~, ~] = fileparts(mfilename('fullpath'));
if exist(fullfile(pathstr, 'libiio', 'libiio_if.m'), 'file')
    addpath(fullfile(pathstr, 'libiio'));
    
    % Initialize the libiio_if object
    handles.libiio_ctrl_dev = libiio_if();
    [ret, err_msg, msg_log] = init(handles.libiio_ctrl_dev, ip_address, ...
        'ad9361-phy', '', 0, 0);
    fprintf('%s', msg_log);
else
    err_msg = ['The libiio_if object was not found. ' ...
        'Note that if the application is being run from a ' ...
        'zip snapshot the ''libiio'' bindings will be missing ' ...
        'from the related directory. Either recursively clone the ' ...
        'repository using git or use the latest official installer ' ...
        'for the filter wizard.'];
    ret = -1;
end

if (ret < 0)
    set(handles.target_get_clock, 'Enable', 'off');
    set(handles.connect2target, 'Enable', 'on');
    set(handles.connect2target, 'String', 'Connect to Target');
    if (~isempty(handles.libiio_ctrl_dev))
        delete(handles.libiio_ctrl_dev);
    end
    handles.libiio_ctrl_dev = {};
    msgbox(err_msg, 'Error','error');
else
    set(handles.connect2target, 'String', 'Connected to Target');
    set(handles.target_get_clock, 'Enable', 'on');
    
    if isfield(handles, 'rfirtaps') && isfield(handles, 'tfirtaps')
        set(handles.save2target, 'Enable', 'on');
    end
    
    % save IP address to restore on next startup
    [pathstr, ~, ~] = fileparts(mfilename('fullpath'));
    cached_ip_file = fullfile(pathstr, '.previous_ip_addr');
    fd = fopen(cached_ip_file, 'wt');
    fprintf(fd, get(handles.IP_num, 'String'));
    fclose(fd);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in design_filter.
function design_filter_Callback(hObject, eventdata, handles)
% hObject    handle to design_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
create_filter(hObject, handles);


function Port_num_Callback(hObject, eventdata, handles)
% hObject    handle to Port_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Port_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Port_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function magnitude_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magnitude_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate magnitude_plot

% --- Executes during object creation, after setting all properties.
function FIR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DAC_by2.
function DAC_by2_Callback(hObject, eventdata, handles)
% hObject    handle to DAC_by2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

HB = cellstr(get(hObject,'String'));
HB = HB{get(hObject,'Value')};
HB = str2double(HB(1:2));

handles.tx.DAC_div = HB;

% reset caldiv so the advanced mode setting is respected
caldiv = get_caldiv(handles);
if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

data2gui(hObject, handles);
handles = guidata(hObject);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function DAC_by2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DAC_by2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FIR.
function FIR_Callback(hObject, eventdata, handles)
% hObject    handle to FIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

HB = cellstr(get(hObject,'String'));
HB = HB{get(hObject,'Value')};
HB = str2double(HB(1:2));

if get(handles.filter_type, 'Value') == 1
    handles.rx.FIR = HB;
else
    handles.tx.FIR = HB;
end

% reset caldiv so the advanced mode setting is respected
caldiv = get_caldiv(handles);
if get(handles.filter_type, 'Value') == 1
    handles.rx.caldiv = caldiv;
else
    handles.tx.caldiv = caldiv;
end

data2gui(hObject, handles);
handles = guidata(hObject);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function results_taps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_taps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Pll_rate.
function Pll_rate_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Pll_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close AD9361_Filter_app.
function AD9361_Filter_app_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to AD9361_Filter_app (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on selection change in which_device.
function which_device_Callback(hObject, eventdata, handles)
% hObject    handle to which_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.design_filter, 'Enable', 'off');
if get(handles.which_device,'Value')==1 || get(handles.which_device,'Value')==2
    set(handles.use_FPGAfilter, 'Value', 0);
    set(handles.FPGA_label, 'Visible', 'off');
    set(handles.use_FPGAfilter, 'Visible', 'off');
    set(handles.FPGA_rate, 'Visible', 'off');
    freq_on(handles);
end


% --- Executes during object creation, after setting all properties.
function which_device_CreateFcn(hObject, eventdata, handles)
% hObject    handle to which_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Saved_Filters.
function Saved_Filters_Callback(hObject, eventdata, handles)
% hObject    handle to Saved_Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_settings(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Saved_Filters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Saved_Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in store_filter.
function store_filter_Callback(hObject, eventdata, handles)
% hObject    handle to store_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = 'ad9361_settings.mat';

if ~ exist(filename, 'file')
    errordlg('I can not find the required files, must be some sort of installation error', ...
        'File Error');
    set(handles.Saved_Filters, 'Visible', 'off');
    return;
end

answer = char(inputdlg('Save As?'));
% exit if the cancel button is hit
if isempty(answer)
    return;
end

options = load(filename);

button = 'Replace';
if isfield(options.ad9361_settings.rx, answer)
    button = questdlg(strcat('Rx setting "', answer, '" exists, replace?'),...
        'Replace', 'Replace', 'Cancel', 'Cancel');
end
if strcmp(button, 'Replace')
    options.ad9361_settings.rx.(answer) = handles.rx;
end

button = 'Replace';
if isfield(options.ad9361_settings.tx, answer)
    button = questdlg(strcat('Tx setting "', answer, '" exists, replace?'),...
        'Replace', 'Replace', 'Cancel', 'Cancel');end
if strcmp(button, 'Replace')
    options.ad9361_settings.tx.(answer) = handles.tx;
end

ad9361_settings = options.ad9361_settings;
save(filename, 'ad9361_settings');
reset_input(hObject, handles);

str = get(handles.Saved_Filters, 'String');
for i = 1:length(str);
    if strcmp(strcat(answer, ' (Rx & Tx)'), str(i))
        set(handles.Saved_Filters, 'Value', i);
        break;
    end
end
load_settings(hObject, handles);


% --- Executes on button press in LockRxTx.
function LockRxTx_Callback(hObject, eventdata, handles)
% hObject    handle to LockRxTx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in HB2.
function HB2_Callback(hObject, eventdata, handles)
% hObject    handle to HB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

HB = cellstr(get(hObject,'String'));
HB = HB{get(hObject,'Value')};
HB = str2double(HB(1:2));

if get(handles.filter_type, 'Value') == 1
    handles.rx.HB2 = HB;
    handles.rx.HBs = handles.rx.HB1 * HB * handles.rx.HB3;
else
    handles.tx.HB2 = HB;
    handles.tx.HBs = handles.tx.HB1 * HB * handles.tx.HB3;
end

guidata(hObject, handles);
data2gui(hObject, handles);


% --- Executes during object creation, after setting all properties.
function HB2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in HB3.
function HB3_Callback(hObject, eventdata, handles)
% hObject    handle to HB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirty(hObject, handles);
handles = guidata(hObject);

HB = cellstr(get(hObject,'String'));
HB = HB{get(hObject,'Value')};
HB = str2double(HB(1:2));

if get(handles.filter_type, 'Value') == 1
    handles.rx.HB3 = HB;
    handles.rx.HBs = handles.rx.HB1 * handles.rx.HB2 * HB;
else
    handles.tx.HB3 = HB;
    handles.tx.HBs = handles.tx.HB1 * handles.tx.HB2 * HB;
end

data2gui(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function HB3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in Response_Type.
function Response_Type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Response_Type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
h = get(eventdata.OldValue,'String');
switch h
    case 'Lowpass'
        set(handles.LP_MagSpecs, 'Visible', 'off');
    case 'Root Raised Cosine'
        set(handles.RRC_MagSpecs, 'Visible', 'off');
    case 'Bandpass'
        set(handles.LP_MagSpecs, 'Visible', 'off');
        set(handles.Fcenter_label, 'Visible', 'off');
        set(handles.Fcenter, 'Visible', 'off');
    case 'Equalize'
        set(handles.Freq_Specs, 'Visible', 'on');
end

h = get(eventdata.NewValue,'String');
switch h
    case 'Lowpass'
        set(handles.LP_MagSpecs, 'Visible', 'on');
    case 'Root Raised Cosine'
        set(handles.RRC_MagSpecs, 'Visible', 'on');
    case 'Bandpass'
        set(handles.LP_MagSpecs, 'Visible', 'on');
        set(handles.Fcenter_label, 'Visible', 'on');
        set(handles.Fcenter, 'Visible', 'on');
    case 'Equalize'
        set(handles.Freq_Specs, 'Visible', 'off');
end

display_default_image(hObject);

function Fcenter_Callback(hObject, eventdata, handles)
% hObject    handle to Fcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Fcenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function results_Apass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_Apass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% vim: set et sw=4 ts=4 ft=matlab:


% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15


% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pluto.
function pluto_Callback(hObject, eventdata, handles)
% hObject    handle to pluto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pluto contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pluto
set(handles.design_filter, 'Enable', 'on');

% --- Executes during object creation, after setting all properties.
function pluto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pluto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function FIR_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FIR_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on FVTool_deeper and none of its controls.
function FVTool_deeper_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to FVTool_deeper (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
