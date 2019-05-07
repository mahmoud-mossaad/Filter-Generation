function varargout = untitled(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 04-May-2019 17:11:26

% Begin initialization code - DO NOT EDIT
global zh ph z p pt
global ax1 surface_display_opts
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
    action = 'init';
    z = [0 -1]' ;
    p = [1/2+1/2*j 1/2-1/2*j]' ;
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

switch action
    case 'init'
        z_surface_CameraPos=[3.7719  -15.7111  275.3541];
        z_surface_CameraUpVec=[-0.1237    0.5153   23.1286];
        surface_display_opts = 0;
        set(0,'defaultaxesfontsize',10)
  
        [zh,ph,cruff]=zplaneplot(z,p);
        ax1 = gca;
        ylabel('Imaginary Part')
        xlabel('Real Part')
 
 
        set(zh,'buttondownfcn','zpgui(''zeroclick'')',...
            'markersize',8,'linewidth',1)
        set(ph,'buttondownfcn','zpgui(''poleclick'')',...
            'markersize',8,'linewidth',1)
        
 
    case 'zeroclick'
 
        set(gcf,'userdata','')
        set(gcf,'windowbuttonmotionfcn','set(gcf,''userdata'',''motion'')')
        set(gcf,'windowbuttonupfcn','set(gcf,''userdata'',''up'')')
 
        ind = find(zh==gco);
        %set(zh(ind),'erasemode','xor')
      %  set(Lresp,'erasemode','xor')
        pair = length(get(zh(ind),'xdata'))==2;
        done = 0;
 
        pt = get(ax1,'currentpoint');
        pt = pt(1,1:2);
        title(['selected position: ' num2str(pt) 'j'])
        while ~done
            waitfor(gcf,'userdata')
            switch get(gcf,'userdata')
                case 'motion'
                    pt = get(ax1,'currentpoint');
                    pt = pt(1,1:2);
                    %title(['selected position: ' num2str(pt) 'j'])
                    %if pair
                    if 1
                        set(zh(ind),'xdata',[pt(1) pt(1)],'ydata',[pt(2) -pt(2)])
                    else
                        set(zh(ind),'xdata',pt(1),'ydata',pt(2))
                    end
                   
 
                    zpgui('recompute')
                case 'up'
                    done = 1;
            end
            set(gcf,'userdata','')
        end
        set(gcf,'windowbuttonmotionfcn','')
        set(gcf,'windowbuttonupfcn','')
        %set(zh(ind),'erasemode','normal')
        %set(Lresp,'erasemode','normal')
        %set(ax2,'ylimmode','auto')
        %ylim = get(ax2,'ylim');
        %Y = get(Lresp,'ydata');
        %if ylim(1)>min(Y)-3 | ylim(2)<max(Y)+3
         %   set(ax2,'ylim',[min(Y)-3 max(Y)+3])
        %end
 
    case 'poleclick'
        
        set(gcf,'userdata','')
        set(gcf,'windowbuttonmotionfcn','set(gcf,''userdata'',''motion'')')
        set(gcf,'windowbuttonupfcn','set(gcf,''userdata'',''up'')')
 
        ind = find(ph==gco);
       % set(ph(ind),'erasemode','xor')
        %set(Lresp,'erasemode','xor')
        pair = length(get(ph(ind),'xdata'))==2;
        done = 0;
 
        pt = get(ax1,'currentpoint');
        pt = pt(1,1:2);
        title(['selected position: ' num2str(pt) 'j'])
        while ~done
            waitfor(gcf,'userdata')
            switch get(gcf,'userdata')
                case 'motion'
                    pt = get(ax1,'currentpoint');
                    pt = pt(1,1:2);
                    title(['selected position: ' num2str(pt) 'j'])
                    if 1
                        set(ph(ind),'xdata',[pt(1) pt(1)],'ydata',[pt(2) -pt(2)])
                    else
                        set(ph(ind),'xdata',pt(1),'ydata',pt(2))
                    end
 
                    zpgui('recompute')
                case 'up'
                    done = 1;
            end
            set(gcf,'userdata','')
        end

    
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
