function varargout = task5(varargin)
% TASK5 MATLAB code for task5.fig
%      TASK5, by itself, creates a new TASK5 or raises the existing
%      singleton*.
%
%      H = TASK5 returns the handle to a new TASK5 or the handle to
%      the existing singleton*.
%
%      TASK5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TASK5.M with the given input arguments.
%
%      TASK5('Property','Value',...) creates a new TASK5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before task5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to task5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help task5

% Last Modified by GUIDE v2.5 09-May-2019 10:12:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @task5_OpeningFcn, ...
                   'gui_OutputFcn',  @task5_OutputFcn, ...
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

%Code for Signal

% --- Executes just before task5 is made visible.
function task5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to task5 (see VARARGIN)

%% More Signal Code
global zh ph Lresp Nfft b a z h Y Oresp p fmax signal outResp online matResp conjugate
global ax1 ax2 
 
fmax = 44100;
signal = [];
outResp = plot([]);
online = 1;
conjugate=0;
    
 z = [0]' ;
        p = [1/2+1/2*j 1/2-1/2*j]' ;
       
        set(0,'defaultaxesfontsize',10)
 
        axes(handles.ax_zp)
 
        [zh,ph,cruff]=zplaneplot(z,p);
        %set(cruff,'hittest','off')
        ax1 = gca;
        ylabel('Imaginary Part')
        xlabel('Real Part')
 
        
 
        axes(handles.ax_response)
 
        [b,a]=zp2tf(z,p,1);
 
        Nfft = 512;
 
        Y = fft(b,Nfft)./fft(a,Nfft);
         
        Lresp = plot((0:Nfft-1)/Nfft*2-1, 20*log10(fftshift(abs(Y))),'linewidth',2,'color',[0 0 0]);
        %Lresp = plot((0:Nfft-1)/Nfft*2-1, fftshift(abs(Y)),'linewidth',2);
        ax2 = gca;
        
        axes(handles.ax_our_response);
        gain = Get_gain_manual(z,p,1);
        Oresp = plot(linspace(0,fmax,2001),gain,'linewidth',2,'color',[0,0,0]);
                
        set(ax2,'xlim',[0 1])
        %set(ax2,'ylim',[max(min(Lresp),0.1) max(Lresp)])
         
        set(ax2,'xtick', [0:1/4:1])
        set(ax2,'xticklabel', {'0','\pi/4','\pi/2','2\pi/3','\pi'})
        %get(ax2,'xticklabel')
        grid on
        xlabel('Frequency (\Omega)')
        ylabel('Magnitude (dB) ')
         
        %set(Lresp,'erasemode','xor')
 
        set(zh,'buttondownfcn',@zeroClick,...
            'markersize',8,'linewidth',1)
        set(ph,'buttondownfcn',@poleClick,...
            'markersize',8,'linewidth',1)
        
        axes(handles.ax_matlab_resp);
        gain = Rania_Gain(z,p);
        matResp = plot(linspace(0,fmax,2001),gain);
        
       
        
        
% Choose default command line output for task5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes task5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% More Signal Code
       

% --- Outputs from this function are returned to the command line.
function varargout = task5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cb_conjugate.
function cb_conjugate_Callback(hObject, eventdata, handles)
% hObject    handle to cb_conjugate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_conjugate
global conjugate
conjugate = get(handles.cb_conjugate,'value');


% --- Executes on mouse press over axes background.
function ax_zp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ax_zp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function zeroClick(object, eventdata)
        global zh ax1 ax2 Lresp conjugate
        set(gcf,'userdata','')
        set(gcf,'windowbuttonmotionfcn','set(gcf,''userdata'',''motion'')') %drag
        set(gcf,'windowbuttonupfcn','set(gcf,''userdata'',''up'')')%mouse is left
 
        ind = find(zh==gco);
        %set(zh(ind),'erasemode','xor')
      %  set(Lresp,'erasemode','xor')
        pair = conjugate;
        done = 0;
        
        %remove zero 
        if eventdata.Button == 3%right click
           delete(zh(ind));
           zh(ind) = [];
           recompute(1);
           return;
        end
        
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
                    if pair
                        set(zh(ind),'xdata',[pt(1) pt(1)],'ydata',[pt(2) -pt(2)])
                    else
                        set(zh(ind),'xdata',pt(1),'ydata',pt(2))
                    end
 
                    recompute(1)
                case 'up'
                    done = 1;
            end
            set(gcf,'userdata','')
        end
        set(gcf,'windowbuttonmotionfcn','')
        set(gcf,'windowbuttonupfcn','')
        set(ax2,'ylimmode','auto')
        ylim = get(ax2,'ylim');
        Y = get(Lresp,'ydata');
        if ylim(1)>min(Y)-3 | ylim(2)<max(Y)+3
            set(ax2,'ylim',[min(Y)-3 max(Y)+3])
        end
        
    function poleClick(object,eventdata)
        global zh ph Lresp Nfft b a z h Y conjugate
        global ax1 ax2 ax3 ax4 z_surface_CameraPos z_surface_CameraUpVec surface_display_opts
        set(gcf,'userdata','')
        set(gcf,'windowbuttonmotionfcn','set(gcf,''userdata'',''motion'')')
        set(gcf,'windowbuttonupfcn','set(gcf,''userdata'',''up'')')
 
        ind = find(ph==gco);
        if eventdata.Button == 3
           delete(ph(ind));
           ph(ind) = [];
           recompute(1);
           return;
        end
       % set(ph(ind),'erasemode','xor')
        %set(Lresp,'erasemode','xor')
        pair = conjugate;
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
                    if pair
                        set(ph(ind),'xdata',[pt(1) pt(1)],'ydata',[pt(2) -pt(2)])
                    else
                        set(ph(ind),'xdata',pt(1),'ydata',pt(2))
                    end
 
                    recompute(1)
                case 'up'
                    done = 1;
            end
            set(gcf,'userdata','')
        end
        set(gcf,'windowbuttonmotionfcn','')
        set(gcf,'windowbuttonupfcn','')
        %set(ph(ind),'erasemode','normal')
        %set(Lresp,'erasemode','normal')
        set(ax2,'ylimmode','auto')
        Y = get(Lresp,'ydata');
        ylim = get(ax2,'ylim');
        if ylim(1)>min(Y)-3 | ylim(2)<max(Y)+3
            set(ax2,'ylim',[min(Y)-3 max(Y)+3])
        end
        
   function recompute(args)
        global zh ph Lresp Nfft b a z h Y p Oresp signal inpResp outResp conjugate online matResp
        global ax1 ax2 ax3 ax4 z_surface_CameraPos z_surface_CameraUpVec surface_display_opts
        z = [];
        p = [];
        b = 1;
        a = 1;
        for i=1:length(zh)
            zx = get(zh(i),'xdata');
            zy = get(zh(i),'ydata');
            if length(zx)==1
                b = conv(b,[1 -(zx+sqrt(-1)*zy)]);
            else
                b = conv(b,[1 -2*zx(1) zx(1).^2+zy(1).^2]);
            end
            z = [z zx+sqrt(-1)*zy];
 
        end
        for i=1:length(ph)
            px = get(ph(i),'xdata');
            py = get(ph(i),'ydata');
            if length(px)==1
                a = conv(a,[1 -(px+sqrt(-1)*py)]);
            else
                a = conv(a,[1 -2*px(1) px(1).^2+py(1).^2]);
            end
            p = [p px+sqrt(-1)*py];
        end
 
        Y = fft(b,Nfft)./fft(a,Nfft);
        %Y = Y/max(abs(Y));
        set(Lresp,'ydata',20*log10(fftshift(abs(Y))))
        gain = Get_gain_manual(transp(z),transp(p),1);
        set(Oresp,'ydata',gain);
        gain = Rania_Gain(transp(z),transp(p));
        gain = real(gain);
        set(matResp,'ydata',gain);
        
        if ~online
        filteredSignal = Digital_filter(z,p,signal);
        set(outResp,'ydata',filteredSignal);
        end
       


% --- Executes on button press in btn_addzero.
function btn_addzero_Callback(hObject, eventdata, handles)
% hObject    handle to btn_addzero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zh ph Lresp Nfft b a z h Y
        global ax1 
if length(zh)>0
            zh(end+1) =  line(.5,0,'parent',ax1,'buttondownfcn',@zeroClick,...
                'markersize',8,'linewidth',1,'marker','o','linestyle','none');
        else
            zh = line(.5,0,'parent',ax1,'buttondownfcn',@zeroClick,...
                'markersize',8,'linewidth',1,'marker','o','linestyle','none');
            set(findobj('string','Remove Zeros'),'enable','on')
        end
        recompute(1);

% --- Executes on button press in btn_addpole.
function btn_addpole_Callback(hObject, eventdata, handles)
% hObject    handle to btn_addpole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zh ph Lresp Nfft b a z h Y
        global ax1 
if length(ph)>0
            ph(end+1) =  line(.5,0,'parent',ax1,'buttondownfcn',@poleClick,...
                'markersize',8,'linewidth',1,'marker','x','linestyle','none');
        else
            ph =  line(.5,0,'parent',ax1,'buttondownfcn',@poleClick,...
                'markersize',8,'linewidth',1,'marker','x','linestyle','none');
            set(findobj('string','Remove Poles'),'enable','on')
        end

        recompute(1);


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global signal z p inpResp outResp online
    [FileName,FilePath]= uigetfile();
    ExPath = fullfile(FilePath, FileName);
    var = load(ExPath);   
    signal = var.val;
    axes(handles.ax_input);
    inpResp = plot(signal);
    axes(handles.ax_output);
    filtered_signal = Digital_filter(transp(z),transp(p),signal);
    outResp = plot(filtered_signal);
    if online
        online_processing();
    end
    

% --- Executes on button press in offline.
function offline_Callback(hObject, eventdata, handles)
% hObject    handle to offline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of offline
global online
online = get(handles.online,'value');


% --- Executes on button press in online.
function online_Callback(hObject, eventdata, handles)
% hObject    handle to online (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of online
global online
online = get(handles.online,'value');

function online_processing()
    global p z signal outResp inpResp online
    inputBuffer = zeros(20,1);
    outputBuffer = zeros(20,1);
    pz = poly(z);
    pp = poly(p);
   
    for i = 1:length(signal)
        out = 0;
        %update z/p if changed in z plane
        pz = poly(z);
        pp = poly(p);
        
        inputBuffer = [signal(1,i); inputBuffer];
        %diff casual equation
        for j = 1:length(pz)
           out = out + pz(1,j)*inputBuffer(j); 
        end
        
        for j = 2:length(pp)
           out = out - pp(1,j)*outputBuffer(j-1); 
        end
        out = out * 1/pp(1,1);%divide a0
       
        set(outResp,'ydata',flipud(outputBuffer));%flipping 3shan y7sl display s7
        set(inpResp,'ydata',flipud(inputBuffer));
        drawnow
        if ~online
            return;
        end
        outputBuffer = [out; outputBuffer];
        
    end



function et_fs_Callback(hObject, eventdata, handles)
% hObject    handle to et_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_fs as text
%        str2double(get(hObject,'String')) returns contents of et_fs as a double
global fmax  Oresp matResp
fmax = str2double(get(handles.et_fs, 'String'))/2
set(Oresp,'xdata',linspace(0,fmax,2001));
set(matResp,'xdata',linspace(0,fmax,2001));

% --- Executes during object creation, after setting all properties.
function et_fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
