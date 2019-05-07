function varargout = zpgui(varargin)
% ZPGUI  Zero Pole dragging Graphic User Interface
% Allows you to add, remove, and move the zeros and
% poles of a filter interactively.
% To begin, type:
%    zpgui
 
%   Author: Tom Krauss, 9/1/98
%   Adapted by David Dorran (2011) to force poles and zeros to be conjugate
%   pairs (maybe not the most flexible feature but it saves me time). 
%   Updated in 2012 to plot the z-surface.

%   J.W. Peltenburg, 18 dec 2005
%   Fixed some warnings and errors concerning older version of MATLAB, this
%   script is now compatible with r2014b.
%   Added impulse response to the lower right corner
%   Changed 3d mesh to 3d surface plot
%   Removed pushbutton to remove 3d view area fill
%   Other minor cosmetic changes were made

%   J.W. Peltenburg, 14 december 2015
%   Changed "normalised" frequency from [0...1] to [0...2pi]

 
global zh ph Lresp Nfft b a z h Y
global ax1 ax2 ax3 ax4 z_surface_CameraPos z_surface_CameraUpVec surface_display_opts
 
if nargin == 0
    action = 'init';
    z = [0 -1]' ;
    p = [1/2+1/2*j 1/2-1/2*j]';
elseif nargin >= 3
    p = varargin{1}';
    z = varargin{2}';
    action = 'init';
else
        action = varargin{1};
end
if nargin == 3
    jpg_filename = varargin{3};
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
        
    case 'addzero'
        if length(zh)>0
            zh(end+1) =  line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''zeroclick'')',...
                'markersize',8,'linewidth',1,'marker','o','linestyle','none');
        else
            zh = line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''zeroclick'')',...
                'markersize',8,'linewidth',1,'marker','o','linestyle','none');
            set(findobj('string','Remove Zeros'),'enable','on')
        end
        zpgui('recompute')
    case 'removezero'
        delete(zh(end))
        zh(end)=[];
        if length(zh)==0
            set(gco,'enable','off')
        end
        zpgui('recompute')
    case 'addpole'
        if length(ph)>0
            ph(end+1) =  line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''poleclick'')',...
                'markersize',8,'linewidth',1,'marker','x','linestyle','none');
        else
            ph =  line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''poleclick'')',...
                'markersize',8,'linewidth',1,'marker','x','linestyle','none');
            set(findobj('string','Remove Poles'),'enable','on')
        end
        zpgui('recompute')
    case 'removepole'
        delete(ph(end))
        ph(end)=[];
        if length(ph)==0
            set(gco,'enable','off')
        end
        zpgui('recompute')
end