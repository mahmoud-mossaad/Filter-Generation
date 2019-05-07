function varargout = zpgui(varargin)
% ZPGUI  Zero Pole dragging Graphic User Interface
% Allows you to add, remove, and move the zeros and
% poles of a filter interactively.
% To begin, type:
%    zpgui
 
%   Author: Tom Krauss, 9/1/98
%   Adapted by David Dorran (2011) to force poles and zeros to be conjugate
%   pairs (maybe not the most flexible feature but it saves me time). Updated in 2012 to plot the z-surface.
 
global zh ph Lresp Nfft b a z h Y
global ax1 ax2 ax3 z_surface_CameraPos z_surface_CameraUpVec surface_display_opts
 
if nargin == 0
    action = 'init';
    z = [0.29-0.41j 0.29+0.41j]' ;
    p = [-0.6+0.73j -0.6-0.73j]' ;
elseif nargin &amp;gt;= 2
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
        surface_display_opts = 0;;
        set(0,'defaultaxesfontsize',10)
 
        subplot(2,2,1)
 
        [zh,ph,cruff]=zplaneplot(z,p);
        %set(cruff,'hittest','off')
        ax1 = gca;
        ylabel('Imaginary Part')
        xlabel('Real Part')
 
        uicontrol('style','pushbutton',...
            'string','Add Zeros',...
            'fontsize',8,...
            'units','normalized',...
            'position',[.15 .46 .18 .04],...
            'callback','zpgui(''addzero'')');
        uicontrol('style','pushbutton',...
            'string','Remove Zeros',...
            'fontsize',8,...
            'units','normalized',...
            'position',[.55 .46 .18 .04],...
            'callback','zpgui(''removezero'')');
        uicontrol('style','pushbutton',...
            'units','normalized',...
            'fontsize',8,...
            'position',[.35 .46 .18 .04],'string','Add Poles',...
            'callback','zpgui(''addpole'')');
        uicontrol('style','pushbutton',...
            'units','normalized',...
            'fontsize',8,...
            'position',[.75 .46 .18 .04],'string','Remove Poles',...
            'callback','zpgui(''removepole'')');
        uicontrol('style','pushbutton',...
            'units','normalized',...
            'fontsize',8,...
            'position',[.5 .96 .18 .04],'string','Toggle surface display',...
            'callback','zpgui(''toggle_surface_display'')');
 
        subplot(2,1,2)
 
        [b,a]=zp2tf(z,p,1);
 
        Nfft = 512;
 
        Y = fft(b,Nfft)./fft(a,Nfft);
 
        Lresp = plot((0:Nfft-1)/Nfft*2-1, 20*log10(fftshift(abs(Y))),'linewidth',2);
        ax2 = gca;
        set(ax2,'xlim',[0 1])
 
        %set(ax2,'xtick', [0:1/16:0.5])
        %set(ax2,'xticklabel', '0 \pi/8  \pi/8  \pi/8  \pi/8  \pi/8  \pi/8  \pi/8 ')
        %get(ax2,'xticklabel')
        grid on
        xlabel('Normalised Frequency (x\pi rads/sample)')
        ylabel('Magnitude (db) ')
 
        set(Lresp,'erasemode','xor')
 
        set(zh,'buttondownfcn','zpgui(''zeroclick'')',...
            'markersize',8,'linewidth',1)
        set(ph,'buttondownfcn','zpgui(''poleclick'')',...
            'markersize',8,'linewidth',1)
        subplot(2,2,2)
        plot_z_surface(p,z);
        ax3 = gca;
 
        set(gcf, 'ToolBar', 'figure');
 
        if(nargin == 3) % create some figures ans save them to a file
            figure
            set(gcf,'position',  [ 232   387   417   279]);
            plot_z_surface(p,z);
            saveas(gcf, strcat('z_surface_', jpg_filename, '.bmp'), 'bmp');
            figure
            set(gcf,'position',  [ 232   387   417   279]);
            plot((0:Nfft-1)/Nfft - .5, 20*log10(fftshift(abs(Y))),'linewidth',2);
            set(gca,'xlim',[0 .5])
            grid on
            xlabel('Normalised Frequency')
            ylabel('Magnitude (db) ')
            saveas(gcf,  strcat('z_response_', jpg_filename, '.bmp'), 'bmp');
            figure
            set(gcf,'position',  [ 232   387   417   279]);
            zplaneplot(z,p)
            ylabel('Imaginary Part')
            xlabel('Real Part')
            saveas(gcf, strcat('z_pz_', jpg_filename, '.bmp'),  'bmp');
 
            figure
            set(gcf,'position',  [232   605   392    61])
            b_coeff_vals = poly(z)
            a_coeff_vals = poly(p)
 
            text(0.05,0.1 ,strcat('a = [',num2str(a_coeff_vals) ,']'))
            text(0.05,0.8 ,strcat('b = [',num2str(b_coeff_vals) ,']'))
            set(gca, 'visible','off')
            pause(0.1)
            saveas(gcf, strcat('z_ba_', jpg_filename, '.bmp'),'bmp');
            close all
        end
    case 'toggle_surface_display'
        surface_display_opts = surface_display_opts + 1;
 
        if(surface_display_opts &amp;gt;= 2)
            surface_display_opts = 0;
        end
        zpgui('recompute')
    case 'addzero'
        if length(zh)&amp;gt;0
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
        if length(ph)&amp;gt;0
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
                    title(['selected position: ' num2str(pt) 'j'])
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
        set(zh(ind),'erasemode','normal')
        set(Lresp,'erasemode','normal')
        set(ax2,'ylimmode','auto')
        ylim = get(ax2,'ylim');
        Y = get(Lresp,'ydata');
        if ylim(1)&amp;gt;min(Y)-3 | ylim(2)&amp;lt;max(Y)+3
            set(ax2,'ylim',[min(Y)-3 max(Y)+3])
        end
 
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
        set(gcf,'windowbuttonmotionfcn','')
        set(gcf,'windowbuttonupfcn','')
        set(ph(ind),'erasemode','normal')
        set(Lresp,'erasemode','normal')
        set(ax2,'ylimmode','auto')
        Y = get(Lresp,'ydata');
        ylim = get(ax2,'ylim');
        if ylim(1)&amp;gt;min(Y)-3 | ylim(2)&amp;lt;max(Y)+3
            set(ax2,'ylim',[min(Y)-3 max(Y)+3])
        end
 
    case 'recompute'
 
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
        subplot(2,2,2)
 
        z_surface_CameraPos = get(ax3, 'CameraPosition');
        z_surface_CameraUpVec = get(ax3, 'CameraUpVector');
        plot_z_surface(p,z, z_surface_CameraPos, z_surface_CameraUpVec, surface_display_opts);
 
end