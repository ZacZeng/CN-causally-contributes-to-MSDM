function [] = GUI_4()
% Demonstrate how to make a multiline editbox.
% Produces a GUI with an editbox on the left and a listbox on the right.
% The user is invited to enter text into the editbox, either hitting return
% at the end of each line or letting it wrap automatically.  When the
% button is pushed, each line of text from the editbox is placed as an
% entry into the listbox.  Notice the difference between how a wrapped line
% is treated and a returned line is treated in the lisbox.
%
%
% Author:  Matt Fig
% Date:  7/15/2009

S.fh = figure('units','pixels',...
              'position',[450 450 400 200],...
              'menubar','none',...
              'name','Verify Password.',...
              'resize','off',...
              'numbertitle','off',...
              'name','GUI_4');
S.ed = uicontrol('style','edit',...
                 'units','pix',...
                 'position',[10 60 190 120],...
                 'min',0,'max',2,...  % This is the key to multiline edits.
                 'string',{'Enter text here'; 'then push the button.'},...                 
                 'fontweight','bold',...
                 'horizontalalign','center',...
                 'fontsize',11);
S.ls = uicontrol('style','list',...
                 'units','pix',...
                 'position',[210 60 180 120],...
                 'backgroundcolor','w',...
                 'HorizontalAlign','left');
S.pb = uicontrol('style','push',...
                 'units','pix',...
                 'position',[10 10 380 40],...
                 'HorizontalAlign','left',...
                 'string','Transfer',...
                 'fontsize',14,'fontweight','bold',...
                 'callback',{@pb_call,S});
uicontrol(S.ed)   % Give the editbox control. 


function [] = pb_call(varargin)
% Callback for edit.
S = varargin{3};
% Get the string from the edit box.  Note that since the editbox is a
% multiline editbox (max-min>2), the string returned is a cell array.
E = get(S.ed,'string');  
set(S.ls,'string',E)  % Now set the listbox string to the value in E.