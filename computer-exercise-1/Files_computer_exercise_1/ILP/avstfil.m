avstandavb = get(gcbo,'String'); 
normavbrott = 0;
avb = avstandavb;
whitecomnotmy
avstandavbrott = str2num(avstandavb);
if avstandavbrott <= 0
  men1= uicontrol('Parent',men, ...
	'Backgroundcolor',[1 1 1], ...
	'Foreground',[1 0 0], ...
	'Units','points', ...
	'Position',[10 40 280 20], ...	
	'Style','text', ...
 	'Fontsize',12, ...
	'string','The distance must be grater than 0', ...
	'Tag','StaticText1');
  return
else
men1= uicontrol('Parent',men, ...
	'Backgroundcolor',[1 1 1], ...
	'Foreground',[0 0 0], ...
	'Units','points', ...
	'Position',[10 122 140 20], ...	
	'Style','text', ...
 	'Fontsize',12, ...
	'string',avb, ...
	'Tag','StaticText1');
end
if yyy == 1
  whitelosilp
  setoffilp;
  rita
  setonilp;
end
yyy = 0;

