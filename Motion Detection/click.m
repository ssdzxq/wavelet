function click( src, event )
%CLICK Summary of this function goes here
%   Detailed explanation goes here
global hfig
xy = get(hfig, 'CurrentPoint');
hpos = get(hfig,'Position');
apos = get(gca,'Position');
x = (xy(1)-apos(1)*hpos(3))/(apos(3)*hpos(3));
y = (xy(2)-apos(2)*hpos(4))/(apos(4)*hpos(4));
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
x = x*(xlim(2)-xlim(1)) + xlim(1);
y = y*(ylim(2)-ylim(1)) + ylim(1);
sprintf('x = %f:\t',x);
sprintf('y = %f:\n',y);
end

