function SuperTitle(titleText)
% puts one title at the top of a subplot

set(gcf,'NextPlot','add');
axes;
h = title(titleText);
set(gca,'Visible','off');
set(h,'Visible','on');
end