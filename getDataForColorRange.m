function [x,y,idx] = getDataForColorRange(axesHandle,c)

childHandle = axesHandle.Children;
idxR = childHandle.CData > c(1);
idxL = childHandle.CData < c(2);
idx = logical(idxR.*idxL);
x = childHandle.XData(idx);
y = childHandle.YData(idx);
