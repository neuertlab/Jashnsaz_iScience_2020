function subplotHJ = subplotHJ(nrows, ncolmns, kk, dy, dx)
n=nrows; m=ncolmns; k=kk(1); 
W=(1-(m+1)*dx)/m; % width of each subplot                                                          
H=(1-(n+1)*dy)/n; % Height of each subplot
[j,i]=ind2sub([m,n],k);
% get position for left bottom corner of each subplot wrt the left bottom
% corner of the figure which is at (0, 0) and the whole figure has width=1
% and height=1. 
leftspace = (j-1)*W + j*dx; 
bottomspace = 1-i*(H+dy);
subplotHJ=subplot('Position',[leftspace bottomspace W*length(kk)+(length(kk)-1)*dx H]);





