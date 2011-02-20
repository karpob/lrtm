function [goodtimes]=get_diff(dos,mac,col,m)

for mm=1:m
    goodtimes(mm)=100*(dos(mm,col)-mac(mm,col))/dos(mm,col);
end