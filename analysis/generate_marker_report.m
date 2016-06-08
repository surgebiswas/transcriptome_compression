function  generate_marker_report( model, tids, Y )

s = model.S;
p = diff(-[1, model.punexp]);
fid = fopen('marker_report.txt', 'w+');

for i =1  : length(s)
    pp = prctile(Y(:, s(i)), [5 10 25 50 75 90 95]);
    fprintf( fid, '%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n',...
        tids{s(i)}, (100*p(i)), pp(1), pp(2), pp(3), pp(4), pp(5), pp(6), pp(7));
end
fclose(fid);


end

