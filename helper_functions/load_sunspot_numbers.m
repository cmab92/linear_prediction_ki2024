function [t, x] = load_sunspot_numbers()
    url = 'http://www.sidc.be/silso/INFO/snmtotcsv.php';
    filename = 'sunspot_data.csv';
    if not(exist(filename, 'file'))
        disp("download sunspot numbers")
        websave(filename,url);
    end
    data = readmatrix(filename);  
    t = data(:,3);
    x = data(:,4);
end

