function dat = readAtt(filename)

	fid = fopen(filename);
	dims = fread(fid,[1 3],'uint32');
	dat = fread(fid,[dims(1) dims(2)*dims(3)],'uint8');
	dat = reshape(dat,dims);
	fclose(fid);
	%Add code to investigate data here...
        return 
end
