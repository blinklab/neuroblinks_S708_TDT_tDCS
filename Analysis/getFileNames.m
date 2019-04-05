function names=getFileNames(dirstruct)

m=length(dirstruct);
names=cell(m,1);

for i=1:m
    names{i}=dirstruct(i).name;
end

end