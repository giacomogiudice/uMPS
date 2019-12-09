function defaults = mergestruct(defaults,mods)
labels = fieldnames(mods);
for ind = 1:length(labels)
    l = labels{ind};
    if isstruct(mods.(l))
        defaults.(l) = mergestruct(defaults.(l),mods.(l));
    else
        defaults.(l) = mods.(l);
    end
end
end
