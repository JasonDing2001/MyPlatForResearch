function y = NormalizeObj(y)
    ymin = min(y,[],1);
    ymax = max(y,[],1);
    y = (y - ymin) ./ (ymax - ymin);
end