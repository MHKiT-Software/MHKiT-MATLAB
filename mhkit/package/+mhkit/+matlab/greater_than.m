function result = greater_than(maximum_version)
    disp(maximum_version);
    result = isMATLABReleaseOlderThan(maximum_version);
end
