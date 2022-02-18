function reloadPy()
%RELOADPY Reloads Python for use after Python code changes
    warning('off','MATLAB:ClassInstanceExists')
    clear classes
    mod = py.importlib.import_module('mhkit');
    py.importlib.reload(mod);
end
