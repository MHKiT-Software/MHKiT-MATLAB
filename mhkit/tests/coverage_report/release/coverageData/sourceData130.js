var sourceData130 = {"FileName":"/Users/asimms/Desktop/Programming/mhkit_matlab_simms_dev/MHKiT-MATLAB-2/mhkit/utils/reloadPy.m","RawFileContents":["function reloadPy()","%RELOADPY Reloads Python for use after Python code changes","    warning('off','MATLAB:ClassInstanceExists')","    clear classes","    mod = py.importlib.import_module('mhkit');","    py.importlib.reload(mod);","end","",""],"CoverageDisplayDataPerLine":{"Function":{"LineNumber":1,"Hits":0,"StartColumnNumbers":0,"EndColumnNumbers":19,"ContinuedLine":false},"Statement":[{"LineNumber":3,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":47,"ContinuedLine":false},{"LineNumber":4,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":18,"ContinuedLine":false},{"LineNumber":5,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":46,"ContinuedLine":false},{"LineNumber":6,"Hits":0,"StartColumnNumbers":4,"EndColumnNumbers":29,"ContinuedLine":false}]}}