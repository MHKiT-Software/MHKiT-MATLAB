Return to MHKiT `documentation <https://mhkit-code-hub.github.io/MHKiT/index.html>`_

QC module
===================

The QC module includes quality control functions from Pecos, see https://pecos.readthedocs.io for more details.

===========================================  =========================
Functions                                    Description
===========================================  =========================
qc_corrupt                                   Check for corrupt data 
qc_delta                                     Check for stagant data and/or abrupt changes in the data using the difference between max and min values within a rolling window
qc_increment                                 
qc_missing                                   Check for missing data
qc_outlier                                   Check for outliers using normalized data within a rolling window
qc_range                                     Check for data outside the expected range
qc_timestamp                                 Check time series for missing, non-monotonic, and duplicate timestamps
===========================================  ========================= 


.. Note::
    The names of the functions below are of the convention path.path.path.function. Only the function name is used when calling the function in MATLAB. For example, to call on "mhkit.qc.check_timestamp" simply 
    use check_timestamp(). 

.. automodule:: mhkit.qc
    :members:
    :no-undoc-members:
    :show-inheritance:
    
.. autofunction:: mhkit.qc.check_timestamp

.. autofunction:: mhkit.qc.check_missing

.. autofunction:: mhkit.qc.check_corrupt

.. autofunction:: mhkit.qc.check_range

.. autofunction:: mhkit.qc.check_delta

.. autofunction:: mhkit.qc.check_outlier