Installation
========================================

MHKit Matlab consists of Matlab code which runs MHKit Python, therefore installation of both packages is necessary. 

MHKit Python
-------------

MHKiT requires Python (tested on 3.6, and 3.7) along with several Python 
package dependencies.  Information on installing and using Python can be found at 
https://www.python.org/.  Python distributions, such as Anaconda,
are recommended to manage the Python interface.  
Anaconda Python distributions include the Python packages needed to run MHKiT.


MHKiT can be installed using pip, git, or a downloaded zip file.  

**pip:** To install MHKiT using pip::

	pip install mhkit
	
**git**: To install MHKiT using git::

	git clone https://code.primre.org/mhkit/mhkit-python 
	cd mhkit-python
	python setup.py install

**zip file**: To install MHKiT using a downloaded zip file, go to https://code.primre.org/mhkit/mhkit-python, 
select the "Download" button and then select "Download zip".
This downloads a zip file called mhkit-python-master.zip.
The software can then be installed by unzipping the file and running setup.py::

	unzip mhkit-python-master.zip
	cd mhkit-python-master
	python setup.py install	
	
Required Python package dependencies include:

* **Pandas**: used for data storage and analysis, http://pandas.pydata.org
* **Numpy**: used for data storage and analysis, http://www.numpy.org
* **Scipy**: used for numerical methods, statistics, and signal processing, https://docs.scipy.org
* **Matplotlib**: used to produce figures, http://matplotlib.org
* **Requests**: used to get data from websites, https://requests.readthedocs.io/
* **Pecos**: used for quality control analysis, https://pecos.readthedocs.io/

MHKiT-Matlab
--------------

mhkit_python_utils package
--------------------------
mhkit_python_utils is a helper package for running MHKiT- Matlab. From https://code.primre.org/mhkit/matlab, download setup.py and mhkit_python_utils. 
Run setup.py on your machine by running::

	python3 setup.py install

Matlab Requirements
--------------------
Matlab 2013b or later is required to run Python from Matlab. MHKit Matlab has been tested on versions 2019b and 2018b.

MHKit Matlab Installation
--------------------------
Download mhkit.mltbx from https://code.primre.org/mhkit/matlab. 
In Matlab, navigate to the folder where you downloaded mhkit.mltbx to, double click on it, and the toolbox will install automatically. 

Setup Matlab Environment
--------------------------

Open Matlab and in the terminal type::

    pyversion

If the resulting Python version is not 3.6, or 3.7 open a Window or Mac terminal window and type::

    python3 -c "import sys; print(sys.executable)"

If the resulting path_to_exe indicates Python 3.6, or 3.7, copy the path and in the Matlab terminal run::

    pyversion('<path_to_exe>')

Note: Mac computers come with Python 2.7 pre-installed. MHKit does not work with Python 2.7.  A second version of 
Python (3.6, or 3.7) will need to be installed on your machine. DO NOT DELETE Python 2.7. Use the above steps to assure 
Matlab is running the proper version of Python. 

Test Install Example
---------------------
To test that your install of MHKit worked correctly, run the following in your Matlab terminal:

	[x,y]=circular(30)

The results should be: 

	x = 30

	y = 1.1310e+04





