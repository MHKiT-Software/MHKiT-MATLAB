"""
Type conversion utilities for MHKiT-Python to MATLAB data transfer.

This module provides utilities for converting Python data structures into MATLAB-compatible
formats while preserving type information and handling index data explicitly. It supports
conversion of basic numeric types, NumPy arrays, Pandas Series, and DataFrames.

The module is designed to work with the MHKiT-MATLAB typecast_from_mhkit_python function
to ensure reliable data transfer between the two environments.
"""

import array
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Union, Any, Optional

import pandas as pd
import xarray as xr
import numpy as np


# Type definitions
NumericType = Union[int, float, np.number]
ArrayType = Union[np.ndarray]
DataFrameType = Union[pd.DataFrame]
SeriesType = Union[pd.Series]


class DataType(Enum):
    """Enumeration of supported data types for MATLAB conversion.

    These types are used to indicate the format of data being sent to MATLAB,
    allowing proper type casting on the MATLAB side.
    """

    INT = "int"
    DOUBLE = "double"
    SCALAR = "scalar"
    ARRAY = "array"
    MIXED = "mixed"
    CAPTURE_LENGTH_MATRIX = "capture_length_matrix"


@dataclass
class ConversionResult:
    """Container for converted data and associated metadata.

    Attributes
    ----------
    type : str
        String identifier of the data type for MATLAB parsing
    data : Any
        The converted data in a MATLAB-compatible format
    index : Optional[Dict[str, np.ndarray]]
        Index information containing 'name' and 'data' keys if present,
        None otherwise
    """

    type: str
    data: Any
    index: Optional[Dict[str, np.ndarray]] = None


def get_numeric_type(value: NumericType) -> str:
    """Determine the appropriate MATLAB numeric type.

    Parameters
    ----------
    value : NumericType
        Python numeric value to check

    Returns
    -------
    str
        Either 'int' or 'double' depending on the input type
    """
    if isinstance(value, (int, np.integer)):
        return DataType.INT.value
    return DataType.DOUBLE.value


def convert_scalar(data: NumericType) -> ConversionResult:
    """Convert scalar values to MATLAB-compatible format.

    Parameters
    ----------
    data : NumericType
        Single numeric value to convert

    Returns
    -------
    ConversionResult
        Converted scalar with appropriate type information

    Notes
    -----
    - Integers are preserved as integers in MATLAB
    - Floating point numbers are converted to MATLAB doubles
    """
    return ConversionResult(
        type=get_numeric_type(data),
        data=float(data) if isinstance(data, (float, np.floating)) else int(data),
        index=None,
    )


def convert_numpy_array(data: Union[np.ndarray, array.array]) -> ConversionResult:
    """Convert numpy arrays to MATLAB-compatible format.

    Parameters
    ----------
    data : np.ndarray
        NumPy array to convert

    Returns
    -------
    ConversionResult
        Array data with type information

    Notes
    -----
    - Arrays are converted to column-major order for MATLAB compatibility
    - Complex arrays are preserved
    """

    # Cast array.array to numpy array
    data = np.array(data)

    return ConversionResult(type=str(type(data)), data=data, index=None)


def convert_index(index: pd.Index) -> np.ndarray:
    """Convert pandas Index to MATLAB-compatible format.

    Parameters
    ----------
    index : pd.Index
        Pandas index to convert

    Returns
    -------
    np.ndarray
        Converted index values

    Notes
    -----
    - DatetimeIndex is converted to Unix timestamps (seconds since epoch)
    - Other index types are converted to their array values
    """
    if isinstance(index, pd.DatetimeIndex):
        return (index.astype("int64") // 10**9).values
    return index.values


def convert_series(data: pd.Series) -> ConversionResult:
    """Convert pandas Series to MATLAB-compatible format.

    Parameters
    ----------
    data : pd.Series
        Pandas Series to convert

    Returns
    -------
    ConversionResult
        Converted data with index information

    Raises
    ------
    TypeError
        If the Series has a MultiIndex

    Notes
    -----
    - Single value Series are converted to scalars
    - Index information is preserved in the index field
    - The Series name is used as the index name if available
    """
    if isinstance(data.index, pd.MultiIndex):
        raise TypeError(
            "MHKiT-MATLAB-Utils: MultiIndex Series are not supported. Please flatten or restructure your data."
        )

    values = np.array(data.values)
    series_name = data.name
    this_type = get_numeric_type(values[0])

    # Handle single value series
    if values.size == 1:
        values = (
            float(values[0])
            if isinstance(values[0], (float, np.floating))
            else int(values[0])
        )
        this_type = get_numeric_type(values)

    # Always include index information
    index_name = data.index.name if data.index.name is not None else "index"
    index_data = convert_index(data.index)

    return ConversionResult(
        type=f"array_{this_type}",
        data={series_name: values},
        index={
            "name": index_name,
            "data": index_data,
        },
    )


def convert_capture_length_matrix(data: pd.DataFrame) -> ConversionResult:
    """Convert capture length matrix format to MATLAB-compatible format.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame in capture length matrix format with x_centers index
        and y_centers columns

    Returns
    -------
    ConversionResult
        Structured data containing x_centers, y_centers, and the matrix values

    Raises
    ------
    ValueError
        If conversion to float arrays fails or data is corrupted

    Notes
    -----
    - Input must be a DataFrame with 'x_centers' index and 'y_centers' columns
    - All values are converted to floating point
    - NaN values are not allowed in center coordinates
    """
    try:
        # Convert x_centers to float array
        x_centers_raw = data.index.get_level_values("x_centers").unique().values
        x_centers = np.array(x_centers_raw, dtype=np.float64)

        # Convert y_centers to float array
        y_centers_raw = data.columns.values
        y_centers = np.array(y_centers_raw, dtype=np.float64)

        # Validation checks
        if len(x_centers) != len(x_centers_raw) or len(y_centers) != len(y_centers_raw):
            raise ValueError(
                "MHKiT-MATLAB-Utils: Capture length matrix label data was lost during float conversion"
            )

        if np.any(np.isnan(x_centers)) or np.any(np.isnan(y_centers)):
            raise ValueError(
                "MHKiT-MATLAB-Utils: Capture length matrix label values could not be converted to float"
            )

        matrix = np.array(data.values)

        return ConversionResult(
            type=DataType.CAPTURE_LENGTH_MATRIX.value,
            data={
                "x_centers": x_centers,
                "y_centers": y_centers,
                "capture_length_matrix": matrix,
            },
            index=None,
        )

    except (ValueError, TypeError) as e:
        raise ValueError(
            f"Failed to convert x_centers or y_centers to floating point array: {str(e)}"
        ) from e


def convert_dataframe(data: pd.DataFrame) -> ConversionResult:
    """Convert pandas DataFrame to MATLAB-compatible format.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame to convert

    Returns
    -------
    ConversionResult
        Converted data with preserved column structure and index information

    Raises
    ------
    TypeError
        If the DataFrame has a MultiIndex

    Notes
    -----
    - Checks for capture length matrix format first
    - Single value columns are converted to scalars
    - Column names are preserved as struct field names in MATLAB
    - Index information is preserved in the index field
    """
    if data.index.name == "x_centers" and data.columns.name == "y_centers":
        return convert_capture_length_matrix(data)

    if isinstance(data.index, pd.MultiIndex):
        raise TypeError(
            "MHKiT-MATLAB-Utils: MultiIndex DataFrames are not supported. Please flatten or restructure your data."
        )

    result = {}
    types = []

    # Always include index information
    index_name = data.index.name if data.index.name is not None else "index"
    index_data = convert_index(data.index)

    for column in data.columns:
        this_value = np.array(data[column].values)

        if this_value.size == 1:
            this_value = (
                float(this_value[0])
                if isinstance(this_value[0], (float, np.floating))
                else int(this_value[0])
            )
            types.append(f"scalar_{get_numeric_type(this_value)}")
        else:
            types.append(f"array_{get_numeric_type(this_value[0])}")

        result[column] = this_value

    output_type = types[0] if len(set(types)) == 1 else DataType.MIXED.value

    return ConversionResult(
        type=output_type,
        data=result,
        index={
            "name": index_name,
            "data": index_data,
        },
    )


def convert_to_matlab_compatible(
    data: Union[NumericType, ArrayType, DataFrameType, SeriesType],
) -> Dict[str, Any]:
    """Convert Python data structures to MATLAB-compatible format.

    This is the main entry point for converting Python data structures to a format
    that can be easily loaded into MATLAB. It handles type conversion, preserves
    index information, and provides explicit type hints for MATLAB parsing.

    Parameters
    ----------
    data : Union[NumericType, ArrayType, DataFrameType, SeriesType]
        The Python data to convert. Supported types include:
        - Basic numeric types (int, float)
        - NumPy arrays
        - Pandas Series
        - Pandas DataFrames
        - Capture length matrices (special DataFrame format)

    Returns
    -------
    Dict[str, Any]
        A dictionary containing:
        - type: string indicating the data type for MATLAB parsing
        - data: the converted data in MATLAB-compatible format
        - index: index information if present (contains 'name' and 'data' keys),
                None otherwise

    Raises
    ------
    TypeError
        - If the input type is not supported (e.g., xarray)
        - If a DataFrame/Series has a MultiIndex

    ValueError
        - If capture length matrix conversion fails
        - If numeric conversion fails

    Notes
    -----
    The output structure is designed to work with the MATLAB function
    typecast_from_mhkit_python for seamless data transfer between environments.
    """
    if isinstance(data, (xr.Dataset, xr.DataArray)):
        raise TypeError(
            "MHKiT-MATLAB-Utils: xarray conversion is not supported. Please convert data to pandas or numpy format first."
        )

    if np.isscalar(data):
        result = convert_scalar(data)
    elif isinstance(data, (np.ndarray, array.array)):
        result = convert_numpy_array(data)
    elif isinstance(data, pd.Series):
        result = convert_series(data)
    elif isinstance(data, pd.DataFrame):
        result = convert_dataframe(data)
    else:
        raise TypeError(f"MHKiT-MATLAB-Utils: Unsupported data type: {type(data)}")

    return {"type": result.type, "data": result.data, "index": result.index}
