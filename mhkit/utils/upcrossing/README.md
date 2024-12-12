# Upcrossing Analysis Functions

This module contains a collection of functions that facilitate **upcrossing analysis** for time-series data. It provides tools to find zero upcrossings, peaks, troughs, heights, periods, and allows for custom user-defined calculations between zero crossings.

## Key Functions

### `upcrossing(t, data)`
Finds the zero upcrossing points in the given time-series data.

**Parameters:**
- `t` (array): Time array.
- `data` (array): Signal time-series data.

**Returns:**
- `inds` (array): Indices of zero upcrossing points.

---

### `peaks(t, data, inds)`
Finds the peaks between zero upcrossings.

**Parameters:**
- `t` (array): Time array.
- `data` (array): Signal time-series data.
- `inds` (array, optional): Indices of the upcrossing points.

**Returns:**
- `peaks` (array): Peak values between the zero upcrossings.

---

### `troughs(t, data, inds)`
Finds the troughs between zero upcrossings.

**Parameters:**
- `t` (array): Time array.
- `data` (array): Signal time-series data.
- `inds` (array, optional): Indices of the upcrossing points.

**Returns:**
- `troughs` (array): Trough values between the zero upcrossings.

---

### `heights(t, data, inds)`
Calculates the height between zero upcrossings. The height is defined as the difference between the maximum and minimum values between each pair of zero crossings.

**Parameters:**
- `t` (array): Time array.
- `data` (array): Signal time-series data.
- `inds` (array, optional): Indices of the upcrossing points.

**Returns:**
- `heights` (array): Height values between the zero upcrossings.

---

### `periods(t, data, inds)`
Calculates the period between zero upcrossings. The period is the difference in time between each pair of consecutive upcrossings.

**Parameters:**
- `t` (array): Time array.
- `data` (array): Signal time-series data.
- `inds` (array, optional): Indices of the upcrossing points.

**Returns:**
- `periods` (array): Period values between the zero upcrossings.

---

### `custom(t, data, func, inds)`
Applies a custom user-defined function between zero upcrossing points.

**Parameters:**
- `t` (array): Time array.
- `data` (array): Signal time-series data.
- `func` (function handle): A custom function that will be applied between the zero crossing periods.
- `inds` (array, optional): Indices of the upcrossing points.

**Returns:**
- `custom_vals` (array): Custom values calculated between the zero crossings.

---

## Author(s)
- **mshabara**

## Date
- 12/12/2024
