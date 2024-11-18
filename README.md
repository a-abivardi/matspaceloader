
# `MatSpaceLoader`

A collection of functions to load a very specific matlab file.

### 1. Initialize Loader

`space = MatSpaceLoader('path/to/workspace.mat')`

### 2. List Variables in Workspace

`space.get_varnames()`


### 3. Search for a Variable by Name or Code

`pizza_column = space.search("Liking for pizza")`\
`found varsVARS: 20698-0.0, varsHeader: Liking for pizza (0.0)`

or

`pizza_column = space.search(20698)`\
`found varsVARS: 20698-0.0, varsHeader: Liking for pizza (0.0)`


### 4. Quick Load from Workspace

`vars = space.load('vars')  # load entire array`\
`vars_subset = space.load('vars', cols=[1, 3], rows=[10, 20])`


### 5. Load Data into Pandas DataFrame (preferred one)


`df_pizza_coffee = thespace.load_pandas('vars', headers='varsHeader',cols=[space.search('Liking for pizza'), space.search('Liking for coffee')], sub_list='path/subjects.txt', store=True)`

- store will store df in space.df

### 6. Load DataFrame from multiple variables at once

This works for variables, which match the length of either 'subject_IDs' (=session) or 'subject_IDs_unique'. In cases of mixed session and unique data, the unique data will be duplicated for subjects with 2 available sessions.

`df_multi_variables = space.load_pandas(['birth_date', 'vars', 'IDPs'], headers=['_birth_date', 'varsHeader', 'IDP_names'], cols=[[0], space.search('Liking for pizza'), [3,4,6,20]], return_index=True)`

- return_index will add the original row indexes for each variable as column (matlab and python convention separately)

### 7. Clean and reassemble dataframe functions

- If you need the dataframe without index and subject columns (e.g., for further processing), do this:

`df.clean() # removes subject and index columns`

- Need them back, do this:

`df.reassemble() # reinserts the columns`


### Dependencies:
- h5py, pandas and numpy

### for further documentation see code / docstrings

