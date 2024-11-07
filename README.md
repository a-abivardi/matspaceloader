
# `MatSpaceLoader`

### 1. Initialize Loader


`thespace = MatSpaceLoader('path/to/workspace.mat')`

### 2. List Variables in Workspace

`thespace.get_varnames()`


### 3. Search for a Variable by Name or Code

`pizza_column = thespace.search("Liking for pizza")`\
`found varsVARS: 20698-0.0, varsHeader: Liking for pizza (0.0)`

or

`pizza_column = thespace.search(20698)`\
`found varsVARS: 20698-0.0, varsHeader: Liking for pizza (0.0)`


### 4. Quick Load a Variable from Workspace

`vars = thespace.load('vars')  # load entire array`
`vars_subset = thespace.load('vars', cols=[1, 3], rows=[10, 20])`


### 5. Load Data into Pandas DataFrame (use this one)


`df_pizza_coffee = thespace.load_pandas('vars', headers='varsHeader',cols=[thespace.search('Liking for pizza'), thespace.search('Liking for coffee')], sub_list='path/subjects.txt', return_index=True, store=True)`


### 6. Load  DataFrame from multiple variables at once 

`df_multi_variables = thespace.load_pandas(['birth_date', 'vars', 'IDPs'], headers=['_birth_date', 'varsHeader', 'IDP_names'], cols=[[0], thespace.search('Liking for pizza'), [3,4,6,20]])`