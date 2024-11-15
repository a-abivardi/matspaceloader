import h5py
import numpy as np
import pandas as pd
from itertools import chain
import os

class MatSpaceLoader():
    def __init__(self, workspace):
        
        """ 
        Initialize with the workspace file path.

        Parameters
        ----------
        workspace (str): path to MATLAB workspace
        """
        if not os.path.exists(workspace):
            raise  FileNotFoundError('workspace not found, try again')
        self.workspace = workspace
        self.df = None

    def load(self, variable_name, cols=None, rows=None, is_matlab_index=False, string_index=slice(None)):

        """ 
        Simple load of a variable from workspace mat file as numpy array.
        
        Parameters
        ----------
        variable_name (str): name of variable
        cols, rows (list or slice, optional): columns/rows to return
        is_matlab_index (bool, optional): index data is converted to python convention (done automatically if "index" in variable name)
        string_index (bool, optional): use this for indexing 1d string array

        Returns:
            np.array: loaded data in array
        """

        # check which loader to use
        mode = self._checks(variable_name)

        if mode == 'string':
            print('retrieving data from dataset with h5py references and converting to strings')
            if rows or cols:
                print('rows and cols are ignored when loading string array, use string_index')
            return self._load_strings(variable_name, index=string_index)

        if mode == 'matrix':
            print('loading matlab matrix')
            return self._load_matrix(variable_name, cols=cols, rows=rows, is_matlab_index=is_matlab_index)

    def load_pandas(self, variable_names, headers=None, cols=None, rows=None, header_index=None, header_index_pos=None,
                    subid_n = (67470, 62823), store=False, return_index=False, sub_list = None):
        
        """ 
        Load of a variable from workspace mat file as pandas dataframe with subject_ids.
        
        Parameters
        ----------
        variable_names (str or list[str]): name/-s of variable/-s to load. 
        headers (str or list[str], optional): individual header or list of headers to extract from workspace. 
                                                custom headers can be passed as one list per matlab variable.
        header_index (str, list or slice, optional): index for header. can be index or string for extraction from workspace.
        header_index_pos (int): if passing multiple headers, position of index has to be defined. currently implemented for single header only.
        cols, rows (list or slice, optional): columns/rows to return. 
        subid_n (tuple(int,int)): number of subjects with 2 scans or 1 scan, respectively (may need to be adjusted to workspace).
        store (bool, default: False): if true dataframe will be stored in MatSpaceLoader object as object.df
        return_index (bool, default: False): if true original python and matlab index of extracted variable(/-s) will be added as column to dataframe.
        sub_list (path or list or array): path to subject ID list text file or list or array of subject_IDs. full subject_IDs will be used where available
                                                and if provided, otherwise reduction to or use of subject_ID_orig. 


        Returns:
            pd.DataFrame with custom functions: the requested dataframe with subject IDs, when possible. multi-variable frames only available for subject ID compatible 
            variables. when variables with mixed subject indexes are used (i.e., imaging session index 'subject_IDs' + original 'subject_IDs_orig/_unique'), 
            the original index data will be duplicated in cases of double entries for a subject (i.e., both imaging sessions available). in multiframe
            mode the length of the variable_names list must match cols/rows lists and headers list, when provided.
            (Workspace indexes (using python convention) are added when returning a subset of subject rows.)

            Custom functions:
                DataFrame.clean(): removes all indexes and subject_ids from DataFrame but keeps them in a safe location
                DataFrame.reassemble(): reinserts all indexes and subject_ids  

        """

        # lots of input checking
        variable_names, headers, cols, rows, sub_list = self._process_load_pandas_inputs(variable_names, headers=headers, cols=cols, rows=rows, 
                                                                                         header_index=header_index, header_index_pos=header_index_pos, 
                                                                                         sub_list = sub_list)

        # create dataframe
        df = self._make_multiframe(variable_names, headers=headers, cols=cols, rows=rows, return_origshape=True,
                        subid_n = subid_n, return_index=return_index, sub_list = sub_list)

        if store:
            self.df = df
        return df
 
    def get_varnames(self):

        """ 
        This returns all variable names contained in the workspace as a list
        """

        with h5py.File(self.workspace, 'r') as f:
            variable_names = [name for name in f.keys() if name not in ['I', '#refs#', 'ans']]
        return variable_names
    
    def search(self, search_infer = None, search_code=None, search_string=None):        
        
        """ 
        Search for UKBB column indexes by code or name in vars and Svars 
        """

        if search_infer and (search_code is None) and (search_string is None):
            if isinstance(search_infer, int):
                search_code = str(search_infer) + '-'

            elif isinstance(search_infer, str):
                try:
                    if '-' not in search_infer:
                        search_code = str(int(search_infer)) + '-' 
                    else:
                        search_code = str(int(search_infer))
                except:

                    search_string = search_infer 

        if (search_code is None) == (search_string is None):
            raise ValueError('Please provide either search code or string')
        
        self._load_headers()

        if search_code:
            if isinstance(search_code, int):
                search_code = str(search_code) + '-'
            elif isinstance(search_code, str):
                if '-' not in search_code:
                    search_code = search_code + '-'
            else:
                raise ValueError('search code needs to be string or int')
            
            found_codes = [self._find_code(VARS, search_code) for VARS in [self.varsVARS, self.SvarsVARS]]
        
        if search_string:
            found_codes = [self._find_code(header, search_string) for header in [self.varsHeader, self.SvarsHeader]]

        if found_codes[0]:
            [print(f'found varsVARS: {self.varsVARS[idx]}, varsHeader: {self.varsHeader[idx]}') for idx in found_codes[0]]
        
        if found_codes[1]:
            [print(f'found SvarsVARS: {self.SvarsVARS[idx]}, SvarsHeader: {self.SvarsHeader[idx]}') for idx in found_codes[1]]

        if found_codes[0] and found_codes[1]:
            print('returning dictionary with vars and Svars indexes')
            return dict(vars_index = found_codes[0], Svars_index = found_codes[1])
        elif found_codes[0]:
            print('returning vars indexes')
            return found_codes[0]
        elif found_codes[1]:
            print('returning Svars indexes')
            return found_codes[1]    
        else:
            print('found nothing, try again')
            return None
        
    @staticmethod
    def code_loader(codes_path):
        """
        A function to load UKBB codes from the first column of a textfile, as it would appear 
        if you were to copy-paste the data-fields from the UKBB browser (e.g., from  https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100060)
        into it.

        Codes are returned as a list of strings with '-' appended to them (subject to change)

        Example input text:

        Field ID    Description
        20126	Bipolar and major depression status
        20122	Bipolar disorder status
        20127	Neuroticism score
        20124	Probable recurrent major depression (moderate)
        20125	Probable recurrent major depression (severe)
        20123	Single episode of probable major depression  
        """
        codes = []
        with open(codes_path, 'r') as file:
            for line in file:
                line_cols = line.split('\t')
                if len(line_cols) < 2:
                    continue
                if ('(pilot)' in line_cols[1]) or ('Field' in line_cols[0]):
                    continue
                else:
                    code = line_cols[0] + '-'
                    codes.append(code)
        return codes
    

    def idxs_from_codes(self, codes):

        """ 
        This gets column indexes for vars and Svars using list of UKBB codes (e.g., from code_loader).
        
        Returns: List of columns for vars or Svars with a location specifier ('vars', 'Svars'). 
        if columns are found in both variables, a tuplev(vars_indexes, Svars_indexes) and location 'both' are returned.
        """
 
        self._load_headers()

        vars_indexes = []
        Svars_indexes = []
        for code in codes:

            vars_indexes.append(self._find_code(self.varsVARS, code))
            Svars_indexes.append(self._find_code(self.SvarsVARS, code))

        # flatten lists if not already flat
        try:
            vars_indexes = list(chain.from_iterable(vars_indexes))
        except TypeError:
            pass

        try:
            Svars_indexes = list(chain.from_iterable(Svars_indexes))
        except TypeError:
            pass

        if vars_indexes and Svars_indexes:
            return (vars_indexes, Svars_indexes), 'both'
        
        if vars_indexes:
            return vars_indexes, 'vars'
        
        if Svars_indexes:
            return Svars_indexes, 'Svars'
        
        else:
            return None, 'nowhere'
        
    class CustomPandas(pd.DataFrame):
        """
        custom functionlity for pandas dataframe.
        clean(): removes IDs and index columns from dataframe but keeps them safe
        reassemble(): puts everything back together (albeit not exactly in the same order)
        """
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.cleaned = False 
  
        def clean(self):
            the_columns = ['subject_IDs','index_python', 'index_matlab']
            the_columns_loc = self.columns.str.contains('|'.join(the_columns))
            self.attrs['the_columns'] = self.loc[:, the_columns_loc]
            self.drop(columns=self.columns[self.columns.str.contains('|'.join(the_columns))], inplace=True)
            self.cleaned = True 
 
        def reassemble(self):
            if self.cleaned:
                for i, col in enumerate(self.attrs['the_columns'].columns):
                    self.insert(loc=i, column=col, value=self.attrs['the_columns'][col])
                self.cleaned = False
                print('dataframe reassembled. Caution: IDs and indexes have been reinserted at start of dataframe not in original location.')
            else:
                print('nothing to reassemble')
        

    def _process_load_pandas_inputs(self, variable_names, headers=None, cols=None, rows=None, header_index=None, header_index_pos=None, sub_list = None):

        """
        This checks inputs of load pandas and puts into lists where needed. Probably still incomplete.
        """

        # put single variable and header into list
        if not isinstance(variable_names, list):
            variable_names = [variable_names]

            if headers and (not isinstance(headers, list) or (isinstance(headers, list) and len(headers) > 1)):
                headers = [headers]

        
        if isinstance(cols, list):
            if not isinstance(cols[0], list):
                cols = [cols]
                
        elif isinstance(cols, int):
            cols = [[cols]]


        # header lists into arrays
        if headers:
            headers = [np.array(header) if isinstance(header, list) else header for header in headers]

        # apply header index if requested
        if header_index:
            if isinstance(header_index, str):
                header_index = self._load_matrix(header_index)

            if len(headers) > 1:
                if header_index_pos is None:
                    raise ValueError('position of header to be indexed required for more than one header')
            else:
                header_index_pos = 0

            if isinstance(headers[header_index_pos], str):
                headers[header_index_pos] = self._load_strings(headers[header_index_pos])[header_index]
            else:
                headers[header_index_pos] = headers[header_index_pos][header_index]

        #load sub list if provided
        if sub_list is not None:
            if isinstance(sub_list, str):
                try:
                    sub_list = np.loadtxt(sub_list)
                except (OSError, ValueError) as e:
                    raise ValueError(f"Failed to load subject list from '{sub_list}': {e}")
                
        if isinstance(rows, pd.core.series.Series):
            rows = rows.to_numpy().tolist()
    
        # check variable contents (subject to removal or change to warning)
        if any(self._checks(var) == 'string' for var in variable_names):
            raise TypeError('at least one variable contains strings only, please use numerical data')

        return variable_names, headers, cols, rows, sub_list


    def _find_code(self, string_list, code):
        """ 
        helper for search
        """
        return [index for index, string in enumerate(string_list) if string.startswith(code)]

    def _load_headers(self, headerlist = ['varsVARS', 'SvarsVARS', 'varsHeader', 'SvarsHeader']):
        """ 
        load headers for search
        """
        for attr in headerlist:
            if not hasattr(self, attr):
                setattr(self, attr, self._load_strings(attr))

    def _checks(self, variable_name):

        """
        check if variable in workspace and for kind of data (matrix or string)
        """

        with h5py.File(self.workspace, 'r') as f:
            if variable_name not in f:
                raise KeyError(f"Variable '{variable_name}' not found in workspace.")
            dataset = f[variable_name]
            if dataset.dtype.kind == 'O':
                if isinstance(dataset[()].flatten()[0], h5py.Reference):
                    return('string')
                else:
                    print('dataset of object type but no h5py reference in first entry - check output')
                    return('matrix')
            else:
                return('matrix')

    def _load_matrix(self, variable_name, cols=None, rows=None, is_matlab_index = False, return_origshape=False):

        """
        loader for non-string data
        """

        with h5py.File(self.workspace, 'r') as f:
            data = f[variable_name][:].T
            if return_origshape:
                origshape = data.shape
            # get specified cols and rows:
            data = self._slice_matrix(data, rows, cols)
        if is_matlab_index or ('index' in variable_name.lower()):
            print('converting index from matlab to python convention')
            data = (data - 1).astype(int)
        if return_origshape:
            return data.squeeze(), origshape
        else:
            return data.squeeze()

    def _load_strings(self, variable_name, index=slice(None)):

        """
        loader for string data with conversion
        """

        with h5py.File(self.workspace, 'r') as f:
            references =f[variable_name][()].flatten()
            string_list = np.array(["".join(chr(c.item()) for c in f[r][:]) for r in references])[index]
        return string_list
    
    def _make_dataframe(self, variable_name, header, header_name=None, cols=None, rows=None, return_origshape=True,
                        return_index=False):
        # load data and combine with header into dataframe

        """
        helper for making the pandas dataframe
        """

        if header is not None:
            if any(header):
                if not isinstance(header, list) and not isinstance(header, np.ndarray):
                    header = [header]
        if cols is not None and header is not None:
            if len(header) > 1:
                header = header[cols]

        
        data, data_origshape = self._load_matrix(variable_name, cols=cols, rows=rows, return_origshape=return_origshape)

        df = pd.DataFrame(data=data, columns=header)
        
        if header_name:
            df.columns.name = header_name

        # add workspace index if slicing rows
        if rows:
            if isinstance(rows, slice):
                workspace_index = pd.DataFrame(data=[i for i in range(*rows.indices(rows.stop))], columns=[f'{variable_name}_index_python'])
            else:
                workspace_index = pd.DataFrame(data=rows, columns=[f'{variable_name}_index_python'])
            workspace_index[f'{variable_name}_index_matlab'] = workspace_index[f'{variable_name}_index_python'] + 1
            df = pd.concat((workspace_index, df), axis=1)
        elif not rows and return_index:

            workspace_index = pd.DataFrame(data=df.index, columns=[f'{variable_name}_index_python'])
            workspace_index[f'{variable_name}_index_matlab'] = workspace_index[f'{variable_name}_index_python'] + 1
            df = pd.concat((workspace_index, df), axis=1)

        return df, data_origshape
    
    def _make_multiframe(self, variable_names, headers=None, cols=None, rows=None, return_origshape=True,
                         subid_n = (67470, 62823), return_index=False, sub_list =None):

        """
        function for assembly of multiframe
        """

        if cols == None:
            cols = [None for x in variable_names]
        if headers == None:
            headers = [None for x in variable_names]

        
        if not (len(variable_names) == len(headers) == len(cols)):
            raise ValueError('variable names, headers and cols list lengths must match when provided')

        if not (rows is None or sub_list is None):
            raise ValueError('provide either rows or sub_list not both. in multiframe sublists are preferable.')

        # load or set headers
        for cnt, header in enumerate(headers):
            if isinstance(header, str):
                if header.startswith('_'):
                    headers[cnt] = header[1:]
                else:
                    headers[cnt] = self._load_strings(header)
            
            else:
                pass


        for cnt, (var, head, col) in enumerate(zip(variable_names, headers, cols)):

            if cnt == 0:
                df_main, origshape_main = self._make_dataframe(var, head, cols=col, rows=rows, 
                                                               header_name=None, return_origshape=return_origshape,
                                                               return_index=return_index)

                if not origshape_main[0] in subid_n:
                    raise ValueError('multiframe df needs subject ID compatible matrices')

                df_main = self._add_subject_ids(df_main, origshape_main, rows=rows, subid_n=subid_n)

                if sub_list is not None:
                    df_main, python_index = self._df_sublister(df_main, sub_list)
                    df_main[f'{var}_index_python'] = python_index
                    df_main[f'{var}_index_matlab'] = python_index + 1
            else:
                df_next, origshape = self._make_dataframe(var, head, cols=col, rows=rows, header_name=None,
                                                          return_origshape=return_origshape, return_index=return_index)
                
                if not origshape[0] in subid_n:
                    raise ValueError('multiframe df needs subject ID compatible matrices')
                
                if rows:
                    if origshape_main[0] != origshape[0]:
                        raise NotImplementedError("rows implemented for same length df's only. use subject lists")

                df_next = self._add_subject_ids(df_next, origshape, rows=rows)
                
                if sub_list is not None:
                    df_next, python_index = self._df_sublister(df_next, sub_list)
                    df_next[f'{var}_index_python'] = python_index
                    df_next[f'{var}_index_matlab'] = python_index + 1
                
                if 'subject_IDs' in df_main and 'subject_IDs' in df_next:
                    df_main = pd.merge(df_main, df_next, on=['subject_IDs', 'subject_IDs_orig'], how='outer')
                else:
                    df_main = pd.merge(df_main, df_next, on='subject_IDs_orig', how='outer')
                
        # IDs to int after merging
        if 'subject_IDs_orig' in df_main:
            df_main['subject_IDs_orig'] = df_main['subject_IDs_orig'].astype('int')
        if 'subject_IDs' in df_main:
            df_main['subject_IDs'] = df_main['subject_IDs'].astype('int')

        return self.CustomPandas(df_main)

    def _df_sublister(self, df, sub_list):

        """
        Given a sub_list and a dataframe, this returns the dataframe for the subjects in the order of the sub_list.
        For details see _df_sublister_core function below
        """

 
        if sub_list[0] > 20000000:
            
            if 'subject_IDs' in df:

                df_with_subs, orig_index = self._df_sublister_core(df, sub_list, 'subject_IDs')
        
            elif 'subject_IDs' not in df:
                sub_list = np.mod(sub_list, 10000000)
                df_with_subs, orig_index = self._df_sublister_core(df, sub_list, 'subject_IDs_orig')
        
        elif sub_list[0] < 10000000 and 'subject_IDs_orig' in df:
            df_with_subs, orig_index = self._df_sublister_core(df, sub_list, 'subject_IDs_orig')
        
        else:
            raise RuntimeError('Something is wrong with the sublisting')


        return df_with_subs, orig_index
    

    def _df_sublister_core(self, df, sub_list, sub_col):

        """
        Actual sublisting is done here.
        """

        # put sub_list into dataframe
        sub_list = pd.DataFrame(data=sub_list, columns=[sub_col])

        # extract subjects contained in sub_list from input dataframe
        df = df.loc[df[sub_col].isin(sub_list[sub_col])]

        # store original index of these subjects
        df.loc[:, 'orig_index']=df.index

        # merge input dataframe on sub_list dataframe to ensure same ordering as in the input list
        df = pd.merge(sub_list, df, on=sub_col, how='left')

        return df.drop(columns=['orig_index']), df['orig_index']


    def _slice_matrix(self, data, rows, cols):

        """
        helper for slicing
        """

        if rows is None and cols is None:
            return data
        if rows is None:
            return data[:, cols]
        if cols is None:
            return data[rows, :]
        if any(isinstance(x, slice) for x in (rows, cols)):
            return data[rows, cols]
        else:
            return data[np.ix_(rows, cols)]

    def _load_subject_ids(self, data_origshape, rows=None, subid_n = (67470, 62823)):

        """
        helper for loading the IDs
        """

        if data_origshape[0]== max(subid_n):
            return self._load_matrix('subject_IDs',rows=rows).astype(int), self._load_matrix('subject_IDs_orig',rows=rows).astype(int), 
        if data_origshape[0]==min(subid_n):
            return self._load_matrix('subject_IDs_unique', rows=rows).astype(int)
        
    def _add_subject_ids(self, df, data_origshape, rows=None, subid_n = (67470, 62823)):
        """
        helper for adding sub_ids to df
        """
        if data_origshape[0]== max(subid_n):
            subs, subs_orig = self._load_subject_ids(data_origshape, rows, subid_n=subid_n)
            sub_df = pd.DataFrame({'subject_IDs': subs, 'subject_IDs_orig': subs_orig})
        elif data_origshape[0]== min(subid_n):
            subs_unique = self._load_subject_ids(data_origshape, rows, subid_n=subid_n)
            sub_df = pd.DataFrame({'subject_IDs_orig': subs_unique})
        else:
            print("data doesn't fit any subject dimension, frame returned without subject IDs")
            return df
        return pd.concat((sub_df, df), axis=1)
    

    @staticmethod
    def code_loader(codes_path):
        """
        A function to load UKBB codes from the first column of a textfile, as it would appear 
        if you were to copy-paste the data-fields from the UKBB browser (e.g., from  https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100060)
        into it.

        Codes are returned as a list of strings with '-' appended to (subject to change)

        Example txt:

        Field ID    Description
        20126	Bipolar and major depression status
        20122	Bipolar disorder status
        20127	Neuroticism score
        20124	Probable recurrent major depression (moderate)
        20125	Probable recurrent major depression (severe)
        20123	Single episode of probable major depression  
        """
        codes = []
        with open(codes_path, 'r') as file:
            for line in file:
                line_cols = line.split('\t')
                if len(line_cols) < 2:
                    continue
                if ('(pilot)' in line_cols[1]) or ('Field' in line_cols[0]):
                    continue
                else:
                    code = line_cols[0] + '-'
                    codes.append(code)
        return codes
    

if __name__ == '__main__':
    pass