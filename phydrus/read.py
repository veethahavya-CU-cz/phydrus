"""
The read module contains methods that can be used to read in- and output
files. The methods can be used stand-alone but are also available from the
Model object. All the methods return the data as Pandas DataFrames.

Examples
--------
>>> import phydrus as ps
>>> ps.read.read_obs_node()

or

>>> ml.read_obs_node()

"""

from pandas import read_csv, DataFrame, to_numeric

from .decorators import check_file_path


def read_profile(path="PROFILE.OUT"):
    """
    Method to read the PROFILE.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the profile out file. default is "PROFILE.OUT".

    Returns
    -------
    data: pandas.DataFrame
        Pandas with the profile data

    """
    # Detect if the file uses 'depth' or 'x' as the depth header
    start_str = "depth"
    idx_col = "n"
    try:
        with open(path, "r") as f:
            content = f.read(2000) # Read beginning of file
            if " x " in content or " x\n" in content or " x\r" in content:
                start_str = "x"
                idx_col = None # Use positional index mapping or numeric
            elif "depth" in content.lower():
                start_str = "depth"
                idx_col = "n"
    except Exception:
        pass

    data = _read_file(path=path, start=start_str, idx_col=idx_col)
    return data


def read_run_inf(path="RUN_INF.OUT", usecols=None):
    """
    Method to read the RUN_INF.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the run_inf out file. default is "RUN_INF.OUT".
    usecols: list of str optional
        List with the names of the columns to import. Default is all columns.

    Returns
    -------
    data: pandas.DataFrame
        Pandas with the run_inf data

    """
    data = _read_file(path=path, start="TLevel", idx_col="TLevel",
                      usecols=usecols)
    return data


@check_file_path
def read_i_check(path="I_CHECK.OUT"):
    """
    Method to read the I_CHECK.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the I_Check out file. default is "I_Check.OUT".

    Returns
    -------
    data: dict
        Dictionary with the node as a key and a Pandas DataFrame as a value.

    """
    names = ["theta", "h", "log_h", "C", "K", "log_K", "S", "Kv"]

    end = []
    start = 0

    with open(path) as file:
        # Find the starting line
        for i, line in enumerate(file.readlines()):
            if "theta" in line:
                start = i
            elif "end" in line and start > 0:
                end.append(i)

        data = {}

        for i, e in enumerate(end):
            file.seek(0)  # Go back to start of file

            # Read data into a Pandas DataFrame
            nrows = e - start - 2
            data[i] = read_csv(file, skiprows=start + 1, nrows=nrows,
                               skipinitialspace=True, sep=r"\s+",
                               names=names, dtype=float)
            start = e
        return data


def read_tlevel(path="T_LEVEL.OUT", usecols=None):
    """
    Method to read the T_LEVEL.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the t_level out file. default is "T_LEVEL.OUT".
    usecols: list of str optional
        List with the names of the columns to import. By default
        only the real fluxes are imported and not the cumulative
        fluxes. Options are: "rTop", "rRoot", "vTop", "vRoot", "vBot",
        "sum(rTop)", "sum(rRoot)", "sum(vTop)", "sum(vRoot)", "sum(vBot)",
        "hTop", "hRoot", "hBot", "RunOff", "sum(RunOff)", "Volume",
        "sum(Infil)", "sum(Evap)", "TLevel", "Cum(WTrans)", "SnowLayer".

    Returns
    -------
    data: pandas.DataFrame
        Pandas with the t_level data

    """
    data = _read_file(path=path, start="rTop", idx_col="Time",
                      remove_first_row=True, usecols=usecols)
    return data


def read_alevel(path="A_LEVEL.OUT", usecols=None):
    """
    Method to read the A_LEVEL.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the t_level out file. default is "A_LEVEL.OUT".
    usecols: list of str optional
        List with the names of the columns to import.

    Returns
    -------
    data: pandas.DataFrame
        Pandas with the a_level data

    """
    data = _read_file(path=path, start="Time", idx_col="Time",
                      remove_first_row=True, usecols=usecols)
    return data


def read_solute(path="SOLUTE1.OUT"):
    """
    Method to read the SOLUTE.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the solute out file. default is "SOLUTE1.OUT".

    Returns
    -------
    data: pandas.DataFrame
        Pandas with the a_level data

    """
    data = _read_file(path=path, start="Time", idx_col="Time",
                      remove_first_row=True)
    return data


@check_file_path
def _read_file(path, start, end="end", usecols=None, idx_col=None,
               remove_first_row=False):
    """
    Internal method to read Hydrus output files.

    Parameters
    ----------
    path: str
        String with the filepath.
    start: str
        String that determines the start of the data to be imported.
    end: str, optional
        String that determines the end of the data to be imported.
    usecols: list, optional
        List with the names of the columns to import. Default is all columns.
    idx_col: str, optional
        String with the name used for the index column.
    remove_first_row: bool, optional
        Remove the first row if True. Default is False.

    Returns
    -------
    data: pandas.DataFrame
        Pandas DataFrame with the imported data.

    """
    with open(path) as file:
        lines = file.readlines()
        s = None
        e = len(lines)
        
        # Find the starting line
        for i, line in enumerate(lines):
            if s is None and start in line:
                s = i
            elif s is not None and end in line:
                e = i
                break
                
        if s is None:
            raise ValueError(f"Could not find start string '{start}' in {path}")
        file.seek(0)  # Go back to start of file

        # Read data into a Pandas DataFrame
        data = read_csv(file, skiprows=s, nrows=e - s - 1, usecols=usecols,
                        index_col=idx_col, skipinitialspace=True,
                        sep=r"\s+")

        def convert_to_numeric(df):
            # Convert columns
            df = df.apply(lambda col: to_numeric(col, errors="coerce"))
            # Convert index
            df.index = to_numeric(df.index, errors="coerce")
            return df

        if remove_first_row:
            # Check if first row is a units row (usually contains strings like [L], [T], etc.)
            # We check the first row's index and values for brackets.
            try:
                # Convert the entire first row to a list of strings and check for brackets
                check_strs = [str(data.index[0])] + [str(v) for v in data.iloc[0].values]
                if any("[" in s or "]" in s for s in check_strs):
                    data = data.drop(index=data.index[0])
            except (IndexError, AttributeError):
                pass
            
        data = convert_to_numeric(data)

    return data


@check_file_path
def read_obs_node(path="OBS_NODE.OUT", nodes=None, conc=False, cols=None):
    """
    Method to read the OBS_NODE.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the OBS_NODE out file. default is
        "OBS_NODE.OUT".
    nodes: list of ints, optional
        nodes to imoport
    conc: boolean, optional
    cols: list of strs, optional
        List of strings with the columns to read.

    Returns
    -------
    data: dict
        Dictionary with the node as a key and a Pandas.DataFrame as a value.

    """
    data = {}
    start = None
    end = None

    with open(path) as file:
        for i, line in enumerate(file.readlines()):
            # Match the actual data header: "time  h  theta  Temp..."
            # Case-insensitive match for robust parsing
            l_lower = line.lower()
            if "time" in l_lower and "theta" in l_lower:
                start = i      # this line becomes the pandas column header
                data_start = i + 2  # skip header + units row
            elif "end" in line:
                end = i
                break

    if start is None or end is None:
        return data

    # Read the full table (header + units + data), using the header line as columns
    df1 = read_csv(path, skiprows=start, index_col=0,
                   nrows=end - start - 1,
                   skipinitialspace=True, sep=r"\s+", engine="c")

    # Check if first row is a units row (usually contains strings like [L], [T], etc.)
    try:
        # Convert the entire first row to a list of strings and check for brackets
        check_strs = [str(df1.index[0])] + [str(v) for v in df1.iloc[0].values]
        if any("[" in s or "]" in s for s in check_strs):
            df1 = df1.drop(index=df1.index[0])
    except (IndexError, AttributeError):
        pass

    # Convert data and index to numeric
    df1 = df1.apply(lambda col: to_numeric(col, errors="coerce"))
    df1.index = to_numeric(df1.index, errors="coerce")

    if cols is None:
        cols = ["h", "theta", "Temp"]
    if conc:
        cols.append("Conc")

    for i, node in enumerate(nodes):
        if i > 0:
            usecols = [f"{c}.{i}" for c in cols]
        else:
            usecols = cols

        df = df1.loc[:, usecols]
        df.columns = cols
        data[node] = df

    return data


@check_file_path
def read_nod_inf(path="NOD_INF.OUT", times=None):
    """
    Method to read the NOD_INF.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the NOD_INF out file. default is "NOD_INF.OUT".
    times: int, optional
        Create a DataFrame with nodal values of the pressure head, the water
        content, the solution and sorbed concentrations, and temperature,
        etc, at the time "times". default is None.

    Returns
    -------
    data: dict
        Dictionary with the time as a key and a Pandas DataFrame as a value.

    """
    # Fixed column names from NOD_INF.OUT header
    _cols = ["Node", "Depth", "Head", "Moisture", "K", "C",
             "Flux", "Sink", "Kappa", "v/KsTop", "Temp"]

    use_times = []
    time_line_indices = []   # Line numbers of each "Time:" line

    with open(path) as file:
        lines = file.readlines()

    # Pass 1: find all "Time:" line positions and their values
    for i, line in enumerate(lines):
        if " Time:" in line and "Date" not in line:
            try:
                t = float(line.split(":")[1].strip())
                use_times.append(t)
                time_line_indices.append(i)
            except ValueError:
                pass

    if times is None:
        times = use_times

    # Build block boundaries: each block goes from its "Time:" line
    # to just before the next "Time:" line (or end of file)
    block_ends = time_line_indices[1:] + [len(lines)]

    def safe_to_numeric(col):
        try:
            return to_numeric(col)
        except ValueError:
            return col

    data = {}
    for t_idx, t_end, time in zip(time_line_indices, block_ends, use_times):
        if time not in times:
            continue
        # Within this block, find the "Node  Depth  Head..." header line
        node_header = None
        for j in range(t_idx, t_end):
            if "Node" in lines[j] and "Depth" in lines[j]:
                node_header = j
                break
        if node_header is None:
            continue
        # Data starts 3 lines after the Node header:
        #   node_header+1: units [L][L][-]...
        #   node_header+2: blank
        #   node_header+3: first data row
        data_start = node_header + 3
        block_lines = [l for l in lines[data_start:t_end] if l.strip()]
        if not block_lines:
            continue
        from io import StringIO
        df = read_csv(StringIO("".join(block_lines)), sep=r"\s+",
                      skipinitialspace=True, header=None, names=_cols)
        df = df.apply(safe_to_numeric)
        df.index = range(1, len(df) + 1)
        data[time] = df

    if len(data) == 1:
        return next(iter(data.values()))
    else:
        return data


@check_file_path
def read_balance(path="BALANCE.OUT", usecols=None):
    """
    Method to read the BALANCE.OUT output file.

    Parameters
    ----------
    path: str, optional
        String with the name of the run_inf out file. default is "BALANCE.OUT".
    usecols: list of str optional
        List with the names of the columns to import. By default:
        ['Area','W-volume','In-flow','h Mean','Top Flux', 'Bot Flux',
        'WatBalT','WatBalR'].

    Returns
    -------
    data: pandas.DataFrame
        Pandas with the balance data

    """
    if usecols is None:
        usecols = ["Area", "W-volume", "In-flow", "h Mean", "Top Flux",
                   "Bot Flux", "WatBalT", "WatBalR"]

    lines = open(path).readlines()
    use_times = []
    start = []
    end = [16]
    for i, line in enumerate(lines):
        for x in usecols:
            if x in line:
                line1 = x
                line2 = line.replace("\n", "").split(" ")[-1]
                line3 = line.replace("  ", " ").split(" ")[-2]
                lines[i] = [line1, line2, line3]

        if "Time" in line and "Date" not in line:
            time = float(
                line.replace(" ", "").split("]")[1].replace("\n", ""))
            use_times.append(time)
        if "Area" in line:
            start.append(i)
        if "WatBalR" in line:
            end.append(i + 1)
        if "Sub-region" in line:
            subreg = line.replace("  ", " ").replace("\n", "").split(" ")[-1]

    data = {}
    for s, e, time in zip(start, end, use_times):
        df = DataFrame(lines[s:e]).set_index(0).T
        index = {}
        for x in range(int(subreg) + 1):
            index[x + 1] = x
            df = df.rename(index=index)
        data[time] = df

    return data
