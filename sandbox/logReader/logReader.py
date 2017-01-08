def getLogNumbers(path):
    """
    Get the simulation numbers from the log file.
    
    Parameters
    ----------
    path : str
        The path of the log file
        
    Retruns
    -------
    logNumbers : dict
        The numbers stored in a dict
    """

    with open(path,"r") as f:
        for line in f:
            if "Sim Time  |  RHS evals  | Wall Time |" in line:
               line = line.replace("\n", "")
               header = line.split("|")
               # Last part of list should now be split by " "
               headerWithSpaces = header.pop().split(" ")
               header = [*header, *headerWithSpaces]
               header = [head.replace(" ","") for head in header if head != ""]
               data   = OrderedDict({header[0]:[]})
               for head in header[0:]:
                   data[head]=[]
           
               break
           
        newlineCount = 0
        for line in f:
            if line == "\n":
                if newlineCount > 0:
                    break
                else:
                    newlineCount += 1
                    continue
        
            line = line.replace("\n", "")
            line = line.split(" ")
            line = [column for column in line if column != ""]
            for key, value in zip(data.keys(), line):
                data[key].append(float(value))
            
            
    for key in data.keys():
        data[key] = tuple(data[key])
        
    return dict(data)

