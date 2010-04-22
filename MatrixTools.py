import sys, re, string

class Matrix:

    def __init__(self):
        pass
    
def addOptions( parser ):
    """add matrices to option parser."""

    parser.add_option("-f", "--format", dest="format", type="string",
                      help="format."  )

    parser.add_option( "--row-names", dest="row_names", type="string",
                      help="list of row names separated by ','."  )
    
    parser.add_option( "--col-names", dest="col_names", type="string",
                      help="list of col names separated by ','."  )

    parser.add_option( "--file-row-names", dest="file_row_names", type="string",
                      help="filename with row names."  )
    
    parser.add_option( "--file-col-names", dest="file_col_names", type="string",
                      help="filename with col names."  )

    parser.add_option( "--numeric", dest="numeric", action="store_true",
                      help="row and column titles are numeric."  )

    parser.add_option( "--asymmetric", dest="asymmetric", action="store_true",
                      help="matrix is asymmetric."  )

    parser.add_option( "--default", dest="default", type="string",
                      help="default value for missing values."  )

    parser.add_option( "--default-diagonal", dest="default_diagonal", type="string",
                      help="default value for missing values on diagonal."  )

    parser.add_option( "--input-format", dest="input_format", type="choice",
                      choices=("row-col-weight", "", "row-col-weight-replicates"),
                      help="input format."  )

    parser.set_defaults(
        default = "0",
        default_diagonal = "0",
        format = "string",
        asymmetric = False,
        is_numeric = False,
        row_names = None,
        col_names = None,
        file_row_names = None,
        file_col_names = None,
        input_format = "row-col-weight",
        )

    
def getMapTokens( options ):

    map_token2row, map_token2col = {}, {}
    
    if options.file_row_names:
        row_tokens = map( lambda x: string.split(x[:-1], "\t")[0], open(options.file_row_names, "r").readlines())
        for row_token in row_tokens:
            map_token2row[row_token] = len(map_token2row)

    if options.row_names:
        for x in options.row_names.split(","):
            map_token2row[x] = len(map_token2row)

    if options.file_col_names:
        col_tokens = map( lambda x: string.split(x[:-1], "\t")[0], open(options.file_col_names, "r").readlines())
        for col_token in col_tokens:
            map_token2col[col_token] = len(map_token2col)

    if options.col_names:
        for x in options.col_names.split(","):
            map_token2col[x] = len(map_token2col)

    if not options.asymmetric:
        if map_token2row and not map_token2col:
            map_token2col = map_token2row
        elif map_token2col and not map_token2row:
            map_token2row = map_token2col

    return map_token2row, map_token2col

def getMatrixFromEdges( lines, options, in_map_token2row = {}, in_map_token2col = {}):
    """read matrix from lines
    """

    # remove comments
    lines = filter(lambda x:x[0] != "#" and len(x[:-1]) > 0, lines)

    if in_map_token2row:
        map_token2row = in_map_token2row
    else:
        map_token2row = {}

    if in_map_token2col:
        map_token2col = in_map_token2col
    else:
        map_token2col = {}

    if options.format == "string":

        has_row_names = len(map_token2row) > 0
        has_col_names = len(map_token2col) > 0
        
        ## if either row/column names are not given:
        if not map_token2row or not map_token2col:

            row_tokens = map(lambda x: string.split(x[:-1], "\t")[0], lines )
            col_tokens = map(lambda x: string.split(x[:-1], "\t")[1], lines )

            if options.is_numeric:
                row_tokens = map(float, row_tokens)
                col_tokens = map(float, col_tokens)                    
                row_tokens.sort()
                col_tokens.sort()
                row_tokens = map(str, row_tokens)
                col_tokens = map(str, col_tokens)                    
            else:
                row_tokens.sort()
                col_tokens.sort()

            if not has_row_names:
                for row_token in row_tokens:
                    if row_token not in map_token2row:
                        map_token2row[row_token] = len(map_token2row)
            if not has_col_names:
                for col_token in col_tokens:
                    if col_token not in map_token2col:
                        map_token2col[col_token] = len(map_token2col)

        if not options.asymmetric:
            for col_token in map_token2col.keys():
                if col_token not in map_token2row:
                    map_token2row[col_token] = len(map_token2row)            
            map_token2col = map_token2row

        matrix = [ [ options.default for j in range(len(map_token2col))] for i in range(len(map_token2row)) ]

        if len(map_token2col) == len(map_token2row):
            for j in range(len(map_token2col)):
                matrix[j][j] = options.default_diagonal

        ## return matrix
        m = Matrix()

        if options.input_format == "row-col-weight":
            for line in lines:
                row_token, col_token, weight = string.split(line[:-1], "\t")[:3]
                matrix[map_token2row[row_token]][map_token2col[col_token]] = weight
                if not options.asymmetric:
                    matrix[map_token2col[col_token]][map_token2row[row_token]] = weight
            
        elif options.input_format == "row-col-weight-replicates":
            replicates = [ [ 0 for j in range(len(map_token2col))] for i in range(len(map_token2row)) ]
            for line in lines:
                row_token, col_token, weight, n = string.split(line[:-1], "\t")[:4]
                matrix[map_token2row[row_token]][map_token2col[col_token]] = weight
                replicates[map_token2row[row_token]][map_token2col[col_token]] = int(n)
                if not options.asymmetric:
                    matrix[map_token2col[col_token]][map_token2row[row_token]] = weight
                    replicates[map_token2col[col_token]][map_token2row[row_token]] = int(n)                     
            m.mReplicates = replicates
            
        col_tokens = map_token2col.items()
        col_tokens.sort( lambda x,y: cmp(x[1],y[1]))
        row_tokens = map_token2row.items()
        row_tokens.sort( lambda x,y: cmp(x[1],y[1]))

        m.mMatrix = matrix
        m.mMapRow2Token = row_tokens
        m.mMapCol2Token = col_tokens
        m.mMapToken2Row = map_token2row
        m.mMapToken2Col = map_token2col        

    return matrix