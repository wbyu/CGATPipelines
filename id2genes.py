import sys, re, string, os, getopt, time

USAGE = """python %s < stdin > stdout

substitute tokens

OPTIONS:
-c, --create            create substitution list
-a, --apply=            apply substitution list given by file.
-i, --invert            pairs of substitution patterns is to,from
-r, --regex-token       regular expression for tokens (has to create one pair of brackets)
-p, --pattern-sub       pattern for substitution
-m, --multiple          do multiple substitutions per row
-o, --columns-token     substitute tokens in columns
-e, --echo              echo susbstituted column
-f, --filter            remove lines not matching anything
-n, --inplace           do inplace subsitutions of all files on command line
-b, --backup            keep backup (with ending .bak)
-s, --skip              skip on error
--error=           output file for errors
""" % sys.argv[0]

param_long_options = ["create=", "regex-token=", "pattern-sub=",
                      "apply=", "invert", "multiple", "columns-token=", "echo",
                      "filter", "inplace", "backup", "skip", "error=" ]
param_short_options = "c:r:p:a:imo:fnbs"

param_create = None
param_regex_token = None
param_pattern_sub = "%s"
param_apply = None
param_invert = None
param_multiple = None
param_columns_token = None
param_echo = False
param_filter = False
param_filename_error = None

param_inplace = False
param_backup = False
param_separator = "|"
param_skip = False

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)
                                      

    except getopt.error, msg:
        print USAGE, msg
        sys.exit(1)

    for o,a in optlist:
        if o in ("-c", "--create"):
            param_create = a
        elif o in ("-r", "--regex-token"):
            param_regex_token = re.compile( a )
        elif o in ("-p", "--pattern-sub"):
            param_pattern_sub = a
        elif o in ("-a", "--apply"):
            param_apply = a
        elif o in ("-i", "--invert"):
            param_invert = 1
        elif o in ("-m", "--multiple"):
            param_multiple = 1
        elif o in ("-o", "--columns-token"):
            param_columns_token = map(lambda x: int(x) - 1,string.split(a, ","))
        elif o in ("-e", "--echo"):
            param_echo = True
        elif o in ("-f", "--filter"):
            param_filter = True
        elif o in ("-n", "--inplace"):
            param_inplace = True
        elif o in ("-b", "--backup"):
            param_backup = True
        elif o in ("-s", "--skip"):
            param_skip = True
        elif o == "--error":
            param_filename_error = a

    if param_inplace and len(args) == 0:
        raise "please supply file name(s) for inplace substitution."

    if param_create:
        replaced = {}

    if param_filename_error:
        err = open(param_filename_error, "a")
    else:
        err = sys.stderr
        
    file_id = 0

    keys = {}
    
    if not param_apply:
        raise "please specify ids."

    infile = open(param_apply, "r")
    for line in infile:
        if line[0] == "#": continue
        a = line[:-1].split("\t")[0]
        k = param_separator.join(a.split(param_separator)[:2])
        keys[k] = a
            
    files = args
    
    if not param_inplace and len(args) == 0:
        files = ["-"]

    rx = re.compile("([^%s\s>]+[%s][^%s\s]+[%s][^%s\s]+[%s][^%s\s]+)" % ((param_separator,) * 7) )
    
    for file in files:

        close_infile = False
        close_outfile =False
        
        if file == "-":
            infile = sys.stdin
            outfile = sys.stdout
        else:
            if param_inplace:
                os.rename( file, file + ".bak" )
                infile = open( file + ".bak", "r")
                outfile = open( file, "w")
                close_infile = True
                close_outfile = True
            else:
                infile = open( file, "r")
                outfile = sys.stdout
                close_infile = True
        
        for line in infile:
            if line[0] == "#":
                outfile.write(line)
                continue

            error = False
            
            if param_columns_token:
                data = line[:-1].split("\t")
                keep = True
                for c in param_columns_token:
                    k = param_separator.join(data[c].split(param_separator)[:2])
                    if k in keys:
                        if param_create:
                            replaced[data[c]] = keys[k]
                        data[c] = keys[k]
                    else:
                        err.write( "warning: did not find gene for prediction %s in line %s\n" % (k, line[:-1]))
                        error = True

                if error and param_skip: continue
                    
                outfile.write(string.join(data, "\t") + "\n")
            else:
                frags = []
                last_from = 0
                for m in rx.finditer( line[:-1]):
                    frags.append( line[last_from:m.start(1)] )
                    kk = m.groups(1)[0]
                    k = param_separator.join(kk.split(param_separator)[:2])
                    if k in keys:
                        frags.append( keys[k] )
                        if param_create:
                            replaced[kk] = keys[k]
                    else:
                        err.write( "warning: did not find gene for prediction %s in line %s\n" % (k, line[:-1]))
                        error = True
                        
                    last_from = m.end(1)
                    
                if error and param_skip: continue

                frags.append(line[last_from:])
                line = "".join(frags)

                outfile.write(line)
                
        if close_outfile: outfile.close()
        if close_infile: infile.close()

        if param_create:
            outfile_replaced = open( param_create, "w")
            for a,b in replaced.items():
                outfile_replaced.write("%s\t%s\n" % (a,b))
            outfile_replaced.close()
        
        if param_inplace and not param_backup:
            os.remove( file + ".bak" )

    if param_filename_error:
        err.close()