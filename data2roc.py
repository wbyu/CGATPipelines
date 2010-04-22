USAGE=""" calculate a roc curve (true positive rate versus false positive rate)

python data2roc.py [OTPIONS] < in > out

Input is from stdin, first field is +/- or 1/0 for a true
or false positive, respectively. The input needs to be
sorted appropriately. 

If the options --false-negative is set, the input is +/- or 1/0 for a 
true positive or false negative, respectively.

Output on stdout is a tab-separated table with the following
columns:

TP: true positives
FP: false positives
TPR: true positive rate  = true_positives /  predicted
P: predicted
FPR: false positive rate = false positives  / predicted
value: value
"""

import string, re, getopt, sys, optparse

import Experiment

def doMultiple( options ):

    first = True
    last_val = None
    headers, hists = [], []
    
    for line in sys.stdin:

        if line[0] == "#": continue

        data = line[:-1].split("\t")

        if first:
            headers = data[1:]
            ncols = len(headers)
            first = False
        else:
            val = data[0]
            if val != last_val:
                if last_val != None:
                    hists.append( (last_val, this_hist) )
                this_hist = [0] * ncols
                last_val = val

            for x in range(ncols):
                this_hist[x] += int(data[x+1])
                
    if last_val != None: hists.append( (last_val, this_hist) )

    if not headers: return
    options.stdout.write( "value\ttotal\t%s\t%s\n" % ("\t".join(headers),
                                                      "\t".join( ["p%s" % x for x in headers] ) ) )
    
    ntotal = 0
    hist = [0] * ncols
    for val, this_hist in hists:
        for x in range(ncols): hist[x] += this_hist[x]
        ntotal += sum( this_hist )
        options.stdout.write( "%s\t%i\t%s\t%s\n" % (val, ntotal,
                                                    "\t".join(["%i" % x for x in hist]),
                                                    "\t".join(["%f" % (float(x) / ntotal) for x in hist]) ) )


def main( argv = None ):

    if not argv: argv = sys.argv
 
    parser = optparse.OptionParser( version = "%prog version: $Id: data2roc.py 2782 2009-09-10 11:40:29Z andreas $", usage = USAGE)

    parser.add_option("-p", "--positives", dest="positives", type="float",
                      help="total number of true positives. If not set, take all true positives encountered [default=%default]." )

    parser.add_option("-m", "--monotonous", dest="monotonous", action="store_true",
                      help="bin data, so that sensitivity decreases monotonously [default=%default]." )

    parser.add_option("-b", "--bin", dest="bin", action="store_true",
                      help="bin data [default=%default]." )

    parser.add_option("-u", "--multiple", dest="multiple", action="store_true",
                      help="perform multiple flag analysis [default=%default]." )

    parser.add_option("-f", "--false-negatives", dest="false_negatives", action="store_true",
                      help="a negative flag indicates a false negative (and not a false positive) [default=%default]." )

    parser.set_defaults(
        positives = None,
        false_negatives = False,
        skip_redundant = True,
        monotonous = False,
        bin = False,
        raw = False,
        multiple = False,
        bin_by_score = True )

    options, args = Experiment.Start( parser, argv = sys.argv )

    if options.multiple:
        doMultiple( options )
        Experiment.Stop()
        sys.exit(0)

    true_positives = 0
    predicted = 0

    values = []

    last_value = None

    ninput, noutput, nskipped = 0, 0, 0

    for line in sys.stdin:
        if line[0] == "#": continue
        
        data = string.split(line[:-1], "\t")

        ninput += 1
        if data[0] not in ("+", "-", "1", "0"): 
            nskipped += 1
            continue
        
        if not options.bin_by_score:
            if not options.bin or last_value != data[1]:
                values.append( (true_positives, predicted, data[1]) )
        else:
            if last_value != None and last_value != data[1]:
                values.append( (true_positives, predicted, last_value) )

        predicted += 1

        if data[0] in ("+", "1"):
            true_positives += 1

        last_value = data[1]

    values.append( (true_positives, predicted, last_value) )
    values.append( (true_positives, predicted, data[1]) )

    if true_positives == 0:
        raise "# no true positives!"
    
    
    if options.positives == None:
        if options.false_negatives:
            positives = float(predicted)
        else:
            positives = float(true_positives)
    else:
        positives = float(options.positives)

    options.stdout.write( "value\tpred\tTP\tFP\tTN\tFN\tTPR\tFPR\tTNR\tFNR\tRTPR\tRFNR\n" )

    last_positives = None
    last_tpr = None

    for true_positives, predicted, value in values:

        if (predicted == 0):
            predicted = 1

        if options.false_negatives:
            false_negatives = predicted - true_positives
            false_positives = 0
            true_negatives = 0
        else:
            true_negatives = 0
            false_negatives = positives - predicted
            false_positives = predicted - true_positives

        tpr = float(true_positives) / predicted
        fpr = float(false_positives) / (true_positives + false_negatives )
        fnr = float(false_negatives) / positives
        tnr = 0

        # relative rates 
        rfpr = float(false_positives) / predicted
        rfnr = float(false_negatives) / predicted

        if options.monotonous and last_tpr and last_tpr < tpr:
            continue

        if options.skip_redundant and true_positives == last_positives:
            continue
        
        if (predicted > 0):
            options.stdout.write( "%s\t%i\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n" %\
                                      (value,
                                       predicted,
                                       true_positives,
                                       false_positives,
                                       true_negatives,
                                       false_negatives,
                                       tpr, fpr, tnr, fnr,
                                       rfpr, rfnr ) )


            noutput ++ 1

        last_positives = true_positives
        last_tpr = tpr

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))

    Experiment.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )