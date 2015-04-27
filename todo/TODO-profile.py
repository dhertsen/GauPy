#! /usr/bin/python

import sys
import itertools
import copy
import gaussian
import numpy as np
from collections import OrderedDict
import pickle
import subprocess

# Check if the script is run locally.
echo_vsc = subprocess.Popen(['echo $VSC_INSTITUTE_LOCAL'], stdout=subprocess.PIPE, shell=True)
location, stderr = echo_vsc.communicate()
run_locally = location.strip() != 'gent'

from matplotlib import rc
import matplotlib as mpl

# Use latex fonts and text rendering if run locally.
if run_locally:
    rc('font',**{'family':'serif','serif':['Computer Modern']})
    rc('text', usetex=True)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import Subplot

def energyProfile(energies,positions,colors,legend,toplabels,bottomlabels,width,space,barlinewidth , connectionlinewidth , barlinestyle , connectionlinestyle , xmargin , ymargin , yaxis , toplabellinespace,bottomlabellinespace , toplabelsize , bottomlabelsize , lines, texts , ylabel, ylabelsize , yticksize , title , titlesize,toplabelcolors, bottomlabelcolors, barcolors, connectioncolors, legendposition, legendalignment, axes,**kwargs):

    # Open plot
    fig = plt.figure(1,(3,3))
    ax = Subplot(fig,111)
    fig.add_subplot(ax)

    # X values of positions
    minposition = min([min(position) for position in positions])
    maxposition = max([max(position) for position in positions])
    possiblepositions = range(minposition,maxposition+1)

    # x is a dictionary x[i] = true [xbegin,xend] of position i.
    x = dict()

    # xmin = 0 = x value of the left side of the first position
    # xmax = x value of the right side of the last position
    xmax = 0
    for position in possiblepositions:
        relativeposition = position - minposition
        begin = xmargin + relativeposition * (width + space)
        end = begin + width
        if end > xmax:
            xmax = end
        x[position] = [begin,end]

    # Minimal and maximal energy for layout purposes.
    energymax = max_incomplete_matrix(energies)
    energymin = min_incomplete_matrix(energies)

    # Iterate over profiles.
    for i,profile in enumerate(energies):

        # Iterate over horizontal bars in the same profile.
        for j,energy in enumerate(profile):

            # Select standard color of specified color.
            toplabelcolor = toplabelcolors[i][j] or colors[i]
            bottomlabelcolor = bottomlabelcolors[i][j] or colors[i]
            barcolor = barcolors[i][j] or colors[i]

            # Plot horizontal bars
            try:

                xj = x[positions[i][j]]
                yj = [energy,energy]
                ax.plot(xj,yj,color=barcolor,linewidth=float(barlinewidth[i][j]),linestyle=barlinestyle[i][j],label=legend[i])
                xtext = (xj[0] + xj[1]) * 0.5

                # Plot labels
                if toplabels[i][j]:
                    plt.text(xtext,energy+0.01*(energymax-energymin)*toplabellinespace[i][j],toplabels[i][j],size=toplabelsize[i][j],horizontalalignment='center',verticalalignment='bottom',color=toplabelcolor)

                if bottomlabels[i][j]:
                    plt.text(xtext,energy-0.01*(energymax-energymin)*bottomlabellinespace[i][j],bottomlabels[i][j],size=bottomlabelsize[i][j],horizontalalignment='center',verticalalignment='top',color=bottomlabelcolor)

            except:
                pass

            # Plot connection lines
            connectioncolor = connectioncolors[i][j] or colors[i]

            try:
                xbegin = x[positions[i][j]][1]
                xend = x[positions[i][j+1]][0]
                xj = [xbegin,xend]
                ybegin = energy
                yend = energies[i][j+1]
                yj = [ybegin,yend]
                ax.plot(xj,yj,color=connectioncolor,linewidth=float(connectionlinewidth[i][j]),linestyle=connectionlinestyle[i][j])
            except:
                pass

    # Plot extra lines
    if lines:
        for line in lines:

            # Relative x coordinates can be used. '1begin'->beginning of bar 1,
            # same for '1middle' and '1end'.
            for xdelimiter in 'x1','x2':
                line[xdelimiter] = relativeX(xdictionary=x,coordinate=line[xdelimiter])

            # Aliases and default values.
            x1 = line['x1']
            x2 = line['x2']
            y1 = line['y1']
            y2 = line['y2']
            linecolor = line.setdefault('color','black')
            linestyle = line.setdefault('style','-')
            linewidth = int(line.setdefault('width',1))

            # Plot extra line.
            ax.plot([x1,x2],[y1,y2],color=linecolor,linestyle=linestyle,linewidth=float(linewidth))

    # Plot extra texts
    if texts:
        for text in texts:

            # Defaults.
            textcolor = text.setdefault('color','black')
            textsize = text.setdefault('size','medium')
            # Relative x coordinates can be used.
            text['x'] = relativeX(xdictionary=x,coordinate=text['x'])
            plt.text(text['x'],text['y'],text['text'],size=text['size'],horizontalalignment='center',verticalalignment='center',color=text['color'])

    # Axes dimensions
    if not axes:
        absoluteYMargin = (ymargin*0.01) * (energymax-energymin)
        plt.axis([0, xmax+xmargin, energymin - absoluteYMargin, energymax + absoluteYMargin])
    else:
        plt.axis([0, xmax+xmargin, axes['ymin'],axes['ymax']])

    # Title
    plt.title(title,fontsize=titlesize)

    # y axis

    show_y_axis = {}

    if yaxis == 'both':
        show_y_axis['left'] = show_y_axis['right'] = True
    elif yaxis == 'right' or yaxis =='left':
        other_side = (set(['right','left']) - set([yaxis])).pop()
        show_y_axis[yaxis] = True
        show_y_axis[other_side] = False
    else:
        show_y_axis['left'] = show_y_axis['right'] = False

    for side in show_y_axis:
        if show_y_axis[side]:
            ax.axis[side].set_label(ylabel)
            ax.axis[side].label.set_size(ylabelsize)
            ax.axis[side].major_ticklabels.set_size(yticksize)
            ax.axis[side].toggle(ticklabels=True,label=True)
        else:
            ax.axis[side].set_visible(False)

    # Hide x axes
    ax.axis['top'].set_visible(False)
    ax.axis['bottom'].set_visible(False)

    # Legend

    # Empty legend will be provided as [None ...]
    showLegend = sum([bool(i) for i in legend])
    if showLegend:
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        if legendposition:
            legendx = legendposition['x']
            legendy = legendposition['y']
            plt.legend(by_label.values(), by_label.keys(),loc = legendalignment, bbox_to_anchor=[legendx, legendy])
        else:
            plt.legend(by_label.values(), by_label.keys(),loc = legendalignment)

    # Show plot
    return plt


def add_incomplete_matrices(matrix1,matrix2,multiplication=1):

    a = copy.deepcopy(matrix1)
    b = copy.deepcopy(matrix2)

    for i,row in enumerate(a):
        for j,element in enumerate(row):
            try:
                # First, try to add the corresponding (i,j) element,
                # i.e. matrix1 and matrix2 have the same dimensions.
                new = float(a[i][j]) + multiplication * float(b[i][j])
                a[i][j] = new
            except:
                try:
                    # Then, handle the case in which matrix2 is a list of one-element
                    # lists, each containing a single value.
                    new = float(a[i][j]) + multiplication * float(b[i][0])
                    a[i][j] = new
                except:
                    # As a last resort, interpret b as a list of values.
                    new = float(a[i][j]) + multiplication * float(b[i])
                    a[i][j] = new
    return a

def multiply_incomplete_matrices(matrix1,matrix2):

    a = copy.deepcopy(matrix1)
    b = copy.deepcopy(matrix2)

    for i,row in enumerate(a):
        for j,element in enumerate(row):
            try:
                # First, try to multiply the corresponding (i,j) elements,
                # i.e. matrix1 and matrix2 have the same dimensions.
                new = float(a[i][j]) * float(b[i][j])
                a[i][j] = new
            except:
                try:
                    # Then, handle the case in which matrix2 is a list of one-element
                    # lists, each containing a single value.
                    new = float(a[i][j]) * float(b[i][0])
                    a[i][j] = new
                except:
                    # As a last resort, interpret b as a list of values.
                    new = float(a[i][j]) * float(b[i])
                    a[i][j] = new
    return a

def print_incomplete_matrix(matrix,decimals=1):
    for i,row in enumerate(matrix):
        if type(row) == list:
            print(('%-15.4f' * len(row)) % tuple(row))
        else:
            print('%-15.4f' % row)

def max_incomplete_matrix(matrix):

    maximum = 0

    if type(matrix[0]) == list:
        maxima = []
        for row in matrix:
            maxima.append(max(row))
        maximum = max(maxima)
    else:
        maximum = max(matrix)

    return maximum

def min_incomplete_matrix(matrix):

    minimum = 0

    if type(matrix[0]) == list:
        minima = []
        for row in matrix:
            minima.append(min(row))
        minimum = min(minima)
    else:
        minimum = min(matrix)

    return minimum


def parseLine(line,replace=None):

    """
     Split a single line into a list of strings according to this scheme:
       - whitespaces, equality signs and double quotes are used to split lines
       - no splitting within double quotes
       - only empty double quotes can be used to include an empty string in the list
       - equal signs are included in the list
     Example:
       '    option  arg1=bla1  arg2="bla bla 2" arg3=""   '
       --> ['option','arg1','=','bla1','arg2','=','bla bla 2','arg3','=','']
    List items will be substituted according to the replace dictionary.
    Comment lines starting with '#" will be ignored.
    """

    parsedLine = []

    # Parse non-comment lines
    if not line.strip().startswith('#'):

        # Split according to quotes first, in this way the evenly numbered fragments
        # constitute the text inside the quotes, et vice versa.
        quoteSplit = line.split('"')

        # Iterate over all quote splitted fragments.
        for i,quoteSplittedFragment in enumerate(quoteSplit):
            # Text within the quotes may be added without further modification.
            if i%2:
                parsedLine.append(quoteSplittedFragment)
            # Text within the quotes will be processed further.
            else:
                fullySplittedOutsideQuotes = [a.partition('=') for a in quoteSplittedFragment.split()]
                fullySplittedOutsideQuotes = itertools.chain(*fullySplittedOutsideQuotes)
                fullySplittedOutsideQuotes = [item for item in fullySplittedOutsideQuotes if item]
                parsedLine.extend(fullySplittedOutsideQuotes)

        # Substitute strings in the list according to replace dictionary
        if replace:
            parsedLine = [a if not a in replace else replace[a] for a in parsedLine]

    return parsedLine

def parseBlockOption(blockName,config,dimensions=None,replace={'none':''},stripBlock=True,stripBlockType=list):

    # Discuss that it's nested list of strings
    # Dimension expansion: way of horizontal, vertical expanding. Ignoring items possible as well.
    # non-existent blocks ==> [] output
    # block delignation, case insensitive, starts with &, may include other words, end must include end


    # Parse the config file line-wise.
    if type(config) == str:
        configLines = config.splitlines()
    elif type(config) == list:
        configLines = copy.copy(config)

    # Gather all errors in this list.
    errors = []

    # Lengths of the lines, i.e. before dimension expansion
    lengths = []

    # Block that will be returned
    block = []

    # Default values for the beginning and the end of the block.
    startBlock = 0
    endBlock = 0

    # Look for the start and the end of the block.
    for i,configLine in enumerate(configLines):

        # Dummy manipulations of the line.
        lowStripLine = configLine.lower().strip()
        lowStripSplitLine = lowStripLine.split()

        # Checks of the line.
        startingSign = lowStripLine.startswith('&')
        blockNameOnLine = blockName.lower() in lowStripLine
        endKeywordOnLine = 'end' in lowStripSplitLine or '&end' in lowStripSplitLine

        # Use checks to mark blocks.
        if startingSign and blockNameOnLine:
            if endKeywordOnLine:
                endBlock = i
            else:
                startBlock = i

    # If both startBlock and endBlock == 0, the block was not found.
    if not (startBlock or endBlock):
        errors.append(blockName+' block not found.')

    else:
        blockLength = endBlock - startBlock - 1

        # Start after end line, or end line missing.
        if blockLength < 0:
            errors.append(blockName+' block: incorrect placement of start and/or end lines.')

        # Empty block is not valid.
        elif blockLength == 0:
            errors.append(blockName+' block is empty.')

        # Now, in non erroneous case.
        else:

            # List of all lines in the block.
            blockLines = configLines[startBlock+1:endBlock]

            # Remove empty lines.
            blockLines =  [line for line in blockLines if line.strip()]

            # Flag for special type of block containing '=' fragments
            # in their parsed lines. If '=' occurs in the block between
            # double quotes this is of course not true. See further.
            specialType = False

            for i,blockLine in enumerate(blockLines):

                # Split the line in the correct way.
                parsedLine = parseLine(line=blockLine,replace=replace)

                # Block options with lines  that contain '=' in their parsed versions,
                # will be parsed in a special way. They will be parsed as a list of
                # dictionaries (1 dictionary per line):
                #
                # 'x=1 y=4 text="test" \n x=2 y=4 text="ijsje"'
                # --> [ {'x':'1','y':'4','text':'test'} , {'x':'2','y':'4','text':'ijsje'} ]
                #
                # For these blocks, all dimensions will be ignored.

                specialType = specialType or ('=' in parsedLine)
                if specialType:
                    keys = parsedLine[::3]
                    values = parsedLine[2::3]
                    parsedLine = dict(zip(keys,values))
                    dimensions = None

                # Expand the dimensions of the line if necessary and possible
                # i.e. if non-empty line i in the block has less 'items' than
                # dimensions[i], and dimensions[i] % number of 'items' == 0, then
                # the items will be duplicated until dimensions[i] items are reached.
                #
                # Example:
                # parsedLine = ['a','b']
                # dimensions[i] = 6
                # --> parsedLine = ['a','b','a','b','a','b']
                #
                # The above is true if dimensions is a list of integers. If dimensions
                # is an integer itself, a column will be created from the block. In that
                # case, only the first element on the line is selected into the parsed block.

                lineLength = len(parsedLine)
                lengths.append(lineLength)

                if type(dimensions) == list:
                    lineDimension = dimensions[i]
                    try:
                        parsedLine = [copy.copy(parsedLine[i%lineLength]) for i in range(lineDimension)]
                    except:
                        errors.append(blockName+' block: error in expanding dimensions.')
                elif type(dimensions) == int:
                    parsedLine = parsedLine[0]
                block.append(parsedLine)

            # Check for superfluous dimensions. The shortage of lines will be taken care of
            # in a similar way as the number of words on a line:
            #
            # If len(dimensions) % len(block) == 0, then the rows will be duplicated until the number
            # of dimensions if fulfilled.

            if type(dimensions) == list and len(block) < len(dimensions):

                try:
                    # Number of lines
                    blockLength = len(block)
                    # Wanted number of lines
                    blockDimension = len(dimensions)

                    # Iterate over the 'nonspecified' lines.
                    for i in range(blockLength,blockDimension):

                        # Every extra line refers to a reference line, i.e. if 3 lines are
                        # in the block and 5 are asked for, the 4th line refers to the 1st
                        # and the 5th to the 2nd.  # Length of the (unexpanded original line)
                        originalLength = lengths[i%blockLength]
                        # (Unexpanded) original line
                        originalLine = block[i%blockLength][:originalLength]
                        # Length of the extra line
                        extraLineDimension = dimensions[i]

                        extraLine = [copy.copy(originalLine[j%originalLength]) for j in range(extraLineDimension)]
                        block.append(extraLine)
                except:
                    errors.append(blockName+' block: error in expanding dimensions.')

            # If dimensions is an integer, a column matrix is created with the proper dimensions.
            # The expasions happens in the usual way. Superfluous lines are ignored, missing lines
            # are filled in from the lines above.
            elif type(dimensions) == int and len(block) < dimensions:

                try:
                    # Number of lines
                    blockLength = len(block)
                    block = [block[i%blockLength] for i in range(dimensions)]
                except:
                    errors.append(blockName+' block: error in expanding dimensions.')


    # If stripBlock, the block will be stripped from the original config and this new config
    # will be returned. stripBlockType (str,list) determines the type of the return.
    if stripBlock and (stripBlockType in [str,list]):

        # If the (not startBlock == endBlock) would not be included, the first line of the config
        # would be deleted if parseBlockOption() is called for a non-existent block.
        if not startBlock == endBlock:
            configLines[startBlock:endBlock+1] = []

        if stripBlockType == str:
            configLines = '\n'.join(configLines)
        return block,configLines,errors
    else:
        return block,errors


# only works for nested lists
def evaluateBlock(block,parameter='gibbs',returnFilenames=False):

    # Store file-evaluated values in a dictionary for efficiency
    evaluated = {}

    # Do not modify the original block.
    evaluatedBlock = copy.deepcopy(block)

    # In this nested list the only cell with a value!=None
    # are the cells in the corresponding block that were actually
    # filenames.
    foundFilenames = copy.deepcopy(block)

    # Iterate over all cells
    for i,row in enumerate(block):

        # Two-dimensional matrices.
        if type(row) == list:
            for j,cell in enumerate(row):
                try:
                    # First check whether they are numerical values and convert
                    # them to floats.
                    evaluatedBlock[i][j] = float(evaluatedBlock[i][j])
                    foundFilenames[i][j] = None
                except:
                    try:
                        # Then, try to evaluate them gaussian output file names. It can be
                        # either a single filename or a group op filenames that was enclosed
                        # by double quotes.

                        # A list of the single filename or the quote-enclosed filenames.
                        filenames = evaluatedBlock[i][j].split()

                        # Evaluated values for the filenames.
                        # Check if the files have been evaluated before.
                        # If one of the values within the double quotes is an integer, it is stored
                        # as an int in the values list. Such a value will be multiplied with filename
                        # that follow:
                        # 'file1 2 file2' --> value(file1) + 2 * value(file2)

                        # Keep the integers as strings here.
                        values = [ float(evaluated.get(filename,gaussian.analyze_log(filename)[parameter])) if not filename.isdigit() else filename for filename in filenames ]
                        # Substitute the integers by (integer-1)*next element.
                        values = [(int(fragment)-1) * values[k+1] if type(fragment) == str else fragment for k,fragment in enumerate(values) ]

                        # Set stored values
                        for k,filename in enumerate(filenames):
                            evaluated[filename] = values[k]

                        # Evaluate this cell as the sum of all.
                        evaluatedBlock[i][j] = sum(values)

                    # If none of these conversions work, the original string is maintained.
                    except:
                        foundFilenames[i][j] = None
                        print(evaluatedBlock[i][j]+' could not be parsed.')

        # Column matrices, same strategy.
        else:
            try:
                evaluatedBlock[i] = float(row)
                foundFilenames[i] = None
            except:
                try:
                    filenames = evaluatedBlock[i].split()
                    values = [ float(evaluated.get(filename,gaussian.analyze_log(filename)[parameter])) if not filename.isdigit() else filename for filename in filenames ]
                    values = [ (int(fragment)-1) * values[k+1] if type(fragment) == str else fragment for k,fragment in enumerate(values) ]
                    for k,filename in enumerate(filenames):
                        evaluated[filename] = values[k]
                    evaluatedBlock[i] = sum(values)
                except:
                    foundFilenames = None

    # Print all the evaluated files.
    try:
        length_filenames = max([len(key) for key in evaluated])
        print
        print('%s of the following files was evaluated:' % parameter)
        for evaluated_file in evaluated:
            print('    %-*s %s' % (length_filenames+5,evaluated_file,evaluated[evaluated_file]))
    except:
        pass

    # Return the evaluated block.
    if returnFilenames:
        return evaluatedBlock,foundFilenames
    else:
        return evaluatedBlock

def evaluateEnergies(energies,corrections=0,reference=0,conversion=1,parameter={'energies':'gibbs','reference':'gibbs','corrections':'gibbs'},store_energies=None):

    # Evaluate the energies
    E,foundFilenames = evaluateBlock(energies,parameter=parameter.get('energies','gibbs'),returnFilenames=True)

    print
    print('energies:')
    print_incomplete_matrix(E)

    # Evaluate the energy corrections
    Ecorr = evaluateBlock(corrections,parameter=parameter.get('corrections','gibbs'))

    print
    print('corrections:')
    print_incomplete_matrix(Ecorr)

    # Add the corrections to the energy.
    Ecorrected = add_incomplete_matrices(matrix1 = E, matrix2 = Ecorr)

    print
    print('corrected energies:')
    print_incomplete_matrix(Ecorrected)

    # Evaluate the conversion factors. In this case 'evaluating' only implicates
    # the conversion to a matrix of floats.
    conversion = evaluateBlock(conversion)

    print
    print('conversions:')
    print_incomplete_matrix(conversion)

    # Parse the reference.
    parsedReference = copy.deepcopy(reference)

    # Construct a ndarray with numerical values that can be substracted
    # from the energies.

    # If reference is passed as a single keyword it can be either a numerical value,
    # 'value-#', 'value-#-#' or a filename. Convert this to a column matrix of the correct dimensions.

    if not type(reference) == list:
        parsedReference = [[copy.deepcopy(parsedReference)] for i in energies]

    # If reference is passed as a list (column matrix), it describes the reference value for each
    # row. Such a value can be either a numerical value, a filename or a string of
    # the form 'value-#' to refer to the numerical value of position# in that profile
    # or a string of the form 'value-#-#' to refer to entry row #, column #.
    #
    # The reference value will be taken from the corrected energy (E=E+Ecorr). In this way,
    # gibbs free energies will be referenced to gibss free energies and not to electronic
    # energies.
    #
    # If the reference was passed as a single keyword, we will use the column matrix
    # above for further manipulations.

    # Just subtitute the 'value-#' and 'value-#-#' entries, and pass reference to evaluateBlock().


    for i,row in enumerate(parsedReference):
        if str(row).lower().startswith('value'):
            try:
                refindices = [int(a) for a in str(row).split('-')[1:]]
                if len(refindices) == 1:
                    parsedReference[i] = Ecorrected[i][refindices[0]-1]
                elif len(refindices) == 2:
                    parsedReference[i] = Ecorrected[refindices[0]-1][refindices[1]-1]
            except:
                pass


    Eref = evaluateBlock(parsedReference,parameter=parameter.get('reference','gibbs'))
    print
    print('references:')
    print_incomplete_matrix(Eref)

    Erelative = add_incomplete_matrices(matrix1 = Ecorrected, matrix2 = Eref, multiplication = -1)
    print
    print('relative, corrected energies (not converted):')
    print_incomplete_matrix(Erelative)

    Efinal = multiply_incomplete_matrices(matrix1 = Erelative, matrix2 = conversion)
    print
    print('relative, corrected energies (converted):')
    print_incomplete_matrix(Efinal)
    print

    if store_energies:
        try:
            energy_overview = {'energies':E, 'corrections':Ecorr, 'references':Eref, 'profile':Efinal}
            pickle.dump(energy_overview,open(store_energies,'wb'))
        except:
            pass

    return Efinal,foundFilenames

def parseLineOption(line,replace={'none':''}):

    # Line is passed as a string, convert into its fragments
    # with parseLine()
    parsedLine = parseLine(line,replace=replace)

    # Evaluate extended line options of the type
    # keyword param1=value1 param2=value2
    if '=' in parsedLine:
        keyword = parsedLine[0]
        parameters =  parsedLine[1::3]
        values = parsedLine[3::3]

        # Return a dictionary {keyword:{param1:value1,param2:value2}}
        return {keyword:dict(zip(parameters,values))}

    # Empty lines.
    elif not parsedLine:
        return dict()

    # Evaluated line options of the type
    # keyword value.
    else:
        keyword = parsedLine[0]
        value = parsedLine[1]

        # Return dictionary {keyword:value}
        return {keyword:value}

def formatOption(optionValue,optionFormat):
    if optionFormat == bool:
        if optionValue == 'yes':
            optionValue = True
        elif optionValue == 'no':
            optionValue = False
    elif optionFormat == float:
        try:
            optionValue = float(optionValue)
        except:
            pass
    elif optionFormat == int:
        try:
            optionValue = int(optionValue)
        except:
            pass
    return optionValue

# already formatted may be inside
def formatOptions(optionValues,optionFormats):
    for option in optionFormats:
        try:
            optionFormat = optionFormats[option]
            optionValue = optionValues[option]

            if optionFormat in [bool,float,int]:
                if type(optionValue) == list:
                    for i,row in enumerate(optionValue):
                        for j,column in enumerate(row):
                            optionValues[option][i][j] = formatOption(optionValue[i][j],optionFormat)
                else:
                    optionValues[option] = formatOption(optionValue,optionFormat)
            elif type(optionFormat) == dict:
                for parameter in optionFormat:
                    try:
                        parameterValue = optionValue[parameter]
                        parameterFormat = optionFormat[parameter]
                        optionValues[option][parameter] = formatOption(parameterValue,parameterFormat)
                    except:
                        pass
        except:
            pass

# normally, both return and set
def expandLineOption(lineOption,dimensions):

    # Keyword and value of the line option.
    key = lineOption.keys()[0]
    value = lineOption.values()[0]

    # Request for a nested list, i.e. dimensions
    # is a nested list of integers.
    if type(dimensions) == list:
        # This is the normal, i.e. the value of the line option is given as a single string or number.
        if not type(value) == list:
            expandedValue = [[copy.copy(value) for j in range(dimensions[i])] for i in range(len(dimensions))]
        # If the value of the line option is given as one dimensional list and a higher dimension is required,
        # this will also be taken care of. This was added for the default value of 'references': it is provided
        # as list(range(1,max(dimensions)+1)). This has to be expanded into a 2D matrix.
        else:
            expandedValue = [[copy.copy(value[j%len(value)]) for j in range(dimensions[i])] for i in range(len(dimensions))]
    elif type(dimensions) == int:
        expandedValue = [copy.copy(value) for i in range(dimensions)]
    lineOption[key] = expandedValue
    return expandedValue

def parseLabels(toplabels,bottomlabels,energies,foundFilenames,roundenergies=1):


    # Iterate over every entry in toplabels and bottomlabels.
    for labels in toplabels,bottomlabels:
        for i,labelrow in enumerate(labels):
            for j,label in enumerate(labelrow):

                # Substitute 'energy' by the energies in enegries.
                # Round according to roundenergies.
                if label == 'energy':
                    labels[i][j] = str(round(energies[i][j],roundenergies))

                # Substitute 'filename' by the filename without extension,
                # if the filename was present in the original energies block.
                # Otherwise, '' will be substituted.
                elif label == 'filename':
                    if foundFilenames[i][j]:
                        filenameWithoutExtension = foundFilenames[i][j].replace('.log','')
                        labels[i][j] = filenameWithoutExtension
                    else:
                        labels[i][j] = ''

                # Substitute the unique part of the filename, i.e. the part that
                # is different from every other filename in the block. Filenames
                # are compared as '-' seperated fragments. If no filename was present
                # in foundFilenames, substitute ''.
                elif label == 'uniquefilename':
                    filej = foundFilenames[i][j]
                    if filej == None:
                        labels[i][j] = ''
                    else:
                        longestUniqueString = ''
                        for k,filelistk in enumerate(foundFilenames):
                            for l,filel in enumerate(filelistk):
                                if (filej != filel) and filel:
                                    filejroot = filej.replace('.log','')
                                    filelroot = filel.replace('.log','')
                                    uniqueParts = set(filejroot.split('-')) - set(filelroot.split('-'))
                                    uniqueString = '-'.join(uniqueParts)
                                    if len(uniqueString) > len(longestUniqueString):
                                        longestUniqueString = uniqueString
                        labels[i][j] = longestUniqueString

                # The same thing as 'uniquefilename', but per profile (row).
                elif label == 'uniqueprofilefilename':
                    filej = foundFilenames[i][j]
                    if filej == None:
                        labels[i][j] = ''
                    else:
                        longestUniqueString = ''
                        for k,filek in enumerate(foundFilenames[i]):
                            if (filej != filek) and filek:
                                filejroot = filej.replace('.log','')
                                filekroot = filek.replace('.log','')
                                uniqueParts = set(filejroot.split('-')) - set(filekroot.split('-'))
                                uniqueString = '-'.join(uniqueParts)
                                if len(uniqueString) > len(longestUniqueString):
                                    longestUniqueString = uniqueString
                        labels[i][j] = longestUniqueString

def relativeX(xdictionary,coordinate):
    try:
        if 'begin' in coordinate:
            bar = int(coordinate.replace('begin',''))
            barposition = 0
            relx = xdictionary[bar][barposition]
        elif 'end' in coordinate:
            bar = int(coordinate.replace('end',''))
            barposition = 1
            relx = xdictionary[bar][barposition]
        elif 'middle' in coordinate:
            bar = int(coordinate.replace('middle',''))
            relx = float(sum(xdictionary[bar]))/2
        return relx
    except:
        return coordinate


# Obtaining the files
inputFiles = sys.argv[1:]
cfgFiles = [inputFile for inputFile in inputFiles if '.cfg' in inputFile]
logFiles = [inputFile for inputFile in inputFiles if '.log' in inputFile]

# Iterate over every config file
for cfgFile in cfgFiles:

    # Options
    options = dict()

    # Store errors
    errors = []

    # Lines in the config file
    cfgLines = open(cfgFile).readlines()

    # Parse the energies. This block will determine
    # to which dimensions the other blocks will be expanded.
    energyBlock,cfgLines,energyErrors = parseBlockOption(blockName='energies',config=cfgLines)

    # Store the errors and the parsed values.
    # The values are stored in options['originalEnergies'], because
    # the filenames in this block may be used for labels afterwards.
    options['originalEnergies'] = energyBlock
    errors.append(energyErrors)

    # Set the dimensions
    # dimensions is a list of integer with dimensions[i] = number of items on row i
    dimensions = [len(row) for row in energyBlock]
    # nrows = number of rows = dimension of column matrix
    nrows = len(dimensions)

    # Set the default options with the correct dimensions.
    # Be aware entries like 'none' and 'value#' can only be used in the
    # config file and not in this dictionary of default values. This dictionary.keys()
    # provides an overview of possible options.
    options.update({'showplot':False,'saveplot':False,'conversions':1,'corrections':0,'references':0,'roundenergies':1,'colors':'black','positions':list(range(1,max(dimensions)+1)),'barlinewidth':3,'connectionlinewidth':1.5,'barlinestyle':'-','connectionlinestyle':'--','xmargin':3,'ymargin':10,'yaxis':'both','toplabellinespace':1,'bottomlabellinespace':1,'toplabelsize':'large','bottomlabelsize':'large','ylabelsize':'medium','yticksize':'medium','titlesize':'large','legend':'','toplabels':None,'bottomlabels':None,'lines':None,'texts':None,'title':'','width':7,'space':2,'ylabel':'','toplabelcolors':None,'bottomlabelcolors':None,'barcolors':None,'axes':None,'size':None,'connectioncolors':None,'energytype':{'energies':'gibbs','reference':'gibbs','corrections':'gibbs'},'legendalignment':'upper right','legendposition':None, 'margins':None})

    # Most block options, besides energies, can be provided either as block options or
    # as single line options. In both cases, some of them will be parsed as column matrices,
    # i.e. a single property for every profile. Some of them will be parsed as full matrices
    # with the same dimensions as the 'energies' block, i.e. a single property for every
    # horizontal energy bar.
    #
    # A special type of block option (see optionsListsOfDictionaries) can only be provided as
    # a block.

    # So the possible block options are ['energies'] + optionsColumnMatrix + optionsFullyExpanded + optionsListsOfDictionaries.

    # Options that will be parsed as column matrices.
    optionsColumnMatrix = ['colors','references','legend','conversions']
    # Options that will be parsed as matrices with the same dimensions as energies.
    optionsFullyExpanded = ['toplabels','bottomlabels','positions','corrections','toplabelcolors','bottomlabelcolors','barcolors','barlinewidth','barlinestyle','toplabellinespace','bottomlabellinespace','toplabelsize','bottomlabelsize','connectioncolors','connectionlinewidth','connectionlinestyle']
    # Options that will be parsed as non-expanded lists of dictionaries.
    optionsListsOfDictionaries = ['texts','lines']


    # Parse the other block options, strip them from cfgLines,
    # store the errors and values.
    for blockOption in optionsFullyExpanded + optionsColumnMatrix + optionsListsOfDictionaries:

        # Set the correct dimensions.
        if blockOption in optionsColumnMatrix:
            blockDimensions = nrows
        elif blockOption in optionsFullyExpanded:
            blockDimensions = dimensions
        else:
            blockDimensions = None

        # Parse the block.
        blockValues,cfgLines,blockErrors = parseBlockOption(blockName=blockOption,config=cfgLines,dimensions=blockDimensions)

        # Only overwrite default options values if the blocks contains information.
        if blockValues:
            options[blockOption] = blockValues
        # If a standard value was given, no new value was found in the config file:
        # expand the default value to the right dimensions. If the option is in
        # optionsListsOfDictionaries do not expand.
        elif blockOption in options and not blockOption in optionsListsOfDictionaries:
            options[blockOption] = expandLineOption({blockOption:options[blockOption]},blockDimensions)

        errors.append(blockErrors)

    # Parse all the other lines as line options. Comment lines
    # will be ignored automatically. Only store the options if the corresponding
    # block option has not been found and parsed. Line options
    # have priority over block options: a line option will always
    # overwrite given block and default values.

    # Counter for unsigned blocks
    blocksignCounter = 0

    # Parse the lines.
    for i,line in enumerate(cfgLines):

        # Block sign encounter = unparsed block.
        blockDelimiter = line.strip().startswith('&')
        if blockDelimiter:
            blocksignCounter += 1
        inBlock = blocksignCounter % 2

        # True line options, i.e. not in an unparsed block.
        if not inBlock and not blockDelimiter:

            # Parse the line.
            lineOption = parseLineOption(line)

            # Expand to the correct dimensions.
            if set(lineOption) & set(optionsColumnMatrix):
                lineDimension  = nrows
                expandLineOption(lineOption,lineDimension)
            elif set(lineOption) & set(optionsFullyExpanded):
                lineDimension = dimensions
                expandLineOption(lineOption,lineDimension)

            # Line options that correspond to blocks that must be parsed as lists
            # of dictionaries.
            elif set(lineOption) & set(optionsListsOfDictionaries):
                errors.append(lineOption.keys()[0]+' cannot be a line option.')
                lineOption = dict()

            # Non-empty lines with unknown line options.
            elif lineOption and not set(lineOption) & set(options):
                errors.append(lineOption.keys()[0]+' is not a valid line option.')
                lineOption = dict()

            options.update(lineOption)

    # Report unparsed blocks
    if blocksignCounter:
        errors.append('%i unparsed blocks' % (blocksignCounter/2))

    # Convert the option values to the correct type. Non-specified options will keep their string values.
    optionFormats = {'conversions':float, 'width':float, 'space':float, 'barlinewidth':float, 'connectionlinewidth':float, 'barlinestyle':float, 'connectionlinestyle':float, 'xmargin':float, 'ymargin':float, 'yaxis':str, 'toplabellinespace':float,'bottomlabellinespace':float, 'toplabelsize':float, 'bottomlabelsize':float, 'ylabelsize':float, 'yticksize':float, 'saveplot':bool, 'showplot':bool,'positions':int,'axes':{'ymin':float,'ymax':float},'size':{'width':float,'height':float},'roundenergies':int, 'legendposition':{'x':float,'y':float}, 'margins':{'left':float,'right':float,'bottom':float,'top':float,'wspace':float,'hspace':float}}
    formatOptions(optionValues=options,optionFormats=optionFormats)

    # Convert the energies.
    if not options['saveplot']:
        store_energies = options['saveplot']
    else:
        store_energies = options['saveplot'][:options['saveplot'].find('.')] + '.pickle'

    options['energies'],options['foundFilenames'] = evaluateEnergies(energies=options['originalEnergies'],corrections=options['corrections'],reference=options['references'],conversion=options['conversions'],parameter=options['energytype'],store_energies=store_energies)

    # Parse the labels
    parseLabels(toplabels=options['toplabels'],bottomlabels=options['bottomlabels'],foundFilenames=options['foundFilenames'],energies=options['energies'],roundenergies=options['roundenergies'])

    # Make the scheme
    profiles = energyProfile(**options)

    # Adjust size
    if options.get('size'):
        profiles.gcf().set_size_inches(options['size']['width'],options['size']['height'])

    # Adjust margins
    if options['margins']:
        plt.subplots_adjust(**options['margins'])

    # Save and show.
    if options['saveplot']:
        profiles.savefig(options['saveplot'],bbox_inches='tight',transparent=True, format='pdf')
        print(options['saveplot']+' written.')
        print
    reallyshowplot = not options['saveplot'] or options['showplot']
    if reallyshowplot:
        profiles.show()

    # Close the figure after each config file.
    profiles.close()

    # document py files
    # put everything in one py file

    # make scheme of how options can be included,
    # general info how they are expanded

    # line options have priority
    # legend style and color when they are changed?
