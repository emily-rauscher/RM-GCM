import sys

### How to use ——
### This python script is intended to replace the sed commands in the compile_nopg_restart script
### Parameters are added in order of : krun, begday, tspd, kits, lrstrt, tfrc
### Must have all 6 of these parameters in order to correctly assign values

def replacelines(seek, new):
    print(seek + ' is', new)
    with open('fort.7', 'r') as file:
        data = file.readlines()

    #Creates an array that returns -1 for each line that does not contain 'seek'
    # a value of != -1 means the location that 'seek' is found in that line
    line = [s.find(seek) for s in data]
    for notfound in line:
        if notfound != -1:
            found = notfound
    location = line.index(found)
    spaces = ''
    for i in range(0, found):
        spaces += ' '
    data[location] = spaces + seek + '     = '+ new + ',\n'
    with open('fort.7', 'w') as file:
        file.writelines(data)
    return


if len(sys.argv) !=7:
    print('\n                 ********** WARNING ************')
    print('****** system arguments did not read in correctly')
    print('****** make sure system argument order is krun, begday, tspd, kits, lrstrt, tfrc')
    print('****** and make sure that you include all of these params \n')
    raise Exception('read above')

#sys.argv are arguments given after the restart_fort7.py when typing in terminal
unused, krun, begday, tspd, kits, lrstrt, tfrc = sys.argv

replacelines('KRUN', krun)
replacelines('BEGDAY', begday)
replacelines('TSPD', tspd)
replacelines('KITS', kits)
replacelines('LRSTRT', lrstrt)
replacelines('TFRC', tfrc)
