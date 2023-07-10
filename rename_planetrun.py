from createfort7dictionary import fort7dict
import os
import glob

name = []
layers = []


molef = fort7dict['MOLEF']
cloudtype = 'ALLCLOUDS_'
for x in range (0,len(molef)):
    if molef[x] == 0:
        cloudtype = 'NUCCLOUDS_'
if molef == [0,0,0,0,0,0,0,0,0,0,0,0,0]:
    cloudtype = ''



for x in fort7dict['THECOMMENT']:
    if x == ',':
         break
    name.append(x)

start = len(name) + 1
cutlayers = fort7dict['THECOMMENT'][start:]
for y in cutlayers:
    if y == ' ':
        break
    layers.append(y)

hazetype = fort7dict['HAZETYPE']
density = fort7dict['MTLX']
DENSE = ''
if density == 1:
    DENSE = '_DENSE'

rename_planetrun = ''.join(name) + '_' + str(cloudtype) + ''.join(layers) + 'LAYERS_' + str(hazetype).upper()  + DENSE
os.chdir('..')
os.system('pwd')
os.system('ls')
files = glob.glob('*')
if len(files) == 1:
    if files[0] == rename_planetrun:
        print('folder already named ' + rename_planetrun)
    else:
        os.system('mv -f * ' + rename_planetrun)
        print('renamed ' + files[0] + ' folder to ' + rename_planetrun)
        os.chdir(rename_planetrun)
else:
    raise Exception('renaming code cannot run if there is more than one file in source folder!!!')

