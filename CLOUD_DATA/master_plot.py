import numpy as np
import matplotlib.pylab as plt
import shutil
import os
plt.style.use('science.mplstyle')
species = ['Al2O3', 'CaTiO3', 'Cr', 'Fe', 'Ni', 'KCl', 'Mg2SiO4', 'MnS', 'Na2S', 'SiO2', 'VO', 'ZnS']
types = ['qext', 'pi0', 'gg']
os.system('pwd')
kind = input('Rosseland or Planck or spectral? ')
if kind[0].upper() == 'R':
    name = 'rosselandMean'
elif kind[0].upper() == 'P':
    name = 'PlanckMean'
elif kind[0].upper() == 'S':
    name = 'wav'
else:
    raise Exception('Please type rosseland or planck or spectral')

for sp in species:
    # read in k_ext data
    fname = sp+'_wav_qext.txt'

    f = open(fname, 'r')
    line = f.readline().split()
    na = int(line[0])
    nT = int(line[1])
    sigma = float(line[2])
    line = f.readline().split()
    a = np.zeros(na)
    a[:] = line[:]
    line = f.readline().split()
    T= np.zeros(nT)
    T[:] = line[:]

    #print(na, nT)
    #print(a[:])
    #print(T[:])

    kext = np.loadtxt(fname,skiprows=3)


    fig = plt.figure()

    kext_lev = np.linspace(np.min(kext),np.max(kext),10)
    con_kext = plt.contourf(T, a, kext, kext_lev)
    plt.colorbar(con_kext, label=r' $Q_{\rm ext}$ [cm$^{2}$ particle$^{-1}$]')
    plt.ylabel('median radius [um]')

    if name != 'wav':
        plt.xlabel('temperature [K]')
    else:
        plt.xlabel(r'wavelength [$\mu$m]')

    plt.yscale('log')


    plt.title('sigma: ' + str(sigma))

    plt.savefig(sp+'_qext.png',dpi=144,bbox_inches='tight')

    # read in a data
    fname = sp+'_wav_pi0.txt'

    f = open(fname, 'r')
    line = f.readline().split()
    na = int(line[0])
    #print(na, '\n')
    nT = int(line[1])
    #print(nT, '\n')
    sigma = float(line[2])
    #print(sigma, '\n')
    line = f.readline().split()
    #print(line, '\n')
    a = np.zeros(na)
    a[:] = line[:]
    #print(a)
    line = f.readline().split()
    T= np.zeros(nT)
    T[:] = line[:]

    #print(na, nT)
    #print(a[:])
    #print(T[:])

    w = np.loadtxt(fname,skiprows=3)

    fig = plt.figure()

    w_lev = np.linspace(np.min(w),np.max(w),10)
    con_w = plt.contourf(T, a, w, w_lev)
    plt.colorbar(con_w, label=r'$\omega$ [-]')
    plt.yscale('log')
    plt.ylabel('median radius [um]')

    if name != 'wav':
        plt.xlabel('temperature [K]')
    else:
        plt.xlabel(r'wavelength [$\mu$m]')

    plt.title('sigma: ' + str(sigma))

    plt.savefig(sp+'_a.png',dpi=144,bbox_inches='tight')

    # read in g data
    fname = sp+'_wav_gg.txt'

    f = open(fname, 'r')
    line = f.readline().split()
    na = int(line[0])
    nT = int(line[1])
    sigma = float(line[2])
    line = f.readline().split()
    a = np.zeros(na)
    a[:] = line[:]
    line = f.readline().split()
    T= np.zeros(nT)
    T[:] = line[:]

    #print(na, nT)
    #print(a[:])
    #print(T[:])

    g = np.loadtxt(fname,skiprows=3)

    fig = plt.figure()

    g_lev = np.linspace(np.min(g),np.max(g),10)

    con_g = plt.contourf(T, a, g, g_lev)
    plt.colorbar(con_g, label=r'$g$ [-]')
    plt.yscale('log')
    plt.ylabel('median radius [um]')

    if name != 'wav':
        plt.xlabel('temperature [K]')
    else:
        plt.xlabel(r'wavelength [$\mu$m]')


    plt.title('sigma: ' + str(sigma))

    plt.savefig(sp+'_g.png',dpi=144,bbox_inches='tight')
    print('plots made for ' + sp)
    plt.close()

#raise Exception(';ofsidj;sdf')
header = False
for spec in species:
    for typ in types:
        fname = spec + '_wav_' + typ + '.txt'
        print('editing ' + fname)
        # if typ == 'kext':
        #    typ = 'qext'

        f = open(fname, 'r')

        line = f.readline().split()
        # print(line)

        na = int(line[0])
        nT = int(line[1])
        sigma = float(line[2])
        line = f.readline().split()
        #print('\n', line, '\n')
        a = np.zeros(na)
        a[:] = line[:]
        line = f.readline().split()
        #print('THIRD PRINT HERE\n', line, '\n')
        T = np.zeros(nT)
        T[:] = line[:]

        # print(na, nT)
        # print(a[:])
        # print(T[:])
        wavelengths = np.zeros(nT)
        nwl = na

        # min and max wavelengths [um]
        wlmin = 0.1
        wlmax = 20

        # Set up wavelength grid (logspaced recommended)
        wl = np.logspace(np.log10(wlmin), np.log10(wlmax), nwl)

        # For GCM spectral models at certain bands, apply the wavelengths manually
        # e.g. below we have the 11 bands of Tiffany Kataria of MITgcm fame.
        # nwl = 12
        # wl = [0.260, 0.420, 0.610, 0.850, 1.320, 2.020, 2.500, 3.500, 4.400, 8.70, 20.00 ,324.68]

        # Output wavelengths to file


        value = np.loadtxt(fname, skiprows=3)

        # kext = np.log10(np.loadtxt(fname))

        if not header:
            rad = open('radius_array_for_cloud_scattering_data_in_microns.txt', 'w')
            wav = open('wavelength_array_for_cloud_scattering_data_in_microns.txt', 'w')
            for r in a:
                rad.write(str(r))
                rad.write('\n')
            rad.close()
            wav = open('wavelength_array_for_cloud_scattering_data_in_microns.txt', 'w')
            for i in range(nwl):
                wav.write(str(wl[i]) + '\n')
            wav.close()

            #shutil.copyfile('radius_array_for_cloud_scattering_data_in_microns.txt',
                            #'finished_data/radius_array_for_cloud_scattering_data_in_microns.txt')
            #shutil.copyfile('wavelength_array_for_cloud_scattering_data_in_microns.txt',
                           # 'finished_data/wavelength_array_for_cloud_scattering_data_in_microns.txt')
            header = True
        value_file = open(spec + '_' + name + '_' + typ + '.txt', 'w')
        for val in value:
            val = str(np.flip(val)).strip('[]').split()
            vals = []
            for v in val:
                v = float(v)
                vals.append('{:.3e}'.format(v))
            vals = '  '.join(vals)
            value_file.write(str(vals).strip('[]').replace("'", " "))
            value_file.write('\n')
        #shutil.copyfile(spec + '_wav_' + typ + '.txt', 'finished_data/' + spec + '_wav_' + typ + '.txt')

        f.close()
        # for q in qext:
        '''
        fig = plt.figure()

        kext_lev = np.linspace(np.min(kext),np.max(kext),10)

        con_kext = plt.contourf(T,a,kext,kext_lev)

        plt.colorbar(con_kext,label=r'$Qc_{\rm ext}$')
        plt.yscale('log')

        plt.ylabel('median radius [um]')
        plt.xlabel('temperature [K]')

        plt.title('sigma: ' + str(sigma))

        plt.savefig(sp+'_kext.png',dpi=144,bbox_inches='tight')
        '''

        f.close()

try:
    os.system('mkdir ' + name + '_png')
except:
    print('folder already made')
os.system('cp *.png ' + name + '_png')
os.system('rm *.png')

if name != 'wav':
    os.system('rm *_wav_*')
else:
    os.system('rm *lognorm*')
