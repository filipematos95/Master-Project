'''
Script to perfrom cross-match between two catalogs. The input should be catalog1 catalog2 radius radius_unit [flag0] [flag1]

'''

import numpy as np 
import pandas as pd 
from astropy import units as u 
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from astropy.io import fits
import warnings
import sys
try:
 from tabulate import tabulate
 tabulate_module = True
except ImportError, e:
 warnings.warn("The module tabulate is not installed. The code will run whitout it but the output might be limited")
 tabulate_module = False


def catalog_read(catalog_name):
	if catalog_name[-4:] == 'fits':
		#fits Catalog
		hdulist = fits.open('./Catalogs/'+catalog_name)
		tbdata = hdulist[1].data
		name = tbdata.field('NAME')
		ra = tbdata.field('RA') #degree
		dec = tbdata.field('DEC') #degree
		#if tbdata.field('JMAG') in tbdata:
		jmag = tbdata.field('JMAG')
		return name,ra,dec,jmag
	elif catalog_name[-3:] == 'csv':
		#csv catalog
		data = pd.read_csv('./Catalogs/'+catalog_name)

		name = data['NAME'].values
		ra = data['RA'].values #Hour
		dec = data['DEC'].values #Degree
		if 'MAG' in data.columns:
			mag = data['MAG'].values #Days
			if 'PER' in data.columns:
				per = data['PER'].values
				return name,ra,dec,mag,per
			return name,ra,dec,mag
		if 'PER' in data.columns:
			per = data['PER'].values
			if 'MAG' in data.columns:
				mag = data['MAG'].values
				return name,ra,dec,per,mag
			return name,ra,dec,per
	else:
		print 'Unrecognized format, please use csv or fits format'
		print catalog_name[-3:]
		name, ra,dec = ['nan',0,0]
	return name,ra,dec

def get_catalog_name(catalog_name):
 
	if catalog_name[-4:] == 'fits':
		#fits Catalog
		catalog_name = catalog_name[:-5]
	elif catalog_name[-3:] == 'csv':
		#csv catalog
		catalog_name = catalog_name[:-4]
	else:
		print 'Unrecognized format, please use csv or fits format'
		print catalog_name[-3:]
		catalog_name = 'UNREC'
	return catalog_name

def Catalogs(catalog):
    return {
        'HST': 'hst.csv',
        'HST_EXOPLANETS' : 'hst_exoplanets_25.csv',
        'FUSE': 'fuse.csv',
        'FUSE_EXOPLANETS': 'fuse_exoplanets_magper_25.csv',
        'EXOPLANETS_PER': 'exoplanets_magper.csv',
        'EXOPLANETS': 'exoplanets.csv',
        'MDWARF': 'mdwarf.fits',
        'MDWARF_MAG':'mdwarf_mag.fits',
        'TESS' : 'tess_cvz.csv',
        'LOFAR': 'lofar.csv',
        'TESS_HST': 'tess_cvz_hst_90.csv',
        'LOFAR_HST': 'lofar_hst_5.csv',
        'CHANDRA': 'chandra.csv',
        'CHANDRA_EXOPLANETS_21': 'chandra_exoplanets_21.csv',
        'CHANDRA_EXOPLANETS': 'chandra_exoplanets_25.csv',
        'LOFAR_EXOPLANETS_5': 'lofar_exoplanets_5.csv',
        'LOFAR_EXOPLANETS_2_5': 'lofar_exoplanets_2_5.csv',
        'XMM': 'xmm2.csv',
        'XMM_EXOPLANETS': 'xmm2_exoplanets_magper_25.csv',
        'TESS_FIELD1': 'tess_field1.csv'
    }.get(catalog, catalog)

#User input
flags = ['NONE','NONE'] 

if len(sys.argv) >=  5:
	catalog1  = Catalogs(sys.argv[1])
	catalog2 = Catalogs(sys.argv[2])
	radius = sys.argv[3]
	radius_unit = sys.argv[4]
if len(sys.argv) >=  6:
	flags[0] = sys.argv[5]
if len(sys.argv) ==  7:
	flags[1] = sys.argv[6]

else: 
	"Please put all the necessary arguments"

if radius_unit == 'degree':
	radius_unit = u.degree
elif radius_unit == 'arcminute':
	radius_unit = u.arcminute
elif radius_unit == 'arcsecond':
	radius_unit = u.arcsecond
else:
	print "Please use 'degree' or 'arcsecond' for the radius unit"


#Catalog reading 

name1,ra1,dec1 = catalog_read(catalog1)
print(type(name1[0]))
#name1 = [str(x)[:-2] for x in name1] #We are only interested in the host star (otherwise comment it)
name1,unique_index = np.unique(name1,return_index = True)
ra1 = ra1[unique_index]
dec1 = dec1[unique_index]

if flags[0] == 'MAG':
	name2,ra2,dec2,mag2 = catalog_read(catalog2)
	if flags[1] == 'PER':
		name2,ra2,dec2,mag2,per2 = catalog_read(catalog2)
elif flags[0] == 'PER':
	name2,ra2,dec2,per2 = catalog_read(catalog2)
	if flags[1] == 'MAG':
		name2,ra2,dec2,per2,mag2 = catalog_read(catalog2)
else: 
	name2,ra2,dec2 = catalog_read(catalog2)

name2 = [str(x)[:-1] for x in name2]
name2,unique_index = np.unique(name2,return_index = True)
ra2 = ra2[unique_index]
dec2 = dec2[unique_index]

#Cross Match
# Unit Selection
unit1 = u.hour
unit2 = u.hour

if np.max(ra1) > 25.0:
	unit1 = u.degree
if np.max(ra2) > 25.0:
	unit2 = u.degree 

catalog1_coord = SkyCoord(ra=ra1*unit1,dec=dec1*u.degree)
catalog2_coord = SkyCoord(ra=ra2*unit2, dec=dec2*u.degree)

catalog1 = get_catalog_name(catalog1)
catalog2 = get_catalog_name(catalog2)

file_name = catalog1 + '_' + catalog2 + '_' + str(radius)
file_name2_ = file_name.replace('.', '_') + '.csv'
file_name_ = file_name.replace('.', '_') + '.txt'

if not os.path.exists('./Results/'):
    os.makedirs('./Results/')
if not os.path.exists('./Results/Table/'):
    os.makedirs('./Results/')
if not os.path.exists('./Results/LATEX/'):
    os.makedirs('./Results/')

file_name = './Results/' + file_name_
file_name2 = './Results/' + file_name2_
file_namet = './Results/Table/' + file_name_
file_name3 = './Results/MulCross/' + file_name2_
file_name_latex  = "./Results/LATEX/" + file_name_

f = open(file_name,'w')
f3 = open(file_name3,'w')
fl = open(file_name_latex,'w')

f3.write('NAME1' + ',' + 'NAME2' + '\n')
print "Simple cross-match for " + catalog1  + ' and ' + catalog2
print " " 

print "Simple cross-match with a radius of " + str(radius)

f.write("Simple cross-match for " + catalog1  + ' and ' + catalog2 + '\n')
f.write('\n' + "Simple cross-match with a radius of " + str(radius) + ' '  + '\n')

f2 = open(file_name2,'w')

if flags[0] == 'MAG':
	f.write('NAME1 -> NAME2 -> MAG\n')
	f2.write('NAME,RA,DEC,MAG \n')
	if flags[1] == 'PER':
		f.write('NAME -> RA -> DEC -> MAG -> PER\n')
		f2.write('NAME,RA,DEC,MAG,PER \n')
elif flags[0] == 'PER':
	f.write('NAME1 -> NAME2 -> PER\n')
	f2.write('NAME,RA,DEC,PER\n')
	if flags[1] == 'MAG':
		f.write('NAME -> RA -> DEC -> PER -> MAG\n')
		f2.write('NAME,RA,DEC,PER,MAG\n')
else: 
	f.write('NAME -> RA -> DEC\n')
	f2.write('NAME,RA,DEC\n')

idxc, idxcatalog, d2d, d3d = catalog2_coord.search_around_sky(catalog1_coord, float(radius)*radius_unit)
rows = []
for i in range(0,len(idxcatalog)):
	if ((name2[idxcatalog[i]] != name2[idxcatalog[i-1]])): #and (name1[idxc[i]] != name1[idxc[i-1]]):
		print name1[idxc[i]], ' -> ', name2[idxcatalog[i]]
		if flags[0] == 'MAG':
			rows.append([str(name1[idxc[i]]),str(name2[idxcatalog[i]]),str(mag2[idxcatalog[i]])])
			f.write(str(name1[idxc[i]]) + ' -> '+ str(name2[idxcatalog[i]]) + ' -> ' + str(mag2[idxcatalog[i]]) +'\n')
			if flags[1] == 'PER':
				f.write(str(name1[idxc[i]]) + ' -> '+ str(name2[idxcatalog[i]]) + ' -> ' + str(mag2[idxcatalog[i]]) + ' -> ' +  str(per2[idxcatalog[i]])  +'\n')
				rows.append([str(name1[idxc[i]]),str(name2[idxcatalog[i]]),str(mag2[idxcatalog[i]]), str(per2[idxcatalog[i]])])
		if flags[0] == 'PER':
			f.write(str(name1[idxc[i]]) + ' -> '+ str(name2[idxcatalog[i]]) + ' -> ' + str(per2[idxcatalog[i]]) +'\n')
			rows.append([str(name1[idxc[i]]),str(name2[idxcatalog[i]]),str(per2[idxcatalog[i]])])
			fl.write(str(name1[idxc[i]]) + '&' + str(name2[idxcatalog[i]]) + '&' + str(per2[idxcatalog[i]]) + '\\\\ \n')
			if flags[1] == 'MAG':
				f.write(str(name1[idxc[i]]) + ' -> '+ str(name2[idxcatalog[i]]) + ' -> ' + str(per2[idxcatalog[i]]) + ' -> ' +  str(mag2[idxcatalog[i]])  +'\n')		
				rows.append([str(name1[idxc[i]]),str(name2[idxcatalog[i]]),str(per2[idxcatalog[i]]),str(mag2[idxcatalog[i]])])
		else:
			f.write(str(name1[idxc[i]]) + ' -> '+ str(name2[idxcatalog[i]]) +'\n')		
			rows.append([str(name1[idxc[i]]),str(name2[idxcatalog[i]])])
		#f.write(str(name1[idxc[i]]) + '&'+ name2[idxcatalog[i]] + '&' + str(per2[idxcatalog[i]]) +'\n')
		#f.write(str(name1[idxc[i]]) + '&'+ name2[idxcatalog[i]] + '&' + '\n')
			f2.write(str(name2[idxcatalog[i]]) + ',' + str(ra2[idxcatalog[i]]) + ',' + str(dec2[idxcatalog[i]]) + '\n')
			f3.write(str(name1[idxc[i]]) + ',' + str(name2[idxcatalog[i]]) + '\n')
			fl.write(str(name1[idxc[i]]) + '&' + str(name2[idxcatalog[i]]) + '\\\\ \n')


if tabulate_module:
	ft = open(file_namet,'w')
	if flags[0] == 'MAG':
		ft.write(tabulate(rows,['NAME1','NAME2','MAG']))
		if flags[1] == 'PER':
			ft.write(tabulate(rows,['NAME1','NAME2','MAG','PER']))
	elif flags[0] == 'PER':
		ft.write(tabulate(rows,['NAME1','NAME2','PER']))
		if flags[1] == 'MAG':	
			ft.write(tabulate(rows,['NAME1','NAME2','PER','MAG']))
	else: 
		ft.write(tabulate(rows,['NAME1','NAME2']))
	ft.close()


f.close()
f2.close()
f3.close()
fl.close()

name2,name2_index = np.unique(name2[idxcatalog],return_index = True)
#per2 = per2[idxcatalog]
#per2 = per2[name2_index]
print(" \n \n LaTeX table output: \n \n ")
for i in range(len(name2)):
	print(str(name2[i]) + ' \\\\ ')

print(' \n \n There are ' + str(len(name2)) + ' exoplanets in both catalogs')
