import numpy as np
import struct
import string
import glob
import os
import re

def parse_filename(filename):
	"""
	#   PARSE_FILENAME    Break up a full-path filename into its component
	#   parts to check the extension, make it more readable, and extract the step
	#   number.  
	#
	#   PATH,BASENAME,STEP,EXT = PARSE_FILENAME(FILENAME)
	#
	#   E.g. If FILENAME='/home/Blast.0000.bin', then PATH='/home',
	#   BASENAME='Blast', STEP='0000', and EXT='bin'.
	#
	"""


	path=os.path.dirname(filename)
	if path[-3:] == 'id0': 
		path=path[:-3]
		mpi_mode=True
	else:
		path=path+os.path.sep
		mpi_mode=False

	base=os.path.basename(filename)
	id=base[:-9]
	step, ext=base[-8:].split('.')

	return path,id,step,ext,mpi_mode

def parse_line(line, grid):
	sp = line.strip().split()

	if "vtk" in sp:
		grid['vtk_version'] = sp[-1]
	elif "time=" in sp:
		time_index = sp.index("time=")
        	grid['time'] = float(sp[time_index+1].rstrip(','))
        	if 'level' in sp: grid['level'] = int(sp[time_index+3].rstrip(','))
        	if 'domain' in sp: grid['domain'] = int(sp[time_index+5].rstrip(','))  
		if sp[0] == "PRIMITIVE": 
			grid['prim_var_type']=True
	elif "DIMENSIONS" in sp:
		grid['Nx'] = np.array(sp[-3:]).astype('int')
	elif "ORIGIN" in sp:
		grid['left_edge'] = np.array(sp[-3:]).astype('float64')
	elif "SPACING" in sp:
		grid['dx'] = np.array(sp[-3:]).astype('float64')
	elif "CELL_DATA" in sp:
		grid['ncells'] = int(sp[-1])
	elif "SCALARS" in sp:
		grid['read_field'] = sp[1]
		grid['read_type'] = 'scalar'
	elif "VECTORS" in sp:
		grid['read_field'] = sp[1]
		grid['read_type'] = 'vector'
	elif "NSTARS" in sp:
		grid['nstar'] = eval(sp[1])
	elif "POINTS" in sp:
		grid['nstar'] = eval(sp[1])
		grid['ncells'] = eval(sp[1])
        elif "Athena++" in sp:
                grid['dx'] = np.array([1., 1., 1.])
                grid['left_edge'] = np.array([0., 0., 0])
                grid['time'] = 0.


class AthenaDomain(object):
	def __init__(self,filename,grids=None,setgrid=True, noread_mpi=False):
		self.flist = glob.glob(filename)
		if len(self.flist) == 0: 
			print 'no such file: %s' % filename
		dir, id, step, ext, mpi = parse_filename(filename)
                if noread_mpi:
                        mpi = False
		self.dir = dir
		self.id = id
		self.step = step
		self.ext = ext
		self.starfile = os.path.join(dir+'id0/','%s.%s.%s.%s' % (id,step,'starpar',ext))
		self.mpi = mpi
		if mpi: self.flist += glob.glob(os.path.join(dir,'id*/%s-id*.%s.%s' % (id, step, ext)))
		self.flist.sort()
		self.ngrids = len(self.flist)
		if setgrid:
			if grids==None: self.grids=self._setup_grid()
			else: 
				if grids[0]['filename'] != self.flist[0]:
					for g,f in zip(grids,self.flist): g['filename']=f
				self.grids=grids
			self.domain=self._setup_domain(self.grids)
			self._setup()
	
	def _setup(self):
		self.domain['data']={}

	def _setup_domain(self,grids):
		domain = {}
		left_edges = np.empty((self.ngrids,3), dtype='float32')
		dxs = np.empty((self.ngrids,3), dtype='float32')
		Nxs = np.ones_like(dxs)
		for nproc,g in enumerate(grids):
			left_edges[nproc,:] = g['left_edge']
			Nxs[nproc,:] = g['Nx']
			dxs[nproc,:] = g['dx']

		right_edges = left_edges + Nxs*dxs

                try:
		    left_edge = left_edges.min(0)
		    right_edge = right_edges.max(0)
                except ValueError:
                    left_edge = 0.
                    right_edge = 1.

		gis = np.round((left_edges - left_edge)/dxs)
		
		for nproc,g in enumerate(grids):
			g['is']=gis[nproc,:]

		domain['left_edge'] = left_edge
		domain['right_edge'] = right_edge
		domain['dx'] = dxs[0,:]
		domain['Lx'] = right_edge - left_edge
		domain['center'] = 0.5*(right_edge + left_edge)
		domain['Nx'] = np.round(domain['Lx']/domain['dx']).astype('int')
		domain['ndim'] = 3 # should be revised
		file = open(self.flist[0],'rb')
		tmpgrid = {}
		tmpgrid['time']=None
		while tmpgrid['time'] is None:
			line = file.readline()
			parse_line(line,tmpgrid)
		file.close()
		domain['time'] = tmpgrid['time']
		domain['data'] = {}
		domain['field_map']=None

		return domain


	def _setup_grid(self):
		grids=[]
		for nproc in range(self.ngrids):
			file = open(self.flist[nproc],'rb')
			grid = {}
			grid['filename']=self.flist[nproc]
			grid['read_field'] = None
			grid['read_type'] = None
			while grid['read_field'] is None:
				grid['data_offset']=file.tell()
				line = file.readline()
				parse_line(line, grid)
			file.close()
			grid['Nx'] -= 1
			grid['Nx'][grid['Nx'] == 0] = 1
			grid['dx'][grid['Nx'] == 1] = 1.
			grid['right_edge'] = grid['left_edge'] + grid['Nx']*grid['dx']
			#grid['field_map']=None

			grids.append(grid)

		return grids

class AthenaDataSet(AthenaDomain):
	def _setup(self):
		for i,g in enumerate(self.grids):
#			if g['field_map']==None: 
#				print "Setting %d-th grid" % (i)
#				self._set_field_map(g)
#			else: 
			g['data']={}
		if self.domain['field_map']==None:
			self.domain['field_map'] = self._set_field_map(self.grids[0])
		fm = self.domain['field_map']
		if 'cell_centered_B' in fm.keys():
			fm['magnetic_field']=fm['cell_centered_B']
			fm.pop('cell_centered_B')
		elif 'face_centered_B' in fm.keys():
			fm['magnetic_field']=fm['face_centered_B']
			fm.pop('face_centered_B')
		nscal=0
                if 'specific_scalar[0]' in fm.keys():
			keys=fm.keys()
			for k in keys:
				if k.startswith('specific_scalar'):
					newkey=re.sub("\s|\W","",k)
					fm[newkey] = fm.pop(k)
					nscal += 1
		self.field_list=fm.keys()

		
		derived_field_list=[]
		derived_field_list.append('number_density')
		if 'magnetic_field' in self.field_list:
			derived_field_list.append('magnetic_field1')
			derived_field_list.append('magnetic_field2')
			derived_field_list.append('magnetic_field3')
			derived_field_list.append('magnetic_energy1')
			derived_field_list.append('magnetic_energy2')
			derived_field_list.append('magnetic_energy3')
			derived_field_list.append('magnetic_pressure')
			derived_field_list.append('plasma_beta')
			derived_field_list.append('alfven_velocity1')
			derived_field_list.append('alfven_velocity2')
			derived_field_list.append('alfven_velocity3')
			derived_field_list.append('magnetic_stress')
			derived_field_list.append('magnetic_stress1')
			derived_field_list.append('magnetic_stress2')
			derived_field_list.append('magnetic_stress3')
		if 'velocity' in self.field_list:
			derived_field_list.append('velocity1')
			derived_field_list.append('velocity2')
			derived_field_list.append('velocity3')
			derived_field_list.append('kinetic_energy1')
			derived_field_list.append('kinetic_energy2')
			derived_field_list.append('kinetic_energy3')
			derived_field_list.append('momentum1')
			derived_field_list.append('momentum2')
			derived_field_list.append('momentum3')
			derived_field_list.append('reynold_stress')
			derived_field_list.append('reynold_stress1')
			derived_field_list.append('reynold_stress2')
			derived_field_list.append('reynold_stress3')
		if 'pressure' in self.field_list:
			derived_field_list.append('sound_speed')
			derived_field_list.append('temperature')
			derived_field_list.append('T1')
		if 'gravitational_potential' in self.field_list:
			derived_field_list.append('potential_energy')
			derived_field_list.append('gravity_stress')
			derived_field_list.append('gravity_stress1')
			derived_field_list.append('gravity_stress2')
			derived_field_list.append('gravity_stress3')
		if nscal > 0:
			for n in range(nscal):
				derived_field_list.append('scalar%d' % n)
		self.domain['nscal']=nscal
		self.derived_field_list=derived_field_list

	def _set_field_map(self,grid):
		file=open(grid['filename'],'rb')
		file.seek(0,2)
		eof = file.tell()
		offset = grid['data_offset']
		file.seek(offset)

		field_map={}

		if grid.has_key('Nx'): Nx=grid['Nx']

		while offset < eof:

			line=file.readline()
			sp = line.strip().split()
			#print line,sp
				
			field_map[sp[1]] = {}
			field_map[sp[1]]['read_table']=False

			if "SCALARS" in line:
				tmp=file.readline()
				field_map[sp[1]]['read_table']=True
				field_map[sp[1]]['nvar'] = 1
			elif "VECTORS" in line:
				field_map[sp[1]]['nvar'] = 3
			else:
				print 'Error: '+sp[0] + ' is unknown type'
				raise TypeError

			field_map[sp[1]]['offset']=offset
			field_map[sp[1]]['ndata']=field_map[sp[1]]['nvar']*grid['ncells']
			if sp[1] == 'face_centered_B1':
				field_map[sp[1]]['ndata']=(Nx[0]+1)*Nx[1]*Nx[2]
			elif sp[1] == 'face_centered_B2':
				field_map[sp[1]]['ndata']=Nx[0]*(Nx[1]+1)*Nx[2]
			elif sp[1] == 'face_centered_B3':
				field_map[sp[1]]['ndata']=Nx[0]*Nx[1]*(Nx[2]+1)
				
			if sp[2]=='int': dtype='i'
			elif sp[2]=='float': dtype='f'
			elif sp[2]=='double': dtype='d'
			field_map[sp[1]]['dtype']=dtype
			field_map[sp[1]]['dsize']=field_map[sp[1]]['ndata']*struct.calcsize(dtype)
			file.seek(field_map[sp[1]]['dsize'],1)
			offset = file.tell()
			tmp=file.readline()
			if len(tmp)>1: file.seek(offset)
                        else: offset = file.tell()

		#grid['field_map'] = field_map
		#grid['data']={}
		return field_map


	def _read_field(self,file_pointer,field_map):
		ndata=field_map['ndata']
		dtype=field_map['dtype']
		file_pointer.seek(field_map['offset'])
		file_pointer.readline() # HEADER
		if field_map['read_table']: file_pointer.readline()
		data = file_pointer.read(field_map['dsize'])
		var = np.asarray(struct.unpack('>'+ndata*dtype,data))

		return var

	def _read_grid_data(self,grid,field):
		gd=grid['data']
		if gd.has_key(field):
			return

		file=open(grid['filename'],'rb')
		#fm=grid['field_map']
		fm=self.domain['field_map']
		nx1=grid['Nx'][0]
		nx2=grid['Nx'][1]
		nx3=grid['Nx'][2]

		if field == 'face_centered_B1': nx1=nx1+1
		if field == 'face_centered_B2': nx2=nx2+1
		if field == 'face_centered_B3': nx3=nx3+1

		nvar=fm[field]['nvar']
		var = self._read_field(file,fm[field])
		if nvar == 1: 
			var.shape = (nx3, nx2, nx1)
		else: 
			var.shape = (nx3, nx2, nx1, nvar)
		file.close()
		grid['data'][field]=var
		if nvar == 3: self._set_vector_field(grid,field)


	def _get_grid_data(self,grid,field):
		gd=grid['data']
		if gd.has_key(field):
			return gd[field]
		elif field in self.field_list:
			self._read_grid_data(grid,field)
			return gd[field]
		elif field in self.derived_field_list:
			data=self._get_derived_field(grid,field)
			return data
		else:
			print field,' is not supported'

	def _set_vector_field(self,grid,vfield):
		gd=grid['data']
		gd[vfield+'1'] = gd[vfield][:,:,:,0]
		gd[vfield+'2'] = gd[vfield][:,:,:,1]
		gd[vfield+'3'] = gd[vfield][:,:,:,2]

	def _get_derived_field(self,grid,field):
		import astropy.constants as c
		gd=grid['data']
		if gd.has_key(field):
			return gd[field]
		elif field.startswith('velocity'):
			self._read_grid_data(grid,'velocity')
			return gd[field]
		elif field.startswith('magnetic_field'):
			self._read_grid_data(grid,'magnetic_field')
			return gd[field]
		elif field.startswith('number_density'):
			self._read_grid_data(grid,'density')
			return gd['density']
		elif field.startswith('kinetic_energy'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'velocity')
			den=gd['density']
			v1=gd['velocity1']
			v2=gd['velocity2']
			v3=gd['velocity3']
			if field is 'kinetic_energy1': return 0.5*den*v1**2
			if field is 'kinetic_energy2': return 0.5*den*v2**2
			if field is 'kinetic_energy3': return 0.5*den*v3**2
		elif field.startswith('momentum'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'velocity')
			den=gd['density']
			v1=gd['velocity1']
			v2=gd['velocity2']
			v3=gd['velocity3']
			if field is 'momentum1': return den*v1
			if field is 'momentum2': return den*v2
			if field is 'momentum3': return den*v3
		elif field.startswith('magnetic_energy'):
			self._read_grid_data(grid,'magnetic_field')
			B1=gd['magnetic_field1']
			B2=gd['magnetic_field2']
			B3=gd['magnetic_field3']
			if field is 'magnetic_energy1': return B1**2/(8*np.pi)
			if field is 'magnetic_energy2': return B2**2/(8*np.pi)
			if field is 'magnetic_energy3': return B3**2/(8*np.pi)
		elif field.startswith('magnetic_pressure'):
			self._read_grid_data(grid,'magnetic_field')
			B1=gd['magnetic_field1']
			B2=gd['magnetic_field2']
			B3=gd['magnetic_field3']
			if field is 'magnetic_pressure': return (B1**2+B2**2+B3**2)/(8*np.pi)
		elif field.startswith('plasma_beta'):
			vfield='magnetic_field'
			self._read_grid_data(grid,'pressure')
			self._read_grid_data(grid,vfield)
			B1=gd[vfield+'1']
			B2=gd[vfield+'2']
			B3=gd[vfield+'3']
			press=gd['pressure']
			if field is 'plasma_beta': return press*(8.0*np.pi)/(B1**2+B2**2+B3**2)
		elif field.startswith('alfven_velocity'):
			vfield='magnetic_field'
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,vfield)
			den=gd['density']
			B1=gd[vfield+'1']
			B2=gd[vfield+'2']
			B3=gd[vfield+'3']
			if field is 'alfven_velocity1': return np.sqrt(B1**2/(4*np.pi*den))
			if field is 'alfven_velocity2': return np.sqrt(B2**2/(4*np.pi*den))
			if field is 'alfven_velocity3': return np.sqrt(B3**2/(4*np.pi*den))
		elif field.startswith('sound_speed'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'pressure')
			den=gd['density']
			press=gd['pressure']
			return np.sqrt(press/den)
		elif field.startswith('temperature'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'pressure')
			den=gd['density']*u['density']
			press=gd['pressure']*u['pressure']
			T1=(press/den*c.m_p/c.k_B).cgs
			return T1*1.4/1.1
		elif field.startswith('T1'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'pressure')
			den=gd['density']*u['density']
			press=gd['pressure']*u['pressure']
			return (press/den*c.m_p/c.k_B).cgs
		elif field.startswith('potential'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'gravitational_potential')
			den=gd['density']
			pot=gd['gravitational_potential']
			return -den*pot
		elif field.startswith('magnetic_stress'):
			vfield='magnetic_field'
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,vfield)
			B1=gd[vfield+'1']
			B2=gd[vfield+'2']
			B3=gd[vfield+'3']
			if field is 'magnetic_stress1': return B2*B3/(4*np.pi) 
			if field is 'magnetic_stress2': return B1*B3/(4*np.pi) 
			if field is 'magnetic_stress3': return B1*B2/(4*np.pi) 
			return B1*B2/(4*np.pi)
		elif field.startswith('reynold_stress'):
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'velocity')
			den=gd['density']
			v1=gd['velocity1']
			v2=gd['velocity2']
			v3=gd['velocity3']
			if field is 'reynold_stress1': return den*v2*v3
			if field is 'reynold_stress2': return den*v1*v3
			if field is 'reynold_stress3': return den*v1*v2
			return den*v1*v2
		elif field.startswith('gravity_stress'):
			self._read_grid_data(grid,'gravitational_potential')
			phi=gd['gravitational_potential']
			dx=grid['dx']
			g1,g2,g3=gradient(phi,dx)
			if field is 'gravity_stress1': return g2*g3/4/np.pi
			if field is 'gravity_stress2': return g1*g3/4/np.pi
			if field is 'gravity_stress3': return g1*g2/4/np.pi
			return  g1*g2/4/np.pi
		elif field.startswith('scalar'):
			scal = field[6:]
			self._read_grid_data(grid,'density')
			self._read_grid_data(grid,'specific_scalar'+scal)
			den=gd['density']
			sscal=gd['specific_scalar'+scal]
			return sscal*den
		
	def read_all_data(self,field):
		#fm=self.grids[0]['field_map']
		fm=self.domain['field_map']
			
		dnx=self.domain['Nx']
		if field in self.field_list:
		  if fm[field]['nvar']==3:
			data=np.empty((dnx[2],dnx[1],dnx[0],3),dtype=fm[field]['dtype'])
		  else:
			data=np.empty((dnx[2],dnx[1],dnx[0]),dtype=fm[field]['dtype'])
		  if field is 'face_centered_B1':
			data=np.empty((dnx[2],dnx[1],dnx[0]+1),dtype=fm[field]['dtype'])
		  if field is 'face_centered_B2':
			data=np.empty((dnx[2],dnx[1]+1,dnx[0]),dtype=fm[field]['dtype'])
		  if field is 'face_centered_B3':
			data=np.empty((dnx[2]+1,dnx[1],dnx[0]),dtype=fm[field]['dtype'])
                elif field in self.derived_field_list:
		  data=np.empty((dnx[2],dnx[1],dnx[0]),dtype=fm['density']['dtype'])

		for g in self.grids:
			gis=g['is']
			gnx=g['Nx']
			gie=gis+gnx
			gd=self._get_grid_data(g,field)
			if field in self.field_list and fm[field]['nvar']==3:
				data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0],:]=gd
			else:
				if gie[0] == dnx[0] and field is 'face_centered_B1':
				    data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]+1]=gd
				elif gie[1] == dnx[1] and field is 'face_centered_B2':
				    data[gis[2]:gie[2],gis[1]:gie[1]+1,gis[0]:gie[0]]=gd
				elif gie[2] == dnx[2] and field is 'face_centered_B3':
				    data[gis[2]:gie[2]+1,gis[1]:gie[1],gis[0]:gie[0]]=gd
				else:
				    gd=gd[0:gnx[2],0:gnx[1],0:gnx[0]]
				    data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=gd
			
		return data

	def read_starvtk(self):
		if not os.path.isfile(self.starfile): return {}
		file=open(self.starfile,'rb')
		star = {}
		star['filename']=self.starfile
		star['read_field'] = None
		star['read_type'] = None
		while star['read_field'] is None:
			star['data_offset']=file.tell()
			line = file.readline()
			parse_line(line, star)

		nstar=star['nstar']
		#print nstar
		fm=self._set_field_map(star)
		id=self._read_field(file,fm['star_particle_id'])
		mass=self._read_field(file,fm['star_particle_mass'])
		age=self._read_field(file,fm['star_particle_age'])
		pos=self._read_field(file,fm['star_particle_position']).reshape(nstar,3)
		vel=self._read_field(file,fm['star_particle_velocity']).reshape(nstar,3)
		file.close()
		star=[]
		for i in range(nstar):
			star.append({})

		for i in range(nstar):
			star_dict = star[i]
			star_dict['id']=id[i]
			star_dict['mass']=mass[i]
			star_dict['age']=age[i]
			star_dict['v1']=vel[i][0]
			star_dict['v2']=vel[i][1]
			star_dict['v3']=vel[i][2]
			star_dict['x1']=pos[i][0]
			star_dict['x2']=pos[i][1]
			star_dict['x3']=pos[i][2]

		return star


