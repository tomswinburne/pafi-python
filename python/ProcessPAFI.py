import glob,os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import simps,cumtrapz
from scipy.interpolate import interp1d

kB = 8.617e-5


class ProcessPAFI:
	def __init__(self, csv_files, validity_fn = lambda x : x['MaxJump']<0.5) :
		"""
			Take a PAFI csv data dump and produce ensemble average and standard deviation
			`csv_files` : file name or list of file names
			`data` : pandas DataFrame in PAFI format
			`validity_criteria` : a function acting on a Dataframe which returns valid data
		"""
		self.validity_fn = validity_fn

		csv_data = self.load_csv_data(csv_files)
		if len(csv_data)==0:
			print("No csv files could  be found!")
			return

		self.csv_data = pd.concat(csv_data)

		self.establish()


	def establish(self):
		"""
			Initial processing of csv data
		"""
		self.nWorkers = max(set(self.csv_data['WorkerID'])) + 1

		self.validate()

		self.prepare()

		self.ensemble_data = self.ensemble_average()

	def load_csv_data(self,csv_files):
		"""
			Return list of DataFrames
		"""
		existant_csv_files = []
		if isinstance(csv_files,str):
			csv_files = [csv_files]

		for csv_file in csv_files:

			if os.path.exists(csv_file):
				existant_csv_files += [csv_file]
			else:
				print("csv file {} could not be found!".format(csv_file))

		csv_data = []
		for csv_file in existant_csv_files:
			csv_data += [pd.read_csv(csv_file)]

		return csv_data

	def validate(self):
		"""
			Apply validity criteria to all data
		"""
		self.valid = np.ones(self.csv_data.shape[0],bool)
		self.nValid = self.valid.sum()
		if not self.validity_fn is None:
			 self.valid = self.validity_fn(self.csv_data)

	def prepare(self):
		"""
			Separate csv data into parameters (i.e. metadata) and fields (i.e. data)
		"""
		data = self.csv_data[self.valid]


		# SYNTAX: every field before and including "DevFile" is a parameter == metadata
		self.nParams = np.where(data.columns=='DevFile')[0][0]
		self.param_names = list(data.columns[:self.nParams])
		self.field_names = list(data.columns[self.nParams+1:])

		# SYNTAX: every field after is simulation data
		self.data = data.drop('DevFile', axis=1).to_numpy()

		# parameters
		self.param_data = pd.DataFrame(data = self.data[:,:self.nParams],
									   columns = self.param_names).drop_duplicates()

	def ensemble_average(self):
		"""
			Perform an ensemble average of all simulation data
			returns:
				ensemble_data
				DataFrame with a columns for each parameter,
				each field (ensemble average) and field_std (ensemble std dev)
		"""
		# ensemble data
		ensemble_columns = self.param_names + self.field_names + [f+"_std" for f in self.field_names]
		ensemble_data = []
		for p in self.param_data.to_numpy():
			# take mean and standard deviation across ensemble
			sel = (self.data[:,:self.nParams] == p).min(axis=1)
			ensemble_data += [list(self.data[sel].mean(axis=0)) + list(self.data[sel].std(axis=0)[self.nParams:])]
		return pd.DataFrame(data = np.r_[ensemble_data], columns = ensemble_columns)


	def remesh(self,data,coord,density = 10):
		"""
			Splined rediscretization
			data: numpy array (n_data,n_fields)
			coordinate: numpy array (n_data)
			first column is coordinate, currently assumed uniform
			density : integer, optional. Default = 10
				Remeshing density
		"""
		assert data.shape[0] == coord.size

		spl_data = np.zeros((density*data.shape[0],data.shape[1]))

		sp_c = interp1d(np.linspace(0,1,coord.size),coord,kind="linear")
		spl_coord = sp_c(np.linspace(0,1,spl_data.shape[0]))

		for ii in range(data.shape[1]):
			spl_data[:,ii] = interp1d(coord, data[:,ii],kind='cubic')(spl_coord)

		return spl_coord,spl_data


	def ensemble_integrate(self,x='ReactionCoordinate',y='aveF',coord=None,remesh=10):
		"""
			Integrate an ensemble averaged field
			x: str field for integration
			coord: remapping of x
			y: data to be integrated


		"""
		param_data = self.param_data.drop(x,axis=1).drop_duplicates()
		param_fields = list(param_data.columns) + [x,y+"_int",y+'_int_std']
		csv_data = []

		for p in param_data.iterrows():
			try:
				ee = self.ensemble_data.copy()
			except:
				# TDOD: die elegantly, raise error....
				print("Cannot find ensemble_data!")

			for key,val in p[1].iteritems():
				ee = ee[ee[key]==val]
			data = np.r_[[ee[w].to_numpy() for w in [y,y+"_std"]]].T

			r_coord = ee[x].to_numpy()

			if (not coord is None) and (coord.size==r_coord.size):
					x_coord = coord.copy()
			else:
				x_coord = r_coord.copy()

			r_coord,_ = self.remesh(data[r_coord.argsort(),:],r_coord,remesh)
			x_coord,data = self.remesh(data[x_coord.argsort(),:],x_coord,remesh)

			chain_scaling = np.gradient(r_coord,x_coord)


			data[:,0] *= chain_scaling
			data[:,1] *= chain_scaling
			y_a = np.append(np.zeros(1),cumtrapz(data[:,0],x_coord))
			y_e = np.append(np.zeros(1),cumtrapz(data[:,1],x_coord))
			csv_data += [[*list(p[1].to_dict().values()),x_coord,y_a,y_e]]


		return pd.DataFrame(data=csv_data,columns=param_fields)

	def add_csv_data(self,csv_files):
		"""
		Add an additional csv data file and reestablish
		"""
		csv_data = self.load_csv_data(csv_files)
		if len(csv_data)==0:
			print("No extra csv files could  be found!")
			return

		new_csv_data = pd.concat(csv_data)

		self.csv_data = pd.concat([self.csv_data,new_csv_data])


		self.establish()
