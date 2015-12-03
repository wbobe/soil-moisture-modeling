# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:06:09 2015
@author: Evan Coopersmith

Genetic algorithm implementation for soil moisture model calibration.  The
approach is published in HESS (Coopersmith et al, 2014)
"""

"""This is a set of scripts containing functions required for genetic algorithm implementation for parameter fitting"""
import random
from scipy.stats.stats import pearsonr
import SoilMoisture_CalVal as SM
import numpy as np

"""Values for genetic algorithms that effect the range of values examined and the requirements for convergence."""
_initial_pop = 20
_generations = 1000 #This can be made shorter for less important trials to save computational time
_param_ranges_E = [(0.05, 0.05, 365), (0.00001, 0.00001, 0)]
_param_ranges_D = [(80,50,10), (0,0,0)] ##Remember, CRN SM is scaled 0-1, SCAN is 0-100.
_gen_type = "LatinHypercube" #LatinHypercube or MonteCarlo
_rep_type = "Two-Parent" #Two-Parent or Multi-Parent
_r_gen = 4
_num_mates = 12
_mut_prob = 0.5
_mut_mag = 0.03
_death_rate = 0.2
_c2constr = True #Ensures the sinusoids do not fall below zero.
_GA_Eta_Termination = 1#12 #How many generations shall we spend without improving our fits of Eta without terminating?
_GA_Diag_Termination = 1#50 #How many generations shall we spend without improving our fits for the DSME without terminating?

def GetValidSubset(inv_seq, grow_ind, time_series):
    """Given that there are time stamps that are either not within the defined growing season or invalid in
    terms of available data, this function returns a subset of the (time_series) where both criteria are met.
    (inv_seq) contains a vector that is TRUE if the data are invalid and grow_ind is TRUE if we are within 
    the growing season, as defined by parameters above."""      
    return [t for i,g,t in zip(inv_seq, grow_ind, time_series) if (g and not i)]
    
def PrintPopulationAndFitnesses(population, fitnesses):
	for p,f in zip(population, fitnesses):
		print(p,f)    
	return None
 
def CalculateFitnesses(proc, population, t_stamps, z, n, p_series, inseason_sm, fitnesses, b_series, inv_seq, grow_ind):
	if proc == "FittingEta":
		return GetFitnessE(population, t_stamps, z, grow_ind, n, p_series, inseason_sm, fitnesses, inv_seq) 
	elif proc == "FittingDiag":
		return GetFitnessD(population, b_series, n, grow_ind, inseason_sm, fitnesses, inv_seq)	

def UpdateBestCreature(fitnesses, population, g, BEST):
	"""Store the "fittest" creature found thus far."""
	TopFit = max(fitnesses)
	if(TopFit > BEST[len(BEST)-2]):
		best_ind = fitnesses.index(TopFit)
		BEST = population[best_ind] + [TopFit] + [g + 1]
		print("A new 'best' creature has emerged: ", BEST)    
	return BEST  

def FullGA(depth, memory, grow_start, grow_end, inv_seq, proc, sm_seq, t_stamps, p_series, b_series):
	"""Runs the full genetic algorithms, with a few parametric options that are stored globally
	This, more than any other function, could be enhanced by parallel computing"""
	p_ranges = _param_ranges_E if proc == 'FittingEta' else _param_ranges_D 
     #Initialize the population
	BEST = [-999999999] * (len(p_ranges[0]) + 1)
	population, fitnesses = InitializePopulation(_initial_pop, p_ranges, _gen_type, _c2constr)
	BEST.append(-999999999) #Best fitness begins as -9999, last value, generation count, will climb with each iteration
	print("Initialized the population, which contains", len(population), "creatures and", len(fitnesses), "fitness scores")
	grow_ind = [(True if (t >= grow_start and t <= grow_end) else False) for t in t_stamps]
	
	inseason_sm = GetValidSubset(inv_seq, grow_ind, sm_seq)
	print("Length of the remaining valid sequence:", (len(inseason_sm)))
	if len(inseason_sm) == 0: return -999999, -999999, -999999

	#Calculate fitness for each member of that initial population.
	print("Calculating fitnesses for the initial population...")
	fitnesses = CalculateFitnesses(proc, population, t_stamps, depth, memory, p_series, inseason_sm, fitnesses, b_series, inv_seq, grow_ind)
	print("Length of fitnesses:", len(fitnesses), "best living creature:", max(fitnesses) + 1)
	print("Initial Population:"); PrintPopulationAndFitnesses(population, fitnesses)

	for g in range(_generations):
		print("Generation #:", (g + 1)) #Calculate fitness after mutation - this does not need to happen in the very first generation
		if g > 0: 
			fitnesses = CalculateFitnesses(proc, population, t_stamps, depth, memory, p_series, inseason_sm, fitnesses, b_series, inv_seq, grow_ind)
			print("Population, generation:", g); PrintPopulationAndFitnesses(population, fitnesses)
		BEST = UpdateBestCreature(fitnesses, population, g, BEST) #Store the "fittest" creature found thus far.
		
		#If we have not improved in a pre-selected number of generations, immediately return...
		if proc == "FittingEta":
			if BEST[-1] < g - _GA_Eta_Termination:
				print("No improvement since generation", BEST[-1], "terminating..."); return BEST[0], BEST[1], BEST[2]
		elif proc == "FittingDiag":
			if BEST[-1] < g - _GA_Diag_Termination:
				print("No improvement since generation", BEST[-1], "terminating..."); return BEST[0], BEST[1], BEST[2]
				
		#Perform the selection process, in this case, fitness proportional
		mating_population = Selection(population, fitnesses, _num_mates)
		#Create new creatures using MatingPopulation
		new_creatures = Reproduction(mating_population, _rep_type)
		
		#Generate a couple random creatures to avoid sticking at a local max
		if _r_gen > 0:
			random_new, new_fits = InitializePopulation(_r_gen, p_ranges, _gen_type, _c2constr)
			new_creatures += random_new; new_fits = [-999999999] * len(new_creatures)
			
		#Calculate the fitnesses of these new creatures
		new_fits = CalculateFitnesses(proc, new_creatures, t_stamps, depth, memory, p_series, inseason_sm, new_fits, b_series, inv_seq, grow_ind)		
		population = population + new_creatures; fitnesses = fitnesses + new_fits
		BEST = UpdateBestCreature(fitnesses, population, g, BEST) #Re-check, so as not to lose a "BEST" due to death/mutation
				
		#Kill off weaker creatures
		population, fitnesses = Death(population, fitnesses, _death_rate)
		print("Length of fitnesses:", len(fitnesses), "best living creature:", max(fitnesses) + 1)
		print("Best creature produced thus far: ", BEST)
		
         #Mutation - allow random variation amongst creatures
		population, fitnesses = Mutation(population, fitnesses, _mut_prob, _mut_mag, p_ranges, _c2constr)
	return BEST[0], BEST[1], BEST[2]
		
def InitializePopulation(i_pop, ranges, g_type, c2c):
	"""Generate initial creatures within the parameter space, via Monte Carlo or Latin Hypercube methodology.
	'ranges' is a 2 * dim matrix with the maximum acceptable values on the 1st row and the minimum values on the 2nd"""
	population = []
	for i in range(i_pop):
		population.append([random.uniform(low,high) for low,high in zip(ranges[1],ranges[0])])
	if g_type == "MonteCarlo":
		pass
	elif g_type == "LatinHypercube":	#We don't "fill the cube," but rather, ensure no 'box' has multiple creatures within
		diffs = [(r1 - r2) * 1.0 / i_pop for r1, r2 in zip(ranges[0], ranges[1])]
		dim = len(ranges[0])
		for j in range(dim):
			pop_param = []
			for p in population:  #sort all of the values of a given parameter (generates a random order)
				pop_param.append(p[j])
				Order = np.array(pop_param).ravel().argsort()
			NewRand = []
			for p in range(i_pop):  #fill the relevant boxes
				NewRand.append(ranges[1][j] + Order[p]*diffs[j] + random.uniform(0,1)*diffs[j])
			for p in range(i_pop): 
				population[p][j] = NewRand[p]
	else:
		print("method of population generation is not implemented")
		return None
	for p in range(i_pop): #Checking to ensure c2 does not exceed c1...this is to avoid unrealistic pot_evap estimates (i.e. < 0)
		if population[p][1] > population[p][0]: population[p][1] = random.uniform(0,population[p][0])
	return population, [-999999999]*len(population)

def GetFitnessE(pop, t_s, dep, g_ind, mem, precips, snipped_sms, fits, inv_seq):
	"""Using GA, maximizes the correlation bewteen the actual soil moisture time series (sm_vals) and the one generated
	by the various solutions in pop.  """
	fitnesses = []
	for i, (p, f) in enumerate(zip(pop, fits)):
		if f <= -99999: #If we don't already know the fitness...
			print("evaluating creature #", (i + 1))
			sim_seq = list(np.array(SM.AddEtaSeries(p, t_s)))							
			sim_b_series = SM.GetBeta(precips, sim_seq, dep, mem)
			snipped_b = GetValidSubset(inv_seq, g_ind, sim_b_series)
			fitnesses.append(pearsonr(np.array(ConvertToNative(snipped_b)), np.array(ConvertToNative(snipped_sms)))[0] - 1)
		else: #...and if we do.
			fitnesses.append(f) 
	return fitnesses

def ConvertToNative(vector):
	"""Given a (vector) of a non-native data type, i.e. numpy.float64, return a vector of the native type"""
	output = []
	for val in vector:
		if type(val) in [float, int]: #if this is already a native type
			output.append(val)
		else:	
			output.append(np.asscalar(val))
	return output

def GetFitnessD(pop, b_series, mem, g_ind, snipped_sms, fits, inv_seq):
	"""Using GA, minimizes the squared errors between a predicted and empirical soil moisture series"""
	fitnesses = []
	for i, (p,f) in enumerate(zip(pop, fits)):
		if f < -99999: #If we don't already know the fitness
			print("evaluating creature #", (i + 1))
			sms = SM.GetSMSeries(b_series, p)
			sm_estimates = GetValidSubset(inv_seq, g_ind, sms) 
			SSE = np.sum([(est - sm)**2 for est, sm in zip(sm_estimates, snipped_sms)])
			fitnesses.append(SSE * -1)  #This ensures fitnesses are negative and we climb up towards zero.
		else: #...and if we do.
			fitnesses.append(f)
	return fitnesses	
	
def Selection(pop, fits, mates):
	"""Given a population of possible solutions and their fitnesses, select the group that will reproduce.
	The higher the fitness, the more likely the solution is to be chosen for mating."""
	min_fit = min(fits)
	#Because fitnesses are negative, let's make 'em positive and non-zero
	fits = [f - 1.01 * min_fit for f in fits]
	highest_to_lowest = np.argsort(fits)[::-1] 
	fits = np.take(fits, highest_to_lowest)
	sorted_pop = [pop[ind] for ind in highest_to_lowest]
	selected_mates, selected_fits = [],[]
	
	for m in range(mates):
		sum_fit = sum(fits)
		Ps = [f / sum_fit for f in fits] #Normalized: sum of all Ps = 1	
		roulette_wheel = np.cumsum(Ps) #Now we can pick with a random unif from 0-1
		mate_pick = random.random()
		which_ind = (roulette_wheel > mate_pick).nonzero()
		index = min(which_ind[0]) if np.sum(which_ind) > 0 else 0
		selected_mates.append(sorted_pop[index]) #Grab the relevant creature
		selected_fits.append(fits[index]) #...and its fitness
		sorted_pop.remove(sorted_pop[index]); fits = np.delete(fits, index) #now remove for the next selection
	return selected_mates
	
def Reproduction(s_mates, type):
	"""Takes the selected mates as inputs, then returns 1/2 (rounded down) the number of creatures as 'offspring' """
	new_creatures = []
	#Now, generate some new creatures
	if type == "Two-Parent":
		while(len(s_mates)>1):
			#Pick the parents (two random from the population)
			parent_pick = np.argsort([random.random() for i in range(len(s_mates))])[0:2]
			cross_point = int(random.random() * (len(s_mates[0]) - 1)) + 1 #Choose crossover point
			#Create and store the new 'child'
			new_creatures.append(s_mates[parent_pick[0]][:cross_point] + s_mates[parent_pick[1]][cross_point:])
			#Remove the parents from the sample of mates from which we choose
			s_mates.remove(s_mates[max(parent_pick)]); s_mates.remove(s_mates[min(parent_pick)])
				
	elif type == "Multi-Parent":
		for i in range(int(len(s_mates)/2)):
			child = []
			for j in range(len(s_mates[0])):
				which_parent = int(random.random()*len(s_mates))+1
				child.append(s_mates[which_parent])
	else:
		print("This type of reproduction is not enabled")
		return None
			
	return new_creatures		

def Death(pop, fits, d_rate):
	"""Kill a fixed proportion (d_rate) of creatures, rounded_down.  
	Warning: if d_rate > (0.5*NumMates) / (0.5*NumMates + InitialPop), population will reach 0."""
	num_deaths = int(d_rate * len(pop))
	if num_deaths > 0:
		for i in range(num_deaths):
			tot_fit = sum(fits)
			p_deaths = [f * 1.0 / tot_fit for f in fits] #likelihood of death
			roulette_wheel = np.cumsum(p_deaths)
			scythe = min((roulette_wheel >= min(random.random(), .99999)).nonzero()[0]) #which creature buys it
			print("Killed creature: ", pop[scythe], " Fitness: ", fits[scythe])
			fits = np.delete(fits, scythe); pop.remove(pop[scythe])
	return pop, fits
	
def Mutation(pop, fits, mut_prob, mut_mag, ranges, c2c, beta_prime = False):
	"""Allow random mutations for every chromosome of every creature.  Mutations are Gaussian with no bias.
	Mutations cannot exit the allowable ranges. """
	for i in range(len(pop)):						
		for j in range(len(pop[0])):
			mutation_random = random.random()
			if mutation_random < mut_prob:	#Will a mutation occur?
				pop[i][j] = max(pop[i][j] + random.gauss(0, 1) * mut_mag * (ranges[0][j] - ranges[1][j]),0)
				fits[i] = -999999999 #re-assign so that it will be calculated...
				if c2c and not beta_prime and pop[i][0] <= pop[i][1]:
					pop[i][1] = pop[i][0] - 0.0001
	return pop, fits
					

		
	
	
	
	
