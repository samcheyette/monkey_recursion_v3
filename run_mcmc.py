
from mcmc_helpers import *
from numpy.random import exponential, dirichlet, normal
from numpy import mean, std, arange, multiply, array, add
from scipy.stats import expon, pearsonr
from scipy.stats import dirichlet as dirpdf 
import time

def propose(alpha, beta):
	log_alpha_new = math.log(alpha) + normal(0,0.1)
	alpha_new = math.exp(log_alpha_new)
 	#alpha_new = max(1e-3, alpha_new)

	beta_new = []

	#multiply: keeps mean of dirichlet but
	# decreases distance of jumps
	beta_new = multiply(200, beta)

	beta_new = dirichlet(beta_new)
	return alpha_new, beta_new


def flatten_data(all_data):
	flt_data = {}
	sumd = 0.0
	for d in all_data:
		sum_part = float(sum([d[k] for k in d]))
		for k in d:
			if k not in flt_data:
				flt_data[k] = 0.0
			flt_data[k] += d[k]/sum_part
			sumd += d[k]/sum_part

	flt_data_lst = []
	for d in flt_data:
		flt_data_lst.append((d, flt_data[d]/sumd))

	flt_data_lst = sorted(flt_data_lst, key=lambda tup: -tup[1])



	return flt_data_lst

def run_MCMC():

	hypotheses = get_hyps(types, 4)
	out_hyps = get_hyps_gen(hypotheses)
	out_hyps = filter_hyps(out_hyps, 
					thresh=filter_thresh)

	out_hyps = remove_dups(out_hyps)
	hyps_used = copy.copy(out_hyps.keys())

	#need to preprocess the data from steve's files
	all_data_bef = getCountData(file, careAbout, who, subset)
	all_data = []
	for pers in all_data_bef:
		new_dct = {}
		for key in pers:
			new_dct["".join(list(key))] = pers[key]
		all_data.append(copy.deepcopy(new_dct))



	flt_data_lst = flatten_data(all_data)


	#dirichlet hyper-parameters 
	#alpha (int) ~ exponential(lambda)
	#beta (vector) ~ dirichlet(1)
	#theta_i (vector) ~ dirichlet(alpha*beta)
	#strategies (vector) | n_i (int) ~ multinomial(theta_i)
	#responses (vector) | strategies (vector) ~ multinomial(theta_strats)
	
	lmbda = 0.02
	
	beta_ones = [1.0/len(hyps_used) for _ in xrange(len(hyps_used))]
	
	#NEED TO FIT THETAS FOR EACH PARTICIPANT!!!!
	acc_tot = 0
	time_orig = time.time()
	dct_beta = {}
	all_betas = {}
	norm_const = 0.0
	all_alpha = 0.0
	baseline_const = None
	ent = []

	for chain in xrange(MCMC_CHAINS):
		alpha = 1/(0.5 * lmbda)
		beta = dirichlet([1.0 for _ in xrange(len(hyps_used))])
		t_0 = time.time()
		steps = 0
		posterior = None
		prior = None
		likelihood = None

		while steps < MCMC_STEPS:
			prop_ab = propose(alpha, beta)
			prop_alpha = prop_ab[0]
			prop_beta = prop_ab[1]

			norm_bet = [max(prb, 1e-3) 
					for prb in prop_beta]

			pr_alpha = math.log(expon.pdf(prop_alpha, 
										scale=1.0/lmbda))
			
			pr_beta = 0.0
			pr_beta = dirpdf.logpdf(beta_ones, normalise(norm_bet)) 


			#check entropy of pr_beta (test)
			tmp_ent = 0.0
			for b in beta:
				if b > 0.0:
					tmp_ent += -b * math.log(b, 2.0)

			ent.append(tmp_ent)

			#print pr_beta

			alph_bet_prop = multiply(prop_alpha, prop_beta)

			posterior_prop = 0.0
			nsamp = 10
			draw = nsamp * 10
			thetas_prop = dirichlet(alph_bet_prop, draw)
			likepropsum = 0.0
			priorpropsum = 0.0
	
			#sample from dirichlet(alpha * beta) to estimate
			#likelihood 
			for samp in xrange(nsamp):
				#need to sample from dirichlet about alpha * <beta>
				#(can't just use point estimate)
				#to get accurate uncertainty 

				lkhd_prop = 0.0
				prior_prop = 0.0
				prob_d = 0.0

				#for all participants in a population
				#find how likely 

				for data in all_data:

					theta_prop = thetas_prop[random.randint(0,draw-1)]
					theta_resp_prop = flatten_theta(theta_prop,
										out_hyps, hyps_used)
			
					p_d_given_h = prob_d_given_h(theta_resp_prop, 
								data, sm=1e-2)


					lkhd_prop += p_d_given_h
				
					#what else needs to be in the prior?


				###DO THIS####
				#the constant is a very hacky way to protect
				#against overflow.... need to find a better
				#way to do this....
				if baseline_const == None:
					constant=600 #max exp() before overflow

					baseline_const = -(lkhd_prop + constant)
					baseline_const -= pr_beta + pr_alpha


 				posterior_prop += (math.exp(lkhd_prop + pr_beta + pr_alpha +
 								 prior_prop + baseline_const))
 				likepropsum += lkhd_prop 
 				priorpropsum += pr_beta + pr_alpha


			rand_acc =  ACC_TEMPERATURE * random.random()
			#choose whether to accept or reject proposal 
			#(metropolis hastings)
			if ((posterior == None or posterior == 0.0 or
				 posterior_prop > posterior or
				 ((rand_acc) <
				 	min(1.0, 
				 		posterior_prop/posterior))) and 
				 posterior_prop != 0.0):
				alpha = copy.deepcopy(prop_alpha)
				beta = copy.deepcopy(prop_beta)
				posterior = posterior_prop
				prior = priorpropsum
				likelihood = likepropsum
				acc_tot += 1



			#print posterior_prop/posterior


			#####################################
			#######JUST PRINT SOME STUFF#########

			flt_beta = flatten_theta(beta, out_hyps, hyps_used)
			if int(time.time() - t_0) >= 1:
				hs = []
				ps = []
				tot_steps =float(MCMC_STEPS * chain + 
								steps+1)
				print "WHO:      ", who
				print "CHAIN:    ", chain + 1
				print "TIME PER: ", (time.time() - time_orig)/tot_steps
				print "STEPS:    ", steps
				print "ACC_PROP: ", acc_tot/tot_steps
				print "ALPHA:    ", alpha
				print "POSTERIOR:", posterior
				print "LIKELHOOD:", likelihood
				print "PRIOR:    ", prior
				print "ENTROPY:  ", sum(ent)/float(len(ent))

				for b in xrange(len(beta)):
					if beta[b] > 0.005:
						hs.append((hyps_used[b], beta[b]))
				srt_b = sorted(hs, key=lambda tup: -tup[1])
				print "TOP 10 HYPS:"
				for x in srt_b[:10]:
					print x[0], x[1]


				print "..."
				print "DATA v MOD:  "
				for d in flt_data_lst[:10]:
					if d[0] not in flt_beta:
						betad=0.0
					else:
						betad=flt_beta[d[0]]
					print d[0], d[1], betad

				print ">" * 50
				print
				t_0 = time.time()
			#############################################
			############################################

			if steps > PCT_BURNIN * MCMC_STEPS:
				for b in xrange(len(beta)):
					if hyps_used[b] not in dct_beta:
						dct_beta[hyps_used[b]] = 0.0
					dct_beta[hyps_used[b]] += beta[b]
				norm_const += 1
				all_alpha += alpha

				for b in flt_beta:
					if b not in all_betas:
						all_betas[b] = 0.0
					all_betas[b] += flt_beta[b]

			steps += 1

	for b in dct_beta:
		dct_beta[b] = dct_beta[b]/float(norm_const)

	for b in all_betas:
		all_betas[b] = all_betas[b]/float(norm_const)
	all_alpha = all_alpha/float(norm_const)


	#print some more stuff
	lst_beta = []
	for b in dct_beta.keys():
		lst_beta.append((b, dct_beta[b]))
	lst_beta = sorted(lst_beta, key=lambda tup: -tup[1])
	print "TOP ALPHA: ", all_alpha
	print "TOP H FINAL: "
	for b in lst_beta[:10]:
		print b[0], b[1]
	print


	return all_alpha, dct_beta, flt_data_lst, out_hyps, all_betas



if __name__ == "__main__":

	###MCMC PARAMS ####
	MCMC_STEPS = 3000
	MCMC_CHAINS = 5
	PCT_BURNIN= 0.5
	ACC_TEMPERATURE = 1.0


	####NOISE PARAMS###
	#(filter_thresh=0.5 fully filters out all invalid strategies)
	filter_thresh = 0.5 

	out_file = "mixture_output/output2.csv"
	whos = ["Kids", "Monkeys", "Tsimane"]
	monks = ["Monkeys"]
	tsim =["Tsimane"]
	kids = ["Kids"]
	all_mods_traces = []
	all_mods_distr = []
	all_data_who = []
	hyps_used = []


	for who in whos:
		careAbout = "Order pressed"
		if who == "Kids":
			file = "stevesdata/RecursionKids.csv"
			start = 0 
			subset= {}

		elif who == "Monkeys":
			MCMC_STEPS *= 4
			file = "stevesdata/RecursionMonkey.csv"
			start = 100
			subset={"Exposure": "2"}

		else:
			MCMC_STEPS /= 4
			file =  "stevesdata/RecursionTsimane.csv"
			start = 1000
			subset= {}


		model_res = run_MCMC()
		alpha_w = model_res[0]
		beta_w = model_res[1]
		data_w = model_res[2]
		out_hyps = model_res[3]
		thetas_w = model_res[4]
		for k in data_w[:10]:
			thet_p = 0.0
			if k[0] in thetas_w:
				thet_p = thetas_w[k[0]]
			print k[0], k[1], thet_p



		all_data_who.append((who, alpha_w,
					 copy.deepcopy(beta_w),
						copy.deepcopy(data_w),
						 copy.deepcopy(thetas_w)))

	output_full_model(all_data_who, out_file, out_hyps)