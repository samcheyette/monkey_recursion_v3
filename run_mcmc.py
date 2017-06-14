
from mcmc_helpers import *
from numpy.random import exponential, dirichlet, normal
from numpy import mean, std, arange, multiply, array
from scipy.stats import expon
import time

def propose(alpha, beta):
	log_alpha_new = math.log(alpha) + normal(0,0.1)
	alpha_new = math.exp(log_alpha_new)
 	#alpha_new = max(1e-3, alpha_new)

	beta_new = []

	#multiply: keeps mean of dirichlet but
	# decreases distance of jumps
	beta_new = multiply(200., beta)
	beta_new = dirichlet(beta_new)

	return alpha_new, beta_new


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


	#dirichlet hyper-parameters 
	#alpha (int) ~ exponential(lambda)
	#beta (vector) ~ dirichlet(1)
	#theta_i (vector) ~ dirichlet(alpha*beta)
	#strategies (vector) | n_i (int) ~ multinomial(theta_i)
	#responses (vector) | strategies (vector) ~ multinomial(theta_strats)
	
	lmbda = 0.1
	alpha = 100#1/lmbda
	beta = dirichlet([1.0 for _ in xrange(len(hyps_used))])

	
	#NEED TO FIT THETAS FOR EACH PARTICIPANT!!!!
	acc_tot = 0

	dct_beta = {}
	norm_const = 0.0
	all_alpha = 0.0

	for chain in xrange(MCMC_CHAINS):
		t_0 = time.time()
		steps = 0
		posterior = None
		prior = None
		likelihood = None
		while steps < MCMC_STEPS:
			prop_ab = propose(alpha, beta)
			prop_alpha = prop_ab[0]
			prop_beta = prop_ab[1]

			alph_bet_prop = multiply(prop_alpha, prop_beta)
			posterior_prop = 0.0
			lkhd_prop = 0.0
			prior_prop = 0.0
			nsamp = 25

			#sample from dirichlet(alpha * beta) to estimate
			#likelihood 
			for _ in xrange(nsamp):
				#need to sample from dirichlet about alpha * <beta>
				#(can't just use point estimate)
				#to get accurate uncertainty 

				theta_prop = dirichlet(alph_bet_prop)
				theta_resp_prop = flatten_theta(theta_prop,
									out_hyps, hyps_used)

				#for all participants in a population
				#find how likely 
				for data in all_data:

					lkhd_prop += prob_d_given_h(theta_resp_prop, 
								data, sm=1e-3)
				
					#what else needs to be in the prior?
					prior_prop += math.log(expon.pdf(prop_alpha, 
										scale=1.0/lmbda)) 

			posterior_prop = (prior_prop + lkhd_prop - 
							nsamp * len(all_data) * 
								math.log(nsamp * len(all_data)))

			rand_acc =  ACC_TEMPERATURE * math.log(random.random())
			#choose whether to accept or reject proposal 
			#(metropolis hastings)
			if (posterior == None or
				 posterior_prop > posterior or
				 ((rand_acc) < 
				 	min(0.0, 
				 		posterior_prop - posterior))):
				alpha = copy.deepcopy(prop_alpha)
				beta = copy.deepcopy(prop_beta)
				posterior = posterior_prop
				prior = prior_prop
				likelihood = lkhd_prop
				acc_tot += 1


			#just print some stuff every 2 seconds
			hs = []
			ps = []
			if int(time.time() - t_0) >= 2:
				print "WHO:      ", who
				print "CHAIN:    ", chain + 1
				print "STEPS:    ", steps
				print "ACC_PROP: ", acc_tot/float(MCMC_STEPS * chain + 
								steps+1)
				print "ALPHA:    ", alpha
				print "POSTERIOR:", posterior
				print "LIKELHOOD:", likelihood
				print "PRIOR:    ", prior
				for b in xrange(len(beta)):
					if beta[b] > 0.005:
						hs.append((hyps_used[b], beta[b]))
				srt_b = sorted(hs, key=lambda tup: -tup[1])
				print "TOP 8:"
				for x in srt_b[:8]:
					print x[0], x[1]
				print
				t_0 = time.time()


			if steps > PCT_BURNIN * MCMC_STEPS:
				for b in xrange(len(beta)):
					if hyps_used[b] not in dct_beta:
						dct_beta[hyps_used[b]] = 0.0
					dct_beta[hyps_used[b]] += beta[b]
				norm_const += 1
				all_alpha += alpha


			steps += 1

	for b in dct_beta:
		dct_beta[b] = dct_beta[b]/float(norm_const)
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


	return all_alpha, dct_beta, all_data, out_hyps



if __name__ == "__main__":

	###MCMC PARAMS ####
	MCMC_STEPS = 10000
	MCMC_CHAINS = 2
	PCT_BURNIN=0.05
	ACC_TEMPERATURE = 1.0

	####NOISE PARAMS###
	#(filter_thresh=0.5 fully filters out all invalid strategies)
	filter_thresh = 0.5 

	out_file = "mixture_output/output.csv"
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
			file = "stevesdata/RecursionMonkey.csv"
			start = 100
			subset={"Exposure": "2"}

		else:
			file =  "stevesdata/RecursionTsimane.csv"
			start = 1000
			subset= {}


		model_res = run_MCMC()
		alpha_w = model_res[0]
		beta_w = model_res[1]
		data_w = model_res[2]
		out_hyps = model_res[3]

		all_data_who.append((who, alpha_w,
					 copy.deepcopy(beta_w),
						copy.deepcopy(data_w)))

	output_full_model(all_data_who, out_file, out_hyps)